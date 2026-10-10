#!/usr/bin/env python3
"""Block-local path tests for pairwise cell-type contrasts, from read classes.

Tests are defined from block-local molecule counts (collect_local_path_reads.py
--annotated) and subject folds only, without EC caches: every annotated block
and cell-type pair with at least --min-subjects subjects observed in both types
and at least --min-molecules pooled molecules. Within each test, mature paths
whose label-blind pooled share is below --min-share are dropped and read
classes are projected onto the retained paths. Output follows the
run_paired_path_test.py shard layout for the split and long-read assessments.
"""

import argparse
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from tealeaf.sc.junction_benchmark import benjamini_hochberg
from tealeaf.sc.local_path_reads import path_read_opportunities, pooled_path_shares
from tealeaf.sc.local_path_test import PRIMERS, fit_local_block

READ_LENGTH = 151


def project_mask(mask, kept, n_paths, precursor):
    """Old mask over n_paths mature paths plus precursor -> retained-path mask."""
    value = sum(((mask >> old) & 1) << new for new, old in enumerate(kept))
    if precursor:
        value |= ((mask >> n_paths) & 1) << len(kept)
    return value


def declared_tests(counts, blocks, subjects, min_subjects, min_molecules, allowed=None):
    """Sorted (path key, level a, level b) tests passing the label-blind screen.

    allowed, if given, is a function (gene_id, level a, level b) -> bool that
    restricts the contrasts that are evaluated by either endpoint.
    """
    local = counts.loc[counts.subject.isin(subjects)]
    totals = local.groupby(["path_key", "subject", "cell_type"])["count"].sum()
    tests = []
    for key, frame in totals.groupby(level="path_key"):
        table = frame.droplevel("path_key").unstack("cell_type", fill_value=0)
        types = sorted(table.columns)
        for first in range(len(types)):
            for second in range(first + 1, len(types)):
                a, b = types[first], types[second]
                if allowed is not None and not allowed(blocks[key]["gene_id"], a, b):
                    continue
                paired = ((table[a] > 0) & (table[b] > 0)).sum()
                if paired >= min_subjects and table[[a, b]].to_numpy().sum() >= min_molecules:
                    tests.append((key, a, b))
    return tests


def merge(root, shard_count):
    paths = [root / f"shard_{index}" / "paired_path.tsv" for index in range(shard_count)]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"incomplete cohort, missing {missing[:3]}")
    table = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    table["fdr"] = benjamini_hochberg(table.p_value.to_numpy())
    (root / "merged").mkdir(exist_ok=True)
    table.to_csv(root / "merged" / "paired_path.tsv", sep="\t", index=False, na_rep="NA")
    print(f"{root}: {len(table)} tests, {int(table.converged.sum())} converged, {int((table.fdr < .05).sum())} BH", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--subject-folds", type=Path, required=True)
    parser.add_argument("--cohort", choices=("0", "1", "full"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=32)
    parser.add_argument("--min-subjects", type=int, default=4)
    parser.add_argument("--min-molecules", type=float, default=20.)
    parser.add_argument("--min-share", type=float, default=.02)
    parser.add_argument("--no-anchors", action="store_true")
    parser.add_argument("--precursor", action="store_true")
    parser.add_argument("--primer-offset", action="store_true")
    parser.add_argument("--nodes", type=int, default=9)
    parser.add_argument("--count-only", action="store_true")
    parser.add_argument("--universe", type=Path, help="matched split universes (gene_id, pair_id 'A||B'); their contrasts are tested")
    parser.add_argument("--lr-cell-types", help="semicolon-separated cell types with a long-read mapping; all their pairs are tested")
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    if args.merge:
        merge(args.output_dir, args.shard_count)
        return
    folds = pd.read_csv(args.subject_folds, sep="\t", dtype={"subject": str})
    subjects = set(folds.subject if args.cohort == "full" else folds.loc[folds.fold.eq(int(args.cohort)), "subject"])
    counts = pd.read_csv(args.local_reads / "counts.tsv.gz", sep="\t")
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    allowed = None
    if args.universe is not None or args.lr_cell_types:
        pairs = set()
        if args.universe is not None:
            universe = pd.read_csv(args.universe, sep="\t", usecols=["gene_id", "pair_id"]).drop_duplicates()
            pairs = {(gene.split(".")[0], tuple(sorted(pair.split("||")))) for gene, pair in zip(universe.gene_id, universe.pair_id)}
        lr_types = set(args.lr_cell_types.split(";")) if args.lr_cell_types else set()
        allowed = lambda gene, a, b: (gene.split(".")[0], (a, b)) in pairs or (a in lr_types and b in lr_types)
    tests = declared_tests(counts, blocks, subjects, args.min_subjects, args.min_molecules, allowed)
    if args.count_only:
        print(f"cohort {args.cohort}: {len(tests)} tests, {len({key for key, _, _ in tests})} blocks", flush=True)
        return
    tests = tests[args.shard_index::args.shard_count]
    keys = {key for key, _, _ in tests}
    local = counts.loc[counts.path_key.isin(keys) & counts.subject.isin(subjects)]
    lookup = {(row.path_key, row.subject, row.cell_type, row.primer, int(row.mask)): float(row.count) for row in local.itertuples(index=False)}
    grouped = {}
    for (path, subject, cell_type, primer, mask), value in lookup.items():
        grouped.setdefault(path, []).append((subject, cell_type, primer, mask, value))
    rows_out, usage_out, opportunity_cache = [], [], {}
    started = time.monotonic()
    for key, level_a, level_b in tests:
        block = blocks[key]
        n_paths = len(block["signatures"])
        test_id = f"{block['block_id']}|cell_type|{level_a}|{level_b}"
        levels = (level_a, level_b)
        entries = [entry for entry in grouped[key] if entry[1] in levels]
        present = sorted({(subject, levels.index(cell_type)) for subject, cell_type, _, _, value in entries if value > 0})
        subject_array, labels = np.array([subject for subject, _ in present]), np.array([label for _, label in present])
        record = {"test_id": test_id, "block_id": block["block_id"], "gene_id": block["gene_id"], "contrast": "cell_type_pairwise", "level_a": level_a, "level_b": level_b, "method": "block-local path classes, random-subject Beta binomial", "n_annotated_paths": n_paths, "n_paths": 0, "path_signatures": "[]", "n_subjects": 0, "median_gene_umis": np.nan, "statistic": 0., "p_value": 1., "raw_p_value": 1., "converged": False, "mean_difference": "[]", "mean_difference_norm": np.nan, "error": ""}
        begin = time.monotonic()
        try:
            if (key, n_paths) not in opportunity_cache:
                opportunity_cache[(key, n_paths)] = path_read_opportunities(block["paths"], READ_LENGTH)
            full_opportunities = opportunity_cache[(key, n_paths)]
            pooled = {mask: sum(lookup.get((key, subject, levels[label], primer, mask), 0.) for subject, label in present for primer in PRIMERS) for mask in full_opportunities}
            shares = pooled_path_shares(pooled, full_opportunities)[:n_paths]
            shares = shares / shares.sum()
            kept = [index for index in range(n_paths) if shares[index] >= args.min_share]
            if len(kept) < 2:
                raise ValueError("fewer than two expressed paths")
            paths = [block["paths"][index] for index in kept] + ([block["paths"][n_paths]] if args.precursor else [])
            cache_key = (key, tuple(kept), args.precursor)
            if cache_key not in opportunity_cache:
                opportunity_cache[cache_key] = path_read_opportunities(paths, READ_LENGTH)
            opportunities = opportunity_cache[cache_key]
            projected = {}
            row_totals = {}
            for subject, cell_type, primer, mask, value in entries:
                new = project_mask(mask, kept, n_paths, args.precursor)
                if new:
                    projected[(key, subject, cell_type, primer, new)] = projected.get((key, subject, cell_type, primer, new), 0.) + value
                    row_totals[(subject, cell_type)] = row_totals.get((subject, cell_type), 0.) + value
            totals = [row_totals.get((subject, levels[label]), 0.) for subject, label in present]
            signatures = [block["signatures"][index] for index in kept]
            record.update(n_paths=len(kept), path_signatures=json.dumps(signatures), median_gene_umis=float(np.median(totals)), kept_paths=json.dumps(kept), annotated_path_shares=json.dumps(shares.tolist()))
            fields, usage = fit_local_block(projected, key, opportunities, len(kept), subject_array, labels, levels, anchors=not args.no_anchors, primer_offset=args.primer_offset, nodes=args.nodes, precursor=args.precursor)
            record.update({name: json.dumps(value) if isinstance(value, list) else value for name, value in fields.items()})
            record["raw_p_value"] = record["p_value"]
            for target, level, proportion in usage:
                usage_out.append({"test_id": test_id, "block_id": block["block_id"], "gene_id": block["gene_id"], "subject": "model", "cell_type": levels[level], "path": f"Path {target + 1}", "path_number": target + 1, "path_signature": json.dumps(signatures[target]), "proportion": proportion})
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exception:
            record["error"] = str(exception)
        record["runtime_seconds"] = time.monotonic() - begin
        rows_out.append(record)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows_out).to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(usage_out, columns=["test_id", "block_id", "gene_id", "subject", "cell_type", "path", "path_number", "path_signature", "proportion"]).to_csv(args.output_dir / "path_usage.tsv", sep="\t", index=False, na_rep="NA")
    settings = {"subject_fold": None if args.cohort == "full" else int(args.cohort), "test_effect": "cell_type_pairwise", "min_subjects": args.min_subjects, "min_molecules": args.min_molecules, "min_share": args.min_share, "anchors": not args.no_anchors, "precursor": args.precursor, "primer_offset": args.primer_offset, "local_reads": str(args.local_reads), "read_length": READ_LENGTH}
    summary = {"candidates": len(rows_out), "converged": int(sum(row["converged"] for row in rows_out)), "elapsed_seconds": time.monotonic() - started, "candidate_settings": settings}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "candidate_settings"}), flush=True)


if __name__ == "__main__":
    main()
