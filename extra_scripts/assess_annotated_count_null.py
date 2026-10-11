#!/usr/bin/env python3
"""Count-level null calibration of block-local tests over annotated blocks.

Random full-data tests from the run_local_path_tests.py family (label-blind
order, no significance or LR selection) keep their retained paths, the
precursor target and their observed subject, cell-type and primer molecule
totals. Class counts are redrawn with no cell-type effect on mature path
usage (tealeaf.sc.local_path_reads.simulate_null_class_counts) and each draw
is fitted with and without the precursor component. The "precursor shift"
scenario multiplies the precursor fraction of level-b rows by --precursor-shift,
a cell-type difference in unspliced RNA that leaves mature ratios unchanged.
Failed fits stay as p = 1.
"""

import argparse
import gzip
import json
from pathlib import Path
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_local_path_tests import declared_tests
from tealeaf.sc.local_path_reads import path_read_opportunities, pooled_path_shares, project_mask, project_opportunities, simulate_null_class_counts
from tealeaf.sc.local_path_test import PRIMERS, fit_local_block, prepare_block_test

READ_LENGTH = 151
SCENARIOS = {"subject 20, row 20": (20., 20., False), "subject 20, no row variation": (20., None, False), "subject 20, row 5": (20., 5., False), "precursor shift": (20., 20., True)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--subject-folds", type=Path, required=True)
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--lr-cell-types", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hypotheses", type=int, default=200)
    parser.add_argument("--draws", type=int, default=2)
    parser.add_argument("--precursor-shift", type=float, default=4.)
    parser.add_argument("--stored-opportunities", action="store_true")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--collect", action="store_true")
    args = parser.parse_args()
    if args.collect:
        table = pd.concat([pd.read_csv(path, sep="\t") for path in sorted(args.output_dir.glob("shard_*/null.tsv"))], ignore_index=True)
        if len(table) != 2 * args.hypotheses * args.draws * len(SCENARIOS):
            raise ValueError(f"incomplete null family: {len(table)} fits")
        table["multi_path"] = table.n_paths > 2
        summary = table.groupby(["scenario", "model", "multi_path"]).apply(lambda local: pd.Series({"trials": len(local), "converged": int(local.converged.sum()), **{f"reject_{level:g}": float((local.p_value < level).mean()) for level in (.05, .01, .001)}}), include_groups=False).reset_index()
        summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
        table.to_csv(args.output_dir / "trials.tsv.gz", sep="\t", index=False, na_rep="NA")
        print(summary.to_string(index=False), flush=True)
        return
    folds = pd.read_csv(args.subject_folds, sep="\t", dtype={"subject": str})
    counts = pd.read_csv(args.local_reads / "counts.tsv.gz", sep="\t")
    counts = counts.loc[counts.subject.astype(str).isin(set(folds.subject))]
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    stored = json.load(gzip.open(args.local_reads / "opportunities.json.gz", "rt")) if args.stored_opportunities else None
    universe = pd.read_csv(args.universe, sep="\t", usecols=["gene_id", "pair_id"]).drop_duplicates()
    pairs = {(gene.split(".")[0], tuple(sorted(pair.split("||")))) for gene, pair in zip(universe.gene_id, universe.pair_id)}
    lr_types = set(args.lr_cell_types.split(";"))
    tests = declared_tests(counts, blocks, set(folds.subject), 4, 20., lambda gene, a, b: (gene.split(".")[0], (a, b)) in pairs or (a in lr_types and b in lr_types))
    grouped = {}
    for row in counts.itertuples(index=False):
        grouped.setdefault(row.path_key, []).append((str(row.subject), row.cell_type, row.primer, int(row.mask), float(row.count)))
    chosen = []
    for index in np.random.default_rng(20261010).permutation(len(tests)):
        key, level_a, level_b = tests[index]
        levels = (level_a, level_b)
        n_paths = len(blocks[key]["signatures"])
        full = {int(mask): np.asarray(vector, dtype=float) for mask, vector in (stored[key] if stored is not None else path_read_opportunities(blocks[key]["paths"], READ_LENGTH)).items()}
        try:
            prepared = prepare_block_test(key, [entry for entry in grouped[key] if entry[1] in levels], full, n_paths, levels, precursor=True)
        except ValueError:
            continue
        chosen.append((key, levels, prepared))
        if len(chosen) == args.hypotheses:
            break
    records = []
    for key, levels, prepared in chosen[args.shard_index::args.shard_count]:
        size = len(prepared["kept"])
        test_id = f"{blocks[key]['block_id']}|cell_type|{levels[0]}|{levels[1]}"
        opportunities, lookup = prepared["opportunities"], prepared["lookup"]
        rows = list(zip(prepared["subjects"], prepared["labels"].tolist()))
        totals = {(row, primer): sum(lookup.get((key, subject, levels[label], primer, mask), 0.) for mask in opportunities) for row, (subject, label) in enumerate(rows) for primer in PRIMERS}
        pooled = {}
        for (_, _, _, _, mask), value in lookup.items():
            pooled[mask] = pooled.get(mask, 0.) + value
        shares = pooled_path_shares(pooled, opportunities)
        mature_opportunities = project_opportunities(opportunities, list(range(size)), size, False)
        for scenario, (subject_concentration, row_concentration, shift) in SCENARIOS.items():
            scales = {row: np.r_[np.ones(size), args.precursor_shift if label == 1 else 1.] for row, (_, label) in enumerate(rows)} if shift else None
            for draw in range(args.draws):
                rng = np.random.default_rng(zlib.crc32(f"{test_id}|{scenario}|{draw}".encode()))
                simulated = simulate_null_class_counts(totals, opportunities, shares, [subject for subject, _ in rows], rng, subject_concentration=subject_concentration, row_concentration=row_concentration, row_scales=scales)
                with_precursor, mature = {}, {}
                for (row, primer, mask), value in simulated.items():
                    subject, label = rows[row]
                    with_precursor[(key, subject, levels[label], primer, mask)] = float(value)
                    new = project_mask(mask, list(range(size)), size, False)
                    if new:
                        mature[(key, subject, levels[label], primer, new)] = mature.get((key, subject, levels[label], primer, new), 0.) + float(value)
                for model, local_lookup, local_opportunities, precursor in (("anchors", mature, mature_opportunities, False), ("anchors + precursor", with_precursor, opportunities, True)):
                    record = {"test_id": test_id, "n_paths": size, "scenario": scenario, "draw": draw, "model": model, "p_value": 1., "converged": False, "error": ""}
                    try:
                        fields, _ = fit_local_block(local_lookup, key, local_opportunities, size, prepared["subjects"], prepared["labels"], levels, anchors=True, primer_offset=True, precursor=precursor)
                        record.update(p_value=fields["p_value"], converged=fields["converged"], statistic=fields["statistic"], n_converged_path_tests=fields["n_converged_path_tests"])
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exception:
                        record["error"] = str(exception)
                    records.append(record)
    out = args.output_dir / f"shard_{args.shard_index}"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out / "null.tsv", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
