#!/usr/bin/env python3
"""Coverage-only null p-values of the small-p tests of a run_local_path_tests.py run.

Every test of a merged run with p below --threshold is prepared as in that run (settings from shard_0/summary.json:
min_share, structure, junction grouping, anchors, primer offset, precursor), its class counts are redrawn with each
row's observed molecule totals per primer and one set of path shares pooled over both cell types (no subject or row
variation, tealeaf.sc.local_path_reads.simulate_null_class_counts), and it is refitted with the same model. A test
whose median null p-value over the draws is below .05 is coverage-driven (the rule applied to every method in the
evaluation). Failed null fits are p = 1. Writes shard_<i>/null.tsv; --merge concatenates them into null.tsv.gz.
"""

import argparse
import gzip
import json
from pathlib import Path
import zlib

import numpy as np
import pandas as pd

from tealeaf.sc.local_path_reads import pooled_path_shares, same_junction_groups, simulate_null_class_counts
from tealeaf.sc.local_path_test import PRIMERS, fit_local_block, prepare_block_test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-dir", type=Path, required=True, help="run_local_path_tests.py output (merged/paired_path.tsv, shard_0/summary.json)")
    parser.add_argument("--subject-folds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=.01)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--nodes", type=int, default=9)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    if args.merge:
        paths = [args.output_dir / f"shard_{index}/null.tsv" for index in range(args.shard_count)]
        missing = [str(path) for path in paths if not path.exists()]
        if missing:
            raise FileNotFoundError(f"missing shards: {missing}")
        table = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
        table.to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False, na_rep="NA")
        nulls = table[[f"p_null{draw}" for draw in range(args.draws)]]
        print(f"{len(table)} tests, median null p < .05: {int((nulls.median(axis=1) < .05).sum())}, per-draw null p < .05: {float((nulls < .05).mean().mean()):.4f}", flush=True)
        return
    settings = json.loads((args.fit_dir / "shard_0/summary.json").read_text())["candidate_settings"]
    fits = pd.read_csv(args.fit_dir / "merged/paired_path.tsv", sep="\t", usecols=["test_id", "block_id", "level_a", "level_b", "p_value"])
    selected = fits.loc[fits.p_value < args.threshold].sort_values("test_id").reset_index(drop=True).iloc[args.shard_index::args.shard_count]
    folds = pd.read_csv(args.subject_folds, sep="\t", dtype={"subject": str})
    fold = settings["subject_fold"]
    subjects = set(folds.subject if fold is None else folds.loc[folds.fold.eq(int(fold)), "subject"])
    local_reads = Path(settings["local_reads"])
    blocks = json.load(gzip.open(local_reads / "blocks.json.gz", "rt"))
    key_of = {block["block_id"]: key for key, block in blocks.items()}
    keys = {key_of[block] for block in selected.block_id}
    stored = json.load(gzip.open(local_reads / "opportunities.json.gz", "rt"))
    counts = pd.concat(chunk.loc[chunk.path_key.isin(keys) & chunk.subject.isin(subjects)] for chunk in pd.read_csv(local_reads / "counts.tsv.gz", sep="\t", dtype={"subject": str}, chunksize=2_000_000))
    grouped = {}
    for row in counts.itertuples(index=False):
        grouped.setdefault(row.path_key, []).append((row.subject, row.cell_type, row.primer, int(row.mask), float(row.count)))
    records = []
    for test in selected.itertuples(index=False):
        key = key_of[test.block_id]
        block = blocks[key]
        n_paths = len(block["signatures"])
        levels = (test.level_a, test.level_b)
        record = {"test_id": test.test_id, "p_value": test.p_value, "error": ""}
        try:
            full = {int(mask): np.asarray(vector, dtype=float) for mask, vector in stored[key].items()}
            groups = same_junction_groups(block["paths"][:n_paths]) if settings["group_same_junctions"] else None
            entries = [entry for entry in grouped[key] if entry[1] in levels]
            prepared = prepare_block_test(key, entries, full, n_paths, levels, settings["min_share"], settings["precursor"], groups, settings["min_subjects"], settings["min_molecules"], None, settings.get("unpaired", False))
        except (ValueError, KeyError) as exception:
            record["error"] = f"prepare: {exception}"
            records.append(record | {f"p_null{draw}": 1. for draw in range(args.draws)})
            continue
        opportunities, lookup = prepared["opportunities"], prepared["lookup"]
        rows = list(zip(prepared["subjects"], prepared["labels"].tolist()))
        totals = {(row, primer): sum(lookup.get((key, subject, levels[label], primer, mask), 0.) for mask in opportunities) for row, (subject, label) in enumerate(rows) for primer in PRIMERS}
        pooled = {}
        for (_, _, _, _, mask), value in lookup.items():
            pooled[mask] = pooled.get(mask, 0.) + value
        shares = pooled_path_shares(pooled, opportunities)
        for draw in range(args.draws):
            rng = np.random.default_rng(zlib.crc32(f"{test.test_id}|coverage null|{draw}".encode()))
            simulated = simulate_null_class_counts(totals, opportunities, shares, [subject for subject, _ in rows], rng, subject_concentration=None, row_concentration=None)
            null_lookup = {(key, rows[row][0], levels[rows[row][1]], primer, mask): float(value) for (row, primer, mask), value in simulated.items()}
            p_null, converged = 1., False
            try:
                fields, _ = fit_local_block(null_lookup, key, opportunities, len(prepared["kept"]), prepared["subjects"], prepared["labels"], levels, anchors=settings["anchors"], primer_offset=settings["primer_offset"], nodes=args.nodes, precursor=settings["precursor"], primer_opportunities=prepared["primer_opportunities"], unpaired=settings.get("unpaired", False), primer_offset_prior_sd=settings.get("primer_offset_prior_sd"))
                p_null, converged = float(fields["p_value"]), bool(fields["converged"])
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exception:
                record["error"] = str(exception)
            record[f"p_null{draw}"], record[f"converged_null{draw}"] = p_null, converged
        records.append(record)
    out = args.output_dir / f"shard_{args.shard_index}"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out / "null.tsv", sep="\t", index=False, na_rep="NA")
    print(f"shard {args.shard_index}: {len(records)} tests", flush=True)


if __name__ == "__main__":
    main()
