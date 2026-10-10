#!/usr/bin/env python3
"""Count-level null calibration of the block-local path tests.

Random full-data pairwise contrasts with block-local molecules (no
significance or LR selection) keep their observed subject, cell-type and
primer molecule totals; class counts are redrawn with no cell-type effect
(tealeaf.sc.local_path_reads.simulate_null_class_counts) and refitted with the
same fit_local_block call as the real-data driver. Failed fits stay as p = 1.
"""

import argparse
import gzip
import json
from pathlib import Path
import pickle
import zlib

import numpy as np
import pandas as pd

from extra_scripts.collect_local_path_reads import path_key
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from extra_scripts.run_paired_path_test import filtered_inputs
from tealeaf.sc.local_path_reads import path_read_opportunities, pooled_path_shares, simulate_null_class_counts
from tealeaf.sc.local_path_test import PRIMERS, fit_local_block

READ_LENGTH = 151
SCENARIOS = {"subject 20, row 20": (20., 20.), "subject 20, no row variation": (20., None), "subject 20, row 5": (20., 5.)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hypotheses", type=int, default=160)
    parser.add_argument("--draws", type=int, default=2)
    parser.add_argument("--no-anchors", action="store_true")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--collect", action="store_true")
    args = parser.parse_args()
    if args.collect:
        table = pd.concat([pd.read_csv(path, sep="\t") for path in sorted(args.output_dir.glob("shard_*/null.tsv"))], ignore_index=True)
        if len(table) != args.hypotheses * args.draws * len(SCENARIOS):
            raise ValueError(f"incomplete null family: {len(table)} trials")
        table["multi_path"] = table.n_paths > 2
        summary = table.groupby(["scenario", "multi_path"]).apply(lambda local: pd.Series({"trials": len(local), "converged": int(local.converged.sum()), **{f"reject_{level:g}": float((local.p_value < level).mean()) for level in (.05, .01, .001)}}), include_groups=False).reset_index()
        summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
        table.to_csv(args.output_dir / "trials.tsv.gz", sep="\t", index=False, na_rep="NA")
        print(summary.to_string(index=False), flush=True)
        return
    cached = pickle.load(args.candidate_cache.open("rb"))
    counts = pd.read_csv(args.local_reads / "counts.tsv.gz", sep="\t")
    lookup = {(row.path_key, row.subject, row.cell_type, row.primer, int(row.mask)): float(row.count) for row in counts.itertuples(index=False)}
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    observed_keys = set(counts.path_key)
    eligible = [candidate for candidate in cached["candidates"] if path_key(candidate[1], candidate[6]) in observed_keys]
    chosen = np.sort(np.random.default_rng(20261011).choice(len(eligible), args.hypotheses, replace=False))
    candidates = partition_candidates([eligible[index] for index in chosen], args.shard_count)[args.shard_index]
    metadata = filtered_inputs(args.data_cache, cached["settings"])[0]
    records = []
    for candidate in candidates:
        test_id, block_id, _, _, _, _, signatures, rows, _, levels = candidate
        key = path_key(block_id, signatures)
        local_metadata, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
        subjects = local_metadata.mouse.astype(str).to_numpy()
        opportunities = path_read_opportunities(blocks[key]["paths"], READ_LENGTH)
        pairs = sorted(set(zip(subjects, labels.tolist())))
        totals = {(index, primer): sum(lookup.get((key, subject, levels[label], primer, mask), 0.) for mask in opportunities) for index, (subject, label) in enumerate(pairs) for primer in PRIMERS}
        pooled = {mask: sum(lookup.get((key, subject, levels[label], primer, mask), 0.) for subject, label in pairs for primer in PRIMERS) for mask in opportunities}
        try:
            shares = pooled_path_shares(pooled, opportunities)
        except ValueError:
            shares = np.full(len(signatures), 1 / len(signatures))
        for scenario, (subject_concentration, row_concentration) in SCENARIOS.items():
            for draw in range(args.draws):
                record = {"test_id": test_id, "n_paths": len(signatures), "scenario": scenario, "draw": draw, "p_value": 1., "converged": False, "error": ""}
                rng = np.random.default_rng(zlib.crc32(f"{test_id}|{scenario}|{draw}".encode()))
                simulated = simulate_null_class_counts(totals, opportunities, shares, [subject for subject, _ in pairs], rng, subject_concentration=subject_concentration, row_concentration=row_concentration)
                simulated_lookup = {(key, pairs[row][0], levels[pairs[row][1]], primer, mask): float(value) for (row, primer, mask), value in simulated.items()}
                try:
                    fields, _ = fit_local_block(simulated_lookup, key, opportunities, len(signatures), subjects, labels, levels, anchors=not args.no_anchors, primer_offset=True)
                    record.update(p_value=fields["p_value"], converged=fields["converged"], statistic=fields["statistic"], n_converged_path_tests=fields["n_converged_path_tests"])
                except (ValueError, np.linalg.LinAlgError) as exception:
                    record["error"] = str(exception)
                records.append(record)
    out = args.output_dir / f"shard_{args.shard_index}"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out / "null.tsv", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
