#!/usr/bin/env python3
"""Count-level null calibration of the binary EC random-subject grid test.

Random two-path full-data hypotheses (no significance or LR selection) keep
their observed depths, primer totals and compatibility; counts are redrawn
with no true cell-type path effect. Failed fits stay in the denominator as p=1.
"""

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_marginal_grid import binary_grid_test, prepare_grid_likelihood
from tealeaf.sc.path_simulation import simulate_counts

SCENARIOS = {"common composition": {}, "biological precision 20": {"residual_concentration": 20.}, "event mass shift": {"residual_concentration": 20., "event_mass_type_scale": .5}, "within-path tilt": {"residual_concentration": 20., "within_path_type_scale": 1.}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hypotheses", type=int, default=128)
    parser.add_argument("--draws", type=int, default=2)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--collect", action="store_true")
    args = parser.parse_args()
    if args.collect:
        table = pd.concat([pd.read_csv(path, sep="\t") for path in sorted(args.output_dir.glob("shard_*/null.tsv"))], ignore_index=True)
        expected = args.hypotheses * args.draws * len(SCENARIOS)
        if len(table) != expected:
            raise ValueError(f"expected {expected} trials, found {len(table)}")
        rows = []
        for scenario, local in table.groupby("scenario"):
            rows.append({"scenario": scenario, "trials": len(local), "converged": int(local.converged.sum()), **{f"reject_{level:g}": float((local.p_value < level).mean()) for level in (.05, .01, .001)}, **{f"reject_{level:g}_converged": float((local.loc[local.converged].p_value < level).mean()) for level in (.05, .01)}})
        summary = pd.DataFrame(rows)
        summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
        table.to_csv(args.output_dir / "trials.tsv.gz", sep="\t", index=False, na_rep="NA")
        print(summary.to_string(index=False), flush=True)
        return
    cached = pickle.load(args.candidate_cache.open("rb"))
    binary = [candidate for candidate in cached["candidates"] if len(candidate[6]) == 2]
    chosen = np.sort(np.random.default_rng(20261010).choice(len(binary), args.hypotheses, replace=False))
    candidates = partition_candidates([binary[index] for index in chosen], args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    records = []
    for candidate in candidates:
        test_id, _, _, gene, transcripts, path_index, _, rows, _, levels = candidate
        local_metadata, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
        subjects = local_metadata.mouse.astype(str).to_numpy()
        base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
        baseline = pooled_isoform_weights(base, max_iter=250)
        for scenario, options in SCENARIOS.items():
            for draw in range(args.draws):
                start = time.monotonic()
                record = {"test_id": test_id, "scenario": scenario, "draw": draw, "p_value": 1., "converged": False, "error": ""}
                try:
                    rng = np.random.default_rng(zlib.crc32(f"{test_id}|{scenario}|{draw}".encode()))
                    generated = simulate_counts(base, baseline, subjects, rng, .5, labels=labels, path_index=path_index, **options)
                    nuisance = pooled_isoform_weights(generated, max_iter=250)
                    result = binary_grid_test(prepare_grid_likelihood(generated, path_index, labels, subjects, nuisance))
                    record.update(p_value=result["p_value"], converged=result["converged"], statistic=result["statistic"], kappa=result["alternative_concentration"], sigma=result["alternative_subject_sd"], quadrature_error=result["quadrature_error"])
                except (ValueError, np.linalg.LinAlgError) as exception:
                    record["error"] = str(exception)
                record["runtime_seconds"] = time.monotonic() - start
                records.append(record)
    out = args.output_dir / f"shard_{args.shard_index}"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out / "null.tsv", sep="\t", index=False, na_rep="NA")
    (out / "settings.json").write_text(json.dumps({"hypotheses": args.hypotheses, "draws": args.draws, "scenarios": SCENARIOS, "subject_scale": .5, "selection": "uniform random two-path full-data candidates, seed 20261010"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
