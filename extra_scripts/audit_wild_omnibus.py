#!/usr/bin/env python3
"""Subject-level wild residual bootstrap of quantified omnibus compositions."""

import argparse
import json
from pathlib import Path
import pickle
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc import differential, ec_block_glmm
from tealeaf.sc.omnibus import cluster_max_f


def wild_tests(values, labels, subjects, rng, replicates=64):
    design, tested, _, groups = ec_block_glmm.blocked_multilevel_design(labels, subjects)
    null_design = np.delete(design, tested, axis=1)
    projected = null_design @ np.linalg.lstsq(null_design, values, rcond=None)[0]
    residuals = values - projected
    encoded = np.searchsorted(groups, subjects)
    observed = cluster_max_f(values, design, tested, subjects)
    nulls = []
    for replicate in range(replicates):
        signs = rng.choice([-1., 1.], len(groups))
        synthetic = projected + residuals * signs[encoded, None]
        nulls.append({"replicate": replicate, **cluster_max_f(synthetic, design, tested, subjects)})
    return observed, nulls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--test-concentration", type=float, default=32.)
    args = parser.parse_args()
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    outputs, nulls, failures, baselines = [], [], [], {}
    for candidate in partition_candidates(cached["candidates"], args.shard_count)[args.shard_index]:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "path_signatures": json.dumps(signatures), "n_paths": len(signatures), "converged": True}
        try:
            local, _, labels = local_test_design(metadata, rows, tested_levels, "cell_type")
            subjects = local.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local), 1)), subjects, drop_zero=False)
            key = (gene, tuple(rows), tuple(transcripts))
            if key not in baselines:
                baselines[key] = ec_block_glmm.pooled_isoform_weights(base)
            fitted = ec_block_glmm.blocked_multilevel_path_test(base, path_index, labels, subjects, baseline=baselines[key], path_pseudocount=args.test_concentration, path_pseudocount_scaling="total")
            observed, null = wild_tests(fitted["values"], fitted["observation_labels"], fitted["observation_subjects"], np.random.default_rng(zlib.crc32(test_id.encode()) + 712), 64)
            design, tested, levels, _ = ec_block_glmm.blocked_multilevel_design(fitted["observation_labels"], fitted["observation_subjects"])
            proportions = np.asarray([fit.path_proportions for fit in fitted["path_fits"]])
            coefficients = np.linalg.lstsq(design, proportions, rcond=None)[0][tested]
            details = {"n_subjects": fitted["n_subjects"], "strategy": f"CR2 maximum-coordinate wild ILR A{int(args.test_concentration)}", "levels": json.dumps([tested_levels[int(level)] for level in levels]), "adjusted_effects": json.dumps(np.vstack([np.zeros(len(signatures)), coefficients]).tolist())}
            outputs.append({**header, **details, **observed})
            nulls.extend({**header, "n_subjects": fitted["n_subjects"], "strategy": details["strategy"], **record} for record in null)
        except (ValueError, np.linalg.LinAlgError) as exception:
            failures.append({"test_id": test_id, "error": str(exception)})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output_dir / "observed.tsv", sep="\t", index=False)
    pd.DataFrame(nulls).to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")


if __name__ == "__main__":
    main()
