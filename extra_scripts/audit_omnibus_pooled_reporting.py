#!/usr/bin/env python3
"""Pool full-data omnibus reporting without selecting contrasts using LR signs."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc import differential, ec_block_glmm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--profile-mass", action="store_true")
    args = parser.parse_args()
    if args.profile_mass:
        differential.fit_path_perturbation = differential.fit_profiled_path_perturbation
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    output, baselines = [], {}
    for candidate in partition_candidates(cached["candidates"], args.shard_count)[args.shard_index]:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, levels = candidate
        local, _, labels = local_test_design(metadata, rows, levels, "cell_type")
        base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local), 1)), local.mouse.astype(str).to_numpy(), drop_zero=False)
        key = (gene, tuple(rows), tuple(transcripts))
        if key not in baselines:
            baselines[key] = ec_block_glmm.pooled_isoform_weights(base)
        for strategy, balanced in (("pooled local", False), ("primer-balanced pooled local", True)):
            proportions, statuses = [], []
            try:
                for index in range(len(levels)):
                    pooled = [np.asarray(values[labels == index], dtype=float).sum(axis=0) for values in base.counts]
                    positive = [value.sum() for value in pooled if value.sum() > 0]
                    if not positive:
                        raise ValueError("empty reporting level")
                    if balanced:
                        target = np.mean(positive)
                        pooled = [value * target / value.sum() if value.sum() > 0 else value for value in pooled]
                    fit = differential.fit_path_perturbation(pooled, base.compatibility, baselines[key], path_index, path_pseudocount=1., path_pseudocount_scaling="total")
                    proportions.append(fit.path_proportions.tolist())
                    statuses.append(fit.converged)
                output.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "strategy": strategy, "levels": json.dumps(levels), "path_signatures": json.dumps(signatures), "adjusted_effects": json.dumps(proportions), "converged": all(statuses)})
            except (ValueError, np.linalg.LinAlgError):
                output.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "strategy": strategy, "converged": False})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
