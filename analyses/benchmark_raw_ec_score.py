#!/usr/bin/env python3
"""Benchmark a plug-in efficient score screen on real splice blocks."""

from __future__ import annotations

import argparse
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd
import scipy.stats

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    reparameterize_fixed_effects,
    tensor_data,
    treatment_design,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm, ec_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-candidates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=73)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--mode-steps", type=int, default=12)
    return parser.parse_args()


def main():
    args = parse_args()
    with args.data_cache.open("rb") as handle:
        groups, counts, _, _, gene_ecs, designs = pickle.load(handle)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    candidates = cached["candidates"]
    rng = np.random.default_rng(args.seed)
    selected = rng.choice(
        len(candidates),
        size=min(args.max_candidates, len(candidates)),
        replace=False,
    )
    candidates = [candidates[index] for index in sorted(selected)]
    metadata = group_metadata(groups)
    if str(settings["test_effect"]).startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        levels = represented.index[
            represented >= int(settings["min_celltype_mice"])
        ]
        retained = metadata.cell_type.isin(levels).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    objective_cache = {}
    output = []
    for candidate in candidates:
        (
            test_id,
            _,
            _,
            gene,
            transcripts,
            path_index,
            _,
            rows,
            _,
            tested_levels,
        ) = candidate
        started = time.perf_counter()
        local_metadata, nuisance, labels = local_test_design(
            metadata, rows, tested_levels, settings["test_effect"]
        )
        local_counts = tuple(matrix[rows] for matrix in counts)
        base, _, _ = local_gene_data(
            local_counts,
            designs,
            transcripts,
            gene_ecs[gene],
            np.ones((len(local_metadata), 1)),
            local_metadata.mouse.to_numpy(),
            drop_zero=False,
        )
        tested = treatment_design(labels, len(tested_levels))
        null_tensor, alternative_tensor, degrees = (
            ec_block_glmm.block_fixed_effect_tensors(
                nuisance, tested, path_index
            )
        )
        null_data = tensor_data(base, null_tensor)
        null = ec_glmm.fit_laplace(
            null_data,
            max_iter=args.max_iter,
            mode_steps=args.mode_steps,
            mode_gradient="implicit",
            objective_cache=objective_cache,
            objective_cache_key=(
                "global_shape",
                ec_glmm.laplace_objective_shape_key(null_data),
            ),
        )
        alternative_data = tensor_data(base, alternative_tensor)
        initial = reparameterize_fixed_effects(
            null, null_tensor, alternative_tensor
        )
        evaluated = ec_glmm.fit_laplace(
            alternative_data,
            initial=initial,
            max_iter=0,
            mode_steps=args.mode_steps,
            mode_gradient="implicit",
            return_outer_hessian=True,
            objective_cache=objective_cache,
            objective_cache_key=(
                "global_shape",
                ec_glmm.laplace_objective_shape_key(alternative_data),
            ),
        )
        tested_indices = np.arange(
            null_tensor.shape[2], alternative_tensor.shape[2]
        )
        statistic = ec_glmm.efficient_score_statistic(
            evaluated, tested_indices
        )
        fixed_effect_count = int(evaluated["fixed_effect_count"])
        fixed_hessian = evaluated["outer_hessian"][
            :fixed_effect_count, :fixed_effect_count
        ]
        hessian_eigenvalues = np.linalg.eigvalsh(fixed_hessian)
        absolute_eigenvalues = np.abs(hessian_eigenvalues)
        positive_eigenvalues = absolute_eigenvalues[absolute_eigenvalues > 1e-12]
        hessian_condition = (
            np.inf
            if not len(positive_eigenvalues)
            else absolute_eigenvalues.max() / positive_eigenvalues.min()
        )
        output.append(
            {
                "test_id": test_id,
                "degrees": degrees,
                "score_statistic": statistic,
                "score_p_value": scipy.stats.chi2.sf(statistic, degrees),
                "fixed_hessian_minimum_eigenvalue": hessian_eigenvalues.min(),
                "fixed_hessian_maximum_eigenvalue": hessian_eigenvalues.max(),
                "fixed_hessian_condition": hessian_condition,
                "null_converged": null["converged"],
                "null_objective_evaluations": null["objective_evaluations"],
                "elapsed_seconds": time.perf_counter() - started,
            }
        )
        while len(objective_cache) > 16:
            objective_cache.pop(next(iter(objective_cache)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
