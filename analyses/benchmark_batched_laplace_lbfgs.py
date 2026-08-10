#!/usr/bin/env python3
"""Benchmark independent batched L-BFGS on padded Laplace EC-GLMM fits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    tensor_data,
    treatment_design,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm, ec_glmm
from tealeaf.sc.batched_lbfgs import minimize_batched_lbfgs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--batchability", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument(
        "--selection", choices=("largest", "high_cost"), default="largest"
    )
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--mode-steps", type=int, default=12)
    return parser.parse_args()


def filter_metadata_and_counts(groups, counts, settings):
    metadata = group_metadata(groups)
    if str(settings.get("test_effect", "cell_type")).startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        levels = represented.index[
            represented >= int(settings.get("min_celltype_mice", 3))
        ]
        retained = metadata.cell_type.isin(levels).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    return metadata, counts


def prepare_batch(args):
    with args.data_cache.open("rb") as handle:
        groups, counts, _, _, gene_ecs, designs = pickle.load(handle)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached.get("settings", {})
    metadata, counts = filter_metadata_and_counts(groups, counts, settings)
    table = pd.read_csv(args.batchability, sep="\t")
    available = {candidate[0]: candidate for candidate in cached["candidates"]}
    table = table.loc[table.test_id.isin(available)].copy()
    table["fit_pair_key"] = (
        table.isoforms.astype(str)
        + ","
        + table.null_coefficients.astype(str)
        + ","
        + table.alternative_coefficients.astype(str)
    )
    grouped = table.groupby("fit_pair_key")
    if args.selection == "largest":
        key = grouped.size().idxmax()
    else:
        eligible_keys = grouped.size()
        eligible_keys = eligible_keys.index[eligible_keys >= int(args.batch_size)]
        key = (
            table.loc[table.fit_pair_key.isin(eligible_keys)]
            .groupby("fit_pair_key")
            .cost_proxy.median()
            .idxmax()
        )
    eligible = table.loc[table.fit_pair_key.eq(key)].sort_values("cost_proxy")
    if len(eligible) > int(args.batch_size):
        middle = len(eligible) // 2
        half = int(args.batch_size) // 2
        eligible = eligible.iloc[middle - half : middle - half + int(args.batch_size)]
    prepared = []
    for test_id in eligible.test_id:
        candidate = available[test_id]
        _, _, _, gene, transcripts, path_index, _, rows, _, tested_levels = candidate
        local_metadata, nuisance, labels = local_test_design(
            metadata,
            rows,
            tested_levels,
            settings.get("test_effect", "cell_type"),
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
        null_tensor, _, _ = ec_block_glmm.block_fixed_effect_tensors(
            nuisance, tested, path_index
        )
        alternative_tensor = ec_block_glmm.full_fixed_effect_tensor(
            nuisance, tested, len(transcripts) - 1
        )
        prepared.append(
            (
                test_id,
                (gene, tuple(np.asarray(rows, dtype=int))),
                tensor_data(base, null_tensor),
                tensor_data(base, alternative_tensor),
            )
        )
    return key, prepared


def pad_batch(data, initial=None):
    fit_count = len(data)
    isoforms = data[0].n_isoforms
    dimension = isoforms - 1
    coefficient_count = data[0].fixed_effect_tensor.shape[2]
    maximum_samples = max(len(item.clusters) for item in data)
    maximum_clusters = max(len(np.unique(item.clusters)) for item in data)
    maximum_ecs = [
        max(item.counts[primer].shape[1] for item in data)
        for primer in range(len(data[0].counts))
    ]
    tensors = np.zeros((fit_count, maximum_samples, dimension, coefficient_count))
    memberships = np.zeros((fit_count, maximum_samples, maximum_clusters))
    counts = [
        np.zeros((fit_count, maximum_samples, ec_count)) for ec_count in maximum_ecs
    ]
    mappings = [np.zeros((fit_count, ec_count, isoforms)) for ec_count in maximum_ecs]
    for fit, item in enumerate(data):
        samples = len(item.clusters)
        clusters = pd.factorize(item.clusters, sort=True)[0]
        tensors[fit, :samples] = item.fixed_effect_tensor
        memberships[fit, np.arange(samples), clusters] = 1.0
        for primer in range(len(counts)):
            ecs = item.counts[primer].shape[1]
            counts[primer][fit, :samples, :ecs] = item.counts[primer]
            mappings[primer][fit, :ecs] = item.compatibility[primer]
    if initial is None:
        initial = np.stack(
            [ec_glmm._initial_outer(item, "multinomial") for item in data]
        )
    else:
        initial = np.asarray(initial, dtype=float)
        expected = (fit_count, coefficient_count + 1)
        if initial.shape != expected:
            raise ValueError(f"initial has shape {initial.shape}, expected {expected}")
    return tensors, tuple(counts), tuple(mappings), memberships, initial


def batched_laplace_functions(data, mode_steps, initial=None):
    import jax
    import jax.numpy as jnp
    import jax.scipy as jsp

    jax.config.update("jax_enable_x64", True)
    tensors, counts, mappings, memberships, initial = pad_batch(data, initial)
    tensors = jnp.asarray(tensors)
    counts = tuple(jnp.asarray(value) for value in counts)
    mappings = tuple(jnp.asarray(value) for value in mappings)
    memberships = jnp.asarray(memberships)
    dimension = data[0].n_isoforms - 1
    coefficient_count = data[0].fixed_effect_tensor.shape[2]
    cluster_count = memberships.shape[2]

    def observation_log_likelihood(free_logits, observed_rows, local_mappings):
        logits = jnp.concatenate((free_logits, jnp.zeros(1, dtype=free_logits.dtype)))
        abundance = jnp.exp(logits - jnp.max(logits))
        value = jnp.asarray(0.0, dtype=logits.dtype)

        def primer_log_likelihood(arguments):
            observed, mapping = arguments
            mass = mapping @ abundance
            active_ec = jnp.sum(mapping, axis=1) > 0.0
            denominator = jnp.sum(jnp.where(active_ec, mass, 0.0))
            probability = jnp.where(
                active_ec,
                jnp.maximum(mass / denominator, 1e-300),
                1.0,
            )
            total = jnp.sum(observed)
            constant = jsp.special.gammaln(total + 1.0) - jnp.sum(
                jsp.special.gammaln(observed + 1.0)
            )
            return constant + jnp.sum(observed * jnp.log(probability))

        for observed, mapping in zip(observed_rows, local_mappings):
            value += jax.lax.cond(
                jnp.any(mapping > 0.0) & (jnp.sum(observed) > 0.0),
                primer_log_likelihood,
                lambda arguments: jnp.asarray(0.0, dtype=logits.dtype),
                (observed, mapping),
            )
        return value

    row_gradient = jax.vmap(
        jax.grad(observation_log_likelihood, argnums=0),
        in_axes=(0, (0, 0), (None, None)),
    )
    row_hessian = jax.vmap(
        jax.hessian(observation_log_likelihood, argnums=0),
        in_axes=(0, (0, 0), (None, None)),
    )

    def mode_derivatives(
        modes,
        outer,
        tensor,
        observed_1,
        observed_2,
        mapping_1,
        mapping_2,
        membership,
    ):
        coefficients = outer[:coefficient_count]
        variance = jnp.exp(2.0 * outer[coefficient_count])
        fixed = jnp.einsum("ndp,p->nd", tensor, coefficients)
        free_logits = fixed + membership @ modes
        score_rows = row_gradient(
            free_logits,
            (observed_1, observed_2),
            (mapping_1, mapping_2),
        )
        score = membership.T @ score_rows - modes / variance
        curvature_rows = -row_hessian(
            free_logits,
            (observed_1, observed_2),
            (mapping_1, mapping_2),
        )
        hessian = (
            jnp.einsum("nm,nij->mij", membership, curvature_rows)
            + jnp.eye(dimension)[None, :, :] / variance
        )
        return free_logits, score, hessian, variance

    def solve_mode(
        outer,
        tensor,
        observed_1,
        observed_2,
        mapping_1,
        mapping_2,
        membership,
    ):
        current = jnp.zeros((cluster_count, dimension), dtype=jnp.float64)

        def update(_, modes):
            _, score, hessian, _ = mode_derivatives(
                modes,
                outer,
                tensor,
                observed_1,
                observed_2,
                mapping_1,
                mapping_2,
                membership,
            )
            eigen_floor = jnp.maximum(1e-6 - jnp.linalg.eigvalsh(hessian)[:, 0], 0.0)
            stabilized = hessian + eigen_floor[:, None, None] * jnp.eye(dimension)
            step = jnp.linalg.solve(stabilized, score[..., None])[..., 0]
            maximum = jnp.max(jnp.abs(step), axis=1, keepdims=True)
            return modes + jnp.minimum(1.0, 2.0 / (maximum + 1e-12)) * step

        return jax.lax.fori_loop(0, int(mode_steps), update, current)

    @jax.custom_vjp
    def posterior_mode(
        outer,
        tensor,
        observed_1,
        observed_2,
        mapping_1,
        mapping_2,
        membership,
    ):
        return solve_mode(
            outer,
            tensor,
            observed_1,
            observed_2,
            mapping_1,
            mapping_2,
            membership,
        )

    def mode_forward(*arguments):
        modes = solve_mode(*arguments)
        return modes, (modes, *arguments)

    def mode_backward(residuals, cotangent):
        (
            modes,
            outer,
            tensor,
            observed_1,
            observed_2,
            mapping_1,
            mapping_2,
            membership,
        ) = residuals
        _, _, hessian, _ = mode_derivatives(
            modes,
            outer,
            tensor,
            observed_1,
            observed_2,
            mapping_1,
            mapping_2,
            membership,
        )
        solved = jnp.linalg.solve(hessian, cotangent[..., None])[..., 0]
        _, pullback = jax.vjp(
            lambda value: mode_derivatives(
                modes,
                value,
                tensor,
                observed_1,
                observed_2,
                mapping_1,
                mapping_2,
                membership,
            )[1],
            outer,
        )
        return (
            pullback(solved)[0],
            jnp.zeros_like(tensor),
            jnp.zeros_like(observed_1),
            jnp.zeros_like(observed_2),
            jnp.zeros_like(mapping_1),
            jnp.zeros_like(mapping_2),
            jnp.zeros_like(membership),
        )

    posterior_mode.defvjp(mode_forward, mode_backward)

    def objective(
        outer,
        tensor,
        observed_1,
        observed_2,
        mapping_1,
        mapping_2,
        membership,
    ):
        modes = posterior_mode(
            outer,
            tensor,
            observed_1,
            observed_2,
            mapping_1,
            mapping_2,
            membership,
        )
        free_logits, _, hessian, variance = mode_derivatives(
            modes,
            outer,
            tensor,
            observed_1,
            observed_2,
            mapping_1,
            mapping_2,
            membership,
        )
        likelihood = jnp.sum(
            jax.vmap(observation_log_likelihood, in_axes=(0, (0, 0), (None, None)))(
                free_logits,
                (observed_1, observed_2),
                (mapping_1, mapping_2),
            )
        )
        log_sd = outer[coefficient_count]
        negative_joint = -likelihood + 0.5 * jnp.sum(modes**2 / variance)
        negative_joint += (
            cluster_count * dimension * (log_sd + 0.5 * jnp.log(2.0 * jnp.pi))
        )
        sign, log_determinant = jnp.linalg.slogdet(hessian)
        correction = 0.5 * jnp.sum(log_determinant)
        correction -= 0.5 * modes.size * jnp.log(2.0 * jnp.pi)
        return jnp.where(jnp.all(sign > 0), negative_joint + correction, jnp.inf)

    arguments = (tensors, *counts, *mappings, memberships)
    value_and_gradient = jax.jit(jax.vmap(jax.value_and_grad(objective, argnums=0)))
    score_norm = jax.jit(
        jax.vmap(
            lambda outer, tensor, observed_1, observed_2, mapping_1, mapping_2, membership: jnp.max(
                jnp.abs(
                    mode_derivatives(
                        solve_mode(
                            outer,
                            tensor,
                            observed_1,
                            observed_2,
                            mapping_1,
                            mapping_2,
                            membership,
                        ),
                        outer,
                        tensor,
                        observed_1,
                        observed_2,
                        mapping_1,
                        mapping_2,
                        membership,
                    )[1]
                )
            )
        )
    )

    def evaluate(parameters):
        values, gradients = value_and_gradient(parameters, *arguments)
        jax.block_until_ready((values, gradients))
        return np.asarray(values), np.asarray(gradients)

    def evaluate_score(parameters):
        values = score_norm(parameters, *arguments)
        jax.block_until_ready(values)
        return np.asarray(values)

    return initial, evaluate, evaluate_score


def fit_batched(data, max_iter, mode_steps, initial=None):
    initial, evaluate, evaluate_score = batched_laplace_functions(
        data, mode_steps, initial
    )
    coefficient_count = data[0].fixed_effect_tensor.shape[2]
    bounds = [(-20.0, 20.0)] * coefficient_count + [(-8.0, 3.0)]
    started = time.perf_counter()
    result = minimize_batched_lbfgs(
        evaluate,
        initial,
        bounds,
        max_iter=max_iter,
        gradient_tolerance=1e-5,
    )
    score = evaluate_score(result.parameters)
    elapsed = time.perf_counter() - started
    gradient_norm = np.max(np.abs(result.gradients), axis=1) / result.scales
    converged = np.isfinite(result.values) & (gradient_norm <= 1e-4) & (score <= 1e-4)
    return result, elapsed, gradient_norm, score, converged


def fit_scalar(data, max_iter, mode_steps, initial=None):
    fits = []
    started = time.perf_counter()
    if initial is None:
        initial = [None] * len(data)
    for item, starting in zip(data, initial):
        fits.append(
            ec_glmm.fit_laplace(
                item,
                initial=starting,
                max_iter=max_iter,
                mode_steps=mode_steps,
                mode_gradient="implicit",
            )
        )
    return fits, time.perf_counter() - started


def reparameterize_parameters(parameters, source_tensor, target_tensor):
    """Project fitted logits into a target fixed-effect basis."""
    source_count = source_tensor.shape[2]
    fitted_logits = source_tensor.reshape(-1, source_count) @ parameters[:source_count]
    target_coefficients = np.linalg.lstsq(
        target_tensor.reshape(-1, target_tensor.shape[2]),
        fitted_logits,
        rcond=None,
    )[0]
    return np.r_[target_coefficients, parameters[source_count:]]


def main():
    args = parse_args()
    import jax

    key, prepared = prepare_batch(args)
    test_ids = [item[0] for item in prepared]
    null_data = [item[2] for item in prepared]
    alternative_data = []
    alternative_indices = []
    alternative_index_by_key = {}
    for _, alternative_key, _, alternative in prepared:
        if alternative_key not in alternative_index_by_key:
            alternative_index_by_key[alternative_key] = len(alternative_data)
            alternative_data.append(alternative)
        alternative_indices.append(alternative_index_by_key[alternative_key])
    batched_null, batched_null_seconds, null_gradient, null_score, null_converged = fit_batched(
        null_data, args.max_iter, args.mode_steps
    )
    scalar_null, scalar_null_seconds = fit_scalar(
        null_data, args.max_iter, args.mode_steps
    )
    first_test_by_alternative = {}
    for test_index, alternative_index in enumerate(alternative_indices):
        first_test_by_alternative.setdefault(alternative_index, test_index)
    batched_alternative_initial = []
    scalar_alternative_initial = []
    for alternative_index, alternative in enumerate(alternative_data):
        test_index = first_test_by_alternative[alternative_index]
        source = null_data[test_index].fixed_effect_tensor
        batched_alternative_initial.append(
            reparameterize_parameters(
                batched_null.parameters[test_index],
                source,
                alternative.fixed_effect_tensor,
            )
        )
        scalar_alternative_initial.append(
            reparameterize_parameters(
                scalar_null[test_index]["parameters"],
                source,
                alternative.fixed_effect_tensor,
            )
        )
    (
        batched_alternative,
        batched_alternative_seconds,
        alternative_gradient,
        alternative_score,
        alternative_converged,
    ) = fit_batched(
        alternative_data,
        args.max_iter,
        args.mode_steps,
        batched_alternative_initial,
    )
    scalar_alternative, scalar_alternative_seconds = fit_scalar(
        alternative_data,
        args.max_iter,
        args.mode_steps,
        scalar_alternative_initial,
    )
    rows = []
    for index, test_id in enumerate(test_ids):
        alternative_index = alternative_indices[index]
        batched_statistic = max(
            0.0,
            2.0
            * (
                batched_null.values[index]
                - batched_alternative.values[alternative_index]
            ),
        )
        scalar_statistic = max(
            0.0,
            2.0
            * (
                scalar_null[index]["objective"]
                - scalar_alternative[alternative_index]["objective"]
            ),
        )
        rows.append(
            {
                "test_id": test_id,
                "batched_null_converged": null_converged[index],
                "batched_alternative_converged": alternative_converged[
                    alternative_index
                ],
                "scalar_null_converged": scalar_null[index]["converged"],
                "scalar_alternative_converged": scalar_alternative[alternative_index][
                    "converged"
                ],
                "batched_null_objective": batched_null.values[index],
                "scalar_null_objective": scalar_null[index]["objective"],
                "batched_alternative_objective": batched_alternative.values[
                    alternative_index
                ],
                "scalar_alternative_objective": scalar_alternative[alternative_index][
                    "objective"
                ],
                "batched_statistic": batched_statistic,
                "scalar_statistic": scalar_statistic,
                "absolute_statistic_difference": abs(
                    batched_statistic - scalar_statistic
                ),
                "batched_null_gradient_norm": null_gradient[index],
                "batched_alternative_gradient_norm": alternative_gradient[
                    alternative_index
                ],
                "batched_null_mode_score_norm": null_score[index],
                "batched_alternative_mode_score_norm": alternative_score[
                    alternative_index
                ],
                "batched_null_iterations": batched_null.iterations[index],
                "batched_alternative_iterations": batched_alternative.iterations[
                    alternative_index
                ],
                "batched_null_line_search_evaluations": batched_null.line_search_evaluations[
                    index
                ],
                "batched_alternative_line_search_evaluations": batched_alternative.line_search_evaluations[
                    alternative_index
                ],
                "scalar_null_iterations": scalar_null[index]["iterations"],
                "scalar_alternative_iterations": scalar_alternative[alternative_index][
                    "iterations"
                ],
            }
        )
    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    summary = {
        "jax_backend": jax.default_backend(),
        "jax_devices": [str(device) for device in jax.devices()],
        "fit_pair_key": repr(key),
        "tests": len(result),
        "unique_alternatives": len(alternative_data),
        "alternative_initialization": "fitted_null_projection",
        "batched_seconds": batched_null_seconds + batched_alternative_seconds,
        "scalar_seconds": scalar_null_seconds + scalar_alternative_seconds,
        "speedup": (scalar_null_seconds + scalar_alternative_seconds)
        / (batched_null_seconds + batched_alternative_seconds),
        "batched_joint_convergence": int(
            sum(
                null_converged[index]
                and alternative_converged[alternative_indices[index]]
                for index in range(len(test_ids))
            )
        ),
        "scalar_joint_convergence": int(
            sum(
                scalar_null[index]["converged"]
                and scalar_alternative[alternative_indices[index]]["converged"]
                for index in range(len(test_ids))
            )
        ),
        "median_absolute_statistic_difference": float(
            result.absolute_statistic_difference.median()
        ),
        "maximum_absolute_statistic_difference": float(
            result.absolute_statistic_difference.max()
        ),
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
