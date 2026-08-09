#!/usr/bin/env python3
"""Run estimate-once local-path inference on an EC-block manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    partition_candidates,
    treatment_design,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import differential, ec_block_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--covariance",
        choices=("conditional", "profile"),
        default="conditional",
    )
    parser.add_argument("--fit-max-iter", type=int, default=100)
    parser.add_argument("--label-permutations", type=int, default=0)
    parser.add_argument("--wild-null-replicates", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--inference",
        choices=("wald", "gaussian_lrt_bf"),
        default="wald",
    )
    parser.add_argument(
        "--effect-prior-scales",
        default="0.1,0.25,0.5,1.0",
        help="Comma- or semicolon-separated path-ILR prior standard deviations",
    )
    parser.add_argument(
        "--gaussian-solver", choices=("dense", "woodbury"), default="woodbury"
    )
    parser.add_argument(
        "--cluster-adjustment", choices=("cr1", "cr2", "cr3"), default="cr1"
    )
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    return parser.parse_args()


def candidate_rng(seed, test_id, replicate, null_type):
    """Construct a reproducible generator for one test and null draw."""
    test_hash = zlib.crc32(str(test_id).encode("utf-8"))
    type_hash = zlib.crc32(str(null_type).encode("utf-8"))
    return np.random.default_rng(
        np.random.SeedSequence((int(seed), test_hash, int(replicate), type_hash))
    )


def parse_prior_scales(value):
    """Parse positive effect-size prior scales from CLI text."""
    scales = np.asarray([
        float(item) for item in str(value).replace(";", ",").split(",")
        if item.strip()
    ])
    if not len(scales) or not np.isfinite(scales).all() or np.any(scales <= 0):
        raise ValueError("effect prior scales must be positive and finite")
    return tuple(scales)


def reference_contrast_prior_covariance(tested_columns, dimension):
    """Return the invariant covariance induced by iid latent level effects."""
    contrast_covariance = (
        np.eye(int(tested_columns))
        + np.ones((int(tested_columns), int(tested_columns)))
    )
    return np.kron(contrast_covariance, np.eye(int(dimension)))


def permute_labels_within_clusters(labels, clusters, rng):
    """Shuffle observed multi-level labels independently within each cluster."""
    labels = np.asarray(labels, dtype=int)
    clusters = np.asarray(clusters)
    if len(labels) != len(clusters):
        raise ValueError("labels and clusters must align")
    permuted = labels.copy()
    for cluster in np.unique(clusters):
        positions = np.flatnonzero(clusters == cluster)
        permuted[positions] = rng.permutation(permuted[positions])
    return permuted


def effect_summary(result, nuisance_columns, n_paths):
    """Return tested ILR and centered-log path effects with marginal SEs."""
    coefficients = np.asarray(result["coefficients"], dtype=float)
    tested = coefficients[nuisance_columns:]
    dimension = n_paths - 1
    tested_indices = np.concatenate([
        np.arange(column * dimension, (column + 1) * dimension)
        for column in range(nuisance_columns, len(coefficients))
    ])
    covariance = np.asarray(result["coefficient_covariance"], dtype=float)
    tested_covariance = covariance[np.ix_(tested_indices, tested_indices)]
    basis = differential.helmert_basis(n_paths)
    transform = np.kron(np.eye(len(tested)), basis)
    path_effects = tested @ basis.T
    path_variances = np.diag(transform @ tested_covariance @ transform.T)
    path_standard_errors = np.sqrt(np.maximum(path_variances, 0.0)).reshape(
        path_effects.shape
    )
    return tested, path_effects, path_standard_errors


def filtered_inputs(data_cache, settings):
    with data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(
            handle
        )
    metadata = group_metadata(groups)
    subject_folds = settings.get("subject_folds")
    if subject_folds is not None:
        folds = pd.read_csv(subject_folds, sep="\t", dtype={"subject": str})
        selected = set(
            folds.loc[
                folds.fold.eq(settings["subject_fold"]), "subject"
            ].astype(str)
        )
        retained = metadata.mouse.astype(str).isin(selected).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    if str(settings["test_effect"]).startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        levels = represented.index[
            represented >= int(settings["min_celltype_mice"])
        ]
        retained = metadata.cell_type.isin(levels).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    return metadata, counts, genes, gene_transcripts, gene_ecs, designs


def main():
    args = parse_args()
    effect_prior_scales = parse_prior_scales(args.effect_prior_scales)
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    if settings.get("joint_gene_test", False):
        raise ValueError("joint-gene manifests are not supported")
    candidates = cached["candidates"]
    if args.max_candidates is not None:
        candidates = candidates[: args.max_candidates]
    candidates = partition_candidates(candidates, args.shard_count)[
        args.shard_index
    ]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(
        args.data_cache, settings
    )
    test_effect = settings["test_effect"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    observed_rows = []
    null_rows = []
    failures = []
    pooled_cache = {}
    started = time.perf_counter()
    for position, candidate in enumerate(candidates):
        (
            test_id,
            block_id,
            gene_id,
            gene,
            transcripts,
            path_index,
            signatures,
            rows,
            tested_cell_type,
            tested_levels,
        ) = candidate
        try:
            local_metadata, nuisance, labels = local_test_design(
                metadata, rows, tested_levels, test_effect
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata.mouse.to_numpy()
            base, _, totals = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(local_metadata), 1)),
                clusters,
                drop_zero=False,
            )
            cache_key = (gene, tuple(rows), tuple(transcripts))
            if cache_key not in pooled_cache:
                pooled_cache[cache_key] = ec_block_glmm.pooled_isoform_weights(
                    base
                )
            tested = treatment_design(labels, len(tested_levels))
            effect_prior_covariance = reference_contrast_prior_covariance(
                tested.shape[1], len(signatures) - 1
            )
            if args.inference == "wald":
                result = ec_block_glmm.path_wald_test(
                    base,
                    path_index,
                    nuisance,
                    tested,
                    baseline=pooled_cache[cache_key],
                    covariance=args.covariance,
                    max_iter=args.fit_max_iter,
                    cluster_adjustment=args.cluster_adjustment,
                )
            else:
                result = ec_block_glmm.path_gaussian_lrt_test(
                    base,
                    path_index,
                    nuisance,
                    tested,
                    baseline=pooled_cache[cache_key],
                    covariance=args.covariance,
                    max_iter=args.fit_max_iter,
                    effect_prior_scales=effect_prior_scales,
                    effect_prior_covariance=effect_prior_covariance,
                    random_intercept_solver=args.gaussian_solver,
                )
            estimates = result.pop("path_estimates")
            usable_proportions = estimates["proportions"][estimates["usable"]]
            minimum_sample_proportions = usable_proportions.min(axis=1)
            tested_coefficients, path_effects, path_standard_errors = (
                effect_summary(result, nuisance.shape[1], len(signatures))
            )
            observed_rows.append({
                "test_id": test_id,
                "block_id": block_id,
                "gene_id": gene_id,
                "contrast": test_effect,
                "cell_type": tested_cell_type,
                "method": (
                    "path_wald" if args.inference == "wald"
                    else "path_gaussian_lrt_bf"
                ),
                "latent_space": "estimated_path",
                "covariance_method": args.covariance,
                "cluster_adjustment": (
                    args.cluster_adjustment if args.inference == "wald" else ""
                ),
                "n_paths": len(signatures),
                "n_isoforms": base.n_isoforms,
                "n_original_isoforms": base.n_isoforms,
                "n_ecs": len(gene_ecs[gene]),
                "n_samples": len(local_metadata),
                "n_samples_used": int(estimates["usable"].sum()),
                "sample_use_rate": float(estimates["usable"].mean()),
                "n_test_levels": len(tested_levels),
                "degrees_of_freedom": result["degrees_of_freedom"],
                "denominator_degrees_of_freedom": result.get(
                    "denominator_degrees_of_freedom", np.nan
                ),
                "median_gene_umis": float(np.median(totals)),
                "statistic": result["statistic"],
                "lrt_p_value": result["p_value"],
                "p_value": result["p_value"],
                "bic_log_bayes_factor": result.get(
                    "bic_log_bayes_factor", np.nan
                ),
                "mixture_log_bayes_factor": result.get(
                    "mixture_log_bayes_factor", np.nan
                ),
                "effect_prior_scales_json": json.dumps(
                    list(map(float, effect_prior_scales))
                ),
                "log_bayes_factors_json": json.dumps(
                    np.asarray(result.get("log_bayes_factors", [])).tolist()
                ),
                "null_cluster_variance": result.get(
                    "null_cluster_variance", np.nan
                ),
                "null_residual_variance": result.get(
                    "null_residual_variance", np.nan
                ),
                "alternative_cluster_variance": result.get(
                    "alternative_cluster_variance", np.nan
                ),
                "alternative_residual_variance": result.get(
                    "alternative_residual_variance", np.nan
                ),
                "null_converged": True,
                "alternative_converged": True,
                "n_clusters": result["clusters"],
                "path_fit_iterations": int(estimates["iterations"].sum()),
                "max_path_fit_iterations": int(estimates["iterations"].max()),
                "min_path_proportion": float(estimates["proportions"].min()),
                "path_proportion_q01": float(
                    np.quantile(usable_proportions, 0.01)
                ),
                "path_proportion_q05": float(
                    np.quantile(usable_proportions, 0.05)
                ),
                "sample_rate_min_path_below_0.001": float(
                    np.mean(minimum_sample_proportions < 0.001)
                ),
                "sample_rate_min_path_below_0.01": float(
                    np.mean(minimum_sample_proportions < 0.01)
                ),
                "sample_rate_min_path_below_0.05": float(
                    np.mean(minimum_sample_proportions < 0.05)
                ),
                "tested_levels_json": json.dumps(list(map(str, tested_levels))),
                "tested_coefficients_json": json.dumps(
                    tested_coefficients.tolist()
                ),
                "path_effects_json": json.dumps(path_effects.tolist()),
                "path_standard_errors_json": json.dumps(
                    path_standard_errors.tolist()
                ),
            })
            for replicate in range(args.label_permutations):
                try:
                    rng = candidate_rng(
                        args.seed, test_id, replicate, "label_permutation"
                    )
                    permuted_labels = permute_labels_within_clusters(
                        labels, clusters, rng
                    )
                    permuted_tested = treatment_design(
                        permuted_labels, len(tested_levels)
                    )
                    if args.inference == "wald":
                        permuted = ec_block_glmm.path_wald_from_estimates(
                            estimates,
                            nuisance,
                            permuted_tested,
                            clusters,
                            cluster_adjustment=args.cluster_adjustment,
                        )
                    else:
                        permuted = (
                            ec_block_glmm.path_gaussian_lrt_from_estimates(
                                estimates,
                                nuisance,
                                permuted_tested,
                                clusters,
                                effect_prior_scales=effect_prior_scales,
                                effect_prior_covariance=effect_prior_covariance,
                                null_variance_components=(
                                    result["null_cluster_variance"],
                                    result["null_residual_variance"],
                                ),
                                random_intercept_solver=args.gaussian_solver,
                            )
                        )
                    _, permuted_path_effects, _ = effect_summary(
                        permuted, nuisance.shape[1], len(signatures)
                    )
                    null_rows.append({
                        "test_id": test_id,
                        "null_type": "label_permutation",
                        "cluster_adjustment": (
                            args.cluster_adjustment
                            if args.inference == "wald" else ""
                        ),
                        "replicate": replicate,
                        "statistic": permuted["statistic"],
                        "p_value": permuted["p_value"],
                        "degrees_of_freedom": permuted["degrees_of_freedom"],
                        "denominator_degrees_of_freedom": permuted.get(
                            "denominator_degrees_of_freedom", np.nan
                        ),
                        "bic_log_bayes_factor": permuted.get(
                            "bic_log_bayes_factor", np.nan
                        ),
                        "mixture_log_bayes_factor": permuted.get(
                            "mixture_log_bayes_factor", np.nan
                        ),
                        "log_bayes_factors_json": json.dumps(
                            np.asarray(
                                permuted.get("log_bayes_factors", [])
                            ).tolist()
                        ),
                        "tested_levels_json": json.dumps(
                            list(map(str, tested_levels))
                        ),
                        "path_effects_json": json.dumps(
                            permuted_path_effects.tolist()
                        ),
                    })
                except Exception as exc:
                    failures.append({
                        "test_id": test_id,
                        "block_id": block_id,
                        "null_type": "label_permutation",
                        "replicate": replicate,
                        "error": repr(exc),
                    })
            usable = estimates["usable"]
            design = np.column_stack((nuisance, tested))[usable]
            tested_columns = np.arange(nuisance.shape[1], design.shape[1])
            for replicate in range(args.wild_null_replicates):
                try:
                    rng = candidate_rng(
                        args.seed, test_id, replicate, "wild_cluster"
                    )
                    null_values = differential.restricted_wild_cluster_null(
                        estimates["values"][usable],
                        estimates["covariances"][usable],
                        nuisance[usable],
                        clusters[usable],
                        rng,
                    )
                    if args.inference == "wald":
                        null_result = (
                            differential.clustered_multivariate_wald_test(
                                null_values,
                                estimates["covariances"][usable],
                                design,
                                tested_columns,
                                clusters[usable],
                                cluster_adjustment=args.cluster_adjustment,
                            )
                        )
                    else:
                        null_result = differential.gaussian_mixed_model_lrt(
                            null_values,
                            estimates["covariances"][usable],
                            design,
                            tested_columns,
                            clusters[usable],
                            effect_prior_scales=effect_prior_scales,
                            effect_prior_covariance=effect_prior_covariance,
                            random_intercept_solver=args.gaussian_solver,
                        )
                    _, null_path_effects, _ = effect_summary(
                        null_result, nuisance.shape[1], len(signatures)
                    )
                    null_rows.append({
                        "test_id": test_id,
                        "null_type": "wild_cluster",
                        "cluster_adjustment": (
                            args.cluster_adjustment
                            if args.inference == "wald" else ""
                        ),
                        "replicate": replicate,
                        "statistic": null_result["statistic"],
                        "p_value": null_result["p_value"],
                        "degrees_of_freedom": null_result["degrees_of_freedom"],
                        "denominator_degrees_of_freedom": null_result.get(
                            "denominator_degrees_of_freedom", np.nan
                        ),
                        "bic_log_bayes_factor": null_result.get(
                            "bic_log_bayes_factor", np.nan
                        ),
                        "mixture_log_bayes_factor": null_result.get(
                            "mixture_log_bayes_factor", np.nan
                        ),
                        "log_bayes_factors_json": json.dumps(
                            np.asarray(
                                null_result.get("log_bayes_factors", [])
                            ).tolist()
                        ),
                        "tested_levels_json": json.dumps(
                            list(map(str, tested_levels))
                        ),
                        "path_effects_json": json.dumps(
                            null_path_effects.tolist()
                        ),
                    })
                except Exception as exc:
                    failures.append({
                        "test_id": test_id,
                        "block_id": block_id,
                        "null_type": "wild_cluster",
                        "replicate": replicate,
                        "error": repr(exc),
                    })
        except Exception as exc:
            failures.append({
                "test_id": test_id,
                "block_id": block_id,
                "error": repr(exc),
            })
        if (position + 1) % 20 == 0:
            print(json.dumps({
                "tests_complete": position + 1,
                "tests_total": len(candidates),
                "observed_fits": len(observed_rows),
                "failures": len(failures),
                "seconds": time.perf_counter() - started,
            }), flush=True)
    pd.DataFrame(observed_rows).to_csv(
        args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False
    )
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "null.tsv", sep="\t", index=False
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidate_tests": len(candidates),
        "observed_fits": len(observed_rows),
        "failures": len(failures),
        "method": (
            "path_wald" if args.inference == "wald"
            else "path_gaussian_lrt_bf"
        ),
        "inference": args.inference,
        "effect_prior_scales": list(map(float, effect_prior_scales)),
        "gaussian_solver": args.gaussian_solver,
        "covariance_method": args.covariance,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
