#!/usr/bin/env python3
"""Fit subject-paired tests to EC-derived local splice-path compositions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.special import softmax

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    partition_candidates,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import differential, ec_block_glmm
from tealeaf.sc.path_score import paired_path_score_test, paired_subject_centered_test
from tealeaf.sc.path_bias import paired_null_corrected_path_test
from tealeaf.sc.path_score_mixed import mixed_path_score_test, mixed_score_test, MODEL_VERSION as MIXED_SCORE_VERSION


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--path-pseudocount", type=float, default=0.0)
    parser.add_argument(
        "--path-prior-center",
        choices=("uniform", "baseline"),
        default="uniform",
    )
    parser.add_argument(
        "--path-pseudocount-scaling",
        choices=("per_path", "total"),
        default="per_path",
    )
    parser.add_argument("--smoothing-map", type=Path)
    parser.add_argument("--retain-uncertainty", action="store_true")
    parser.add_argument("--uncertainty-scale", type=float, default=1.0)
    parser.add_argument("--uncertainty-scale-map", type=Path)
    parser.add_argument(
        "--uncertainty-scale-grid",
        help="Semicolon-separated scales for an EB profile output.",
    )
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--test-id-file", type=Path)
    parser.add_argument("--export-path-usage", action="store_true")
    parser.add_argument("--paired-inference", choices=("local", "ec-score", "subject-centered", "null-corrected", "mixed-score"), default="local", help="Experimental paired inference alternatives; local remains production.")
    parser.add_argument("--score-null-concentration", type=float, default=1.)
    parser.add_argument("--score-free-null", action="store_true")
    parser.add_argument("--path-fitter", choices=("fixed", "profiled", "free"), default="fixed", help="Local path fit: fixed nuisance (production), profiled block mass, or free transcript weights.")
    return parser.parse_args()


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
    if settings.get("test_effect", "").startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        levels = represented.index[
            represented >= int(settings["min_celltype_mice"])
        ]
        retained = metadata.cell_type.isin(levels).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    return metadata, counts, genes, gene_transcripts, gene_ecs, designs


def signed_null_p_value(
    differences, covariances, rng, retain_uncertainty, uncertainty_scale
):
    signs = rng.choice((-1.0, 1.0), size=(len(differences), 1))
    if retain_uncertainty:
        return differential.paired_measurement_error_test(
            differences * signs,
            covariances,
            uncertainty_scale=uncertainty_scale,
        )["p_value"]
    return differential.paired_mean_test(differences * signs)["p_value"]


def permuted_independent_p_value(
    values,
    covariances,
    labels,
    rng,
    uncertainty_scale,
    biological_variance,
):
    levels = np.unique(labels)
    permuted = rng.permutation(labels)
    design = np.column_stack([
        np.ones(len(values)),
        *[permuted == level for level in levels[1:]],
    ]).astype(float)
    return differential.multivariate_gls_test(
        values,
        (
            np.zeros_like(covariances)
            if uncertainty_scale == 0
            else float(uncertainty_scale) * covariances
        ),
        design,
        np.arange(1, len(levels)),
        biological_variance=biological_variance,
    )["p_value"]


def permuted_blocked_p_value(
    values,
    covariances,
    labels,
    subjects,
    rng,
    uncertainty_scale,
    biological_variance,
):
    permuted = np.asarray(labels).copy()
    subjects = np.asarray(subjects)
    for subject in np.unique(subjects):
        positions = np.flatnonzero(subjects == subject)
        permuted[positions] = rng.permutation(permuted[positions])
    design, tested, _, _ = ec_block_glmm.blocked_multilevel_design(
        permuted, subjects
    )
    return differential.multivariate_gls_test(
        values,
        (
            np.zeros_like(covariances)
            if uncertainty_scale == 0
            else float(uncertainty_scale) * covariances
        ),
        design,
        tested,
        biological_variance=biological_variance,
    )["p_value"]


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    if args.path_fitter != "fixed":
        differential.fit_path_perturbation = {"profiled": differential.fit_profiled_path_perturbation, "free": differential.fit_free_isoform_paths}[args.path_fitter]
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    smoothing = None
    if args.smoothing_map is not None:
        specification = json.loads(args.smoothing_map.read_text())
        map_scaling = specification.get(
            "path_pseudocount_scaling", "per_path"
        )
        if args.path_pseudocount_scaling != map_scaling:
            raise ValueError(
                "smoothing map and requested pseudocount scaling differ"
            )
        if specification.get("selection_scope") == "global":
            smoothing = {
                int(record["gene_fold"]): float(record["alpha"])
                for record in specification["records"]
            }
        else:
            smoothing = {
                (int(record["gene_fold"]), int(record["n_paths"])): float(record["alpha"])
                for record in specification["records"]
            }
        smoothing_folds = int(specification["folds"])
    uncertainty_scaling = None
    if args.uncertainty_scale_map is not None:
        specification = json.loads(args.uncertainty_scale_map.read_text())
        uncertainty_scaling = {
            (int(record["gene_fold"]), int(record["n_paths"])): float(record["uncertainty_scale"])
            for record in specification["records"]
        }
        uncertainty_folds = int(specification["folds"])
    uncertainty_grid = None
    if args.uncertainty_scale_grid:
        uncertainty_grid = np.asarray(
            [float(value) for value in args.uncertainty_scale_grid.split(";")]
        )
        if np.any(uncertainty_grid < 0):
            raise ValueError("uncertainty scale grid must be nonnegative")
    test_effect = settings.get("test_effect")
    if args.paired_inference == "null-corrected" and (test_effect != "cell_type_pairwise" or args.path_prior_center != "uniform" or args.path_pseudocount_scaling != "total" or args.smoothing_map is not None):
        raise ValueError("null-corrected sensitivity requires paired uniform total concentration, without a smoothing map")
    if args.paired_inference == "mixed-score" and (test_effect != "cell_type_pairwise" or args.smoothing_map is not None or args.path_pseudocount != 0):
        raise ValueError("mixed-score sensitivity requires paired contrasts, zero target concentration and no smoothing map")
    if test_effect not in {
        "cell_type",
        "cell_type_pairwise",
        "condition_within_cell_type",
    }:
        raise ValueError("unsupported local path test effect")
    candidates = cached["candidates"]
    if args.test_id_file is not None:
        requested = set(pd.read_csv(args.test_id_file, sep="\t")["test_id"])
        candidates = [candidate for candidate in candidates if candidate[0] in requested]
        missing = sorted(requested - {candidate[0] for candidate in candidates})
        if missing:
            raise ValueError(f"requested test identifiers not found: {missing}")
    if args.max_candidates is not None:
        candidates = candidates[: int(args.max_candidates)]
    candidates = partition_candidates(candidates, args.shard_count)[
        args.shard_index
    ]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(
        args.data_cache, settings
    )
    observed_rows = []
    null_rows = []
    failures = []
    profile_rows = []
    path_usage_rows = []
    pooled_cache = {}
    started = time.perf_counter()
    for candidate in candidates:
        (
            test_id,
            block_id,
            gene_id,
            gene,
            transcripts,
            path_index,
            signatures,
            rows,
            _,
            tested_levels,
        ) = candidate
        try:
            path_pseudocount = args.path_pseudocount
            if smoothing is not None:
                gene_fold = zlib.crc32(str(gene).encode("utf-8")) % smoothing_folds
                if specification.get("selection_scope") == "global":
                    path_pseudocount = smoothing[gene_fold]
                else:
                    path_pseudocount = smoothing[(gene_fold, len(signatures))]
            uncertainty_scale = args.uncertainty_scale
            if uncertainty_scaling is not None:
                gene_fold = zlib.crc32(str(gene).encode("utf-8")) % uncertainty_folds
                uncertainty_scale = uncertainty_scaling[(gene_fold, len(signatures))]
            local_metadata, _, labels = local_test_design(
                metadata, rows, tested_levels, test_effect
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata.mouse.astype(str).to_numpy()
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
            baseline = pooled_cache[cache_key]
            mixed_components = None
            if test_effect == "cell_type_pairwise":
                if args.paired_inference != "local":
                    if args.retain_uncertainty or args.uncertainty_scale_grid or args.uncertainty_scale_map:
                        raise ValueError("experimental paired backends do not support measurement-error testing")
                    if args.paired_inference == "mixed-score":
                        result = mixed_path_score_test(base, path_index, labels, clusters, baseline=baseline, max_iter=args.max_iter, reporting_concentration=1.)
                        mixed_components = result.pop("components")
                        result["path_fits"] = [[SimpleNamespace(path_proportions=proportions) for _, proportions in report] for report in mixed_components.reporting_proportions]
                        result["difference_covariances"] = np.zeros((len(result["differences"]), len(signatures) - 1, len(signatures) - 1))
                    elif args.paired_inference == "null-corrected":
                        result = paired_null_corrected_path_test(base, path_index, labels, clusters, baseline=baseline, concentration=path_pseudocount, max_iter=args.max_iter, reporting_concentration=1.)
                        reported = result["reporting_proportions"].reshape(-1, 2, len(signatures))
                        result["path_fits"] = [[SimpleNamespace(path_proportions=proportion) for proportion in pair] for pair in reported]
                        result["difference_covariances"] = np.zeros((len(result["differences"]), len(signatures) - 1, len(signatures) - 1))
                    elif args.paired_inference == "subject-centered":
                        result = paired_subject_centered_test(base, path_index, labels, clusters, baseline=baseline, concentration=path_pseudocount)
                        result["difference_covariances"] = np.asarray([first.covariance.covariance + second.covariance.covariance for first, second in result["path_fits"]])
                    else:
                        result = paired_path_score_test(base, path_index, labels, clusters, baseline=baseline, denominator_concentration=path_pseudocount, null_concentration=args.score_null_concentration, free_null=args.score_free_null)
                        proxy_fits = []
                        basis = differential.helmert_basis(len(signatures))
                        for fitted, response in zip(result.pop("null_fits"), result["differences"]):
                            anchor = basis.T @ np.log(fitted.path_proportions)
                            proxy_fits.append([SimpleNamespace(path_proportions=softmax(basis @ (anchor - response / 2))), SimpleNamespace(path_proportions=softmax(basis @ (anchor + response / 2)))])
                        result["path_fits"] = proxy_fits
                        # These zeros are unused by paired_mean_test and sign
                        # flips; measurement-error CLI paths are forbidden.
                        result["difference_covariances"] = np.zeros((len(result["differences"]), len(signatures) - 1, len(signatures) - 1))
                else:
                    result = ec_block_glmm.paired_path_test(
                        base,
                        path_index,
                        labels,
                        clusters,
                        baseline=baseline,
                        max_iter=args.max_iter,
                        path_pseudocount=path_pseudocount,
                        path_prior_center=args.path_prior_center,
                        path_pseudocount_scaling=args.path_pseudocount_scaling,
                        retain_uncertainty=args.retain_uncertainty,
                        uncertainty_scale=uncertainty_scale,
                    )
                values = result.pop("differences")
                value_covariances = result.pop("difference_covariances")
                null_labels = None
                null_subjects = None
            elif test_effect == "cell_type":
                result = ec_block_glmm.blocked_multilevel_path_test(
                    base,
                    path_index,
                    labels,
                    clusters,
                    baseline=baseline,
                    max_iter=args.max_iter,
                    path_pseudocount=path_pseudocount,
                    path_prior_center=args.path_prior_center,
                    path_pseudocount_scaling=args.path_pseudocount_scaling,
                    uncertainty_scale=uncertainty_scale,
                )
                values = result.pop("values")
                value_covariances = result.pop("covariances")
                null_labels = result.pop("observation_labels")
                null_subjects = result.pop("observation_subjects")
            else:
                result = ec_block_glmm.independent_path_test(
                    base,
                    path_index,
                    labels,
                    clusters,
                    baseline=baseline,
                    max_iter=args.max_iter,
                    path_pseudocount=path_pseudocount,
                    path_prior_center=args.path_prior_center,
                    path_pseudocount_scaling=args.path_pseudocount_scaling,
                    uncertainty_scale=uncertainty_scale,
                )
                values = result.pop("values")
                value_covariances = result.pop("covariances")
                null_labels = result.pop("subject_labels")
                null_subjects = None
            path_fits = result.pop("path_fits")
            subject_ids = result.pop("subject_ids", None)
            fitted_levels = result.pop("levels")
            if args.export_path_usage and test_effect == "cell_type_pairwise":
                if len(fitted_levels) != len(tested_levels):
                    raise ValueError("fitted and requested cell-type levels differ")
                fitted_level_names = dict(zip(fitted_levels, tested_levels))
                for subject_id, subject_fits in zip(subject_ids, path_fits):
                    for level, fit in zip(fitted_levels, subject_fits):
                        for path_number, (signature, proportion) in enumerate(zip(signatures, fit.path_proportions), start=1):
                            path_usage_rows.append({
                                "test_id": test_id,
                                "block_id": block_id,
                                "gene_id": gene_id,
                                "subject": subject_id,
                                "cell_type": fitted_level_names[level],
                                "path": f"Path {path_number}",
                                "path_number": path_number,
                                "path_signature": json.dumps(signature),
                                "proportion": float(proportion),
                            })
            elif args.export_path_usage and test_effect == "cell_type":
                if len(fitted_levels) != len(tested_levels):
                    raise ValueError("fitted and requested cell-type levels differ")
                fitted_level_names = dict(zip(fitted_levels, tested_levels))
                for subject_id, level, fit in zip(null_subjects, null_labels, path_fits):
                    for path_number, (signature, proportion) in enumerate(zip(signatures, fit.path_proportions), start=1):
                        path_usage_rows.append({
                            "test_id": test_id,
                            "block_id": block_id,
                            "gene_id": gene_id,
                            "subject": subject_id,
                            "cell_type": fitted_level_names[level],
                            "path": f"Path {path_number}",
                            "path_number": path_number,
                            "path_signature": json.dumps(signature),
                            "proportion": float(proportion),
                        })
            result.pop("mean", None)
            result.pop("mean_covariance", None)
            fitted_mean = result.get("mean_difference", values.mean(axis=0) if len(values) else np.array([]))
            observed_rows.append({
                "test_id": test_id,
                "block_id": block_id,
                "gene_id": gene_id,
                "contrast": test_effect,
                "level_a": tested_levels[0] if test_effect == "cell_type_pairwise" else "",
                "level_b": tested_levels[1] if test_effect == "cell_type_pairwise" else "",
                "tested_cell_type": candidate[8] if test_effect == "condition_within_cell_type" else "",
                "method": "local_path",
                "inference_backend": args.paired_inference if test_effect == "cell_type_pairwise" else "local",
                "score_null_concentration": args.score_null_concentration if args.paired_inference == "ec-score" else np.nan,
                "score_free_null": args.score_free_null if args.paired_inference == "ec-score" else False,
                "path_pseudocount": path_pseudocount,
                "path_prior_center": args.path_prior_center,
                "path_pseudocount_scaling": args.path_pseudocount_scaling,
                "retain_uncertainty": args.retain_uncertainty,
                "uncertainty_scale": uncertainty_scale,
                "n_paths": len(signatures),
                "n_isoforms": base.n_isoforms,
                "n_ecs": len(gene_ecs[gene]),
                "n_samples": result.get("n_observations", len(local_metadata)),
                "n_subjects": result["n_subjects"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "median_gene_umis": float(np.median(totals)),
                "statistic": result["statistic"],
                "p_value": result["p_value"],
                "biological_variance": result.get(
                    "biological_variance", np.nan
                ),
                "restricted_objective": result.get(
                    "restricted_objective", np.nan
                ),
                "converged": result["converged"],
                "mean_difference_norm": float(np.linalg.norm(fitted_mean)),
                "mean_difference": json.dumps(np.asarray(fitted_mean).tolist()),
                "denominator_degrees_of_freedom": result.get("denominator_degrees_of_freedom", np.nan),
                "residual_inflation": result.get("residual_inflation", np.nan),
                "chi_square_p_value": result.get("chi_square_p_value", np.nan),
                "path_signatures": json.dumps(signatures),
            })
            if result["converged"]:
                if uncertainty_grid is not None:
                    for scale in uncertainty_grid:
                        if test_effect == "cell_type_pairwise":
                            profiled = differential.paired_measurement_error_test(
                                values,
                                value_covariances,
                                uncertainty_scale=scale,
                            )
                        elif test_effect == "condition_within_cell_type":
                            levels = np.unique(null_labels)
                            profile_design = np.column_stack([
                                np.ones(len(values)),
                                *[null_labels == level for level in levels[1:]],
                            ]).astype(float)
                            profiled = differential.multivariate_gls_test(
                                values,
                                scale * value_covariances,
                                profile_design,
                                np.arange(1, len(levels)),
                            )
                        else:
                            profile_design, profile_tested, _, _ = (
                                ec_block_glmm.blocked_multilevel_design(
                                    null_labels, null_subjects
                                )
                            )
                            profiled = differential.multivariate_gls_test(
                                values,
                                scale * value_covariances,
                                profile_design,
                                profile_tested,
                            )
                        profile_rows.append({
                            "test_id": test_id,
                            "gene": gene,
                            "n_paths": len(signatures),
                            "uncertainty_scale": scale,
                            "restricted_objective": profiled["restricted_objective"],
                        })
                test_hash = zlib.crc32(test_id.encode("utf-8"))
                for replicate in range(args.null_replicates):
                    rng = np.random.default_rng(
                        np.random.SeedSequence((args.seed, test_hash, replicate))
                    )
                    if mixed_components is not None:
                        signs = rng.choice([-1., 1.], size=(len(mixed_components.scores), 1))
                        null_p_value = mixed_score_test(mixed_components.scores * signs, mixed_components.information, mixed_components.biological_shapes)["p_value"]
                    elif test_effect == "cell_type_pairwise":
                        null_p_value = signed_null_p_value(
                            values,
                            value_covariances,
                            rng,
                            args.retain_uncertainty,
                            uncertainty_scale,
                        )
                    elif test_effect == "condition_within_cell_type":
                        null_p_value = permuted_independent_p_value(
                            values,
                            value_covariances,
                            null_labels,
                            rng,
                            uncertainty_scale,
                            result["biological_variance"],
                        )
                    else:
                        null_p_value = permuted_blocked_p_value(
                            values,
                            value_covariances,
                            null_labels,
                            null_subjects,
                            rng,
                            uncertainty_scale,
                            result["biological_variance"],
                        )
                    null_rows.append({
                        "test_id": test_id,
                        "block_id": block_id,
                        "replicate": replicate,
                        "p_value": null_p_value,
                    })
        except Exception as error:
            failures.append({"test_id": test_id, "error": repr(error)})
            if args.paired_inference in ("null-corrected", "mixed-score"):
                # Unlike historical runners, this prototype explicitly retains
                # every failed candidate in the tested/BH/Simes family.
                observed_rows = [row for row in observed_rows if row["test_id"] != test_id]
                null_rows = [row for row in null_rows if row["test_id"] != test_id]
                path_usage_rows = [row for row in path_usage_rows if row["test_id"] != test_id]
                observed_rows.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "contrast": test_effect, "level_a": tested_levels[0], "level_b": tested_levels[1], "method": "local_path", "inference_backend": args.paired_inference, "path_pseudocount": args.path_pseudocount, "path_prior_center": args.path_prior_center, "path_pseudocount_scaling": args.path_pseudocount_scaling, "n_paths": len(signatures), "n_isoforms": len(transcripts), "n_ecs": len(gene_ecs[gene]), "n_samples": len(rows), "n_subjects": 0, "degrees_of_freedom": len(signatures) - 1, "median_gene_umis": np.nan, "statistic": 0., "p_value": 1., "converged": False, "mean_difference_norm": np.nan, "mean_difference": "[]", "path_signatures": json.dumps(signatures)})
                null_rows.extend({"test_id": test_id, "block_id": block_id, "replicate": replicate, "p_value": 1.} for replicate in range(args.null_replicates))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(observed_rows).to_csv(
        args.output_dir / "paired_path.tsv", sep="\t", index=False
    )
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False
    )
    pd.DataFrame(profile_rows).to_csv(
        args.output_dir / "uncertainty_profiles.tsv", sep="\t", index=False
    )
    if args.export_path_usage:
        pd.DataFrame(path_usage_rows).to_csv(
            args.output_dir / "path_usage.tsv", sep="\t", index=False
        )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidates": len(candidates),
        "completed": len(observed_rows),
        "failures": len(failures),
        "elapsed_seconds": time.perf_counter() - started,
        "inference_backend": args.paired_inference,
        "experimental": args.paired_inference != "local",
        "candidate_settings": settings,
        "score_report": "One-step ILR displacement mapped through subject-null softmax, not separately quantified cell-type usage" if args.paired_inference == "ec-score" else None,
        "null_corrected_report": "Independent free-transcript total-A1 arithmetic subject usage; statistics use strong smoothed-minus-expected-null proportion differences" if args.paired_inference == "null-corrected" else None,
        "mixed_score_report": "Independent free-transcript total-A1 arithmetic subject usage; inference uses type-nuisance-profiled EC scores and REML uncertainty-weighted reference ILR contrasts, with experimental small-sample F tails" if args.paired_inference == "mixed-score" else None,
        "mixed_score_version": MIXED_SCORE_VERSION if args.paired_inference == "mixed-score" else None,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
