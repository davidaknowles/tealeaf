#!/usr/bin/env python3
"""Fit and bootstrap-calibrate splice-block EC-count GLMM DS tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd
import scipy.stats

from extra_scripts.run_compositional_splicing import block_equivalence
from extra_scripts.run_differential_splicing import block_mapping, load_blocks
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full


METHODS = (
    "multinomial_full",
    "multinomial_noise_full",
    "dirichlet_multinomial_full",
    "laplace_multinomial",
    "laplace_multinomial_noise",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument(
        "--candidate-cache",
        type=Path,
        help="Cache the screened genome-wide block/test manifest.",
    )
    parser.add_argument(
        "--event-table",
        type=Path,
        help="Optionally restrict tests to inference-eligible rows in this table.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--methods", nargs="+", choices=METHODS, default=METHODS[:3]
    )
    parser.add_argument(
        "--calibration",
        choices=("bootstrap", "lrt_bic"),
        default="bootstrap",
    )
    parser.add_argument(
        "--test-effect",
        choices=(
            "cell_type",
            "cell_type_pairwise",
            "condition_within_cell_type",
        ),
        default="cell_type",
    )
    parser.add_argument("--null-replicates", type=int, default=20)
    parser.add_argument("--min-gene-umis", type=float, default=10.0)
    parser.add_argument("--min-gene-samples", type=int, default=0)
    parser.add_argument("--min-cell-types", type=int, default=2)
    parser.add_argument("--min-celltype-mice", type=int, default=3)
    parser.add_argument("--min-conditions", type=int, default=2)
    parser.add_argument("--min-condition-mice", type=int, default=3)
    parser.add_argument("--max-isoforms", type=int, default=10)
    parser.add_argument("--max-ecs", type=int, default=128)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--null-replicate-retries", type=int)
    parser.add_argument("--vi-samples", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--subject-folds", type=Path)
    parser.add_argument("--subject-fold", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def group_metadata(groups):
    return pd.DataFrame(
        [str(group).rsplit("__", 2) for group in groups],
        columns=["cell_type", "condition", "mouse"],
    )


def covered_celltype_design(
    metadata,
    gene_umis,
    *,
    min_gene_umis,
    min_samples,
    min_cell_types,
    min_celltype_mice,
):
    """Build a gene-specific design from EC-covered pseudobulks."""
    covered = np.asarray(gene_umis) >= float(min_gene_umis)
    if covered.sum() < int(min_samples):
        return None
    represented = metadata.loc[covered].groupby("cell_type")["mouse"].nunique()
    cell_types = sorted(
        represented.index[represented >= int(min_celltype_mice)]
    )
    if len(cell_types) < int(min_cell_types):
        return None
    retained = covered & metadata["cell_type"].isin(cell_types).to_numpy()
    rows = np.flatnonzero(retained)
    metadata = metadata.iloc[rows].reset_index(drop=True)
    nuisance = pd.get_dummies(metadata["condition"], dtype=float).to_numpy()
    cell_type_index = {value: index for index, value in enumerate(cell_types)}
    labels = np.asarray([cell_type_index[value] for value in metadata["cell_type"]])
    subjects = (
        metadata["condition"].astype(str)
        + "__"
        + metadata["mouse"].astype(str)
    ).to_numpy()
    return rows, metadata, nuisance, labels, subjects, cell_types


def covered_celltype_pairwise_designs(
    metadata,
    gene_umis,
    *,
    min_gene_umis,
    min_samples,
    min_celltype_mice,
):
    """Build paired two-cell-type designs from gene-covered subjects."""
    covered = np.asarray(gene_umis) >= float(min_gene_umis)
    cell_types = sorted(metadata.loc[covered, "cell_type"].unique())
    results = []
    for first_index, first in enumerate(cell_types):
        first_subjects = set(
            metadata.loc[covered & metadata.cell_type.eq(first), "mouse"]
        )
        for second in cell_types[first_index + 1 :]:
            second_subjects = set(
                metadata.loc[covered & metadata.cell_type.eq(second), "mouse"]
            )
            shared = first_subjects & second_subjects
            if len(shared) < int(min_celltype_mice):
                continue
            retained = (
                covered
                & metadata.cell_type.isin([first, second]).to_numpy()
                & metadata.mouse.isin(shared).to_numpy()
            )
            if retained.sum() < max(int(min_samples), 2 * len(shared)):
                continue
            rows = np.flatnonzero(retained)
            local = metadata.iloc[rows].reset_index(drop=True)
            subject_counts = local.groupby(["mouse", "cell_type"]).size()
            if len(subject_counts) != 2 * len(shared) or not subject_counts.eq(1).all():
                continue
            nuisance = pd.get_dummies(local["condition"], dtype=float).to_numpy()
            labels = local.cell_type.eq(second).to_numpy(dtype=int)
            subjects = local.mouse.astype(str).to_numpy()
            coverage = (
                rows,
                local,
                nuisance,
                labels,
                subjects,
                (first, second),
            )
            results.append((coverage, None))
    return results


def covered_condition_designs(
    metadata,
    gene_umis,
    *,
    min_gene_umis,
    min_samples,
    min_conditions,
    min_condition_mice,
):
    """Build condition designs separately within each covered cell type."""
    covered = np.asarray(gene_umis) >= float(min_gene_umis)
    results = []
    for cell_type in sorted(metadata["cell_type"].unique()):
        in_cell_type = metadata["cell_type"].eq(cell_type).to_numpy()
        available = covered & in_cell_type
        represented = (
            metadata.loc[available].groupby("condition")["mouse"].nunique()
        )
        conditions = sorted(
            represented.index[represented >= int(min_condition_mice)]
        )
        if len(conditions) < int(min_conditions):
            continue
        retained = (
            available
            & metadata["condition"].isin(conditions).to_numpy()
        )
        if retained.sum() < int(min_samples):
            continue
        rows = np.flatnonzero(retained)
        local = metadata.iloc[rows].reset_index(drop=True)
        condition_index = {
            value: index for index, value in enumerate(conditions)
        }
        labels = np.asarray([
            condition_index[value] for value in local["condition"]
        ])
        nuisance = np.ones((len(local), 1), dtype=float)
        subjects = local["mouse"].astype(str).to_numpy()
        coverage = (rows, local, nuisance, labels, subjects, conditions)
        results.append((coverage, cell_type))
    return results


def treatment_design(labels, n_levels):
    labels = np.asarray(labels, dtype=int)
    result = np.zeros((len(labels), int(n_levels) - 1), dtype=float)
    nonreference = labels > 0
    result[np.flatnonzero(nonreference), labels[nonreference] - 1] = 1.0
    return result


def modeled_gene_umis(counts, designs, ecs, transcripts):
    """Count the primer-specific EC rows retained by the gene likelihood."""
    ecs = np.asarray(ecs)
    total = np.zeros(counts[0].shape[0], dtype=float)
    for observed, mapping in zip(counts, designs):
        supported = np.asarray(
            mapping[ecs][:, transcripts].sum(axis=1)
        ).ravel() > 0
        total += np.asarray(observed[:, ecs[supported]].sum(axis=1)).ravel()
    return total


def candidate_cache_settings(args):
    """Return the screening settings that define a candidate manifest."""
    return {
        "version": 3,
        "data_cache": str(args.data_cache.resolve()),
        "features": str(args.features.resolve()),
        "block_cache": str(args.block_cache.resolve()),
        "event_table": (
            None if args.event_table is None else str(args.event_table.resolve())
        ),
        "test_effect": args.test_effect,
        "min_gene_umis": args.min_gene_umis,
        "min_gene_samples": args.min_gene_samples,
        "min_cell_types": args.min_cell_types,
        "min_celltype_mice": args.min_celltype_mice,
        "min_conditions": args.min_conditions,
        "min_condition_mice": args.min_condition_mice,
        "max_isoforms": args.max_isoforms,
        "max_ecs": args.max_ecs,
        "subject_folds": (
            None if args.subject_folds is None else str(args.subject_folds.resolve())
        ),
        "subject_fold": args.subject_fold,
    }


def supported_partition_key(candidate):
    """Identify tests made equivalent by the EC-supported isoform subset."""
    (
        _, _, gene_id, _, transcripts, path_index, _, rows,
        tested_cell_type, tested_levels,
    ) = candidate
    relabeling = {}
    canonical_paths = []
    for value in np.asarray(path_index, dtype=int):
        if value < 0:
            canonical_paths.append(-1)
        else:
            canonical_paths.append(relabeling.setdefault(value, len(relabeling)))
    return (
        gene_id,
        tuple(np.asarray(transcripts, dtype=int)),
        tuple(canonical_paths),
        tuple(np.asarray(rows, dtype=int)),
        tested_cell_type,
        tuple(tested_levels),
    )


def deduplicate_supported_partitions(candidates):
    """Keep one annotated block for each identifiable tested partition."""
    result = []
    seen = set()
    for candidate in candidates:
        key = supported_partition_key(candidate)
        if key not in seen:
            result.append(candidate)
            seen.add(key)
    return result


def local_test_design(metadata, rows, tested_levels, test_effect):
    """Reconstruct one cached candidate's fixed-effect design."""
    local = metadata.iloc[np.asarray(rows, dtype=int)].reset_index(drop=True)
    tested_column = (
        "cell_type" if test_effect.startswith("cell_type") else "condition"
    )
    level_index = {value: index for index, value in enumerate(tested_levels)}
    labels = np.asarray([level_index[value] for value in local[tested_column]])
    if test_effect.startswith("cell_type"):
        nuisance = pd.get_dummies(local["condition"], dtype=float).to_numpy()
    else:
        nuisance = np.ones((len(local), 1), dtype=float)
    return local, nuisance, labels


def tensor_data(base, tensor):
    return ec_glmm.ECGLMMData(
        base.counts,
        base.compatibility,
        np.ones((len(base.clusters), 1), dtype=float),
        base.clusters,
        fixed_effect_tensor=tensor,
    )


def fit_method(method, data, args, initial=None):
    if method.startswith("laplace_"):
        return ec_glmm.fit_laplace(
            data,
            family="multinomial",
            observation_noise=method == "laplace_multinomial_noise",
            initial=initial,
            max_iter=args.max_iter,
        )
    return ec_glmm_full.fit_variational(
        data,
        family=(
            "dirichlet_multinomial"
            if method == "dirichlet_multinomial_full"
            else "multinomial"
        ),
        objective=(
            "monte_carlo"
            if method == "dirichlet_multinomial_full"
            else "tilted"
        ),
        observation_noise=method == "multinomial_noise_full",
        samples=args.vi_samples,
        seed=args.seed,
        initial=initial,
        max_iter=args.max_iter,
    )


def fit_with_retries(method, data, args, initial=None, retries=None):
    """Continue fits that reach a stopping failure without changing tolerances."""
    fit = fit_method(method, data, args, initial=initial)
    total_iterations = int(fit["iterations"])
    attempts = 1
    retries = args.retries if retries is None else int(retries)
    while not fit["converged"] and attempts <= retries:
        fit = fit_method(method, data, args, initial=fit["parameters"])
        total_iterations += int(fit["iterations"])
        attempts += 1
    fit["total_iterations"] = total_iterations
    fit["attempts"] = attempts
    return fit


def main():
    args = parse_args()
    if args.calibration == "lrt_bic" and any(
        not method.startswith("laplace_") for method in args.methods
    ):
        raise ValueError("lrt_bic calibration requires Laplace methods")
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(
            handle
        )
    features = np.loadtxt(args.features, dtype=str)
    if len(features) != designs[0].shape[1]:
        raise ValueError(
            f"feature file has {len(features)} rows but EC designs have "
            f"{designs[0].shape[1]} transcript columns"
        )
    metadata = group_metadata(groups)
    if args.subject_folds is not None:
        if args.subject_fold is None:
            raise ValueError("--subject-fold is required with --subject-folds")
        folds = pd.read_csv(args.subject_folds, sep="\t", dtype={"subject": str})
        selected = set(
            folds.loc[folds.fold.eq(args.subject_fold), "subject"].astype(str)
        )
        selected_rows = metadata["mouse"].astype(str).isin(selected).to_numpy()
        if not selected_rows.any():
            raise ValueError("subject fold selects no EC-count pseudobulks")
        metadata = metadata.loc[selected_rows].reset_index(drop=True)
        counts = tuple(matrix[selected_rows] for matrix in counts)
    if args.test_effect.startswith("cell_type"):
        globally_represented = metadata.groupby("cell_type")["mouse"].nunique()
        global_cell_types = globally_represented.index[
            globally_represented >= int(args.min_celltype_mice)
        ]
        global_rows = metadata["cell_type"].isin(global_cell_types).to_numpy()
        metadata = metadata.loc[global_rows].reset_index(drop=True)
        counts = tuple(matrix[global_rows] for matrix in counts)
    cache_settings = candidate_cache_settings(args)
    if args.candidate_cache is not None and args.candidate_cache.exists():
        with args.candidate_cache.open("rb") as handle:
            cached = pickle.load(handle)
        if cached.get("settings") != cache_settings:
            raise ValueError("candidate cache does not match screening settings")
        candidates = cached["candidates"]
    else:
        screening_counts = tuple(matrix.tocsc() for matrix in counts)
        blocks = load_blocks(args.block_cache, None)
        block_to_representative, _ = block_equivalence(args.block_cache)
        eligible_blocks = None
        if args.event_table is not None:
            event_table = pd.read_csv(args.event_table, sep="\t")
            eligible_blocks = set(
                event_table.loc[
                    event_table["inference_eligible"], "block_id"
                ].astype(str)
            )
        gene_position = {
            str(gene).split(".", 1)[0]: index
            for index, gene in enumerate(genes)
        }
        gene_information = {}
        candidates = []
        for block in blocks:
            representative = block_to_representative[block.block_id]
            if block.block_id != representative:
                continue
            if eligible_blocks is not None and representative not in eligible_blocks:
                continue
            gene = gene_position.get(block.gene_id.split(".", 1)[0])
            if (
                gene is None
                or not len(gene_ecs[gene])
                or len(gene_ecs[gene]) > args.max_ecs
            ):
                continue
            if gene not in gene_information:
                transcripts = np.asarray(gene_transcripts[gene], dtype=int)
                transcript_support = sum(
                    np.asarray(
                        mapping[gene_ecs[gene]][:, transcripts].sum(axis=0)
                    ).ravel()
                    for mapping in designs
                )
                transcripts = transcripts[transcript_support > 0]
                if not 2 <= len(transcripts) <= args.max_isoforms:
                    gene_information[gene] = None
                    continue
                gene_umis = modeled_gene_umis(
                    screening_counts, designs, gene_ecs[gene], transcripts
                )
                if args.test_effect == "cell_type":
                    coverage = covered_celltype_design(
                        metadata,
                        gene_umis,
                        min_gene_umis=args.min_gene_umis,
                        min_samples=args.min_gene_samples,
                        min_cell_types=args.min_cell_types,
                        min_celltype_mice=args.min_celltype_mice,
                    )
                    coverage_specs = (
                        [(coverage, None)] if coverage is not None else []
                    )
                elif args.test_effect == "cell_type_pairwise":
                    coverage_specs = covered_celltype_pairwise_designs(
                        metadata,
                        gene_umis,
                        min_gene_umis=args.min_gene_umis,
                        min_samples=args.min_gene_samples,
                        min_celltype_mice=args.min_celltype_mice,
                    )
                else:
                    coverage_specs = covered_condition_designs(
                        metadata,
                        gene_umis,
                        min_gene_umis=args.min_gene_umis,
                        min_samples=args.min_gene_samples,
                        min_conditions=args.min_conditions,
                        min_condition_mice=args.min_condition_mice,
                    )
                gene_information[gene] = (transcripts, coverage_specs)
            information = gene_information[gene]
            if information is None:
                continue
            transcripts, coverage_specs = information
            if not coverage_specs:
                continue
            mapped = block_mapping(block, features[transcripts])
            if mapped is None:
                continue
            path_index, signatures = mapped
            for coverage, tested_cell_type in coverage_specs:
                rows, _, _, _, _, tested_levels = coverage
                test_id = block.block_id
                if tested_cell_type is not None:
                    test_id += f"|condition|{tested_cell_type}"
                elif args.test_effect == "cell_type_pairwise":
                    test_id += "|cell_type|" + "|".join(tested_levels)
                candidates.append((
                    test_id,
                    block.block_id,
                    block.gene_id,
                    gene,
                    transcripts,
                    path_index,
                    signatures,
                    rows,
                    tested_cell_type,
                    tuple(tested_levels),
                ))
        candidates = deduplicate_supported_partitions(candidates)
        if args.candidate_cache is not None:
            args.candidate_cache.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.candidate_cache.with_suffix(
                args.candidate_cache.suffix + ".tmp"
            )
            with temporary.open("wb") as handle:
                pickle.dump(
                    {"settings": cache_settings, "candidates": candidates},
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            temporary.replace(args.candidate_cache)
    if not candidates:
        raise ValueError(
            "no candidate block tests passed screening; check annotation ID "
            "compatibility and coverage thresholds"
        )
    candidates = candidates[args.shard_index :: args.shard_count]
    if args.dry_run:
        print(json.dumps({
            "candidate_tests": len(candidates),
            "candidate_blocks": len({row[1] for row in candidates}),
            "shard_index": args.shard_index,
            "shard_count": args.shard_count,
        }, indent=2))
        return

    observed_rows = []
    null_rows = []
    failures = []
    null_cache = {}
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
                metadata, rows, tested_levels, args.test_effect
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata["mouse"].to_numpy()
            base, _, totals = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(local_metadata), 1)),
                clusters,
                drop_zero=False,
            )
            tested = treatment_design(labels, len(tested_levels))
            null_tensor, alternative_tensor, degrees = (
                ec_block_glmm.block_fixed_effect_tensors(
                    nuisance, tested, path_index
                )
            )
            null_data = tensor_data(base, null_tensor)
            completed_methods = {}
            for method in args.methods:
                cache_key = (test_id, method)
                if cache_key not in null_cache:
                    null_initial = None
                    if (
                        method == "laplace_multinomial_noise"
                        and "laplace_multinomial" in completed_methods
                    ):
                        null_initial = np.r_[
                            completed_methods["laplace_multinomial"][
                                "null"
                            ]["parameters"],
                            np.log(0.3),
                        ]
                    null_cache[cache_key] = fit_with_retries(
                        method, null_data, args, initial=null_initial
                    )
                null = null_cache[cache_key]
                alternative_data = tensor_data(base, alternative_tensor)
                if (
                    method == "laplace_multinomial_noise"
                    and "laplace_multinomial" in completed_methods
                ):
                    initial = np.r_[
                        completed_methods["laplace_multinomial"][
                            "alternative"
                        ]["parameters"],
                        np.log(0.3),
                    ]
                else:
                    initial = ec_glmm_full.fixed_effect_warm_start(
                        null, alternative_tensor.shape[2]
                    )
                alternative = fit_with_retries(
                    method, alternative_data, args, initial=initial
                )
                completed_methods[method] = {
                    "null": null,
                    "alternative": alternative,
                }
                if method.startswith("laplace_"):
                    statistic = 2.0 * (
                        null["objective"] - alternative["objective"]
                    )
                else:
                    statistic = 2.0 * (
                        alternative["objective"] - null["objective"]
                    )
                statistic = max(0.0, float(statistic))
                lrt_p_value = (
                    float(scipy.stats.chi2.sf(statistic, degrees))
                    if args.calibration == "lrt_bic"
                    else np.nan
                )
                bic_log_bayes_factor = (
                    0.5 * (
                        statistic - degrees * np.log(len(np.unique(clusters)))
                    )
                    if args.calibration == "lrt_bic"
                    else np.nan
                )
                observed_rows.append({
                    "test_id": test_id,
                    "block_id": block_id,
                    "gene_id": gene_id,
                    "contrast": args.test_effect,
                    "cell_type": tested_cell_type,
                    "level_a": (
                        tested_levels[0]
                        if args.test_effect == "cell_type_pairwise"
                        else None
                    ),
                    "level_b": (
                        tested_levels[1]
                        if args.test_effect == "cell_type_pairwise"
                        else None
                    ),
                    "method": method,
                    "n_paths": len(signatures),
                    "n_isoforms": len(transcripts),
                    "n_ecs": len(gene_ecs[gene]),
                    "n_samples": len(local_metadata),
                    "n_test_levels": len(tested_levels),
                    "degrees_of_freedom": degrees,
                    "median_gene_umis": float(np.median(totals)),
                    "statistic": statistic,
                    "lrt_p_value": lrt_p_value,
                    "bic_log_bayes_factor": bic_log_bayes_factor,
                    "null_converged": null["converged"],
                    "alternative_converged": alternative["converged"],
                    "null_objective": null["objective"],
                    "alternative_objective": alternative["objective"],
                    "alternative_observation_noise_sd": alternative[
                        "observation_noise_sd"
                    ],
                    "alternative_concentration": alternative["concentration"],
                    "null_iterations": null["total_iterations"],
                    "alternative_iterations": alternative["total_iterations"],
                    "null_gradient_norm": null["gradient_norm"],
                    "alternative_gradient_norm": alternative["gradient_norm"],
                    "null_scaled_gradient_norm": null.get(
                        "scaled_gradient_norm", np.nan
                    ),
                    "alternative_scaled_gradient_norm": alternative.get(
                        "scaled_gradient_norm", np.nan
                    ),
                    "null_optimizer_scale": null.get("optimizer_scale", np.nan),
                    "alternative_optimizer_scale": alternative.get(
                        "optimizer_scale", np.nan
                    ),
                    "null_mode_score_norm": null.get("mode_score_norm", np.nan),
                    "alternative_mode_score_norm": alternative.get(
                        "mode_score_norm", np.nan
                    ),
                    "null_attempts": null["attempts"],
                    "alternative_attempts": alternative["attempts"],
                })
                if args.calibration != "bootstrap":
                    continue
                bootstrap_rng = np.random.default_rng(np.random.SeedSequence((
                    args.seed,
                    zlib.crc32(test_id.encode("utf-8")),
                    METHODS.index(method),
                )))
                for replicate in range(args.null_replicates):
                    simulated_counts = ec_block_glmm.simulate_null_counts(
                        null_data,
                        null,
                        bootstrap_rng,
                        family=(
                            "dirichlet_multinomial"
                            if method == "dirichlet_multinomial_full"
                            else "multinomial"
                        ),
                        observation_noise=method in {
                            "multinomial_noise_full",
                            "laplace_multinomial_noise",
                        },
                    )
                    simulated_null_data = ec_glmm.ECGLMMData(
                        simulated_counts,
                        base.compatibility,
                        np.ones((len(base.clusters), 1), dtype=float),
                        base.clusters,
                        fixed_effect_tensor=null_tensor,
                    )
                    simulated_null = fit_with_retries(
                        method,
                        simulated_null_data,
                        args,
                        initial=null["parameters"],
                        retries=args.null_replicate_retries,
                    )
                    simulated_alternative_data = ec_glmm.ECGLMMData(
                        simulated_counts,
                        base.compatibility,
                        np.ones((len(base.clusters), 1), dtype=float),
                        base.clusters,
                        fixed_effect_tensor=alternative_tensor,
                    )
                    simulated_initial = ec_glmm_full.fixed_effect_warm_start(
                        simulated_null, alternative_tensor.shape[2]
                    )
                    simulated_fit = fit_with_retries(
                        method,
                        simulated_alternative_data,
                        args,
                        initial=simulated_initial,
                        retries=args.null_replicate_retries,
                    )
                    null_rows.append({
                        "test_id": test_id,
                        "block_id": block_id,
                        "contrast": args.test_effect,
                        "cell_type": tested_cell_type,
                        "level_a": (
                            tested_levels[0]
                            if args.test_effect == "cell_type_pairwise"
                            else None
                        ),
                        "level_b": (
                            tested_levels[1]
                            if args.test_effect == "cell_type_pairwise"
                            else None
                        ),
                        "method": method,
                        "replicate": replicate,
                        "degrees_of_freedom": degrees,
                        "median_gene_umis": float(np.median(totals)),
                        "statistic": float(
                            2.0
                            * (
                                simulated_fit["objective"]
                                - simulated_null["objective"]
                            )
                        ),
                        "converged": (
                            simulated_null["converged"]
                            and simulated_fit["converged"]
                        ),
                        "null_converged": simulated_null["converged"],
                        "alternative_converged": simulated_fit["converged"],
                        "null_iterations": simulated_null["total_iterations"],
                        "iterations": simulated_fit["total_iterations"],
                        "null_gradient_norm": simulated_null["gradient_norm"],
                        "gradient_norm": simulated_fit["gradient_norm"],
                        "null_attempts": simulated_null["attempts"],
                        "attempts": simulated_fit["attempts"],
                    })
        except Exception as exc:
            failures.append({
                "test_id": test_id,
                "block_id": block_id,
                "error": repr(exc),
            })
        if (position + 1) % 2 == 0:
            print(json.dumps({
                "tests_complete": position + 1,
                "tests_total": len(candidates),
                "observed_fits": len(observed_rows),
                "null_fits": len(null_rows),
                "failures": len(failures),
                "seconds": time.perf_counter() - started,
            }), flush=True)
    pd.DataFrame(observed_rows).to_csv(
        args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False
    )
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "bootstrap_null.tsv.gz", sep="\t", index=False
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidate_tests": len(candidates),
        "candidate_blocks": len({row[1] for row in candidates}),
        "observed_fits": len(observed_rows),
        "null_fits": len(null_rows),
        "failures": len(failures),
        "methods": list(args.methods),
        "calibration": args.calibration,
        "null_replicates": args.null_replicates,
        "min_gene_umis": args.min_gene_umis,
        "min_gene_samples": args.min_gene_samples,
        "min_cell_types": args.min_cell_types,
        "min_celltype_mice": args.min_celltype_mice,
        "min_conditions": args.min_conditions,
        "min_condition_mice": args.min_condition_mice,
        "test_effect": args.test_effect,
        "subject_folds": (
            None if args.subject_folds is None else str(args.subject_folds)
        ),
        "subject_fold": args.subject_fold,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
