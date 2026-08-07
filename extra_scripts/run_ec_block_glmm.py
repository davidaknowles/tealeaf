#!/usr/bin/env python3
"""Fit and bootstrap-calibrate splice-block EC-count GLMM DS tests."""

from __future__ import annotations

import argparse
from collections import defaultdict
import heapq
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
    "laplace_dirichlet_multinomial",
    "laplace_multinomial_noise",
    "laplace_multinomial_slope",
    "laplace_multinomial_slope_collapsed",
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
    parser.add_argument("--pairwise-null-seed", type=int)
    parser.add_argument("--max-candidates", type=int)
    parser.add_argument("--joint-gene-test", action="store_true")
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


def permute_paired_labels(metadata, labels, tested_levels, seed):
    """Swap paired labels within subjects using a reproducible null draw."""
    labels = np.asarray(labels, dtype=int).copy()
    pair_hash = zlib.crc32("||".join(map(str, tested_levels)).encode("utf-8"))
    pairs = []
    swaps = []
    for subject, positions in metadata.groupby("mouse", sort=True).groups.items():
        positions = np.asarray(list(positions), dtype=int)
        if len(positions) != 2 or set(labels[positions]) != {0, 1}:
            raise ValueError("pairwise null permutation requires one row per level")
        subject_hash = zlib.crc32(str(subject).encode("utf-8"))
        rng = np.random.default_rng(
            np.random.SeedSequence((int(seed), pair_hash, subject_hash))
        )
        pairs.append(positions)
        swaps.append(bool(rng.integers(2)))
    if len(swaps) < 2:
        raise ValueError("pairwise null permutation requires two subjects")
    if all(swaps) or not any(swaps):
        rng = np.random.default_rng(
            np.random.SeedSequence((int(seed), pair_hash, 0x51A9))
        )
        position = int(rng.integers(len(swaps)))
        swaps[position] = not swaps[position]
    for positions, swap in zip(pairs, swaps):
        if swap:
            labels[positions] = 1 - labels[positions]
    return labels


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
        "joint_gene_test": args.joint_gene_test,
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


def partition_candidates(candidates, shard_count):
    """Balance shards without splitting candidates that share an alternative."""
    groups = defaultdict(list)
    for candidate in candidates:
        gene = candidate[3]
        rows = tuple(np.asarray(candidate[7], dtype=int))
        groups[(gene, rows)].append(candidate)
    shards = [[] for _ in range(int(shard_count))]
    heap = [(0, index) for index in range(int(shard_count))]
    heapq.heapify(heap)
    ordered = sorted(groups.values(), key=len, reverse=True)
    for group in ordered:
        load, index = heapq.heappop(heap)
        shards[index].extend(group)
        heapq.heappush(heap, (load + len(group), index))
    return shards


def joint_gene_candidates(candidates):
    """Collapse block candidates to one unrestricted test per gene and design."""
    result = []
    seen = set()
    for candidate in candidates:
        (
            _, _, gene_id, gene, transcripts, _, _, rows,
            tested_cell_type, tested_levels,
        ) = candidate
        key = (
            gene,
            tuple(np.asarray(transcripts, dtype=int)),
            tuple(np.asarray(rows, dtype=int)),
            tested_cell_type,
            tuple(tested_levels),
        )
        if key in seen:
            continue
        seen.add(key)
        suffix = "|".join(map(str, tested_levels))
        test_id = f"{gene_id}|joint|{suffix}"
        result.append((
            test_id,
            f"{gene_id}:JOINT",
            gene_id,
            gene,
            transcripts,
            np.arange(len(transcripts), dtype=int),
            tuple((index,) for index in range(len(transcripts))),
            rows,
            tested_cell_type,
            tested_levels,
        ))
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


def tensor_data(base, tensor, random_effect_design=None):
    return ec_glmm.ECGLMMData(
        base.counts,
        base.compatibility,
        np.ones((len(base.clusters), 1), dtype=float),
        base.clusters,
        fixed_effect_tensor=tensor,
        random_effect_design=random_effect_design,
    )


def fit_method(method, data, args, initial=None):
    if method.startswith("laplace_"):
        return ec_glmm.fit_laplace(
            data,
            family=(
                "dirichlet_multinomial"
                if method == "laplace_dirichlet_multinomial"
                else "multinomial"
            ),
            observation_noise=method == "laplace_multinomial_noise",
            random_slopes=method in {
                "laplace_multinomial_slope",
                "laplace_multinomial_slope_collapsed",
            },
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


def reparameterize_fixed_effects(fit, source_tensor, target_tensor):
    """Map fitted logits into an equivalent fixed-effect tensor basis."""
    source_count = source_tensor.shape[2]
    coefficients = np.asarray(fit["parameters"][:source_count], dtype=float)
    fitted_logits = source_tensor.reshape(-1, source_count) @ coefficients
    target_coefficients = np.linalg.lstsq(
        target_tensor.reshape(-1, target_tensor.shape[2]),
        fitted_logits,
        rcond=None,
    )[0]
    return np.r_[target_coefficients, fit["parameters"][source_count:]]


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
        cached_settings = dict(cached.get("settings", {}))
        cached_settings.setdefault("joint_gene_test", False)
        if cached_settings != cache_settings:
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
        if args.joint_gene_test:
            candidates = joint_gene_candidates(candidates)
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
    if args.max_candidates is not None:
        if args.max_candidates < 1:
            raise ValueError("--max-candidates must be positive")
        rng = np.random.default_rng(args.seed)
        selected = rng.choice(
            len(candidates),
            size=min(int(args.max_candidates), len(candidates)),
            replace=False,
        )
        candidates = [candidates[index] for index in sorted(selected)]
    candidates = partition_candidates(candidates, args.shard_count)[
        args.shard_index
    ]
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
    alternative_cache = {}
    pooled_weight_cache = {}
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
            if args.pairwise_null_seed is not None:
                if args.test_effect != "cell_type_pairwise":
                    raise ValueError(
                        "--pairwise-null-seed requires cell_type_pairwise"
                    )
                labels = permute_paired_labels(
                    local_metadata,
                    labels,
                    tested_levels,
                    args.pairwise_null_seed,
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
            random_effect_design = None
            if args.test_effect == "cell_type_pairwise":
                centered_labels = np.asarray(labels, dtype=float)
                centered_labels -= centered_labels.mean()
                random_effect_design = np.column_stack(
                    (np.ones(len(labels), dtype=float), centered_labels)
                )
            null_data = tensor_data(
                base, null_tensor, random_effect_design=random_effect_design
            )
            full_alternative_tensor = ec_block_glmm.full_fixed_effect_tensor(
                nuisance, tested, base.n_isoforms - 1
            )
            full_alternative_data = tensor_data(
                base,
                full_alternative_tensor,
                random_effect_design=random_effect_design,
            )
            collapsed_spec = None
            if "laplace_multinomial_slope_collapsed" in args.methods:
                weight_key = (gene, tuple(rows))
                if weight_key not in pooled_weight_cache:
                    pooled_weight_cache[weight_key] = (
                        ec_block_glmm.pooled_isoform_weights(
                            base.counts, base.compatibility
                        )
                    )
                collapsed_mappings, collapsed_paths, _ = (
                    ec_block_glmm.collapse_within_paths(
                        base.compatibility,
                        pooled_weight_cache[weight_key],
                        path_index,
                    )
                )
                collapsed_base = ec_glmm.ECGLMMData(
                    base.counts,
                    collapsed_mappings,
                    base.design,
                    base.clusters,
                )
                collapsed_null_tensor, collapsed_alternative_tensor, _ = (
                    ec_block_glmm.block_fixed_effect_tensors(
                        nuisance, tested, collapsed_paths
                    )
                )
                collapsed_full_tensor = (
                    ec_block_glmm.full_fixed_effect_tensor(
                        nuisance, tested, collapsed_base.n_isoforms - 1
                    )
                )
                collapsed_spec = (
                    tensor_data(
                        collapsed_base,
                        collapsed_null_tensor,
                        random_effect_design=random_effect_design,
                    ),
                    tensor_data(
                        collapsed_base,
                        collapsed_full_tensor,
                        random_effect_design=random_effect_design,
                    ),
                    collapsed_null_tensor,
                    collapsed_alternative_tensor,
                    collapsed_full_tensor,
                    collapsed_base,
                )
            completed_methods = {}
            for method in args.methods:
                if method == "laplace_multinomial_slope_collapsed":
                    (
                        method_null_data,
                        method_alternative_data,
                        method_null_tensor,
                        method_block_alternative_tensor,
                        method_full_alternative_tensor,
                        method_base,
                    ) = collapsed_spec
                else:
                    method_null_data = null_data
                    method_alternative_data = full_alternative_data
                    method_null_tensor = null_tensor
                    method_block_alternative_tensor = alternative_tensor
                    method_full_alternative_tensor = full_alternative_tensor
                    method_base = base
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
                        method, method_null_data, args, initial=null_initial
                    )
                null = null_cache[cache_key]
                alternative_key = (
                    (gene, tuple(rows), method, block_id)
                    if method == "laplace_multinomial_slope_collapsed"
                    else (gene, tuple(rows), method)
                )
                alternative_reused = alternative_key in alternative_cache
                if not alternative_reused:
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
                        initial = reparameterize_fixed_effects(
                            null,
                            method_null_tensor,
                            method_full_alternative_tensor,
                        )
                    alternative_cache[alternative_key] = fit_with_retries(
                        method,
                        method_alternative_data,
                        args,
                        initial=initial,
                    )
                alternative = alternative_cache[alternative_key]
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
                    "pairwise_null_seed": args.pairwise_null_seed,
                    "n_paths": len(signatures),
                    "n_isoforms": method_base.n_isoforms,
                    "n_original_isoforms": len(transcripts),
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
                    "null_random_slope_sd": null.get("random_slope_sd", 0.0),
                    "alternative_random_slope_sd": alternative.get(
                        "random_slope_sd", 0.0
                    ),
                    "null_random_effect_sd": null["random_effect_sd"],
                    "alternative_random_effect_sd": alternative[
                        "random_effect_sd"
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
                    "alternative_reused": alternative_reused,
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
                        method_null_data,
                        null,
                        bootstrap_rng,
                        family=(
                            "dirichlet_multinomial"
                            if method in {
                                "dirichlet_multinomial_full",
                                "laplace_dirichlet_multinomial",
                            }
                            else "multinomial"
                        ),
                        observation_noise=method in {
                            "multinomial_noise_full",
                            "laplace_multinomial_noise",
                        },
                        random_slopes=method in {
                            "laplace_multinomial_slope",
                            "laplace_multinomial_slope_collapsed",
                        },
                    )
                    simulated_null_data = ec_glmm.ECGLMMData(
                        simulated_counts,
                        method_base.compatibility,
                        np.ones((len(method_base.clusters), 1), dtype=float),
                        method_base.clusters,
                        fixed_effect_tensor=method_null_tensor,
                        random_effect_design=random_effect_design,
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
                        method_base.compatibility,
                        np.ones((len(method_base.clusters), 1), dtype=float),
                        method_base.clusters,
                        fixed_effect_tensor=method_block_alternative_tensor,
                        random_effect_design=random_effect_design,
                    )
                    simulated_initial = ec_glmm_full.fixed_effect_warm_start(
                        simulated_null, method_block_alternative_tensor.shape[2]
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
                        "statistic": float(2.0 * (
                            simulated_null["objective"]
                            - simulated_fit["objective"]
                            if method.startswith("laplace_")
                            else simulated_fit["objective"]
                            - simulated_null["objective"]
                        )),
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
        "pairwise_null_seed": args.pairwise_null_seed,
        "max_candidates": args.max_candidates,
        "joint_gene_test": args.joint_gene_test,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
