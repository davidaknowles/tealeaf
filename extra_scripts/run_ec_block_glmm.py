#!/usr/bin/env python3
"""Fit and bootstrap-calibrate splice-block EC-count GLMM cell-type tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_compositional_splicing import block_equivalence
from extra_scripts.run_differential_splicing import block_mapping, load_blocks
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full


METHODS = (
    "multinomial_full",
    "multinomial_noise_full",
    "dirichlet_multinomial_full",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--block-cache", required=True, type=Path)
    parser.add_argument("--event-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--null-replicates", type=int, default=20)
    parser.add_argument("--min-celltype-mice", type=int, default=3)
    parser.add_argument("--max-isoforms", type=int, default=10)
    parser.add_argument("--max-ecs", type=int, default=128)
    parser.add_argument("--max-iter", type=int, default=300)
    parser.add_argument("--retries", type=int, default=2)
    parser.add_argument("--null-replicate-retries", type=int)
    parser.add_argument("--vi-samples", type=int, default=16)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def celltype_design(groups, minimum_mice):
    metadata = pd.DataFrame(
        [str(group).rsplit("__", 2) for group in groups],
        columns=["cell_type", "condition", "mouse"],
    )
    represented = metadata.groupby("cell_type")["mouse"].nunique()
    cell_types = sorted(represented.index[represented >= int(minimum_mice)])
    retained = metadata["cell_type"].isin(cell_types).to_numpy()
    metadata = metadata.loc[retained].reset_index(drop=True)
    nuisance = pd.get_dummies(metadata["condition"], dtype=float).to_numpy()
    cell_type_index = {value: index for index, value in enumerate(cell_types)}
    labels = np.asarray([cell_type_index[value] for value in metadata["cell_type"]])
    subjects = (
        metadata["condition"].astype(str)
        + "__"
        + metadata["mouse"].astype(str)
    ).to_numpy()
    return retained, metadata, nuisance, labels, subjects, cell_types


def treatment_design(labels, n_levels):
    labels = np.asarray(labels, dtype=int)
    result = np.zeros((len(labels), int(n_levels) - 1), dtype=float)
    nonreference = labels > 0
    result[np.flatnonzero(nonreference), labels[nonreference] - 1] = 1.0
    return result


def tensor_data(base, tensor):
    return ec_glmm.ECGLMMData(
        base.counts,
        base.compatibility,
        np.ones((len(base.clusters), 1), dtype=float),
        base.clusters,
        fixed_effect_tensor=tensor,
    )


def fit_method(method, data, args, initial=None):
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
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(handle)
    features = np.loadtxt(args.features, dtype=str)
    if len(features) != designs[0].shape[1]:
        raise ValueError(
            f"feature file has {len(features)} rows but EC designs have "
            f"{designs[0].shape[1]} transcript columns"
        )
    retained, metadata, nuisance, labels, _, cell_types = celltype_design(
        groups, args.min_celltype_mice
    )
    counts = tuple(matrix[retained] for matrix in counts)
    clusters = metadata["mouse"].to_numpy()
    blocks = load_blocks(args.block_cache, None)
    block_to_representative, _ = block_equivalence(args.block_cache)
    event_table = pd.read_csv(args.event_table, sep="\t")
    eligible_blocks = set(
        event_table.loc[event_table["inference_eligible"], "block_id"].astype(str)
    )
    gene_position = {str(gene): index for index, gene in enumerate(genes)}
    candidates = []
    for block in blocks:
        representative = block_to_representative[block.block_id]
        if block.block_id != representative or representative not in eligible_blocks:
            continue
        gene = gene_position.get(block.gene_id)
        if gene is None or not len(gene_ecs[gene]) or len(gene_ecs[gene]) > args.max_ecs:
            continue
        transcripts = np.asarray(gene_transcripts[gene], dtype=int)
        support = sum(
            np.asarray(mapping[gene_ecs[gene]][:, transcripts].sum(axis=0)).ravel()
            for mapping in designs
        )
        transcripts = transcripts[support > 0]
        if not 2 <= len(transcripts) <= args.max_isoforms:
            continue
        mapped = block_mapping(block, features[transcripts])
        if mapped is None:
            continue
        path_index, signatures = mapped
        candidates.append((block, gene, transcripts, path_index, signatures))
    candidates = candidates[args.shard_index :: args.shard_count]

    observed_rows = []
    null_rows = []
    failures = []
    null_cache = {}
    started = time.perf_counter()
    for position, (block, gene, transcripts, path_index, signatures) in enumerate(candidates):
        try:
            base, _, totals = local_gene_data(
                counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(metadata), 1)),
                clusters,
                drop_zero=False,
            )
            tested = treatment_design(labels, len(cell_types))
            null_tensor, alternative_tensor, degrees = (
                ec_block_glmm.block_fixed_effect_tensors(
                    nuisance, tested, path_index
                )
            )
            null_data = tensor_data(base, null_tensor)
            for method in args.methods:
                cache_key = (block.block_id, method)
                if cache_key not in null_cache:
                    null_cache[cache_key] = fit_with_retries(
                        method, null_data, args
                    )
                null = null_cache[cache_key]
                alternative_data = tensor_data(base, alternative_tensor)
                initial = ec_glmm_full.fixed_effect_warm_start(
                    null, alternative_tensor.shape[2]
                )
                alternative = fit_with_retries(
                    method, alternative_data, args, initial=initial
                )
                statistic = 2.0 * (alternative["objective"] - null["objective"])
                observed_rows.append({
                    "block_id": block.block_id,
                    "gene_id": block.gene_id,
                    "method": method,
                    "n_paths": len(signatures),
                    "n_isoforms": len(transcripts),
                    "n_ecs": len(gene_ecs[gene]),
                    "n_samples": len(metadata),
                    "degrees_of_freedom": degrees,
                    "median_gene_umis": float(np.median(totals)),
                    "statistic": float(statistic),
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
                    "null_attempts": null["attempts"],
                    "alternative_attempts": alternative["attempts"],
                })
                bootstrap_rng = np.random.default_rng(np.random.SeedSequence((
                    args.seed,
                    zlib.crc32(block.block_id.encode("utf-8")),
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
                        observation_noise=method == "multinomial_noise_full",
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
                        "block_id": block.block_id,
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
            failures.append({"block_id": block.block_id, "error": repr(exc)})
        if (position + 1) % 2 == 0:
            print(json.dumps({
                "blocks_complete": position + 1,
                "blocks_total": len(candidates),
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
        "candidate_blocks": len(candidates),
        "observed_fits": len(observed_rows),
        "null_fits": len(null_rows),
        "failures": len(failures),
        "methods": list(args.methods),
        "null_replicates": args.null_replicates,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
