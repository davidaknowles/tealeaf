#!/usr/bin/env python3
"""Benchmark tilted full-covariance L-BFGS schedules on selected real genes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts import run_ec_glmm
from tealeaf.sc import ec_glmm_full


VARIANTS = {
    "baseline25": {
        "tilted_local_steps": 25,
        "differentiate_tilted_local": True,
        "optimizer_maxcor": 10,
    },
    "envelope25": {
        "tilted_local_steps": 25,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
    },
    "envelope8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
    },
    "autodiff8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": True,
        "optimizer_maxcor": 10,
    },
    "envelope8_memory20": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 20,
    },
    "envelope8_ftol7": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "optimizer_ftol": 1e-7,
    },
    "standardized8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "standardize_latent": True,
    },
    "standardized_half8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "standardize_latent": True,
        "latent_standardization_power": 0.5,
    },
    "standardized_quarter8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "standardize_latent": True,
        "latent_standardization_power": 0.25,
    },
    "standardized_threequarter8": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "standardize_latent": True,
        "latent_standardization_power": 0.75,
    },
    "standardized_threequarter8_relaxed": {
        "schedule": "standardized_threequarter8_relaxed"
    },
    "robust_fp32": {"schedule": "robust_fp32"},
    "robust_bfloat16": {"schedule": "robust_bfloat16"},
    "fp32_25_then_fp64": {"schedule": "fp32_25_then_fp64"},
    "standardized8_ftol7": {
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "optimizer_maxcor": 10,
        "optimizer_ftol": 1e-7,
        "standardize_latent": True,
    },
    "centered20_standardized8": {"schedule": "centered20_standardized8"},
    "centered20_standardized8_relaxed": {
        "schedule": "centered20_standardized8_relaxed"
    },
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--genes", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--contrast", choices=("condition", "cell_type"), required=True)
    parser.add_argument("--variants", nargs="+", choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--max-isoforms", type=int, default=6)
    parser.add_argument("--max-ecs", type=int, default=128)
    parser.add_argument("--min-gene-umis", type=float, default=100)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    return parser.parse_args()


def candidate_lookup(dataset, min_gene_umis, max_isoforms, max_ecs):
    groups, counts, genes, gene_transcripts, gene_ecs, designs = dataset
    ec_totals = sum(np.asarray(observed.sum(axis=0)).ravel() for observed in counts)
    candidates = {}
    for gene_index, (transcripts, ecs) in enumerate(zip(gene_transcripts, gene_ecs)):
        if len(transcripts) < 2 or not len(ecs) or len(ecs) > max_ecs:
            continue
        if float(ec_totals[ecs].sum()) < min_gene_umis * len(groups):
            continue
        support = sum(
            np.asarray(mapping[ecs][:, transcripts].sum(axis=0)).ravel()
            for mapping in designs
        )
        supported = transcripts[support > 0]
        if 2 <= len(supported) <= max_isoforms:
            candidates[str(genes[gene_index])] = (gene_index, supported)
    return candidates


def fit_tilted(data, initial, max_iter, variant):
    started = time.perf_counter()
    if variant == "fp32_25_then_fp64":
        warm = ec_glmm_full.fit_variational(
            data,
            objective="tilted",
            observation_noise=True,
            initial=initial,
            max_iter=min(25, max_iter),
            tilted_local_steps=8,
            differentiate_tilted_local=False,
            standardize_latent=True,
            latent_standardization_power=0.75,
            compute_dtype="float32",
        )
        fit = ec_glmm_full.fit_tilted_variational_robust(
            data,
            observation_noise=True,
            initial=warm["parameters"],
            max_iter=max_iter,
            fallback_iter=50,
        )
        fit["iterations"] += warm["iterations"]
        fit["message"] = f"fp32 warm-up: {warm['message']}; fp64: {fit['message']}"
        return fit, time.perf_counter() - started
    if variant in (
        "standardized_threequarter8_relaxed",
        "robust_fp32",
        "robust_bfloat16",
    ):
        compute_dtype = {
            "standardized_threequarter8_relaxed": "float64",
            "robust_fp32": "float32",
            "robust_bfloat16": "bfloat16",
        }[variant]
        fit = ec_glmm_full.fit_tilted_variational_robust(
            data,
            observation_noise=True,
            initial=initial,
            max_iter=max_iter,
            fallback_iter=50,
            compute_dtype=compute_dtype,
        )
        return fit, time.perf_counter() - started
    if variant in (
        "centered20_standardized8",
        "centered20_standardized8_relaxed",
    ):
        first = ec_glmm_full.fit_variational(
            data,
            objective="tilted",
            observation_noise=True,
            initial=initial,
            max_iter=min(20, max_iter),
            tilted_local_steps=8,
            differentiate_tilted_local=False,
        )
        fit = ec_glmm_full.fit_variational(
            data,
            objective="tilted",
            observation_noise=True,
            initial=first["parameters"],
            max_iter=max(0, max_iter - first["iterations"]),
            tilted_local_steps=8,
            differentiate_tilted_local=False,
            standardize_latent=True,
        )
        fit["iterations"] += first["iterations"]
        fit["message"] = f"centered: {first['message']}; standardized: {fit['message']}"
        if variant.endswith("_relaxed") and not fit["converged"]:
            strict = fit
            fit = ec_glmm_full.fit_variational(
                data,
                objective="tilted",
                observation_noise=True,
                initial=strict["parameters"],
                max_iter=50,
                tilted_local_steps=8,
                differentiate_tilted_local=False,
                optimizer_ftol=1e-7,
                standardize_latent=True,
            )
            fit["iterations"] += strict["iterations"]
            fit["message"] = (
                f"{strict['message']}; relaxed continuation: {fit['message']}"
            )
        return fit, time.perf_counter() - started
    fit = ec_glmm_full.fit_variational(
        data,
        objective="tilted",
        observation_noise=True,
        initial=initial,
        max_iter=max_iter,
        **VARIANTS[variant],
    )
    return fit, time.perf_counter() - started


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    with args.data_cache.open("rb") as handle:
        dataset = pickle.load(handle)
    groups, counts, genes, gene_transcripts, gene_ecs, designs = dataset
    candidates = candidate_lookup(
        dataset, args.min_gene_umis, args.max_isoforms, args.max_ecs
    )
    requested = (
        sorted(candidates)
        if args.genes is None
        else [
            value.strip()
            for value in args.genes.read_text().splitlines()
            if value.strip()
        ]
    )
    requested = requested[args.shard_index :: args.shard_count]
    null_design, alternative_design, clusters, _ = run_ec_glmm.fixed_effect_design(
        groups, args.contrast
    )
    rows = []
    failures = []
    for gene in requested:
        try:
            gene_index, transcripts = candidates[gene]
            ecs = gene_ecs[gene_index]
            null_data, retained, totals = run_ec_glmm.local_gene_data(
                counts, designs, transcripts, ecs, null_design, clusters
            )
            alternative_data, _, _ = run_ec_glmm.local_gene_data(
                counts, designs, transcripts, ecs, alternative_design, clusters
            )
            laplace_started = time.perf_counter()
            null_laplace = run_ec_glmm.ec_glmm.fit_laplace(
                null_data,
                family="multinomial",
                observation_noise=True,
                max_iter=args.max_iter,
            )
            alternative_laplace = run_ec_glmm.ec_glmm.fit_laplace(
                alternative_data,
                family="multinomial",
                observation_noise=True,
                max_iter=args.max_iter,
            )
            laplace_seconds = time.perf_counter() - laplace_started
            null_initial = ec_glmm_full.variational_warm_start(
                null_laplace, null_data.design.shape[1]
            )
            alternative_initial = ec_glmm_full.variational_warm_start(
                alternative_laplace, alternative_data.design.shape[1]
            )
            for variant in args.variants:
                null, null_seconds = fit_tilted(
                    null_data, null_initial, args.max_iter, variant
                )
                alternative, alternative_seconds = fit_tilted(
                    alternative_data, alternative_initial, args.max_iter, variant
                )
                rows.append({
                    "gene": gene,
                    "contrast": args.contrast,
                    "variant": variant,
                    "n_isoforms": len(transcripts),
                    "n_ecs": len(ecs),
                    "n_samples": int(retained.sum()),
                    "median_gene_umis": float(np.median(totals[retained])),
                    "null_converged": null["converged"],
                    "alternative_converged": alternative["converged"],
                    "null_objective": null["objective"],
                    "alternative_objective": alternative["objective"],
                    "null_iterations": null["iterations"],
                    "alternative_iterations": alternative["iterations"],
                    "null_gradient_norm": null["gradient_norm"],
                    "alternative_gradient_norm": alternative["gradient_norm"],
                    "null_message": null["message"],
                    "alternative_message": alternative["message"],
                    "null_seconds": null_seconds,
                    "alternative_seconds": alternative_seconds,
                    "fit_seconds": null_seconds + alternative_seconds,
                    "laplace_seconds": laplace_seconds,
                })
        except Exception as error:
            failures.append({"gene": gene, "error": repr(error)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    args.output.with_suffix(".failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
