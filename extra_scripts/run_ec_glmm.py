#!/usr/bin/env python3
"""Compare inference methods for paired-primer EC-count GLMMs genome-wide."""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import (
    aggregate_inputs,
    gene_structures,
    parse_group,
)
from tealeaf.sc import ec_glmm, glm_cv


METHODS = (
    "laplace_multinomial",
    "tilted_elbo",
    "renyi_multinomial",
    "laplace_dirichlet_multinomial",
    "elbo_dirichlet_multinomial",
    "renyi_dirichlet_multinomial",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    parser.add_argument("--barcode-groups", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-cache", type=Path)
    parser.add_argument("--contrast", choices=("condition", "cell_type"), default="condition")
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=METHODS)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-pseudobulk-umis", type=float, default=100_000)
    parser.add_argument("--min-gene-umis", type=float, default=100)
    parser.add_argument("--max-isoforms", type=int, default=6)
    parser.add_argument("--max-ecs", type=int, default=128)
    parser.add_argument("--max-genes", type=int)
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--vi-samples", type=int, default=128)
    parser.add_argument("--renyi-alpha", type=float, default=0.5)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def fixed_effect_design(groups, contrast):
    metadata = pd.DataFrame(
        [parse_group(group) for group in groups],
        columns=["cell_type", "condition", "mouse"],
    )
    nuisance_name = "cell_type" if contrast == "condition" else "condition"
    tested_name = contrast
    nuisance = pd.get_dummies(metadata[nuisance_name], dtype=float).to_numpy()
    tested = pd.get_dummies(metadata[tested_name], dtype=float, drop_first=True).to_numpy()
    if tested.shape[1] == 0:
        raise ValueError(f"{contrast} has fewer than two levels")
    null = nuisance
    alternative = np.column_stack((nuisance, tested))
    return null, alternative, metadata["mouse"].to_numpy(), metadata


def local_gene_data(counts, designs, transcripts, ecs, fixed, clusters):
    local_counts = []
    local_mappings = []
    for observed, mapping in zip(counts, designs):
        local_mapping = mapping[ecs][:, transcripts].tocsr()
        supported = np.asarray(local_mapping.sum(axis=1)).ravel() > 0
        local_counts.append(np.asarray(observed[:, ecs][:, supported].toarray(), dtype=float))
        local_mappings.append(np.asarray(local_mapping[supported].toarray(), dtype=float))
    totals = sum(value.sum(axis=1) for value in local_counts)
    retained = totals > 0
    if not retained.any():
        raise ValueError("gene has no positive observations")
    return ec_glmm.ECGLMMData(
        tuple(value[retained] for value in local_counts),
        tuple(local_mappings),
        fixed[retained],
        clusters[retained],
    ), retained, totals


def fit_one(method, data, args, initial=None):
    if method == "laplace_multinomial":
        return ec_glmm.fit_laplace(
            data, family="multinomial", initial=initial, max_iter=args.max_iter
        )
    if method == "tilted_elbo":
        return ec_glmm.fit_tilted_variational(
            data, initial=initial, max_iter=args.max_iter
        )
    if method == "renyi_multinomial":
        return ec_glmm.fit_variational(
            data,
            family="multinomial",
            alpha=args.renyi_alpha,
            samples=args.vi_samples,
            seed=args.seed,
            initial=initial,
            max_iter=args.max_iter,
        )
    if method == "laplace_dirichlet_multinomial":
        return ec_glmm.fit_laplace(
            data,
            family="dirichlet_multinomial",
            initial=initial,
            max_iter=args.max_iter,
        )
    if method in ("elbo_dirichlet_multinomial", "renyi_dirichlet_multinomial"):
        return ec_glmm.fit_variational(
            data,
            family="dirichlet_multinomial",
            alpha=1.0 if method.startswith("elbo") else args.renyi_alpha,
            samples=args.vi_samples,
            seed=args.seed,
            initial=initial,
            max_iter=args.max_iter,
        )
    raise ValueError(f"unknown method {method}")


def fit_nested(method, null_data, alternative_data, args, starts):
    null_initial = None
    alternative_initial = None
    laplace_method = {
        "tilted_elbo": "laplace_multinomial",
        "elbo_dirichlet_multinomial": "laplace_dirichlet_multinomial",
    }.get(method)
    if laplace_method in starts:
        null_initial = ec_glmm.variational_warm_start(
            starts[laplace_method][0], null_data.design.shape[1]
        )
        alternative_initial = ec_glmm.variational_warm_start(
            starts[laplace_method][1], alternative_data.design.shape[1]
        )
    if method == "renyi_multinomial" and "tilted_elbo" in starts:
        null_initial = ec_glmm.warm_start(
            starts["tilted_elbo"][0], null_data.design.shape[1]
        )
        alternative_initial = ec_glmm.warm_start(
            starts["tilted_elbo"][1], alternative_data.design.shape[1]
        )
    if (
        method == "renyi_dirichlet_multinomial"
        and "elbo_dirichlet_multinomial" in starts
    ):
        null_initial = ec_glmm.warm_start(
            starts["elbo_dirichlet_multinomial"][0],
            null_data.design.shape[1],
        )
        alternative_initial = ec_glmm.warm_start(
            starts["elbo_dirichlet_multinomial"][1],
            alternative_data.design.shape[1],
        )
    null = fit_one(method, null_data, args, initial=null_initial)
    if alternative_initial is None:
        alternative_initial = ec_glmm.warm_start(
            null, alternative_data.design.shape[1]
        )
    alternative = fit_one(method, alternative_data, args, initial=alternative_initial)
    return null, alternative


def prepare_dataset(args):
    """Prepare and aggregate paired EC inputs, optionally under a file lock."""
    def build():
        prepared = glm_cv.prepare_paired_primer_glm_data(
            args.alevin_dir,
            args.salmon_ref,
            args.primer_pairs,
            ec_design="weighted",
            regularization_target="theta",
            min_eq=5,
            min_half_umis=args.min_half_umis,
            primer_sampling_model="oligodt_tpm",
        )
        groups, _, _, counts = aggregate_inputs(args, prepared)
        structures = gene_structures(prepared, args.transcript_to_gene)
        return (groups, counts, *structures)

    if args.data_cache is None:
        return build()
    args.data_cache.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(f"{args.data_cache}.lock")
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.data_cache.is_file():
            with args.data_cache.open("rb") as handle:
                return pickle.load(handle)
        result = build()
        temporary = Path(f"{args.data_cache}.tmp.{os.getpid()}")
        with temporary.open("wb") as handle:
            pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)
        temporary.replace(args.data_cache)
        return result


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must be between zero and shard count")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (
        groups,
        counts,
        genes,
        gene_transcripts,
        gene_ecs,
        designs,
    ) = prepare_dataset(args)
    null_design, alternative_design, clusters, metadata = fixed_effect_design(
        groups, args.contrast
    )
    ec_totals = sum(
        np.asarray(observed.sum(axis=0)).ravel() for observed in counts
    )
    candidates = []
    for gene_index, (transcripts, ecs) in enumerate(zip(gene_transcripts, gene_ecs)):
        if len(transcripts) < 2 or not len(ecs):
            continue
        if len(ecs) > args.max_ecs:
            continue
        total = float(ec_totals[ecs].sum())
        if total < args.min_gene_umis * len(groups):
            continue
        support = sum(
            np.asarray(mapping[ecs][:, transcripts].sum(axis=0)).ravel()
            for mapping in designs
        )
        supported_transcripts = transcripts[support > 0]
        if not 2 <= len(supported_transcripts) <= args.max_isoforms:
            continue
        candidates.append((gene_index, supported_transcripts))
    eligible_genes_total = len(candidates)
    if args.max_genes is not None:
        candidates = candidates[: args.max_genes]
    candidates = candidates[args.shard_index :: args.shard_count]
    rows = []
    failures = []
    started = time.perf_counter()
    for position, (gene_index, transcripts) in enumerate(candidates):
        ecs = gene_ecs[gene_index]
        try:
            null_data, retained, totals = local_gene_data(
                counts, designs, transcripts, ecs, null_design, clusters
            )
            alternative_data, _, _ = local_gene_data(
                counts, designs, transcripts, ecs, alternative_design, clusters
            )
            starts = {}
            ordered_methods = sorted(args.methods, key=METHODS.index)
            for method in ordered_methods:
                method_started = time.perf_counter()
                null, alternative = fit_nested(
                    method, null_data, alternative_data, args, starts
                )
                starts[method] = (null, alternative)
                row = {
                    "gene": genes[gene_index],
                    "method": method,
                    "contrast": args.contrast,
                    "n_isoforms": len(transcripts),
                    "n_ecs": len(ecs),
                    "n_samples": int(retained.sum()),
                    "median_gene_umis": float(np.median(totals[retained])),
                    "evidence_gain": float(
                        2.0
                        * (
                            null["objective"] - alternative["objective"]
                            if method.startswith("laplace")
                            else alternative["objective"] - null["objective"]
                        )
                    ),
                    "null_objective": null["objective"],
                    "alternative_objective": alternative["objective"],
                    "null_converged": null["converged"],
                    "alternative_converged": alternative["converged"],
                    "null_random_effect_sd": null["random_effect_sd"],
                    "alternative_random_effect_sd": alternative["random_effect_sd"],
                    "null_concentration": null["concentration"],
                    "alternative_concentration": alternative["concentration"],
                    "tested_coefficient_norm": float(np.linalg.norm(
                        alternative["coefficients"][null_data.design.shape[1] :]
                    )),
                    "seconds": time.perf_counter() - method_started,
                }
                if "importance_ess" in null:
                    row["null_importance_ess"] = null["importance_ess"]
                    row["alternative_importance_ess"] = alternative[
                        "importance_ess"
                    ]
                if method == "tilted_elbo":
                    null_scores = ec_glmm.evaluate_variational_objectives(
                        null_data,
                        null,
                        alpha=args.renyi_alpha,
                        samples=args.vi_samples,
                        seed=args.seed + 1,
                    )
                    alternative_scores = (
                        ec_glmm.evaluate_variational_objectives(
                            alternative_data,
                            alternative,
                            alpha=args.renyi_alpha,
                            samples=args.vi_samples,
                            seed=args.seed + 1,
                        )
                    )
                    for score_name in null_scores:
                        row[f"null_{score_name}"] = null_scores[score_name]
                        row[f"alternative_{score_name}"] = (
                            alternative_scores[score_name]
                        )
                        row[f"{score_name}_gain"] = 2.0 * (
                            alternative_scores[score_name]
                            - null_scores[score_name]
                        )
                rows.append(row)
        except Exception as exc:
            failures.append({"gene": str(genes[gene_index]), "error": repr(exc)})
        if (position + 1) % 5 == 0:
            print(json.dumps({
                "genes_complete": position + 1,
                "genes_total": len(candidates),
                "fits": len(rows),
                "failures": len(failures),
                "seconds": time.perf_counter() - started,
            }), flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "ec_glmm_fits.tsv", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidate_genes": len(candidates),
        "eligible_genes_total": eligible_genes_total,
        "fits": len(rows),
        "failures": len(failures),
        "contrast": args.contrast,
        "methods": list(args.methods),
        "renyi_alpha": args.renyi_alpha,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
