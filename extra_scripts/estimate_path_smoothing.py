#!/usr/bin/env python3
"""Estimate cross-fitted splice-path smoothing by Laplace empirical Bayes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import zlib

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special

from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_paired_path_test import filtered_inputs
from tealeaf.sc import differential, ec_block_glmm


def parse_grid(value):
    grid = np.asarray([float(item) for item in value.split(";")], dtype=float)
    if not len(grid) or np.any(grid <= 0) or len(np.unique(grid)) != len(grid):
        raise argparse.ArgumentTypeError("alpha grid must be unique and positive")
    return grid


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha-grid", type=parse_grid, default=parse_grid("0.25;0.5;1;2;4;8;16;32;64;128"))
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    return parser.parse_args()


def main():
    args = parse_args()
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    candidates = cached["candidates"]
    if args.max_candidates is not None:
        candidates = sorted(
            candidates,
            key=lambda candidate: zlib.crc32(candidate[0].encode()),
        )[: args.max_candidates]
    candidates = partition_candidates(candidates, args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(
        args.data_cache, cached["settings"]
    )
    rows = []
    for candidate in candidates:
        test_id, _, _, gene, transcripts, path_index, signatures, selected_rows, _, tested_levels = candidate
        try:
            local_metadata, _, labels = local_test_design(
                metadata, selected_rows, tested_levels, cached["settings"]["test_effect"]
            )
            local_counts = tuple(matrix[selected_rows] for matrix in counts)
            clusters = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(local_metadata), 1)),
                clusters,
                drop_zero=False,
            )
            aggregates = []
            subjects = np.unique(clusters)
            ordered = sorted(
                subjects,
                key=lambda value: zlib.crc32(
                    f"{test_id}|{value}".encode()
                ),
            )
            subject_fold = {
                subject: index % 2 for index, subject in enumerate(ordered)
            }
            for held_out in (0, 1):
                training = np.asarray([
                    subject_fold[cluster] != held_out for cluster in clusters
                ])
                if not training.any() or training.all():
                    continue
                training_data = type(base)(
                    tuple(values[training] for values in base.counts),
                    base.compatibility,
                    base.design[training],
                    base.clusters[training],
                )
                baseline = ec_block_glmm.pooled_isoform_weights(training_data)
                for cluster in subjects:
                    if subject_fold[cluster] != held_out:
                        continue
                    for level in np.unique(labels):
                        positions = np.flatnonzero(
                            (clusters == cluster) & (labels == level)
                        )
                        if not len(positions):
                            continue
                        aggregate = tuple(
                            np.asarray(values[positions], dtype=float).sum(axis=0)
                            for values in base.counts
                        )
                        if sum(value.sum() for value in aggregate) > 0:
                            aggregates.append((aggregate, baseline))
            evidence_by_alpha = {float(alpha): [] for alpha in args.alpha_grid}
            for aggregate, baseline in aggregates:
                local_evidence = {}
                for alpha in args.alpha_grid:
                    fit = differential.fit_path_perturbation(
                        aggregate,
                        base.compatibility,
                        baseline,
                        path_index,
                        path_pseudocount=alpha,
                        path_prior_center="baseline",
                    )
                    if fit.converged:
                        value = differential.path_laplace_log_evidence(
                            fit,
                            alpha,
                            prior_center=differential.path_proportions(
                                baseline, path_index
                            ),
                        )
                        if np.isfinite(value):
                            local_evidence[float(alpha)] = value
                if len(local_evidence) == len(args.alpha_grid):
                    for alpha, value in local_evidence.items():
                        evidence_by_alpha[alpha].append(value)
            for alpha, evidence in evidence_by_alpha.items():
                if evidence:
                    rows.append({
                        "test_id": test_id,
                        "gene": gene,
                        "n_paths": len(signatures),
                        "alpha": alpha,
                        "mean_log_evidence": float(np.mean(evidence)),
                        "n_aggregates": len(evidence),
                    })
        except (ValueError, np.linalg.LinAlgError):
            continue
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)


def select_cross_fitted_alpha(
    tables,
    output,
    folds=5,
    minimum_genes=25,
    maximum_slab_alpha=8.0,
):
    """Select a finite EB slab alpha against a high-precision spike."""
    table = pd.concat([pd.read_csv(path, sep="\t") for path in tables], ignore_index=True)
    table["gene_fold"] = table.gene.map(lambda value: zlib.crc32(str(value).encode()) % folds)
    per_gene = table.groupby(["gene", "gene_fold", "n_paths", "alpha"], as_index=False).agg(
        mean_log_evidence=("mean_log_evidence", "mean")
    )
    records = []
    for held_out in range(folds):
        training = per_gene.loc[per_gene.gene_fold.ne(held_out)]
        for n_paths in sorted(per_gene.n_paths.unique()):
            local = training.loc[training.n_paths.eq(n_paths)]
            if local.gene.nunique() < minimum_genes:
                local = training
            local = local.groupby(
                ["gene", "alpha"], as_index=False
            ).mean_log_evidence.mean()
            matrix = local.pivot(
                index="gene", columns="alpha", values="mean_log_evidence"
            ).dropna()
            spike_alpha = float(matrix.columns.max())
            slab_alphas = [
                float(alpha)
                for alpha in matrix.columns
                if alpha < spike_alpha and alpha <= maximum_slab_alpha
            ]
            if not slab_alphas:
                raise ValueError("alpha grid has no finite slab candidates")
            spike = matrix[spike_alpha].to_numpy(dtype=float)
            candidates = []
            for alpha in slab_alphas:
                slab = matrix[alpha].to_numpy(dtype=float)

                def objective(logit_weight):
                    spike_weight = scipy.special.expit(logit_weight)
                    values = np.logaddexp(
                        np.log(spike_weight) + spike,
                        np.log1p(-spike_weight) + slab,
                    )
                    return -float(values.mean())

                optimized = scipy.optimize.minimize_scalar(
                    objective,
                    bounds=(-8.0, 8.0),
                    method="bounded",
                )
                candidates.append((
                    -float(optimized.fun),
                    alpha,
                    float(scipy.special.expit(optimized.x)),
                ))
            mixture_evidence, selected_alpha, spike_weight = max(candidates)
            records.append({
                "gene_fold": held_out,
                "n_paths": int(n_paths),
                "alpha": selected_alpha,
                "training_genes": int(matrix.shape[0]),
                "spike_alpha": spike_alpha,
                "spike_weight": spike_weight,
                "mean_mixture_log_evidence": mixture_evidence,
            })
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"folds": folds, "records": records}, indent=2) + "\n")


if __name__ == "__main__":
    main()
