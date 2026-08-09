#!/usr/bin/env python3
"""Run the estimate-once local-path Wald ablation on an EC-block manifest."""

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
    partition_candidates,
    treatment_design,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm


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
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
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
            result = ec_block_glmm.path_wald_test(
                base,
                path_index,
                nuisance,
                tested,
                baseline=pooled_cache[cache_key],
                covariance=args.covariance,
                max_iter=args.fit_max_iter,
            )
            estimates = result.pop("path_estimates")
            observed_rows.append({
                "test_id": test_id,
                "block_id": block_id,
                "gene_id": gene_id,
                "contrast": test_effect,
                "cell_type": tested_cell_type,
                "method": "path_wald",
                "latent_space": "estimated_path",
                "covariance_method": args.covariance,
                "n_paths": len(signatures),
                "n_isoforms": base.n_isoforms,
                "n_original_isoforms": base.n_isoforms,
                "n_ecs": len(gene_ecs[gene]),
                "n_samples": len(local_metadata),
                "n_samples_used": int(estimates["usable"].sum()),
                "sample_use_rate": float(estimates["usable"].mean()),
                "n_test_levels": len(tested_levels),
                "degrees_of_freedom": result["degrees_of_freedom"],
                "denominator_degrees_of_freedom": result[
                    "denominator_degrees_of_freedom"
                ],
                "median_gene_umis": float(np.median(totals)),
                "statistic": result["statistic"],
                "lrt_p_value": result["p_value"],
                "p_value": result["p_value"],
                "bic_log_bayes_factor": np.nan,
                "null_converged": True,
                "alternative_converged": True,
                "n_clusters": result["clusters"],
                "path_fit_iterations": int(estimates["iterations"].sum()),
                "max_path_fit_iterations": int(estimates["iterations"].max()),
                "min_path_proportion": float(estimates["proportions"].min()),
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
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidate_tests": len(candidates),
        "observed_fits": len(observed_rows),
        "failures": len(failures),
        "method": "path_wald",
        "covariance_method": args.covariance,
        "seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
