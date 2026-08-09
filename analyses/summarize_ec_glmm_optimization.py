#!/usr/bin/env python3
"""Summarize matched EC-GLMM optimizer and concurrency benchmarks."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

STRATEGY_TAGS = {
    "lbfgs": "baseline",
    "shared_alternative_lbfgs": "shared",
    "previous_null_lbfgs": "previous",
    "adam_lr_0.01": "adam01",
    "adam_lr_0.03": "adam03",
    "adam_then_lbfgs": "hybrid",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mouse-root", required=True, type=Path)
    parser.add_argument("--human-root", required=True, type=Path)
    parser.add_argument("--accounting", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_results(root, tag):
    path = root / f"ec_glmm_opt_{tag}_shard0" / "ec_block_glmm.tsv"
    return pd.read_csv(path, sep="\t")


def summarize_dataset(dataset, root, accounting):
    baseline = read_results(root, "baseline")
    records = []
    for strategy, tag in STRATEGY_TAGS.items():
        result = read_results(root, tag)
        matched = baseline.merge(
            result, on=["test_id", "method"], suffixes=("_baseline", "_strategy")
        )
        difference = np.abs(matched.statistic_strategy - matched.statistic_baseline)
        converged = result.null_converged & result.alternative_converged
        row = accounting.loc[
            (accounting.dataset == dataset) & (accounting.strategy == strategy)
        ].iloc[0]
        records.append(
            {
                "dataset": dataset,
                "strategy": strategy,
                "candidate_tests": len(result),
                "matched_baseline_tests": len(matched),
                "jointly_converged_tests": int(converged.sum()),
                "job_seconds": int(row.job_seconds),
                "wall_ratio_to_lbfgs": np.nan,
                "total_cpu_seconds": float(row.total_cpu_seconds),
                "peak_rss_gib": float(row.peak_rss_kib) / 2**20,
                "median_absolute_statistic_difference": float(difference.median()),
                "maximum_absolute_statistic_difference": float(difference.max()),
            }
        )
    result = pd.DataFrame(records)
    baseline_seconds = result.loc[result.strategy.eq("lbfgs"), "job_seconds"].iloc[0]
    result["wall_ratio_to_lbfgs"] = result.job_seconds / baseline_seconds
    parallel = accounting.loc[
        (accounting.dataset == dataset) & accounting.strategy.eq("two_concurrent_lbfgs")
    ].iloc[0]
    parallel_result = pd.concat(
        [
            pd.read_csv(
                root / f"ec_glmm_opt_parallel_shard{index}" / "ec_block_glmm.tsv",
                sep="\t",
            )
            for index in (0, 1)
        ],
        ignore_index=True,
    )
    parallel_matched = baseline.merge(
        parallel_result, on=["test_id", "method"], suffixes=("_baseline", "_strategy")
    )
    parallel_difference = np.abs(
        parallel_matched.statistic_strategy - parallel_matched.statistic_baseline
    )
    records = pd.concat(
        [
            result,
            pd.DataFrame(
                [
                    {
                        "dataset": dataset,
                        "strategy": "two_concurrent_lbfgs",
                        "candidate_tests": int(parallel.candidate_tests),
                        "matched_baseline_tests": int(parallel.candidate_tests),
                        "jointly_converged_tests": int(
                            (
                                parallel_result.null_converged
                                & parallel_result.alternative_converged
                            ).sum()
                        ),
                        "job_seconds": int(parallel.job_seconds),
                        "wall_ratio_to_lbfgs": float(parallel.job_seconds)
                        / baseline_seconds,
                        "total_cpu_seconds": float(parallel.total_cpu_seconds),
                        "peak_rss_gib": float(parallel.peak_rss_kib) / 2**20,
                        "median_absolute_statistic_difference": float(
                            parallel_difference.median()
                        ),
                        "maximum_absolute_statistic_difference": float(
                            parallel_difference.max()
                        ),
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    return records


def main():
    args = parse_args()
    accounting = pd.read_csv(args.accounting, sep="\t")
    summary = pd.concat(
        [
            summarize_dataset("mouse", args.mouse_root, accounting),
            summarize_dataset("human", args.human_root, accounting),
        ],
        ignore_index=True,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
