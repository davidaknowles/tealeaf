#!/usr/bin/env python3
"""Summarize full closed-form EC derivative runs against scalar production."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        nargs=6,
        metavar=(
            "NAME",
            "ANALYTIC_SHARDS",
            "ANALYTIC_RESULTS",
            "PRODUCTION_RESULTS",
            "ANALYTIC_WORKER_HOURS",
            "PRODUCTION_HOURS",
        ),
        required=True,
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def discoveries(table):
    return set(table.loc[table.fdr < 0.05, "test_id"])


def summarize(
    name,
    analytic_shards,
    analytic_results,
    production_results,
    analytic_worker_hours,
    production_hours,
):
    analytic = pd.read_csv(analytic_results, sep="\t")
    production = pd.read_csv(production_results, sep="\t")
    comparison = analytic.merge(
        production,
        on="test_id",
        suffixes=("_analytic", "_production"),
        validate="one_to_one",
    )
    common = (
        comparison.null_converged_analytic
        & comparison.alternative_converged_analytic
        & comparison.null_converged_production
        & comparison.alternative_converged_production
    )
    difference = np.abs(
        comparison.loc[common, "statistic_analytic"]
        - comparison.loc[common, "statistic_production"]
    )
    analytic_calls = discoveries(analytic)
    production_calls = discoveries(production)
    shard_paths = sorted(glob.glob(analytic_shards))
    worker_seconds = sum(
        json.loads((Path(path) / "summary.json").read_text())["seconds"]
        for path in shard_paths
    )
    analytic_worker_hours = float(analytic_worker_hours)
    production_hours = float(production_hours)
    return {
        "dataset": name,
        "tests": len(comparison),
        "analytic_joint_convergence": int(
            (analytic.null_converged & analytic.alternative_converged).sum()
        ),
        "production_joint_convergence": int(
            (production.null_converged & production.alternative_converged).sum()
        ),
        "common_joint_convergence": int(common.sum()),
        "analytic_bh_discoveries": len(analytic_calls),
        "production_bh_discoveries": len(production_calls),
        "shared_bh_discoveries": len(analytic_calls & production_calls),
        "discovery_jaccard": len(analytic_calls & production_calls)
        / len(analytic_calls | production_calls),
        "median_absolute_lrt_difference": float(np.median(difference)),
        "maximum_absolute_lrt_difference": float(np.max(difference)),
        "analytic_fit_hours": worker_seconds / 3600.0,
        "analytic_worker_hours": analytic_worker_hours,
        "production_worker_hours": production_hours,
        "worker_speedup": production_hours / analytic_worker_hours,
    }


def main():
    args = parse_args()
    table = pd.DataFrame([summarize(*dataset) for dataset in args.dataset])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
