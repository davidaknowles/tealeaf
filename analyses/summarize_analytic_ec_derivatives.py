#!/usr/bin/env python3
"""Compare closed-form and autodiff EC multinomial derivative benchmarks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        nargs=3,
        metavar=("NAME", "ANALYTIC_DIR", "AUTODIFF_DIR"),
        required=True,
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load(directory):
    directory = Path(directory)
    table = pd.read_csv(directory / "ec_block_glmm.tsv", sep="\t")
    summary = json.loads((directory / "summary.json").read_text())
    return table, summary


def summarize(name, analytic_directory, autodiff_directory):
    analytic, analytic_summary = load(analytic_directory)
    autodiff, autodiff_summary = load(autodiff_directory)
    comparison = analytic.merge(
        autodiff,
        on="test_id",
        suffixes=("_analytic", "_autodiff"),
        validate="one_to_one",
    )
    common = (
        comparison.null_converged_analytic
        & comparison.alternative_converged_analytic
        & comparison.null_converged_autodiff
        & comparison.alternative_converged_autodiff
    )
    statistic_difference = np.abs(
        comparison.loc[common, "statistic_analytic"]
        - comparison.loc[common, "statistic_autodiff"]
    )
    objective_columns = [
        "null_objective_evaluation_seconds",
        "alternative_objective_evaluation_seconds",
    ]
    return {
        "dataset": name,
        "tests": len(comparison),
        "analytic_joint_convergence": int(
            (analytic.null_converged & analytic.alternative_converged).sum()
        ),
        "autodiff_joint_convergence": int(
            (autodiff.null_converged & autodiff.alternative_converged).sum()
        ),
        "common_joint_convergence": int(common.sum()),
        "median_absolute_lrt_difference": float(np.median(statistic_difference)),
        "maximum_absolute_lrt_difference": float(np.max(statistic_difference)),
        "analytic_elapsed_seconds": float(analytic_summary["seconds"]),
        "autodiff_elapsed_seconds": float(autodiff_summary["seconds"]),
        "elapsed_speedup": float(
            autodiff_summary["seconds"] / analytic_summary["seconds"]
        ),
        "analytic_objective_seconds": float(analytic[objective_columns].sum().sum()),
        "autodiff_objective_seconds": float(autodiff[objective_columns].sum().sum()),
        "objective_speedup": float(
            autodiff[objective_columns].sum().sum()
            / analytic[objective_columns].sum().sum()
        ),
    }


def main():
    args = parse_args()
    rows = [summarize(*dataset) for dataset in args.dataset]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
