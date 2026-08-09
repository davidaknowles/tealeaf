#!/usr/bin/env python3
"""Combine separable batched L-BFGS benchmark outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def summarize(path):
    table = pd.read_csv(path, sep="\t")
    metadata = json.loads(path.with_suffix(".json").read_text())
    batched_joint = table.batched_null_converged & table.batched_alternative_converged
    scalar_joint = table.scalar_null_converged & table.scalar_alternative_converged
    common = batched_joint & scalar_joint
    return {
        "dataset": path.parent.name,
        "benchmark": path.stem,
        "fit_pair_key": metadata["fit_pair_key"].strip("'"),
        "tests": len(table),
        "unique_alternatives": metadata.get("unique_alternatives", len(table)),
        "batched_seconds": metadata["batched_seconds"],
        "scalar_seconds": metadata["scalar_seconds"],
        "speedup": metadata["speedup"],
        "batched_joint_convergence": int(batched_joint.sum()),
        "scalar_joint_convergence": int(scalar_joint.sum()),
        "common_joint_convergence": int(common.sum()),
        "common_median_absolute_statistic_difference": float(
            table.loc[common, "absolute_statistic_difference"].median()
        ),
        "common_maximum_absolute_statistic_difference": float(
            table.loc[common, "absolute_statistic_difference"].max()
        ),
        "median_batched_null_iterations": float(table.batched_null_iterations.median()),
        "median_scalar_null_iterations": float(table.scalar_null_iterations.median()),
        "median_batched_alternative_iterations": float(
            table.batched_alternative_iterations.median()
        ),
        "median_scalar_alternative_iterations": float(
            table.scalar_alternative_iterations.median()
        ),
    }


def main():
    args = parse_args()
    result = pd.DataFrame([summarize(path) for path in args.benchmark])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
