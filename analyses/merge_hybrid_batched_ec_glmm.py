#!/usr/bin/env python3
"""Merge hybrid batched EC-GLMM fits and compare with scalar results."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
import scipy.stats

from analyses.run_hybrid_batched_ec_glmm import alternative_fit_id
from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--null-shards", required=True)
    parser.add_argument("--alternative-shards", required=True)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--batchability", required=True, type=Path)
    parser.add_argument("--scalar-results", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_shards(pattern):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise ValueError(f"no shards matched {pattern}")
    tables = []
    summaries = []
    for directory in paths:
        directory = Path(directory)
        tables.append(pd.read_csv(directory / "fits.tsv", sep="\t"))
        summaries.append(json.loads((directory / "summary.json").read_text()))
    return pd.concat(tables, ignore_index=True), summaries


def main():
    args = parse_args()
    null, null_summaries = read_shards(args.null_shards)
    alternative, alternative_summaries = read_shards(args.alternative_shards)
    with args.candidate_cache.open("rb") as handle:
        candidates = pickle.load(handle)["candidates"]
    dimensions = pd.read_csv(args.batchability, sep="\t").set_index("test_id")
    manifest = pd.DataFrame(
        {
            "test_id": [candidate[0] for candidate in candidates],
            "alternative_fit_id": [
                alternative_fit_id(candidate) for candidate in candidates
            ],
        }
    )
    manifest["degrees_of_freedom"] = manifest.test_id.map(
        dimensions.tested_coefficients
    )
    alternative = alternative.rename(
        columns={
            column: f"alternative_{column}"
            for column in alternative.columns
            if column != "fit_id"
        }
    )
    null = null.rename(
        columns={
            column: f"null_{column}" for column in null.columns if column != "fit_id"
        }
    )
    result = manifest.merge(
        null,
        left_on="test_id",
        right_on="fit_id",
        validate="one_to_one",
    ).merge(
        alternative,
        left_on="alternative_fit_id",
        right_on="fit_id",
        validate="many_to_one",
        suffixes=("", "_alternative_key"),
    )
    result["converged"] = result.null_converged & result.alternative_converged
    result["statistic"] = np.maximum(
        0.0,
        2.0 * (result.null_objective - result.alternative_objective),
    )
    result["lrt_p_value"] = scipy.stats.chi2.sf(
        result.statistic, result.degrees_of_freedom
    )
    result.loc[~result.converged, "lrt_p_value"] = 1.0
    result["bh_q_value"] = benjamini_hochberg(result.lrt_p_value.to_numpy())
    scalar = pd.read_csv(args.scalar_results, sep="\t").drop_duplicates("test_id")
    comparison = result.merge(
        scalar[
            [
                "test_id",
                "null_converged",
                "alternative_converged",
                "statistic",
                "lrt_p_value",
            ]
        ],
        on="test_id",
        suffixes=("_hybrid", "_scalar"),
        validate="one_to_one",
    )
    scalar_p = comparison.lrt_p_value_scalar.where(
        comparison.null_converged_scalar & comparison.alternative_converged_scalar,
        1.0,
    ).fillna(1.0)
    comparison["scalar_bh_q_value"] = benjamini_hochberg(scalar_p.to_numpy())
    common = (
        comparison.converged
        & comparison.null_converged_scalar
        & comparison.alternative_converged_scalar
    )
    difference = np.abs(
        comparison.loc[common, "statistic_hybrid"]
        - comparison.loc[common, "statistic_scalar"]
    )
    hybrid_calls = set(comparison.loc[comparison.bh_q_value <= 0.05, "test_id"])
    scalar_calls = set(comparison.loc[comparison.scalar_bh_q_value <= 0.05, "test_id"])
    summary = {
        "tests": len(comparison),
        "hybrid_joint_convergence": int(comparison.converged.sum()),
        "scalar_joint_convergence": int(
            (
                comparison.null_converged_scalar
                & comparison.alternative_converged_scalar
            ).sum()
        ),
        "common_joint_convergence": int(common.sum()),
        "common_median_absolute_statistic_difference": float(difference.median()),
        "common_maximum_absolute_statistic_difference": float(difference.max()),
        "hybrid_bh_discoveries": len(hybrid_calls),
        "scalar_bh_discoveries": len(scalar_calls),
        "shared_bh_discoveries": len(hybrid_calls & scalar_calls),
        "hybrid_only_bh_discoveries": len(hybrid_calls - scalar_calls),
        "scalar_only_bh_discoveries": len(scalar_calls - hybrid_calls),
        "null_worker_hours": sum(item["seconds"] for item in null_summaries) / 3600.0,
        "alternative_worker_hours": sum(
            item["seconds"] for item in alternative_summaries
        )
        / 3600.0,
        "total_worker_hours": (
            sum(item["seconds"] for item in null_summaries)
            + sum(item["seconds"] for item in alternative_summaries)
        )
        / 3600.0,
        "null_routes": null.null_route.value_counts().to_dict(),
        "alternative_routes": alternative.alternative_route.value_counts().to_dict(),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "hybrid_results.tsv.gz", sep="\t", index=False)
    comparison.to_csv(
        args.output_dir / "scalar_comparison.tsv.gz", sep="\t", index=False
    )
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
