#!/usr/bin/env python3
"""Merge cell-type compositional shards and recompute global FDR."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--min-samples",
        type=int,
        default=0,
        help="Minimum pseudobulk observations for the FDR hypothesis family.",
    )
    parser.add_argument(
        "--pvalue-method",
        choices=("pooled", "per-test"),
        default="pooled",
        help=(
            "Use a leave-block-out pooled permutation null, or retain the "
            "per-test permutation p-values already present in each shard."
        ),
    )
    parser.add_argument(
        "--calibration-bins",
        type=int,
        default=4,
        help=(
            "Quantile bins of median effective count used within each "
            "degrees-of-freedom stratum."
        ),
    )
    return parser.parse_args()


def empirical_null_calibration(table, null, calibration_bins=1):
    table = table.copy()
    null = null.copy()
    if "test_id" not in table:
        table["test_id"] = table["block_id"]
    if "test_id" not in null:
        null["test_id"] = null["block_id"]
    table["calibration_bin"] = 0
    null["_calibration_bin"] = 0
    if calibration_bins > 1 and "median_effective_count" in table:
        edges = np.unique(np.quantile(
            table["median_effective_count"],
            np.linspace(0.0, 1.0, int(calibration_bins) + 1),
        ))
        if len(edges) > 2:
            edges[0] = -np.inf
            edges[-1] = np.inf
            table["calibration_bin"] = pd.cut(
                table["median_effective_count"],
                edges,
                labels=False,
                include_lowest=True,
            ).astype(int)
            test_to_bin = table.set_index("test_id")[
                "calibration_bin"
            ]
            null["_calibration_bin"] = (
                null["test_id"].map(test_to_bin).astype(int)
            )
    raw_table_stratum = (
        table["degrees_of_freedom"].astype(str)
        + ":"
        + table["calibration_bin"].astype(str)
    )
    counts = raw_table_stratum.value_counts()
    table["_stratum"] = [
        value if counts[value] >= 20 else f"rare:{bin_index}"
        for value, bin_index in zip(
            raw_table_stratum, table["calibration_bin"]
        )
    ]
    raw_null_stratum = (
        null["degrees_of_freedom"].astype(str)
        + ":"
        + null["_calibration_bin"].astype(str)
    )
    null["_stratum"] = [
        value
        if counts.get(value, 0) >= 20
        else f"rare:{bin_index}"
        for value, bin_index in zip(
            raw_null_stratum, null["_calibration_bin"]
        )
    ]
    calibrated = np.empty(len(table), dtype=float)
    null_calibrated = np.empty(len(null), dtype=float)
    for stratum, positions in table.groupby("_stratum").groups.items():
        null_positions = null.index[null["_stratum"] == stratum]
        pool = np.sort(null.loc[null_positions, "p_value"].to_numpy())
        test_values = {
            test_id: np.sort(group["p_value"].to_numpy())
            for test_id, group in null.loc[null_positions].groupby(
                "test_id"
            )
        }
        for position in positions:
            test_id = table.at[position, "test_id"]
            own = test_values.get(test_id, np.empty(0))
            value = table.at[position, "asymptotic_p_value"]
            calibrated[position] = (
                1
                + np.searchsorted(pool, value, side="right")
                - np.searchsorted(own, value, side="right")
            ) / (1 + len(pool) - len(own))
        for position in null_positions:
            test_id = null.at[position, "test_id"]
            own = test_values[test_id]
            value = null.at[position, "p_value"]
            null_calibrated[position] = (
                1
                + np.searchsorted(pool, value, side="right")
                - np.searchsorted(own, value, side="right")
            ) / (1 + len(pool) - len(own))
    table["pooled_permutation_p_value"] = calibrated
    table["p_value"] = calibrated
    return table.drop(columns="_stratum"), null_calibrated


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = [
        pd.read_csv(
            shard / "differential_cell_type_compositional.tsv",
            sep="\t",
        )
        for shard in args.shards
    ]
    summaries = [
        json.loads((shard / "validation_summary.json").read_text())
        for shard in args.shards
    ]
    table = pd.concat(tables, ignore_index=True)
    identity = "test_id" if "test_id" in table else "block_id"
    if table[identity].duplicated().any():
        raise ValueError(f"shards contain duplicate {identity} values")
    table.insert(
        1,
        "gene_id",
        table["block_id"].str.rsplit(":B", n=1).str[0],
    )
    if args.pvalue_method == "pooled":
        null = pd.concat([
            pd.read_csv(shard / "permutation_null.tsv.gz", sep="\t")
            for shard in args.shards
        ], ignore_index=True)
        table, null_calibrated = empirical_null_calibration(
            table, null, args.calibration_bins
        )
    else:
        null_calibrated = np.empty(0)
    table["inference_eligible"] = (
        table["n_samples"] >= int(args.min_samples)
    )
    table["fdr"] = np.nan
    eligible = table["inference_eligible"] & table["converged"]
    table.loc[eligible, "fdr"] = benjamini_hochberg(
        table.loc[eligible, "p_value"].to_numpy()
    )
    table = table.sort_values("p_value")
    table.to_csv(
        args.output_dir / "differential_cell_type_compositional.tsv",
        sep="\t",
        index=False,
    )
    significant = table[table["fdr"] <= 0.05].copy()
    significant.to_csv(
        args.output_dir / "significant_cell_type_events.tsv",
        sep="\t",
        index=False,
    )
    permutation_tests = sum(
        summary["permutation_tests"] for summary in summaries
    )
    permutation_small = sum(
        summary["permutation_tests"]
        * summary["permutation_p_lt_0.05"]
        for summary in summaries
        if summary["permutation_p_lt_0.05"] is not None
    )
    summary = {
        "seconds_max_shard": max(
            summary["seconds"] for summary in summaries
        ),
        "min_path_proportion": summaries[0]["min_path_proportion"],
        "max_logratio_variance": summaries[0][
            "max_logratio_variance"
        ],
        "max_paths": summaries[0].get("max_paths"),
        "pvalue_method": args.pvalue_method,
        "calibration_bins": (
            args.calibration_bins
            if args.pvalue_method == "pooled"
            else None
        ),
        "candidate_partitions": max(
            summary["candidate_partitions"] for summary in summaries
        ),
        "tests": int(len(table)),
        "min_samples": args.min_samples,
        "tests_inference_eligible": int(eligible.sum()),
        "converged": int(table["converged"].sum()),
        "nominal_p_lt_0.05": int(
            ((table["p_value"] < 0.05) & eligible).sum()
        ),
        "fdr_0.05": int(len(significant)),
        "fdr_0.05_events": int(significant["block_id"].nunique()),
        "fdr_0.05_genes": int(significant["gene_id"].nunique()),
        "permutation_tests": permutation_tests,
        "asymptotic_permutation_p_lt_0.05": (
            permutation_small / permutation_tests
            if permutation_tests
            else None
        ),
        "pooled_permutation_p_lt_0.05": (
            float(np.mean(null_calibrated < 0.05))
            if len(null_calibrated)
            else None
        ),
    }
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
