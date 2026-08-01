#!/usr/bin/env python3
"""Merge and permutation-calibrate splice-block EC GLMM shards."""

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
    parser.add_argument("--depth-bins", type=int, default=4)
    parser.add_argument("--min-stratum-events", type=int, default=20)
    parser.add_argument("--calibration-lower", type=float, default=0.025)
    parser.add_argument("--calibration-upper", type=float, default=0.075)
    return parser.parse_args()


def add_depth_bins(table, null, bins):
    table = table.copy()
    null = null.copy()
    table["depth_bin"] = 0
    for method, positions in table.groupby("method").groups.items():
        values = table.loc[positions, "median_gene_umis"]
        edges = np.unique(np.quantile(values, np.linspace(0, 1, bins + 1)))
        if len(edges) <= 2:
            continue
        edges[0], edges[-1] = -np.inf, np.inf
        table.loc[positions, "depth_bin"] = pd.cut(
            values, edges, labels=False, include_lowest=True
        ).astype(int)
    lookup = table.set_index(["block_id", "method"])["depth_bin"]
    null["depth_bin"] = [
        lookup.loc[(block, method)]
        for block, method in zip(null["block_id"], null["method"])
    ]
    return table, null


def add_calibration_strata(table, minimum_events):
    """Pool rare df groups without mixing methods or depth bins."""
    table = table.copy()
    counts = table.groupby(
        ["method", "degrees_of_freedom", "depth_bin"]
    )["block_id"].transform("nunique")
    table["df_stratum"] = table["degrees_of_freedom"].astype(str)
    table.loc[counts < int(minimum_events), "df_stratum"] = "rare"
    return table


def calibrate(table, null):
    table = table.copy()
    null = null[null["converged"]].copy()
    stratum_lookup = table.set_index(["block_id", "method"])["df_stratum"]
    null["df_stratum"] = [
        stratum_lookup.loc[(block, method)]
        for block, method in zip(null["block_id"], null["method"])
    ]
    keys = ["method", "df_stratum", "depth_bin"]
    null_groups = {key: group for key, group in null.groupby(keys)}
    p_values = []
    calibrated_null_rows = []
    for key, records in table.groupby(keys):
        pool = null_groups.get(key, pd.DataFrame())
        block_null = {
            block: values["statistic"].to_numpy()
            for block, values in pool.groupby("block_id")
        }
        all_values = pool["statistic"].to_numpy()
        for row in records.itertuples():
            own = block_null.get(row.block_id, np.empty(0))
            exceed = np.sum(all_values >= row.statistic) - np.sum(
                own >= row.statistic
            )
            denominator = len(all_values) - len(own)
            p_values.append((row.Index, (1 + exceed) / (1 + denominator)))
        for row in pool.itertuples():
            comparison = np.delete(
                all_values,
                np.flatnonzero(pool["block_id"].to_numpy() == row.block_id),
            )
            calibrated_null_rows.append({
                "method": row.method,
                "block_id": row.block_id,
                "permutation": row.permutation,
                "p_value": (
                    (1 + np.sum(comparison >= row.statistic))
                    / (1 + len(comparison))
                ),
            })
    table["p_value"] = np.nan
    for position, value in p_values:
        table.at[position, "p_value"] = value
    return table, pd.DataFrame(calibrated_null_rows)


def main():
    args = parse_args()
    tables = []
    nulls = []
    failures = []
    for shard in args.shards:
        table_path = shard / "ec_block_glmm.tsv"
        null_path = shard / "permutation_null.tsv.gz"
        if table_path.is_file() and table_path.stat().st_size:
            tables.append(pd.read_csv(table_path, sep="\t"))
        if null_path.is_file() and null_path.stat().st_size:
            nulls.append(pd.read_csv(null_path, sep="\t"))
        failure_path = shard / "failures.json"
        if failure_path.is_file():
            failures.extend(json.loads(failure_path.read_text()))
    table = pd.concat(tables, ignore_index=True).drop_duplicates(
        ["block_id", "method"], keep="last"
    )
    null = pd.concat(nulls, ignore_index=True)
    table, null = add_depth_bins(table, null, args.depth_bins)
    table = add_calibration_strata(table, args.min_stratum_events)
    table, calibrated_null = calibrate(table, null)
    table["fdr"] = np.nan
    eligible = table["null_converged"] & table["alternative_converged"]
    for method, positions in table.loc[eligible].groupby("method").groups.items():
        rejection = np.mean(
            calibrated_null.loc[
                calibrated_null["method"] == method, "p_value"
            ] <= 0.05
        )
        if args.calibration_lower <= rejection <= args.calibration_upper:
            table.loc[positions, "fdr"] = benjamini_hochberg(
                table.loc[positions, "p_value"].to_numpy()
            )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False)
    table[table["fdr"] <= 0.05].to_csv(
        args.output_dir / "significant_ec_block_glmm.tsv", sep="\t", index=False
    )
    calibrated_null.to_csv(
        args.output_dir / "calibrated_permutation_null.tsv.gz",
        sep="\t",
        index=False,
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    summary = []
    for method, records in table.groupby("method"):
        method_null = calibrated_null.loc[
            calibrated_null["method"] == method, "p_value"
        ].to_numpy()
        rejection = float(np.mean(method_null <= 0.05))
        summary.append({
            "method": method,
            "events": len(records),
            "convergence": float(np.mean(
                records["null_converged"] & records["alternative_converged"]
            )),
            "nominal": int(np.sum(records["p_value"] <= 0.05)),
            "fdr_0.05": int(np.sum(records["fdr"] <= 0.05)),
            "null_replicates": len(method_null),
            "null_rejection_0.05": rejection,
            "calibration_acceptable": (
                args.calibration_lower <= rejection <= args.calibration_upper
            ),
        })
    pd.DataFrame(summary).to_csv(
        args.output_dir / "summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
