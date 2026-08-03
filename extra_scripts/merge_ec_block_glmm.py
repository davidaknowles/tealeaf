#!/usr/bin/env python3
"""Merge and bootstrap-calibrate splice-block EC GLMM shards."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
import zlib

import numpy as np
import pandas as pd
from scipy.stats import genpareto

from extra_scripts.run_differential_splicing import benjamini_hochberg
from extra_scripts.run_ec_block_glmm import supported_partition_key


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--candidate-cache",
        type=Path,
        help="Collapse blocks equivalent on the EC-supported isoform subset.",
    )
    parser.add_argument("--depth-bins", type=int, default=1)
    parser.add_argument("--audit-depth-bins", type=int, default=4)
    parser.add_argument("--audit-sample-bins", type=int, default=4)
    parser.add_argument(
        "--calibration-mode",
        choices=("empirical", "gpd_tail"),
        default="empirical",
    )
    parser.add_argument("--tail-quantile", type=float, default=0.9)
    parser.add_argument("--tail-folds", type=int, default=5)
    parser.add_argument("--min-stratum-events", type=int, default=20)
    parser.add_argument("--calibration-lower", type=float, default=0.025)
    parser.add_argument("--calibration-upper", type=float, default=0.075)
    return parser.parse_args()


def supported_partition_representatives(candidate_cache):
    with candidate_cache.open("rb") as handle:
        candidates = pickle.load(handle)["candidates"]
    representatives = {}
    first = {}
    for candidate in candidates:
        key = supported_partition_key(candidate)
        representative = first.setdefault(key, candidate[0])
        representatives[candidate[0]] = representative
    return representatives


def quantile_bin_values(table, column, bins):
    result = pd.Series(0, index=table.index, dtype=int)
    for method, positions in table.groupby("method").groups.items():
        values = table.loc[positions, column]
        edges = np.unique(np.quantile(values, np.linspace(0, 1, bins + 1)))
        if len(edges) <= 2:
            continue
        edges[0], edges[-1] = -np.inf, np.inf
        result.loc[positions] = pd.cut(
            values, edges, labels=False, include_lowest=True
        ).astype(int)
    return result


def depth_bin_values(table, bins):
    return quantile_bin_values(table, "median_gene_umis", bins)


def add_depth_bins(table, null, bins):
    table = table.copy()
    null = null.copy()
    table["depth_bin"] = depth_bin_values(table, bins)
    lookup = table.set_index(["test_id", "method"])["depth_bin"]
    null["depth_bin"] = [
        lookup.loc[(test_id, method)]
        for test_id, method in zip(null["test_id"], null["method"])
    ]
    return table, null


def add_calibration_strata(table, minimum_events):
    """Pool rare df groups across depth without mixing methods."""
    table = table.copy()
    counts = table.groupby(
        ["method", "degrees_of_freedom", "depth_bin"]
    )["test_id"].transform("nunique")
    table["df_stratum"] = table["degrees_of_freedom"].astype(str)
    table.loc[counts < int(minimum_events), "df_stratum"] = "rare"
    table.loc[table["df_stratum"] == "rare", "depth_bin"] = -1
    return table


def calibrate(table, null):
    table = table.copy()
    null = null[null["converged"]].copy()
    stratum_lookup = table.set_index(["test_id", "method"])["df_stratum"]
    null["df_stratum"] = [
        stratum_lookup.loc[(test_id, method)]
        for test_id, method in zip(null["test_id"], null["method"])
    ]
    null.loc[null["df_stratum"] == "rare", "depth_bin"] = -1
    keys = ["method", "df_stratum", "depth_bin"]
    null_groups = {key: group for key, group in null.groupby(keys)}
    p_values = []
    calibrated_null_rows = []
    for key, records in table.groupby(keys):
        pool = null_groups.get(key)
        if pool is None or pool["test_id"].nunique() < 2:
            p_values.extend((position, np.nan) for position in records.index)
            continue
        block_null = {
            test_id: np.sort(values["statistic"].to_numpy())
            for test_id, values in pool.groupby("test_id")
        }
        all_values = np.sort(pool["statistic"].to_numpy())
        for row in records.itertuples():
            own = block_null.get(row.test_id, np.empty(0))
            exceed = (
                len(all_values)
                - np.searchsorted(all_values, row.statistic, side="left")
                - len(own)
                + np.searchsorted(own, row.statistic, side="left")
            )
            denominator = len(all_values) - len(own)
            p_values.append((row.Index, (1 + exceed) / (1 + denominator)))
        calibrated = pool[["method", "test_id", "block_id"]].copy()
        calibrated["replicate"] = (
            pool["replicate"] if "replicate" in pool else pool["permutation"]
        )
        calibrated["p_value"] = np.nan
        for test_id, positions in pool.groupby("test_id").groups.items():
            own = block_null[test_id]
            statistics = pool.loc[positions, "statistic"].to_numpy()
            exceed = (
                len(all_values)
                - np.searchsorted(all_values, statistics, side="left")
                - len(own)
                + np.searchsorted(own, statistics, side="left")
            )
            calibrated.loc[positions, "p_value"] = (
                (1 + exceed) / (1 + len(all_values) - len(own))
            )
        calibrated_null_rows.extend(calibrated.to_dict("records"))
    table["p_value"] = np.nan
    for position, value in p_values:
        table.at[position, "p_value"] = value
    return table, pd.DataFrame(calibrated_null_rows)


def fit_gpd_tail(values, quantile):
    """Fit a peaks-over-threshold tail while retaining the empirical body."""
    values = np.asarray(values, dtype=float)
    threshold = float(np.quantile(values, quantile))
    excess = values[values > threshold] - threshold
    if len(excess) < 20:
        raise ValueError("too few bootstrap statistics for GPD tail fitting")
    shape, _, scale = genpareto.fit(excess, floc=0.0)
    if not np.isfinite(shape) or not np.isfinite(scale) or scale <= 0:
        raise ValueError("invalid GPD tail fit")
    return {
        "threshold": threshold,
        "tail_probability": len(excess) / len(values),
        "shape": float(shape),
        "scale": float(scale),
        "sorted_values": np.sort(values),
    }


def gpd_tail_p_values(statistics, model):
    statistics = np.asarray(statistics, dtype=float)
    values = model["sorted_values"]
    p_values = (
        1 + len(values) - np.searchsorted(values, statistics, side="left")
    ) / (1 + len(values))
    in_tail = statistics > model["threshold"]
    p_values[in_tail] = model["tail_probability"] * genpareto.sf(
        statistics[in_tail] - model["threshold"],
        model["shape"],
        loc=0.0,
        scale=model["scale"],
    )
    return np.maximum(p_values, np.finfo(float).tiny)


def calibrate_gpd_tail(table, null, quantile, folds):
    """Cross-fit tail scores, then calibrate them globally by null rank."""
    table = table.copy()
    null = null[null["converged"]].copy()
    stratum_lookup = table.set_index(["test_id", "method"])["df_stratum"]
    null["df_stratum"] = [
        stratum_lookup.loc[(test_id, method)]
        for test_id, method in zip(null["test_id"], null["method"])
    ]
    null.loc[null["df_stratum"] == "rare", "depth_bin"] = -1
    keys = ["method", "df_stratum", "depth_bin"]
    table["tail_score"] = np.nan
    null["tail_score"] = np.nan
    null_groups = {key: group for key, group in null.groupby(keys)}
    for key, records in table.groupby(keys):
        pool = null_groups.get(key)
        if pool is None or pool["test_id"].nunique() < 2:
            continue
        null_fold = np.asarray([
            zlib.crc32(test_id.encode("utf-8")) % int(folds)
            for test_id in pool["test_id"]
        ])
        record_fold = np.asarray([
            zlib.crc32(test_id.encode("utf-8")) % int(folds)
            for test_id in records["test_id"]
        ])
        for held_out in range(int(folds)):
            test = null_fold == held_out
            cross_fit = fit_gpd_tail(
                pool.loc[~test, "statistic"], quantile
            )
            null.loc[pool.index[test], "tail_score"] = gpd_tail_p_values(
                pool.loc[test, "statistic"], cross_fit
            )
            record_test = record_fold == held_out
            table.loc[
                records.index[record_test], "tail_score"
            ] = gpd_tail_p_values(
                records.loc[record_test, "statistic"], cross_fit
            )

    table["p_value"] = np.nan
    calibrated_null_rows = []
    for method, records in table.groupby("method"):
        pool = null.loc[
            null["method"].eq(method) & null["tail_score"].notna()
        ]
        all_scores = np.sort(pool["tail_score"].to_numpy())
        test_scores = {
            test_id: np.sort(group["tail_score"].to_numpy())
            for test_id, group in pool.groupby("test_id")
        }
        for row in records.loc[records["tail_score"].notna()].itertuples():
            own = test_scores.get(row.test_id, np.empty(0))
            count = (
                np.searchsorted(all_scores, row.tail_score, side="right")
                - np.searchsorted(own, row.tail_score, side="right")
            )
            table.at[row.Index, "p_value"] = (
                (1 + count) / (1 + len(all_scores) - len(own))
            )
        pool_p_values = np.full(len(pool), np.nan)
        for test_id, positions in pool.groupby("test_id").groups.items():
            own = test_scores[test_id]
            scores = pool.loc[positions, "tail_score"].to_numpy()
            count = (
                np.searchsorted(all_scores, scores, side="right")
                - np.searchsorted(own, scores, side="right")
            )
            pool_p_values[pool.index.get_indexer(positions)] = (
                (1 + count) / (1 + len(all_scores) - len(own))
            )
        calibrated_null_rows.extend(pd.DataFrame({
            "method": pool["method"].to_numpy(),
            "test_id": pool["test_id"].to_numpy(),
            "block_id": pool["block_id"].to_numpy(),
            "replicate": (
                pool["replicate"].to_numpy()
                if "replicate" in pool
                else pool["permutation"].to_numpy()
            ),
            "p_value": pool_p_values,
        }).to_dict("records"))
    return table.drop(columns="tail_score"), pd.DataFrame(calibrated_null_rows)


def main():
    args = parse_args()
    tables = []
    nulls = []
    failures = []
    for shard in args.shards:
        table_path = shard / "ec_block_glmm.tsv"
        null_path = shard / "bootstrap_null.tsv.gz"
        if not null_path.is_file():
            null_path = shard / "permutation_null.tsv.gz"
        if table_path.is_file() and table_path.stat().st_size:
            shard_table = pd.read_csv(table_path, sep="\t")
            if "test_id" not in shard_table:
                shard_table["test_id"] = shard_table["block_id"]
            tables.append(shard_table)
        if null_path.is_file() and null_path.stat().st_size:
            shard_null = pd.read_csv(null_path, sep="\t")
            if "test_id" not in shard_null:
                shard_null["test_id"] = shard_null["block_id"]
            nulls.append(shard_null)
        failure_path = shard / "failures.json"
        if failure_path.is_file():
            failures.extend(json.loads(failure_path.read_text()))
    table = pd.concat(tables, ignore_index=True).drop_duplicates(
        ["test_id", "method"], keep="last"
    )
    null = pd.concat(nulls, ignore_index=True)
    if args.candidate_cache is not None:
        representatives = supported_partition_representatives(
            args.candidate_cache
        )
        keep = table["test_id"].map(representatives).eq(table["test_id"])
        retained_ids = set(table.loc[keep, "test_id"])
        table = table.loc[keep].reset_index(drop=True)
        null = null.loc[
            null["test_id"].isin(retained_ids)
        ].reset_index(drop=True)
    table, null = add_depth_bins(table, null, args.depth_bins)
    min_stratum_events = int(args.min_stratum_events)
    if args.calibration_mode == "gpd_tail":
        min_stratum_events = max(min_stratum_events, 100)
    table = add_calibration_strata(table, min_stratum_events)
    if args.calibration_mode == "gpd_tail":
        table, calibrated_null = calibrate_gpd_tail(
            table, null, args.tail_quantile, args.tail_folds
        )
    else:
        table, calibrated_null = calibrate(table, null)
    table["audit_depth_bin"] = depth_bin_values(
        table, args.audit_depth_bins
    )
    audit_lookup = table.set_index(["test_id", "method"])["audit_depth_bin"]
    calibrated_null["audit_depth_bin"] = [
        audit_lookup.loc[(test_id, method)]
        for test_id, method in zip(
            calibrated_null["test_id"], calibrated_null["method"]
        )
    ]
    table["audit_sample_bin"] = quantile_bin_values(
        table, "n_samples", args.audit_sample_bins
    )
    sample_lookup = table.set_index(["test_id", "method"])["audit_sample_bin"]
    calibrated_null["audit_sample_bin"] = [
        sample_lookup.loc[(test_id, method)]
        for test_id, method in zip(
            calibrated_null["test_id"], calibrated_null["method"]
        )
    ]
    table["fdr"] = np.nan
    eligible = (
        table["null_converged"]
        & table["alternative_converged"]
        & table["p_value"].notna()
    )
    for method, positions in table.loc[eligible].groupby("method").groups.items():
        rejection = np.mean(
            calibrated_null.loc[
                calibrated_null["method"] == method, "p_value"
            ] <= 0.05
        )
        depth_rejection = calibrated_null.loc[
            calibrated_null["method"] == method
        ].groupby("audit_depth_bin")["p_value"].apply(
            lambda values: np.mean(values <= 0.05)
        )
        sample_rejection = calibrated_null.loc[
            calibrated_null["method"] == method
        ].groupby("audit_sample_bin")["p_value"].apply(
            lambda values: np.mean(values <= 0.05)
        )
        calibrated = (
            args.calibration_lower <= rejection <= args.calibration_upper
            and depth_rejection.between(
                args.calibration_lower, args.calibration_upper
            ).all()
            and sample_rejection.between(
                args.calibration_lower, args.calibration_upper
            ).all()
        )
        if calibrated:
            table.loc[positions, "fdr"] = benjamini_hochberg(
                table.loc[positions, "p_value"].to_numpy()
            )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False)
    table[table["fdr"] <= 0.05].to_csv(
        args.output_dir / "significant_ec_block_glmm.tsv", sep="\t", index=False
    )
    calibrated_null.to_csv(
        args.output_dir / "calibrated_bootstrap_null.tsv.gz",
        sep="\t",
        index=False,
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    summary = []
    for method, records in table.groupby("method"):
        method_eligible = (
            records["null_converged"]
            & records["alternative_converged"]
            & records["p_value"].notna()
        )
        method_null = calibrated_null.loc[
            calibrated_null["method"] == method, "p_value"
        ].to_numpy()
        rejection = float(np.mean(method_null <= 0.05))
        depth_rejection = calibrated_null.loc[
            calibrated_null["method"] == method
        ].groupby("audit_depth_bin")["p_value"].apply(
            lambda values: np.mean(values <= 0.05)
        )
        sample_rejection = calibrated_null.loc[
            calibrated_null["method"] == method
        ].groupby("audit_sample_bin")["p_value"].apply(
            lambda values: np.mean(values <= 0.05)
        )
        calibrated = (
            args.calibration_lower <= rejection <= args.calibration_upper
            and depth_rejection.between(
                args.calibration_lower, args.calibration_upper
            ).all()
            and sample_rejection.between(
                args.calibration_lower, args.calibration_upper
            ).all()
        )
        summary.append({
            "method": method,
            "events": len(records),
            "convergence": float(np.mean(
                records["null_converged"] & records["alternative_converged"]
            )),
            "nominal": int(np.sum(
                method_eligible & (records["p_value"] <= 0.05)
            )),
            "fdr_0.05": int(np.sum(records["fdr"] <= 0.05)),
            "null_replicates": len(method_null),
            "null_rejection_0.05": rejection,
            "null_rejection_0.01": float(np.mean(method_null <= 0.01)),
            "null_rejection_0.001": float(np.mean(method_null <= 0.001)),
            "null_rejection_depth_min": float(depth_rejection.min()),
            "null_rejection_depth_max": float(depth_rejection.max()),
            "null_rejection_sample_min": float(sample_rejection.min()),
            "null_rejection_sample_max": float(sample_rejection.max()),
            "calibration_acceptable": calibrated,
            "calibration_mode": args.calibration_mode,
            "min_stratum_events": min_stratum_events,
            "tail_quantile": (
                args.tail_quantile
                if args.calibration_mode == "gpd_tail"
                else np.nan
            ),
        })
    pd.DataFrame(summary).to_csv(
        args.output_dir / "summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
