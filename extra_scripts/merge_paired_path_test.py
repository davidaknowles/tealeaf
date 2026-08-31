#!/usr/bin/env python3
"""Merge paired local-path test shards and audit signed-label nulls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

from extra_scripts.run_differential_splicing import benjamini_hochberg
from tealeaf.sc import differential


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-stratum-tests", type=int, default=100)
    parser.add_argument(
        "--calibration",
        choices=("native", "empirical"),
        default="empirical",
    )
    parser.add_argument(
        "--moderate-variances",
        action="store_true",
        help="Moderate scalar paired-test variances across blocks.",
    )
    return parser.parse_args()


def moderate_scalar_tests(table, null):
    """Moderate two-path tests and refit the prior in every null family."""
    table = table.copy()
    null = null.copy()
    selected = table.converged & table.degrees_of_freedom.eq(1)
    local = table.loc[selected]
    sample_size = local.n_subjects.to_numpy(dtype=float)
    residual_df = sample_size - 1
    t_statistics = np.sqrt(
        np.maximum(local.statistic.to_numpy(dtype=float), 0)
    )
    means = local.mean_difference_norm.to_numpy(dtype=float)
    variances = np.divide(
        sample_size * means**2,
        t_statistics**2,
        out=np.full_like(means, np.nan),
        where=t_statistics > 1e-12,
    )
    valid = np.isfinite(variances) & (variances > 0)
    prior_df, prior_variance = differential.fit_variance_prior(
        variances[valid], residual_df[valid]
    )
    table.loc[local.index[valid], "p_value"] = (
        differential.moderated_t_pvalues(
            t_statistics[valid],
            variances[valid],
            residual_df[valid],
            prior_df,
            prior_variance,
        )
    )
    table["variance_prior_df"] = np.nan
    table["variance_prior"] = np.nan
    table.loc[local.index[valid], "variance_prior_df"] = prior_df
    table.loc[local.index[valid], "variance_prior"] = prior_variance
    if null.empty:
        return table, null

    summary = pd.DataFrame({
        "test_id": local.loc[valid, "test_id"].to_numpy(),
        "residual_df": residual_df[valid],
        "zero_centered_ss": (
            residual_df[valid] * variances[valid]
            + sample_size[valid] * means[valid] ** 2
        ),
    }).set_index("test_id")
    null_selected = null.test_id.isin(summary.index)
    null_local = null.loc[null_selected]
    null_df = null_local.test_id.map(summary.residual_df).to_numpy(dtype=float)
    zero_centered_ss = null_local.test_id.map(
        summary.zero_centered_ss
    ).to_numpy(dtype=float)
    null_t = scipy.stats.t.isf(
        np.clip(
            null_local.p_value.to_numpy(dtype=float) / 2,
            1e-300,
            0.5,
        ),
        null_df,
    )
    null_variance = zero_centered_ss / (null_df + null_t**2)
    moderated_null = np.empty(len(null_local), dtype=float)
    replicates = null_local.replicate.to_numpy()
    for replicate in np.unique(replicates):
        positions = replicates == replicate
        replicate_prior_df, replicate_prior_variance = (
            differential.fit_variance_prior(
                null_variance[positions], null_df[positions]
            )
        )
        moderated_null[positions] = differential.moderated_t_pvalues(
            null_t[positions],
            null_variance[positions],
            null_df[positions],
            replicate_prior_df,
            replicate_prior_variance,
        )
    null.loc[null_local.index, "p_value"] = moderated_null
    return table, null


def add_calibration_strata(table, minimum_tests):
    """Group tests by Hotelling dimension and paired-subject count."""
    table = table.copy()
    exact = table.groupby(
        ["degrees_of_freedom", "n_subjects"]
    )["test_id"].transform("nunique")
    by_dimension = table.groupby("degrees_of_freedom")[
        "test_id"
    ].transform("nunique")
    table["calibration_stratum"] = (
        table.degrees_of_freedom.astype(str)
        + "|"
        + table.n_subjects.astype(str)
    )
    rare_exact = exact < int(minimum_tests)
    table.loc[rare_exact, "calibration_stratum"] = (
        table.loc[rare_exact, "degrees_of_freedom"].astype(str) + "|all"
    )
    table.loc[
        rare_exact & by_dimension.lt(int(minimum_tests)),
        "calibration_stratum",
    ] = "rare"
    return table


def empirical_null_calibration(table, null):
    """Calibrate analytic paired-test tails against signed-label nulls."""
    table = table.copy()
    null = null.copy()
    table["raw_p_value"] = table["p_value"]
    null["raw_p_value"] = null["p_value"]
    lookup = table.set_index("test_id")["calibration_stratum"]
    null["calibration_stratum"] = null.test_id.map(lookup)
    table["p_value"] = np.nan
    null["p_value"] = np.nan
    for stratum, positions in null.groupby("calibration_stratum").groups.items():
        pool = null.loc[positions]
        all_values = np.sort(pool.raw_p_value.to_numpy(dtype=float))
        own_values = {
            test_id: np.sort(records.raw_p_value.to_numpy(dtype=float))
            for test_id, records in pool.groupby("test_id")
        }
        observed_positions = table.index[
            table.calibration_stratum.eq(stratum)
        ]
        for position in observed_positions:
            value = float(table.at[position, "raw_p_value"])
            own = own_values.get(table.at[position, "test_id"], np.empty(0))
            count = (
                np.searchsorted(all_values, value, side="right")
                - np.searchsorted(own, value, side="right")
            )
            table.at[position, "p_value"] = (
                1 + count
            ) / (1 + len(all_values) - len(own))
        for test_id, local_positions in pool.groupby("test_id").groups.items():
            values = null.loc[local_positions, "raw_p_value"].to_numpy(
                dtype=float
            )
            own = own_values[test_id]
            counts = (
                np.searchsorted(all_values, values, side="right")
                - np.searchsorted(own, values, side="right")
            )
            null.loc[local_positions, "p_value"] = (
                1 + counts
            ) / (1 + len(all_values) - len(own))
    return table, null


def main():
    args = parse_args()
    tables = []
    nulls = []
    failures = []
    for shard in args.shards:
        table_path = shard / "paired_path.tsv"
        null_path = shard / "paired_path_null.tsv.gz"
        if table_path.is_file() and table_path.stat().st_size:
            tables.append(pd.read_csv(table_path, sep="\t"))
        if null_path.is_file() and null_path.stat().st_size:
            nulls.append(pd.read_csv(null_path, sep="\t"))
        failure_path = shard / "failures.json"
        if failure_path.is_file():
            failures.extend(json.loads(failure_path.read_text()))
    if not tables:
        raise ValueError("no paired path results found")
    table = pd.concat(tables, ignore_index=True).drop_duplicates("test_id")
    null = pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame()
    if args.moderate_variances:
        table, null = moderate_scalar_tests(table, null)
    table = add_calibration_strata(table, args.min_stratum_tests)
    if not null.empty and args.calibration == "empirical":
        table, null = empirical_null_calibration(table, null)
    else:
        table["raw_p_value"] = table["p_value"]
        if not null.empty:
            null["raw_p_value"] = null["p_value"]
            lookup = table.set_index("test_id")["calibration_stratum"]
            null["calibration_stratum"] = null.test_id.map(lookup)
    table["fdr"] = np.nan
    eligible = table.converged & table.p_value.notna()
    table.loc[eligible, "fdr"] = benjamini_hochberg(
        table.loc[eligible, "p_value"].to_numpy()
    )
    family_rows = []
    if not null.empty:
        for replicate, records in null.groupby("replicate"):
            values = records.p_value.to_numpy()
            fdr = benjamini_hochberg(values)
            family_rows.append({
                "replicate": replicate,
                "tests": len(records),
                "reject_0.05": int(np.sum(values <= 0.05)),
                "reject_0.01": int(np.sum(values <= 0.01)),
                "reject_0.001": int(np.sum(values <= 0.001)),
                "bh_0.05": int(np.sum(fdr <= 0.05)),
            })
    families = pd.DataFrame(family_rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False)
    table.loc[table.fdr <= 0.05].to_csv(
        args.output_dir / "significant_paired_path.tsv", sep="\t", index=False
    )
    null.to_csv(args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False)
    families.to_csv(args.output_dir / "null_families.tsv", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    summary = {
        "events": len(table),
        "eligible": int(eligible.sum()),
        "nominal_0.05": int(np.sum(eligible & table.p_value.le(0.05))),
        "fdr_0.05": int(np.sum(table.fdr.le(0.05))),
        "null_rejection_0.05": float((null.p_value <= 0.05).mean()) if len(null) else np.nan,
        "null_rejection_0.01": float((null.p_value <= 0.01).mean()) if len(null) else np.nan,
        "null_rejection_0.001": float((null.p_value <= 0.001).mean()) if len(null) else np.nan,
        "null_mean_bh_0.05": float(families["bh_0.05"].mean()) if len(families) else np.nan,
        "null_max_bh_0.05": int(families["bh_0.05"].max()) if len(families) else 0,
        "failures": len(failures),
        "calibration": args.calibration,
        "moderate_variances": args.moderate_variances,
        "retain_uncertainty": bool(table.retain_uncertainty.any()) if "retain_uncertainty" in table else False,
        "minimum_path_pseudocount": float(table.path_pseudocount.min()),
        "median_path_pseudocount": float(table.path_pseudocount.median()),
        "maximum_path_pseudocount": float(table.path_pseudocount.max()),
    }
    pd.DataFrame([summary]).to_csv(
        args.output_dir / "summary.tsv", sep="\t", index=False
    )
    print(pd.DataFrame([summary]).to_string(index=False))


if __name__ == "__main__":
    main()
