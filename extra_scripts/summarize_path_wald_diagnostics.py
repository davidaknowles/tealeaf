#!/usr/bin/env python3
"""Summarize path-Wald structured nulls and split path-effect directions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-observed", required=True, type=Path)
    parser.add_argument(
        "--discovery-null-shards", required=True, nargs="+", type=Path
    )
    parser.add_argument("--fold-observed", nargs=2, type=Path)
    parser.add_argument(
        "--fold-null-shards", action="append", nargs="+", type=Path
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_null_shards(shards):
    tables = []
    for shard in shards:
        path = shard / "null.tsv"
        if path.is_file() and path.stat().st_size:
            try:
                table = pd.read_csv(path, sep="\t")
            except pd.errors.EmptyDataError:
                continue
            if not table.empty:
                tables.append(table)
    if not tables:
        raise ValueError("no nonempty path-Wald null shards")
    return pd.concat(tables, ignore_index=True).drop_duplicates(
        ["test_id", "null_type", "replicate"], keep="last"
    )


def read_observed(path):
    table = pd.read_csv(path, sep="\t")
    table["joint_converged"] = (
        table["null_converged"].astype(bool)
        & table["alternative_converged"].astype(bool)
        & table["p_value"].notna()
    )
    return table.loc[table.joint_converged].drop_duplicates("test_id")


def calibration_tables(observed, null):
    eligible_ids = set(observed.test_id)
    null = null.loc[null.test_id.isin(eligible_ids)].copy()
    family_rows = []
    for (null_type, replicate), table in null.groupby(
        ["null_type", "replicate"]
    ):
        p_value = table.p_value.to_numpy()
        q_value = benjamini_hochberg(p_value)
        family_rows.append({
            "null_type": null_type,
            "replicate": replicate,
            "tests": len(table),
            "p_le_0.05": float(np.mean(p_value <= 0.05)),
            "p_le_0.01": float(np.mean(p_value <= 0.01)),
            "p_le_0.001": float(np.mean(p_value <= 0.001)),
            "bh_0.05": int(np.sum(q_value <= 0.05)),
            "bh_rate": float(np.mean(q_value <= 0.05)),
        })
    families = pd.DataFrame(family_rows)
    summary_rows = []
    for null_type, table in null.groupby("null_type"):
        values = table.p_value.to_numpy()
        local_families = families.loc[families.null_type.eq(null_type)]
        summary_rows.append({
            "null_type": null_type,
            "null_tests": len(table),
            "replicates": table.replicate.nunique(),
            "pooled_p_le_0.05": float(np.mean(values <= 0.05)),
            "pooled_p_le_0.01": float(np.mean(values <= 0.01)),
            "pooled_p_le_0.001": float(np.mean(values <= 0.001)),
            "median_p": float(np.median(values)),
            "uniform_ks_statistic": float(scipy.stats.kstest(values, "uniform").statistic),
            "mean_family_bh_calls": float(local_families["bh_0.05"].mean()),
            "median_family_bh_calls": float(local_families["bh_0.05"].median()),
            "max_family_bh_calls": int(local_families["bh_0.05"].max()),
        })
    summary = pd.DataFrame(summary_rows)
    empirical_rows = []
    observed_statistic = observed.set_index("test_id").statistic
    for (test_id, null_type), table in null.groupby(["test_id", "null_type"]):
        if test_id not in observed_statistic:
            continue
        statistic = float(observed_statistic[test_id])
        empirical_rows.append({
            "test_id": test_id,
            "null_type": null_type,
            "replicates": len(table),
            "observed_statistic": statistic,
            "empirical_p_value": (
                1.0 + float(np.sum(table.statistic >= statistic))
            ) / (len(table) + 1.0),
        })
    empirical = pd.DataFrame(empirical_rows)
    return null, families, summary, empirical


def stratified_null_calibration(null):
    """Summarize analytic null rejection by numerator and cluster df."""
    local = null.copy()
    local["denominator_df_bin"] = pd.cut(
        local.denominator_degrees_of_freedom,
        bins=[-np.inf, 5, 10, 20, np.inf],
        labels=["1-5", "6-10", "11-20", ">20"],
    ).astype(str)
    rows = []
    for keys, table in local.groupby(
        ["null_type", "degrees_of_freedom", "denominator_df_bin"]
    ):
        rows.append({
            "null_type": keys[0],
            "numerator_df": keys[1],
            "denominator_df_bin": keys[2],
            "tests": len(table),
            "p_le_0.05": float(np.mean(table.p_value <= 0.05)),
            "p_le_0.01": float(np.mean(table.p_value <= 0.01)),
            "median_p": float(table.p_value.median()),
        })
    return pd.DataFrame(rows)


def stratified_observed_calibration(observed, null):
    """Map analytic observed p-values through matched structured-null CDFs."""
    denominator_bins = [-np.inf, 5, 10, 20, np.inf]
    denominator_labels = ["1-5", "6-10", "11-20", ">20"]
    observed = observed.copy()
    observed["denominator_df_bin"] = pd.cut(
        observed.denominator_degrees_of_freedom,
        bins=denominator_bins,
        labels=denominator_labels,
    ).astype(str)
    local_null = null.copy()
    local_null["denominator_df_bin"] = pd.cut(
        local_null.denominator_degrees_of_freedom,
        bins=denominator_bins,
        labels=denominator_labels,
    ).astype(str)
    rows = []
    for null_type, type_null in local_null.groupby("null_type"):
        grouped = {
            key: np.sort(table.p_value.to_numpy())
            for key, table in type_null.groupby(
                ["degrees_of_freedom", "denominator_df_bin"]
            )
        }
        for row in observed.itertuples(index=False):
            key = (row.degrees_of_freedom, row.denominator_df_bin)
            reference = grouped.get(key)
            if reference is None or not len(reference):
                continue
            rows.append({
                "test_id": row.test_id,
                "null_type": null_type,
                "analytic_p_value": row.p_value,
                "calibrated_p_value": (
                    1.0 + np.searchsorted(reference, row.p_value, side="right")
                ) / (len(reference) + 1.0),
                "reference_nulls": len(reference),
                "degrees_of_freedom": row.degrees_of_freedom,
                "denominator_df_bin": row.denominator_df_bin,
            })
    calibrated = pd.DataFrame(rows)
    summary_rows = []
    calibrated["fdr"] = np.nan
    for null_type, positions in calibrated.groupby("null_type").groups.items():
        calibrated.loc[positions, "fdr"] = benjamini_hochberg(
            calibrated.loc[positions, "calibrated_p_value"].to_numpy()
        )
        local = calibrated.loc[positions]
        summary_rows.append({
            "null_type": null_type,
            "tests": len(local),
            "calibrated_p_le_0.05": int(np.sum(local.calibrated_p_value <= 0.05)),
            "calibrated_bh_0.05": int(np.sum(local.fdr <= 0.05)),
            "median_calibrated_p": float(local.calibrated_p_value.median()),
        })
    return calibrated, pd.DataFrame(summary_rows)


def decoded_effects(table):
    result = {}
    for row in table.itertuples(index=False):
        try:
            levels = tuple(json.loads(row.tested_levels_json))
            effects = np.asarray(json.loads(row.path_effects_json), dtype=float)
            standard_errors = np.asarray(
                json.loads(row.path_standard_errors_json), dtype=float
            )
        except (AttributeError, TypeError, ValueError, json.JSONDecodeError):
            continue
        result[row.test_id] = (levels, effects, standard_errors, row)
    return result


def pairwise_path_effects(levels, effects):
    """Express centered-log path effects as invariant level-pair contrasts."""
    if effects.ndim != 2 or effects.shape[0] != len(levels) - 1:
        return {}
    by_level = {levels[0]: np.zeros(effects.shape[1], dtype=float)}
    by_level.update({level: effects[index] for index, level in enumerate(levels[1:])})
    return {
        (first, second): by_level[second] - by_level[first]
        for first_index, first in enumerate(levels)
        for second in levels[first_index + 1 :]
    }


def reference_pair_standard_errors(levels, standard_errors):
    """Return SEs for contrasts directly parameterized against the reference."""
    if standard_errors.ndim != 2 or standard_errors.shape[0] != len(levels) - 1:
        return {}
    return {
        (levels[0], level): standard_errors[index]
        for index, level in enumerate(levels[1:])
    }


def split_direction_tables(fold0, fold1):
    decoded0 = decoded_effects(fold0)
    decoded1 = decoded_effects(fold1)
    common = sorted(set(decoded0) & set(decoded1))
    p0 = fold0.set_index("test_id").loc[common, "p_value"].to_numpy()
    p1 = fold1.set_index("test_id").loc[common, "p_value"].to_numpy()
    q0 = dict(zip(common, benjamini_hochberg(p0)))
    q1 = dict(zip(common, benjamini_hochberg(p1)))
    rows = []
    block_rows = []
    for test_id in common:
        levels0, effects0, standard_errors0, row0 = decoded0[test_id]
        levels1, effects1, standard_errors1, row1 = decoded1[test_id]
        pairs0 = pairwise_path_effects(levels0, effects0)
        pairs1 = pairwise_path_effects(levels1, effects1)
        common_pairs = sorted(set(pairs0) & set(pairs1))
        if not common_pairs:
            continue
        standard_error_pairs0 = reference_pair_standard_errors(
            levels0, standard_errors0
        )
        standard_error_pairs1 = reference_pair_standard_errors(
            levels1, standard_errors1
        )
        flat0 = np.concatenate([pairs0[pair] for pair in common_pairs])
        flat1 = np.concatenate([pairs1[pair] for pair in common_pairs])
        norm = np.linalg.norm(flat0) * np.linalg.norm(flat1)
        block_rows.append({
            "test_id": test_id,
            "gene_id": row0.gene_id,
            "coordinates": len(flat0),
            "cosine_similarity": float(flat0 @ flat1 / norm) if norm else np.nan,
            "fold0_p_value": row0.p_value,
            "fold1_p_value": row1.p_value,
            "both_nominal": bool(row0.p_value <= 0.05 and row1.p_value <= 0.05),
            "both_bh": bool(q0[test_id] <= 0.05 and q1[test_id] <= 0.05),
        })
        for pair in common_pairs:
            se_values0 = standard_error_pairs0.get(
                pair, np.full_like(pairs0[pair], np.nan)
            )
            se_values1 = standard_error_pairs1.get(
                pair, np.full_like(pairs1[pair], np.nan)
            )
            for path in range(len(pairs0[pair])):
                effect0 = pairs0[pair][path]
                effect1 = pairs1[pair][path]
                se0 = se_values0[path]
                se1 = se_values1[path]
                rows.append({
                    "test_id": test_id,
                    "gene_id": row0.gene_id,
                    "level_a": pair[0],
                    "level_b": pair[1],
                    "path_index": path,
                    "fold0_effect": effect0,
                    "fold1_effect": effect1,
                    "fold0_standard_error": se0,
                    "fold1_standard_error": se1,
                    "fold0_z": effect0 / se0 if np.isfinite(se0) and se0 > 0 else np.nan,
                    "fold1_z": effect1 / se1 if np.isfinite(se1) and se1 > 0 else np.nan,
                    "sign_agreement": bool(effect0 * effect1 > 0),
                    "both_nominal": bool(
                        row0.p_value <= 0.05 and row1.p_value <= 0.05
                    ),
                    "both_bh": bool(q0[test_id] <= 0.05 and q1[test_id] <= 0.05),
                })
    coordinates = pd.DataFrame(rows)
    blocks = pd.DataFrame(block_rows)
    summary_rows = []
    masks = {
        "all coordinates": np.ones(len(coordinates), dtype=bool),
        "both block p<=0.05": coordinates.both_nominal.to_numpy(),
        "both block BH<=0.05": coordinates.both_bh.to_numpy(),
        "both |z|>=1": (
            coordinates.fold0_z.abs().ge(1)
            & coordinates.fold1_z.abs().ge(1)
        ).to_numpy(),
        "both |z|>=1.96": (
            coordinates.fold0_z.abs().ge(1.96)
            & coordinates.fold1_z.abs().ge(1.96)
        ).to_numpy(),
    }
    for subset, mask in masks.items():
        local = coordinates.loc[mask]
        summary_rows.append({
            "subset": subset,
            "coordinates": len(local),
            "blocks": local.test_id.nunique(),
            "sign_agreement": float(local.sign_agreement.mean()) if len(local) else np.nan,
            "effect_spearman": (
                float(scipy.stats.spearmanr(local.fold0_effect, local.fold1_effect).statistic)
                if len(local) > 1 else np.nan
            ),
        })
    return coordinates, blocks, pd.DataFrame(summary_rows)


def decoded_null_effects(table, eligible_ids):
    result = {}
    local = table.loc[table.test_id.isin(eligible_ids)]
    for row in local.itertuples(index=False):
        try:
            levels = tuple(json.loads(row.tested_levels_json))
            effects = np.asarray(json.loads(row.path_effects_json), dtype=float)
        except (AttributeError, TypeError, ValueError, json.JSONDecodeError):
            continue
        result[(row.test_id, row.null_type, row.replicate)] = (levels, effects)
    return result


def null_split_sign_agreement(null0, null1, eligible_ids):
    decoded0 = decoded_null_effects(null0, eligible_ids)
    decoded1 = decoded_null_effects(null1, eligible_ids)
    rows = []
    for key in sorted(set(decoded0) & set(decoded1)):
        levels0, effects0 = decoded0[key]
        levels1, effects1 = decoded1[key]
        pairs0 = pairwise_path_effects(levels0, effects0)
        pairs1 = pairwise_path_effects(levels1, effects1)
        common_pairs = set(pairs0) & set(pairs1)
        if not common_pairs:
            continue
        values0 = np.concatenate([pairs0[pair] for pair in sorted(common_pairs)])
        values1 = np.concatenate([pairs1[pair] for pair in sorted(common_pairs)])
        rows.append({
            "test_id": key[0],
            "null_type": key[1],
            "replicate": key[2],
            "coordinates": len(values0),
            "agreements": int(np.sum(values0 * values1 > 0)),
        })
    table = pd.DataFrame(rows)
    if table.empty:
        return table, table
    replicate = table.groupby(["null_type", "replicate"], as_index=False).agg(
        coordinates=("coordinates", "sum"), agreements=("agreements", "sum")
    )
    replicate["sign_agreement"] = replicate.agreements / replicate.coordinates
    summary = replicate.groupby("null_type", as_index=False).agg(
        replicates=("replicate", "nunique"),
        mean_sign_agreement=("sign_agreement", "mean"),
        sd_sign_agreement=("sign_agreement", "std"),
        min_sign_agreement=("sign_agreement", "min"),
        max_sign_agreement=("sign_agreement", "max"),
    )
    return replicate, summary


def save_plots(null, coordinates, output_dir):
    from plotnine import (
        aes,
        facet_wrap,
        geom_abline,
        geom_histogram,
        geom_point,
        ggplot,
        labs,
        theme_bw,
    )

    histogram = (
        ggplot(null, aes("p_value"))
        + geom_histogram(binwidth=0.025, boundary=0)
        + facet_wrap("null_type", scales="free_y")
        + theme_bw()
        + labs(x="Analytic null p-value", y="Count")
    )
    histogram.save(output_dir / "null_pvalue_histogram.pdf", width=7, height=3.5)
    if len(coordinates):
        sampled = coordinates.sample(min(len(coordinates), 30000), random_state=3)
        scatter = (
            ggplot(sampled, aes("fold0_effect", "fold1_effect", color="both_bh"))
            + geom_abline(slope=1, intercept=0, alpha=0.4)
            + geom_point(alpha=0.2, size=0.6)
            + theme_bw()
            + labs(x="Fold 0 centered-log path effect", y="Fold 1 centered-log path effect")
        )
        scatter.save(output_dir / "split_path_effects.pdf", width=5.5, height=5)


def main():
    args = parse_args()
    if (args.fold_observed is None) != (args.fold_null_shards is None):
        raise ValueError("fold observed and null inputs must be supplied together")
    if args.fold_null_shards is not None and len(args.fold_null_shards) != 2:
        raise ValueError("provide exactly two --fold-null-shards groups")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    discovery = read_observed(args.discovery_observed)
    discovery_null = read_null_shards(args.discovery_null_shards)
    filtered_null, families, calibration, empirical = calibration_tables(
        discovery, discovery_null
    )
    calibration_by_df = stratified_null_calibration(filtered_null)
    calibrated_observed, calibrated_observed_summary = (
        stratified_observed_calibration(discovery, filtered_null)
    )
    coordinates = pd.DataFrame()
    blocks = pd.DataFrame()
    directions = pd.DataFrame()
    null_replicates = pd.DataFrame()
    null_directions = pd.DataFrame()
    if args.fold_observed is not None:
        folds = [read_observed(path) for path in args.fold_observed]
        coordinates, blocks, directions = split_direction_tables(*folds)
        fold_null = [read_null_shards(shards) for shards in args.fold_null_shards]
        null_replicates, null_directions = null_split_sign_agreement(
            fold_null[0], fold_null[1], set(coordinates.test_id)
        )
    outputs = {
        "null_family_calibration.tsv": families,
        "null_calibration_summary.tsv": calibration,
        "observed_empirical_pvalues.tsv": empirical,
        "null_calibration_by_df.tsv": calibration_by_df,
        "observed_stratified_calibration.tsv": calibrated_observed,
        "observed_stratified_calibration_summary.tsv": calibrated_observed_summary,
    }
    if args.fold_observed is not None:
        outputs.update({
            "split_path_coordinates.tsv.gz": coordinates,
            "split_block_directions.tsv": blocks,
            "split_direction_summary.tsv": directions,
            "null_split_direction_replicates.tsv": null_replicates,
            "null_split_direction_summary.tsv": null_directions,
        })
    for name, table in outputs.items():
        table.to_csv(args.output_dir / name, sep="\t", index=False)
    save_plots(filtered_null, coordinates, args.output_dir)
    print(calibration.to_string(index=False))
    if args.fold_observed is not None:
        print(directions.to_string(index=False))
        print(null_directions.to_string(index=False))


if __name__ == "__main__":
    main()
