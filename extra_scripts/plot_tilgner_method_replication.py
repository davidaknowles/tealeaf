#!/usr/bin/env python3
"""Plot independent long-read directional replication by discovery method."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from plotnine import aes, coord_cartesian, element_blank, element_text, facet_wrap, geom_col, geom_errorbar, geom_hline, geom_line, geom_point, geom_text, ggplot, labs, scale_color_manual, scale_fill_manual, scale_x_continuous, scale_x_discrete, theme, theme_bw


METHODS = ["Tealeaf pairwise", "Tealeaf omnibus", "LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR", "rMATS", "SUPPA transcript PSI"]
COLORS = {"Tealeaf": "#0B6666", "Tealeaf pairwise": "#0B6666", "Tealeaf omnibus": "#1B9E77", "Isoform usage ratio": "#E7298A", "LeafCutter": "#8C510A", "MAJIQ Heterogen": "#D8B365", "scQuint": "#5AB4AC", "Paired junction CLR": "#762A83", "rMATS": "#4D4D4D", "SUPPA transcript PSI": "#CC79A7"}


def _rank_table(table, max_rank, method_column="method"):
    """Add tie-aware significance ranks and cumulative agreement."""
    table = table.copy()
    table["p_value"] = pd.to_numeric(table["p_value"], errors="coerce")
    table["raw_p_value"] = pd.to_numeric(table["raw_p_value"], errors="coerce") if "raw_p_value" in table else table["p_value"]
    table["statistic"] = pd.to_numeric(table["statistic"], errors="coerce") if "statistic" in table else pd.Series(np.nan, index=table.index)
    table["_sort_raw_p"] = table["raw_p_value"].fillna(table["p_value"])
    table["_sort_statistic"] = -table["statistic"].fillna(-float("inf"))
    table = table.sort_values([method_column, "p_value", "_sort_raw_p", "_sort_statistic", "feature_id"], kind="stable")
    table["rank"] = table.groupby(method_column).cumcount() + 1
    table["p_tie_size"] = table.groupby([method_column, "p_value"])["feature_id"].transform("size")
    table = table[table["rank"] <= max_rank].copy()
    table["n_evaluable"] = table.groupby(method_column)["rank"].transform("size")
    table["n_agreement"] = table.groupby(method_column)["pooled_replicated"].cumsum()
    table["n_observed"] = table.groupby(method_column).cumcount() + 1
    table["cumulative_agreement"] = table["n_agreement"] / table["n_observed"]
    return table.drop(columns=["_sort_raw_p", "_sort_statistic"])


def _isoform_ratio_table(replication_path, ratio_path):
    """Build a long-read table for the path-log-ratio baseline.

    The ratio baseline is the previously implemented path-composition
    likelihood test. Its fitted cell-type coefficients are in reference
    log-ratio coordinates, so source path proportions are compared after the
    same additive 0.5 count stabilization.
    """
    replication = pd.read_csv(replication_path, sep="\t", low_memory=False)
    ratio = pd.read_csv(ratio_path, sep="\t", low_memory=False)
    ratio = ratio[ratio["fdr"].astype(float) < 0.05].copy()
    if "contrast" in ratio.columns and ratio["contrast"].astype(str).str.contains("_vs_").any():
        replication = replication.set_index(["block_id", "level_a", "level_b"])
        rows = []
        for model in ratio.itertuples(index=False):
            level_a, level_b = str(model.contrast).split("_vs_", 1)
            key = (model.block_id, level_a, level_b)
            if key not in replication.index:
                continue
            record = replication.loc[key]
            if isinstance(record, pd.DataFrame):
                record = record.iloc[0]
            if not (bool(record.mapping_complete) and float(record.minimum_pooled_depth) >= 20 and pd.notna(record.pooled_replicated)):
                continue
            try:
                effect = np.asarray(json.loads(model.celltype_logit_effects), dtype=float).ravel()
                counts_a = np.asarray(json.loads(record.counts_a_rep1), dtype=float) + np.asarray(json.loads(record.counts_a_rep2), dtype=float)
                counts_b = np.asarray(json.loads(record.counts_b_rep1), dtype=float) + np.asarray(json.loads(record.counts_b_rep2), dtype=float)
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if effect.size != counts_a.size - 1 or counts_a.size != counts_b.size:
                continue
            logratio_delta = np.log((counts_b + 0.5) / (counts_b.sum() + 0.5 * counts_b.size)) - np.log((counts_a + 0.5) / (counts_a.sum() + 0.5 * counts_a.size))
            permutation_p = getattr(model, "permutation_p_value", model.p_value)
            rows.append({"method": "Isoform usage ratio", "feature_id": model.test_id, "block_id": model.block_id, "p_value": float(model.p_value), "raw_p_value": float(permutation_p) if pd.notna(permutation_p) else float(model.p_value), "statistic": float(model.statistic), "pooled_replicated": bool(effect @ logratio_delta[1:] > 0), "mapping_complete": True, "minimum_pooled_depth": float(record.minimum_pooled_depth)})
        return pd.DataFrame(rows)
    ratio = ratio.set_index("block_id")
    rows = []
    for block_id, group in replication.groupby("block_id", sort=False):
        if block_id not in ratio.index:
            continue
        model = ratio.loc[block_id]
        if isinstance(model, pd.DataFrame):
            model = model.iloc[0]
        if not (bool(model.get("converged", False)) and int(model.n_paths) >= 2):
            continue
        try:
            cell_types = str(model.cell_types).split(",")
            effects = np.asarray(json.loads(model.celltype_logit_effects), dtype=float)
            n_paths = int(model.n_paths)
        except (TypeError, ValueError, json.JSONDecodeError):
            continue
        if effects.shape != (len(cell_types) - 1, n_paths - 1):
            continue
        coefficients = {cell_types[0]: np.zeros(n_paths - 1)}
        coefficients.update({name: effects[i] for i, name in enumerate(cell_types[1:])})
        candidates = []
        for record in group.itertuples(index=False):
            if (record.level_a not in coefficients or record.level_b not in coefficients or int(record.n_paths) != n_paths or not bool(record.mapping_complete) or float(record.minimum_pooled_depth) < 20 or pd.isna(record.pooled_replicated)):
                continue
            counts_a = np.asarray(json.loads(record.counts_a_rep1), dtype=float) + np.asarray(json.loads(record.counts_a_rep2), dtype=float)
            counts_b = np.asarray(json.loads(record.counts_b_rep1), dtype=float) + np.asarray(json.loads(record.counts_b_rep2), dtype=float)
            if counts_a.size != n_paths or counts_b.size != n_paths:
                continue
            logratio_delta = np.log((counts_b + 0.5) / (counts_b.sum() + 0.5 * n_paths))
            logratio_delta -= np.log((counts_a + 0.5) / (counts_a.sum() + 0.5 * n_paths))
            delta = coefficients[record.level_b] - coefficients[record.level_a]
            candidates.append((float(np.linalg.norm(delta)), record, float(delta @ logratio_delta[1:])))
        if candidates:
            _, record, dot = max(candidates, key=lambda item: (item[0], item[1].test_id))
            rows.append({"method": "Isoform usage ratio", "feature_id": block_id, "block_id": block_id, "p_value": float(model.p_value), "raw_p_value": float(model.get("permutation_p_value", model.p_value)), "statistic": float(model.statistic), "pooled_replicated": dot > 0, "mapping_complete": True, "minimum_pooled_depth": float(record.minimum_pooled_depth)})
    return pd.DataFrame(rows)


def rank_agreement_table(path, tealeaf_replication_path, tealeaf_significant_path, omnibus_path=None, isoform_ratio_path=None, max_rank=200):
    """Return cumulative agreement with continuous tie breaks and optional omnibus ranking."""
    table = pd.read_csv(path, sep="\t", low_memory=False)
    if tealeaf_replication_path is not None:
        tealeaf = pd.read_csv(tealeaf_replication_path, sep="\t", low_memory=False)
        tealeaf = tealeaf.rename(columns={"test_id": "feature_id"})
        tealeaf["method"] = "Tealeaf pairwise"
        if tealeaf_significant_path is not None:
            significant = pd.read_csv(tealeaf_significant_path, sep="\t", usecols=["test_id", "p_value", "raw_p_value", "statistic"])
            tealeaf = tealeaf.merge(significant.rename(columns={"test_id": "feature_id"}), on="feature_id", how="left")
        table = pd.concat([table, tealeaf], ignore_index=True, sort=False)
    has_long_read_direction = table["pooled_replicated"].notna()
    table["pooled_replicated"] = table["pooled_replicated"].astype(str).str.lower().eq("true")
    table = table[table["method"].isin(["LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR", "Tealeaf pairwise", "rMATS", "SUPPA transcript PSI"]) & table["mapping_complete"].astype(str).str.lower().eq("true") & (table["minimum_pooled_depth"] >= 20) & has_long_read_direction & table["p_value"].notna()].copy()
    tables = [_rank_table(table, max_rank)]
    if omnibus_path is not None:
        omnibus = pd.read_csv(omnibus_path, sep="\t", low_memory=False)
        omnibus["p_value"] = pd.to_numeric(omnibus.get("empirical_p_value", omnibus.get("p_value")), errors="coerce")
        omnibus = omnibus[omnibus["fdr"].astype(float) < 0.05]
        eligible = table[table["method"].eq("Tealeaf pairwise")].copy()
        eligible = eligible.sort_values(["block_id", "original_effect_norm", "feature_id"], ascending=[True, False, True], kind="stable").drop_duplicates("block_id")
        eligible = eligible.merge(omnibus[["block_id", "p_value"]], on="block_id", how="inner", suffixes=("_pairwise", ""))
        eligible["method"] = "Tealeaf omnibus"
        eligible["feature_id"] = eligible["block_id"]
        tables.append(_rank_table(eligible, max_rank))
    if isoform_ratio_path is not None:
        ratio = _isoform_ratio_table(tealeaf_replication_path, isoform_ratio_path)
        if not ratio.empty:
            tables.append(_rank_table(ratio, max_rank))
    return pd.concat(tables, ignore_index=True, sort=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replication", type=Path)
    parser.add_argument("--tealeaf-replication", type=Path)
    parser.add_argument("--tealeaf-significant", type=Path)
    parser.add_argument("--tealeaf-omnibus", type=Path)
    parser.add_argument("--isoform-ratio", type=Path)
    parser.add_argument("--rank-table", type=Path)
    parser.add_argument("--rank-output", type=Path)
    parser.add_argument("--max-rank", type=int, default=200)
    args = parser.parse_args()
    table = pd.read_csv(args.summary, sep="\t")
    table = table[table["scope"] == "all selected calls"].copy()
    table["method"] = pd.Categorical(table["method"], ["Tealeaf", "LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR", "rMATS", "SUPPA transcript PSI"], ordered=True)
    table["label"] = table.apply(lambda row: f"{row.replication_rate:.0%}\n({int(row.n_tests)})", axis=1)
    plot = ggplot(table, aes("method", "replication_rate", fill="method"))
    plot += geom_col(width=0.72, show_legend=False)
    plot += geom_errorbar(aes(ymin="ci_low", ymax="ci_high"), width=0.16)
    plot += geom_point(aes(y="conditional_null_rate"), shape="x", size=2.5, color="#222222")
    plot += geom_text(aes(label="label"), va="bottom", size=7, nudge_y=0.025)
    plot += facet_wrap("~ endpoint")
    plot += scale_x_discrete(drop=False)
    plot += scale_fill_manual(values=COLORS)
    plot += coord_cartesian(ylim=(0, 1))
    plot += labs(x=None, y="Directional replication rate", title="Independent long-read support for full-data discoveries", caption="Parentheses give eligible calls; error bars are Wilson 95% intervals; crosses show method-specific conditional sign-flip null rates.")
    plot += theme_bw(base_size=10)
    plot += theme(panel_grid_minor=element_blank(), axis_text_x=element_text(rotation=30, ha="right"), plot_title=element_text(size=11), plot_caption=element_text(size=8))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot.save(args.output, width=9.0, height=4.8, units="in", verbose=False)
    if args.rank_output:
        if args.replication is None:
            parser.error("--rank-output requires --replication")
        if args.tealeaf_replication is None:
            parser.error("--rank-output requires --tealeaf-replication")
        ranked = rank_agreement_table(args.replication, args.tealeaf_replication, args.tealeaf_significant, omnibus_path=args.tealeaf_omnibus, isoform_ratio_path=args.isoform_ratio, max_rank=args.max_rank)
        if args.rank_table:
            args.rank_table.parent.mkdir(parents=True, exist_ok=True)
            ranked.to_csv(args.rank_table, sep="\t", index=False, na_rep="NA")
        ranked["method"] = pd.Categorical(ranked["method"], METHODS, ordered=True)
        rank_plot = ggplot(ranked, aes("rank", "cumulative_agreement", color="method", group="method"))
        rank_plot += geom_hline(yintercept=0.5, linetype="dashed", color="#777777", size=0.4)
        rank_plot += geom_line(size=0.9)
        rank_plot += scale_x_continuous(limits=(1, args.max_rank), breaks=list(range(0, args.max_rank + 1, 25))[1:])
        rank_plot += coord_cartesian(ylim=(0, 1))
        rank_plot += scale_color_manual(values=COLORS, drop=False)
        rank_plot += labs(x=f"Significance rank (top {args.max_rank} calls per method)", y="Cumulative positive sign agreement with long reads", title="Long-read agreement across significance rank", caption="Calls are ordered by calibrated p-value, then continuous raw p-value and test statistic within ties. Only mapped calls with at least 20 pooled long-read UMIs per cell type and a finite direction contribute. Tealeaf omnibus points use the strongest Tealeaf pairwise direction per significant block; the dashed line is the 50% orientation null.")
        rank_plot += theme_bw(base_size=10)
        rank_plot += theme(panel_grid_minor=element_blank(), plot_title=element_text(size=11), plot_caption=element_text(size=8), legend_title=element_blank())
        args.rank_output.parent.mkdir(parents=True, exist_ok=True)
        rank_plot.save(args.rank_output, width=7.5, height=4.8, units="in", verbose=False)


if __name__ == "__main__":
    main()
