#!/usr/bin/env python3
"""Plot independent long-read directional replication by discovery method."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from plotnine import aes, coord_cartesian, element_blank, element_text, facet_wrap, geom_col, geom_errorbar, geom_hline, geom_line, geom_point, geom_text, ggplot, labs, scale_color_manual, scale_fill_manual, scale_x_continuous, scale_x_discrete, theme, theme_bw


METHODS = ["Tealeaf", "LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR"]
COLORS = {"Tealeaf": "#0B6666", "LeafCutter": "#8C510A", "MAJIQ Heterogen": "#D8B365", "scQuint": "#5AB4AC", "Paired junction CLR": "#762A83"}


def rank_agreement_table(path, tealeaf_replication_path, tealeaf_significant_path):
    """Return cumulative long-read agreement ordered by each method's p-value."""
    table = pd.read_csv(path, sep="\t", low_memory=False)
    if tealeaf_replication_path is not None:
        tealeaf = pd.read_csv(tealeaf_replication_path, sep="\t", low_memory=False)
        tealeaf = tealeaf.rename(columns={"test_id": "feature_id"})
        tealeaf["method"] = "Tealeaf"
        if tealeaf_significant_path is not None:
            significant = pd.read_csv(tealeaf_significant_path, sep="\t", usecols=["test_id", "p_value"])
            tealeaf = tealeaf.merge(significant.rename(columns={"test_id": "feature_id"}), on="feature_id", how="left")
        table = pd.concat([table, tealeaf], ignore_index=True, sort=False)
    has_long_read_direction = table["pooled_replicated"].notna()
    table["pooled_replicated"] = table["pooled_replicated"].astype(str).str.lower().eq("true")
    table = table[table["method"].isin(METHODS) & table["mapping_complete"].astype(str).str.lower().eq("true") & (table["minimum_pooled_depth"] >= 20) & has_long_read_direction & table["p_value"].notna()].copy()
    table = table.sort_values(["method", "p_value", "feature_id"], kind="stable")
    table["rank"] = table.groupby("method").cumcount() + 1
    table["n_eligible"] = table.groupby("method")["rank"].transform("max")
    table["rank_fraction"] = table["rank"] / table["n_eligible"]
    table["cumulative_agreement"] = table.groupby("method")["pooled_replicated"].cumsum() / table["rank"]
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--replication", type=Path)
    parser.add_argument("--tealeaf-replication", type=Path)
    parser.add_argument("--tealeaf-significant", type=Path)
    parser.add_argument("--rank-output", type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.summary, sep="\t")
    table = table[table["scope"] == "all selected calls"].copy()
    table["method"] = pd.Categorical(table["method"], METHODS, ordered=True)
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
        ranked = rank_agreement_table(args.replication, args.tealeaf_replication, args.tealeaf_significant)
        ranked["method"] = pd.Categorical(ranked["method"], METHODS, ordered=True)
        rank_plot = ggplot(ranked, aes("rank_fraction", "cumulative_agreement", color="method", group="method"))
        rank_plot += geom_hline(yintercept=0.5, linetype="dashed", color="#777777", size=0.4)
        rank_plot += geom_line(size=0.9)
        rank_plot += scale_x_continuous(limits=(0, 1), labels=lambda values: [f"{value:.0%}" for value in values])
        rank_plot += coord_cartesian(ylim=(0, 1))
        rank_plot += scale_color_manual(values=COLORS, drop=False)
        rank_plot += labs(x="Significance rank (fraction of eligible calls)", y="Cumulative positive sign agreement with long reads", title="Long-read agreement across significance rank", caption="Calls are ordered by each method's short-read discovery p-value; only mapped calls with at least 20 pooled long-read UMIs per cell type are included. The dashed line is the 50% orientation null.")
        rank_plot += theme_bw(base_size=10)
        rank_plot += theme(panel_grid_minor=element_blank(), plot_title=element_text(size=11), plot_caption=element_text(size=8), legend_title=element_blank())
        args.rank_output.parent.mkdir(parents=True, exist_ok=True)
        rank_plot.save(args.rank_output, width=7.5, height=4.8, units="in", verbose=False)


if __name__ == "__main__":
    main()
