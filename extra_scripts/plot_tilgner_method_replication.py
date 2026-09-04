#!/usr/bin/env python3
"""Plot independent long-read directional replication by discovery method."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from plotnine import aes, coord_cartesian, element_blank, element_text, facet_wrap, geom_col, geom_errorbar, geom_point, geom_text, ggplot, labs, scale_fill_manual, scale_x_discrete, theme, theme_bw


METHODS = ["Tealeaf", "LeafCutter", "MAJIQ Heterogen", "scQuint", "Paired junction CLR"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
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
    plot += scale_fill_manual(values={"Tealeaf": "#0B6666", "LeafCutter": "#8C510A", "MAJIQ Heterogen": "#D8B365", "scQuint": "#5AB4AC", "Paired junction CLR": "#762A83"})
    plot += coord_cartesian(ylim=(0, 1))
    plot += labs(x=None, y="Directional replication rate", title="Independent long-read support for full-data discoveries", caption="Parentheses give eligible calls; error bars are Wilson 95% intervals; crosses show method-specific conditional sign-flip null rates.")
    plot += theme_bw(base_size=10)
    plot += theme(panel_grid_minor=element_blank(), axis_text_x=element_text(rotation=30, ha="right"), plot_title=element_text(size=11), plot_caption=element_text(size=8))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot.save(args.output, width=9.0, height=4.8, units="in", verbose=False)


if __name__ == "__main__":
    main()
