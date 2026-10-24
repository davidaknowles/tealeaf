#!/usr/bin/env python3
"""Cumulative long-read direction agreement by significance rank, Tealeaf and comparators.

Tealeaf ranks come from assess_block_structure_filter.py (lr_rank.tsv.gz,
scope "all tested", ranking "continuous raw"); comparator ranks from the
long-read rank table of the comparator audit (sign_agreement_by_rank_*.tsv).
The curve at rank r is the fraction of the r top-ranked evaluable tests or
features whose direction agrees with the long reads. Needs plotnine.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as p9

LABELS = {"LeafCutter": "LeafCutter", "MAJIQ Heterogen": "MAJIQ", "scQuint": "scQuint", "SUPPA2 (full data)": "SUPPA2", "SUPPA2 (primer aware)": "SUPPA2 primer-aware", "rMATS paired JCEC": "rMATS"}
COLORS = {"Tealeaf": "#17324D", "MAJIQ": "#168C8C", "LeafCutter": "#7FB3D5", "scQuint": "#E6813B", "rMATS": "#9C6ADE", "SUPPA2": "#999999", "SUPPA2 primer-aware": "#CCCCCC"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tealeaf", type=Path, required=True, help="lr_rank.tsv.gz")
    parser.add_argument("--comparators", type=Path, required=True, help="sign_agreement_by_rank_*.tsv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-rank", type=int, default=200)
    args = parser.parse_args()
    tealeaf = pd.read_csv(args.tealeaf, sep="\t")
    tealeaf = tealeaf.loc[tealeaf.scope.eq("all tested") & tealeaf.ranking.eq("continuous raw")].sort_values("rank").head(args.max_rank)
    agree = tealeaf.pooled_replicated.astype(str).str.lower().eq("true").to_numpy(float)
    curves = [pd.DataFrame({"method": "Tealeaf", "rank": np.arange(1, len(agree) + 1), "agreement": np.cumsum(agree) / np.arange(1, len(agree) + 1)})]
    comparators = pd.read_csv(args.comparators, sep="\t")
    comparators = comparators.loc[comparators.method.isin(LABELS) & comparators["rank"].le(args.max_rank)]
    curves.append(comparators.assign(method=comparators.method.map(LABELS))[["method", "rank", "cumulative_agreement"]].rename(columns={"cumulative_agreement": "agreement"}))
    data = pd.concat(curves, ignore_index=True)
    data["method"] = pd.Categorical(data.method, categories=list(COLORS))
    plot = (p9.ggplot(data, p9.aes("rank", "agreement", color="method"))
            + p9.geom_hline(yintercept=0.5, linetype="dashed", color="#888888")
            + p9.geom_line(size=0.9)
            + p9.scale_color_manual(values=COLORS, name="")
            + p9.scale_y_continuous(limits=(0.4, 1.0))
            + p9.labs(x="Rank (by p-value)", y="Cumulative agreement with long reads", color="")
            + p9.theme_bw() + p9.theme(figure_size=(6.5, 3.8), legend_position="right"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot.save(args.output, verbose=False)
    print(data.groupby("method", observed=True).apply(lambda frame: frame.set_index("rank").agreement.reindex([40, 100, 200]).round(3).tolist(), include_groups=False).to_string(), flush=True)


if __name__ == "__main__":
    main()
