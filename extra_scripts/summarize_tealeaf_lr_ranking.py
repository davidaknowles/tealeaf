#!/usr/bin/env python3
"""Compare reporting/ranking sensitivities without modifying statistical tests."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.replication_audit import ranked_direction_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replication", type=Path, required=True)
    parser.add_argument("--pairwise-tests", type=Path, required=True)
    parser.add_argument("--omnibus-tests", type=Path, required=True)
    parser.add_argument("--effects", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source = pd.read_csv(args.replication, sep="\t", low_memory=False)
    source = source.loc[source.mapping_complete.astype(str).str.lower().eq("true") & source.minimum_pooled_depth.ge(20) & source.pooled_replicated.notna()].copy()
    tests = pd.read_csv(args.pairwise_tests, sep="\t", usecols=["test_id", "p_value", "raw_p_value", "statistic", "fdr"])
    source = source.merge(tests, on="test_id", validate="one_to_one")
    omnibus = pd.read_csv(args.omnibus_tests, sep="\t", usecols=["block_id", "p_value", "raw_p_value", "statistic", "fdr"])
    omnibus = omnibus.loc[omnibus.fdr.lt(.05)].drop(columns="fdr")
    effects = pd.concat([pd.read_csv(path, sep="\t", low_memory=False) for path in sorted(args.effects.glob("shard_*.tsv"))], ignore_index=True)
    if effects.duplicated(["test_id", "strategy"]).any():
        raise ValueError("duplicated reporting sensitivity effects")
    primer_effects = effects.loc[effects.strategy.isin(["primer 0 pooled local", "primer 1 pooled local"])].pivot(index="test_id", columns="strategy", values="effect")
    concordance = {}
    for test_id, record in primer_effects.iterrows():
        first = np.asarray(json.loads(record["primer 0 pooled local"]), dtype=float)
        second = np.asarray(json.loads(record["primer 1 pooled local"]), dtype=float)
        norm = np.linalg.norm(first) * np.linalg.norm(second)
        concordance[test_id] = float(first @ second / norm) if np.isfinite(norm) and norm > 1e-12 else -1.
    source["primer_cosine"] = source.test_id.map(concordance).fillna(-1.)
    methods = ["subject-mean A_report=1", *sorted(effects.strategy.unique())]
    summaries, details, production = [], [], []
    for strategy in methods:
        local = source.copy()
        local["report_effect_norm"] = local.original_effect_norm
        local["report_fallback"] = False
        if strategy != methods[0]:
            report = effects.loc[effects.strategy.eq(strategy), ["test_id", "pooled_replicated", "effect_norm", "converged"]].rename(columns={"pooled_replicated": "new_agreement"})
            local = local.merge(report, on="test_id", how="left", validate="one_to_one")
            valid = local.converged.astype(str).str.lower().eq("true") & local.new_agreement.notna()
            # Do not improve a score by dropping failed or zero reporting effects.
            # Retain the original reporting direction as an explicit fallback.
            local.loc[valid, "pooled_replicated"] = local.loc[valid, "new_agreement"]
            local.loc[valid, "report_effect_norm"] = local.loc[valid, "effect_norm"]
            local["report_fallback"] = ~valid
        for unit in ("pairwise", "omnibus"):
            if unit == "pairwise":
                eligible = local.copy()
            else:
                eligible = local.sort_values(["block_id", "report_effect_norm", "test_id"], ascending=[True, False, True], kind="stable").drop_duplicates("block_id")
                eligible = eligible.drop(columns=["p_value", "raw_p_value", "statistic"]).merge(omnibus, on="block_id", validate="one_to_one")
            eligible["effect_evidence_priority"] = eligible.report_effect_norm * -np.log10(np.maximum(eligible.raw_p_value, 1e-300))
            eligible["primer_agrees"] = eligible.primer_cosine.gt(0)
            eligible["primer_evidence_priority"] = eligible.primer_cosine.clip(lower=0) * -np.log10(np.maximum(eligible.raw_p_value, 1e-300))
            for ranking, columns, ascending in (("calibrated then raw", ["p_value", "raw_p_value", "statistic", "test_id"], [True, True, False, True]), ("continuous raw", ["raw_p_value", "p_value", "statistic", "test_id"], [True, True, False, True]), ("effect magnitude, not significance", ["report_effect_norm", "raw_p_value", "test_id"], [False, True, True]), ("effect times log evidence, not significance", ["effect_evidence_priority", "raw_p_value", "test_id"], [False, True, True]), ("primer concordance then raw, not significance", ["primer_agrees", "raw_p_value", "test_id"], [False, True, True]), ("primer cosine times log evidence, not significance", ["primer_evidence_priority", "raw_p_value", "test_id"], [False, True, True])):
                ranked = eligible.sort_values(columns, ascending=ascending, kind="stable").copy()
                ranked["rank"] = np.arange(1, len(ranked) + 1)
                label = f"{unit}; {strategy}; {ranking}"
                ranked["method"] = label
                summaries.extend(ranked_direction_summary(ranked))
                for row in summaries[-3:]:
                    row.update({"unit": unit, "strategy": strategy, "ranking": ranking, "n_fallback": int(ranked.head(row["cutoff"]).report_fallback.sum())})
                details.append(ranked.head(200))
                if strategy == methods[0] and ranking == "calibrated then raw":
                    current = ranked.head(200).copy()
                    current["method"] = "Tealeaf pairwise" if unit == "pairwise" else "Tealeaf omnibus"
                    current["feature_id"] = current.test_id if unit == "pairwise" else current.block_id
                    current["contrast_id"] = "cell_type__" + current.level_a + "__" + current.level_b
                    current["n_agreement"] = current.pooled_replicated.astype(bool).cumsum()
                    current["n_observed"] = current["rank"]
                    current["cumulative_agreement"] = current.n_agreement / current["rank"]
                    production.append(current)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    effects.to_csv(args.output_dir / "reporting_effects.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summaries).to_csv(args.output_dir / "rank_sensitivity_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(details, ignore_index=True).to_csv(args.output_dir / "rank_sensitivity.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.concat(production, ignore_index=True).to_csv(args.output_dir / "production_rank.tsv", sep="\t", index=False, na_rep="NA")
    from plotnine import aes, coord_cartesian, element_blank, geom_hline, geom_line, ggplot, guide_legend, guides, labs, scale_color_manual, theme, theme_bw
    selected = {f"pairwise; {methods[0]}; calibrated then raw": "Production pairwise", f"pairwise; {methods[0]}; continuous raw": "Continuous raw ranking", "pairwise; pooled local; continuous raw": "Pooled reporting + raw", "pairwise; primer-balanced pooled local; continuous raw": "Balanced pooled reporting + raw"}
    plot_table = pd.concat(details, ignore_index=True)
    plot_table = plot_table.loc[plot_table.method.isin(selected)].copy()
    plot_table["method"] = plot_table.method.map(selected)
    plot_table["pooled_replicated"] = plot_table.pooled_replicated.astype(str).str.lower().eq("true").astype(int)
    plot_table["cumulative_agreement"] = plot_table.groupby("method").pooled_replicated.cumsum() / plot_table["rank"]
    plot = ggplot(plot_table, aes("rank", "cumulative_agreement", color="method")) + geom_line() + geom_hline(yintercept=.5, linetype="dashed") + coord_cartesian(ylim=(0, 1)) + scale_color_manual(values={"Production pairwise": "#0B6666", "Continuous raw ranking": "#8C510A", "Pooled reporting + raw": "#CC79A7", "Balanced pooled reporting + raw": "#1B9E77"}) + labs(x="Rank among long-read-evaluable production discoveries", y="Cumulative directional agreement", color=None, caption="Exploratory reporting/ranking sensitivities; production test p-values and discovery sets are unchanged.") + theme_bw() + theme(legend_position="bottom")
    plot += guides(color=guide_legend(nrow=2))
    plot += theme(legend_title=element_blank())
    plot.save(args.output_dir / "ranking_sensitivity.pdf", width=8, height=5.2, verbose=False)
    print(pd.DataFrame(summaries).query("cutoff == 100")[["unit", "strategy", "ranking", "agreement", "normalized_auc", "n_fallback"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
