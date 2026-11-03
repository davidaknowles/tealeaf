#!/usr/bin/env python3
"""Split replication and long-read direction of differential intron retention tests.

--tests-dir holds intron_retention_tests.py tables tests_{<fold a>,<fold b>,full}.tsv.gz and, for the calibration,
the label-swapped fits tests_<cohort>_null.tsv.gz (either may be split into _shard<k> files). P-values are chi2_1 on the likelihood-ratio statistic and the
scaled chi2_1 (statistic over the null inflation c of each cohort, intron_retention.null_inflation). A filter
keeps tests whose label-blind pooled retained molecules reach the cutoff. Split replication: intron p-values are
combined by Simes within a gene--cell-type pair and the pairs of a gene by Simes, in each fold, on the genes with a
test in both folds; a gene replicates when its BH-adjusted conjunction p (max over folds) is at most 0.05.

Long-read direction: --lr-counts (counts.tsv.gz of longread_molecule_intron_retention.py or
long_read_intron_retention.py) gives retained and spliced molecules per intron and cell type, matched to the
tests by coordinates. The long-read effect is the difference of retained fractions; a full-data test is evaluable
if both cell types have at least 20 long-read molecules retaining or splicing the intron and either (i) at least 5
retain it or (ii) the retained fractions differ by at least 0.05. Evaluable tests are ranked by p-value.
Writes split_metrics.tsv and lr_curves.tsv.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from extra_scripts.evaluate_suppa2_statistics import grouped_pvalues, metric_row
from tealeaf.sc.intron_retention import calibrated_p_value, null_inflation
from tealeaf.sc.replication_audit import ranked_direction_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tests-dir", type=Path, required=True)
    parser.add_argument("--folds", default="0,1", help="the two fold names")
    parser.add_argument("--introns", type=Path, required=True, help="introns.tsv.gz of intron_retention_counts.py")
    parser.add_argument("--lr-counts", type=Path, required=True)
    parser.add_argument("--filters", default="0,10,20,50")
    parser.add_argument("--exclude-types", default="", help="comma-separated cell types whose pairs are dropped")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    folds = args.folds.split(",")
    excluded = {v for v in args.exclude_types.split(",") if v}

    def kept_pairs(frame):
        return frame[~frame.pair.str.split(r"\|\|", regex=True).apply(lambda types: bool(excluded & set(types)))] if excluded else frame

    introns = pd.read_csv(args.introns, sep="\t")
    introns = introns[introns.clean].drop_duplicates(["chromosome", "start", "end"])
    lr = pd.read_csv(args.lr_counts, sep="\t").drop(columns=["intron", "gene_id"]).merge(introns[["chromosome", "start", "end", "intron"]], on=["chromosome", "start", "end"])
    lr = lr.groupby(["intron", "cell_type"])[["retained", "spliced"]].sum()

    def read_tests(name):
        """tests_<name>.tsv.gz, or its shards tests_<name>_shard<k>.tsv.gz; None if neither exists."""
        paths = [args.tests_dir / f"tests_{name}.tsv.gz"] if (args.tests_dir / f"tests_{name}.tsv.gz").exists() else sorted(args.tests_dir.glob(f"tests_{name}_shard*.tsv.gz"))
        return kept_pairs(pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)) if paths else None

    variants = {}
    for cohort in [*folds, "full"]:
        tests = read_tests(cohort)
        variants.setdefault("chi2", {})[cohort] = tests.assign(p_value=tests.p_chi2)
        null = read_tests(f"{cohort}_null")
        if null is not None:
            null = null[null.converged]
            inflation = null_inflation(null.statistic.to_numpy())
            print(f"cohort {cohort}: {len(null):,} null tests, chi2 rejection at .05 {np.mean(null.p_chi2 < .05):.3f}, .01 {np.mean(null.p_chi2 < .01):.4f}, .001 {np.mean(null.p_chi2 < .001):.5f}; inflation {inflation:.3f}", flush=True)
            variants.setdefault("scaled chi2", {})[cohort] = tests.assign(p_value=calibrated_p_value(tests.statistic, inflation))

    def long_read(full):
        a, b = full.pair.str.split(r"\|\|", expand=True, regex=True).T.to_numpy()
        la = np.nan_to_num(lr.reindex(pd.MultiIndex.from_arrays([full.intron, a])).to_numpy())
        lb = np.nan_to_num(lr.reindex(pd.MultiIndex.from_arrays([full.intron, b])).to_numpy())
        effect = lb[:, 0] / np.maximum(lb.sum(1), 1) - la[:, 0] / np.maximum(la.sum(1), 1)
        depth = np.minimum(la.sum(1), lb.sum(1)) >= 20
        return effect, {"LR retained >= 5": depth & (la[:, 0] + lb[:, 0] >= 5), "abs(LR diff) >= .05": depth & (np.abs(effect) >= .05)}

    metrics, curves = [], []
    for name, tables in variants.items():
        for cutoff in map(float, args.filters.split(",")):
            label = f"{name}, retained >= {cutoff:g}"
            kept = {c: t[t.retained.ge(cutoff) & np.isfinite(t.p_value)] for c, t in tables.items()}
            pairs = [grouped_pvalues(kept[c].assign(gene_id=kept[c].gene_id.str.split(".").str[0], pair_id=kept[c].pair), ["gene_id", "pair_id"], "simes") for c in folds]
            shared = set(map(tuple, pairs[0][["gene_id", "pair_id"]].to_numpy())) & set(map(tuple, pairs[1][["gene_id", "pair_id"]].to_numpy()))
            frame = pd.DataFrame(sorted(shared), columns=["gene_id", "pair_id"])
            metric, _ = metric_row("pairs", label, "simes", "simes", [frame.merge(p, on=["gene_id", "pair_id"], how="left") for p in pairs])
            full = kept["full"]
            effect, criteria = long_read(full)
            for criterion, sel in criteria.items():
                ranked = full.assign(lr_effect=effect)[sel].sort_values(["p_value", "intron"]).copy()
                ranked["rank"] = np.arange(1, len(ranked) + 1)
                ranked["method"] = label
                ranked["pooled_replicated"] = np.sign(ranked.effect) * np.sign(ranked.lr_effect) > 0
                for summary in ranked_direction_summary(ranked):
                    curves.append({**summary, "criterion": criterion})
            metric.update({"tests_full": len(full), "bh_full": int(multipletests(full.p_value, .05, "fdr_bh")[0].sum()), **{f"fold{c}_tests_bh": int(multipletests(kept[c].p_value, .05, "fdr_bh")[0].sum()) for c in folds}})
            metrics.append(metric)
    metrics, curves = pd.DataFrame(metrics), pd.DataFrame(curves)
    metrics.to_csv(args.output_dir / "split_metrics.tsv", sep="\t", index=False)
    curves.to_csv(args.output_dir / "lr_curves.tsv", sep="\t", index=False)
    columns = ["shared_genes", *[f"fold{c}_tests_bh" for c in folds], "fold0_bh", "fold1_bh", "replicated_bh", "tests_full", "bh_full"]
    table = metrics.set_index("method")[[c for c in columns if c in metrics]]
    lr100 = curves.query("cutoff == 100").pivot_table(index="method", columns="criterion", values=["agreement", "n_available"])
    lr100.columns = [f"{a} {'ret5' if 'retained' in b else 'diff05'}" for a, b in lr100.columns]
    print(table.join(lr100).to_string())


if __name__ == "__main__":
    main()
