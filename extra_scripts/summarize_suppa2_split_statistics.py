#!/usr/bin/env python3
"""Score observed and synchronized-null SUPPA2 tests with identical aggregation."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.compare_suppa2_primer_aware_tealeaf import load_tealeaf
from extra_scripts.evaluate_suppa2_statistics import grouped_pvalues, metric_row, normalize_pairs
from tealeaf.sc.ds_benchmark import shared_pair_gene_reproducibility


def score(tables, method, replicate=-1):
    pairs = [grouped_pvalues(table, ["gene_id", "pair_id"], "simes") for table in tables]
    metric, details = metric_row("fixed Table 1 SUPPA2 universe", method, "simes", "simes", pairs)
    metric["replicate"] = replicate
    details["replicate"] = replicate
    return metric, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--primer-aware", action="store_true")
    parser.add_argument("--tealeaf-subdir", default="tealeaf_paired_path_total_a32_production", help="Merged, calibrated split fit selected for manuscript Table 1.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    observed = [normalize_pairs(pd.read_csv(args.audit_dir / f"fold{k}_tests.tsv.gz", sep="\t")) for k in (0, 1)]
    null = [normalize_pairs(pd.read_csv(args.audit_dir / f"fold{k}_null.tsv.gz", sep="\t")) for k in (0, 1)]
    repro = args.run_root / "junction_benchmark/reproducibility"
    production = []
    for k in (0, 1):
        table = pd.read_csv(repro / f"fold{k}" / args.tealeaf_subdir / "paired_path.tsv", sep="\t")
        production.append(normalize_pairs(table.loc[table.method.eq("local_path") & table.converged & table.n_subjects.ge(4)]))
    shared = set.intersection(*(set(zip(table.gene_id, table.pair_id)) for table in observed))
    shared &= set.intersection(*(set(zip(table.gene_id, table.pair_id)) for table in production))
    observed = [table.loc[[(gene, pair) in shared for gene, pair in zip(table.gene_id, table.pair_id)]] for table in observed]
    null = [table.loc[[(gene, pair) in shared for gene, pair in zip(table.gene_id, table.pair_id)]] for table in null]
    for k, table in enumerate(observed):
        table.to_csv(args.output_dir / f"matched_fold{k}_tests.tsv.gz", sep="\t", index=False)
    metrics, details, null_metrics, marginal = [], [], [], []
    for method in sorted(observed[0].method.unique()):
        tables = [table.loc[table.method.eq(method)] for table in observed]
        metric, detail = score(tables, method)
        metrics.append(metric)
        details.append(detail)
        print(f"scored method={method} replicated_genes={metric['replicated_bh']}", flush=True)
        selected_null = [table.loc[table.method.eq(method)] for table in null]
        for k, table in enumerate(selected_null):
            marginal.append({"method": method, "fold": k, "null_tests": len(table), "reject_0_05": np.mean(table.p_value <= .05), "reject_0_01": np.mean(table.p_value <= .01), "reject_0_001": np.mean(table.p_value <= .001)})
        for replicate in sorted(selected_null[0].replicate.unique()):
            local = [table.loc[table.replicate.eq(replicate)] for table in selected_null]
            metric, _ = score(local, method, replicate)
            null_metrics.append(metric)
    for label in ("Tealeaf production merged", "Tealeaf prior cross-fitted EB merged", "Tealeaf legacy raw shards", "SUPPA2 archived approximation"):
        tables = []
        for k in (0, 1):
            if label == "Tealeaf production merged":
                table = production[k]
            elif label == "Tealeaf prior cross-fitted EB merged":
                table = pd.read_csv(repro / f"fold{k}/tealeaf_paired_path_uniform_eb_crossfit_total/paired_path.tsv", sep="\t")
                table = table.loc[table.method.eq("local_path") & table.converged]
            elif label == "Tealeaf legacy raw shards":
                table = load_tealeaf(repro / f"fold{k}/tealeaf_paired_path_uniform_eb_crossfit_total_shards")
            else:
                archive = "suppa2_primer_aware" if args.primer_aware else "suppa2"
                table = pd.read_csv(args.repo_root / f"analyses/comparator_suppa_rmats/{archive}/split_data_fold{k}_tests.tsv.gz", sep="\t")
            table = normalize_pairs(table)
            table = table.loc[[(gene, pair) in shared for gene, pair in zip(table.gene_id, table.pair_id)]]
            tables.append(table)
        metric, detail = score(tables, label)
        metrics.append(metric)
        details.append(detail)
    pd.DataFrame(metrics).to_csv(args.output_dir / "split_reproducibility.tsv", sep="\t", index=False)
    pd.concat(details).to_csv(args.output_dir / "gene_pvalues.tsv.gz", sep="\t", index=False)
    pd.DataFrame(null_metrics).to_csv(args.output_dir / "null_gene_reproducibility.tsv", sep="\t", index=False)
    pd.DataFrame(marginal).to_csv(args.output_dir / "null_event_calibration.tsv", sep="\t", index=False)
    supports = pd.concat([pd.read_csv(args.audit_dir / f"fold{k}_support.tsv.gz", sep="\t") for k in (0, 1)])
    selected_tests = pd.concat(observed)[["fold", "contrast_id", "feature_id"]].drop_duplicates()
    supports = supports.merge(selected_tests, on=["fold", "contrast_id", "feature_id"], how="inner", validate="one_to_one")
    support_summary = supports.groupby("fold")[["n_subjects", "nonzero_subjects", "effect_size"]].describe()
    support_summary.columns = ["_".join(column) for column in support_summary.columns]
    support_summary.to_csv(args.output_dir / "support_summary.tsv", sep="\t", na_rep="NA")
    label = "SUPPA2 (primer aware; matched; exact paired)" if args.primer_aware else "SUPPA2 (matched split; exact paired)"
    folds = []
    for k in (0, 1):
        comparator = observed[k].loc[observed[k].method.eq("wilcoxon_exact")].copy()
        comparator = comparator.merge(supports.loc[supports.fold.eq(k), ["contrast_id", "feature_id", "n_subjects", "effect_size"]], on=["contrast_id", "feature_id"], how="left", validate="one_to_one")
        comparator["method"], comparator["effect"] = label, "cell_type"
        comparator.to_csv(args.output_dir.parent / f"split_data_matched_exact_fold{k}_tests.tsv.gz", sep="\t", index=False)
        reference = production[k].copy()
        reference["method"] = "Tealeaf direct local path"
        folds.append(pd.concat([reference, comparator], ignore_index=True))
    canonical_metrics, topk, genes = shared_pair_gene_reproducibility(folds, reference_method="Tealeaf direct local path")
    for name, table in (("metrics.tsv", canonical_metrics), ("topk.tsv", topk), ("genes.tsv.gz", genes)):
        table.to_csv(args.output_dir.parent / f"split_data_tealeaf_comparison_{name}", sep="\t", index=False)
    print(pd.DataFrame(metrics).to_string(index=False), flush=True)
    print(pd.DataFrame(null_metrics).groupby("method")[["fold0_bh", "fold1_bh", "replicated_bh"]].agg(["mean", "max"]), flush=True)


if __name__ == "__main__":
    main()
