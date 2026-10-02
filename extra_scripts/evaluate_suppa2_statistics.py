#!/usr/bin/env python3
"""Compare SUPPA2 event statistics on the matched Tealeaf split universe."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from extra_scripts.compare_suppa2_primer_aware_tealeaf import load_merged_tealeaf, load_tealeaf
from tealeaf.sc.ds_benchmark import (
    benjamini_hochberg,
    cauchy_pvalue,
    simes_pvalue,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    reference = parser.add_mutually_exclusive_group(required=True)
    reference.add_argument("--tealeaf-shards", action="append", type=Path)
    reference.add_argument("--tealeaf-tests", action="append", type=Path, help="Production merged paired_path.tsv; repeat for both folds.")
    parser.add_argument("--comparison", action="append", nargs=3, metavar=("LABEL", "FOLD0", "FOLD1"), required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def normalize_pairs(table):
    table = table.loc[table.p_value.notna()].copy()
    table["gene_id"] = table.gene_id.astype(str).str.split(".").str[0]
    levels = np.sort(table[["level_a", "level_b"]].astype(str).to_numpy(), axis=1)
    table["pair_id"] = levels[:, 0] + "||" + levels[:, 1]
    return table


def combine(values, method):
    return simes_pvalue(values) if method == "simes" else cauchy_pvalue(values)


def grouped_pvalues(table, keys, method):
    local = table.loc[:, [*keys, "p_value"]]
    return (
        local.groupby(keys, sort=False)["p_value"]
        .agg(lambda values: combine(values, method))
        .rename("p_value")
        .reset_index()
    )


def metric_row(label, method, pair_combination, gene_combination, pair_tables):
    genes = []
    for table in pair_tables:
        genes.append(grouped_pvalues(table, ["gene_id"], gene_combination).set_index("gene_id"))
    shared_genes = sorted(set(genes[0].index) & set(genes[1].index))
    p0 = genes[0].loc[shared_genes, "p_value"].to_numpy(float)
    p1 = genes[1].loc[shared_genes, "p_value"].to_numpy(float)
    q0 = benjamini_hochberg(p0)
    q1 = benjamini_hochberg(p1)
    conjunction_q = benjamini_hochberg(np.maximum(p0, p1))
    selected0, selected1 = q0 <= 0.05, q1 <= 0.05
    held = []
    if selected0.any():
        held.append(float(np.mean(p1[selected0] <= 0.05)))
    if selected1.any():
        held.append(float(np.mean(p0[selected1] <= 0.05)))
    rho = float(spearmanr(-np.log10(np.maximum(p0, 1e-300)), -np.log10(np.maximum(p1, 1e-300))).statistic)
    metric = {
        "comparison": label,
        "method": method,
        "pair_combination": pair_combination,
        "gene_combination": gene_combination,
        "shared_gene_pairs": len(set(map(tuple, pair_tables[0][["gene_id", "pair_id"]].to_numpy()))),
        "shared_genes": len(shared_genes),
        "fold0_bh": int(selected0.sum()),
        "fold1_bh": int(selected1.sum()),
        "replicated_bh": int((conjunction_q <= 0.05).sum()),
        "heldout_nominal_replication": float(np.mean(held)) if held else np.nan,
        "spearman_logp": rho,
    }
    details = pd.DataFrame({
        "comparison": label,
        "method": method,
        "pair_combination": pair_combination,
        "gene_combination": gene_combination,
        "gene_id": shared_genes,
        "fold0_p_value": p0,
        "fold1_p_value": p1,
        "conjunction_q_value": conjunction_q,
    })
    return metric, details


def main():
    args = parse_args()
    paths = args.tealeaf_tests or args.tealeaf_shards
    if len(paths) != 2:
        raise ValueError("provide exactly two Tealeaf fold inputs")
    if args.tealeaf_tests:
        tealeaf = [normalize_pairs(load_merged_tealeaf(path)) for path in paths]
    else:
        tealeaf = [normalize_pairs(load_tealeaf(path)) for path in paths]
    metrics, diagnostics, details = [], [], []
    for label, fold0_path, fold1_path in args.comparison:
        folds = [
            normalize_pairs(pd.read_csv(path, sep="\t", compression="infer", low_memory=False))
            for path in (fold0_path, fold1_path)
        ]
        methods = sorted(set(folds[0].method) & set(folds[1].method))
        for method in methods:
            selected = [table.loc[table.method.eq(method)].copy() for table in folds]
            shared = set(map(tuple, tealeaf[0][["gene_id", "pair_id"]].to_numpy()))
            shared &= set(map(tuple, tealeaf[1][["gene_id", "pair_id"]].to_numpy()))
            shared &= set(map(tuple, selected[0][["gene_id", "pair_id"]].to_numpy()))
            shared &= set(map(tuple, selected[1][["gene_id", "pair_id"]].to_numpy()))
            selected = [
                table.loc[[tuple(value) in shared for value in table[["gene_id", "pair_id"]].to_numpy()]]
                for table in selected
            ]
            for fold, table in enumerate(selected):
                values = table.p_value.to_numpy(float)
                counts = table.groupby(["gene_id", "pair_id"]).size()
                diagnostics.append({
                    "comparison": label,
                    "method": method,
                    "fold": fold,
                    "event_tests": len(table),
                    "gene_pairs": len(counts),
                    "median_events_per_pair": float(counts.median()),
                    "maximum_events_per_pair": int(counts.max()),
                    "unique_p_values": int(pd.Series(values).nunique()),
                    "minimum_p_value": float(np.min(values)),
                    "nominal_0_05": int(np.sum(values < 0.05)),
                    "nominal_0_001": int(np.sum(values < 0.001)),
                })
            for pair_combination in ("simes", "cauchy"):
                pair_tables = [
                    grouped_pvalues(table, ["gene_id", "pair_id"], pair_combination)
                    for table in selected
                ]
                matched_tealeaf = [
                    table.loc[
                        [
                            tuple(value) in shared
                            for value in table[["gene_id", "pair_id"]].to_numpy()
                        ]
                    ]
                    for table in tealeaf
                ]
                tealeaf_pair_tables = [
                    grouped_pvalues(table, ["gene_id", "pair_id"], pair_combination)
                    for table in matched_tealeaf
                ]
                for gene_combination in ("simes", "cauchy"):
                    metric, local_details = metric_row(
                        label,
                        method,
                        pair_combination,
                        gene_combination,
                        pair_tables,
                    )
                    metrics.append(metric)
                    details.append(local_details)
                    reference_metric, reference_details = metric_row(
                        f"Tealeaf matched to {label}",
                        "Tealeaf direct local path",
                        pair_combination,
                        gene_combination,
                        tealeaf_pair_tables,
                    )
                    diagnostics.append({
                        "comparison": reference_metric["comparison"],
                        "method": reference_metric["method"],
                        "fold": "matched",
                        "event_tests": int(sum(len(table) for table in matched_tealeaf)),
                        "gene_pairs": reference_metric["shared_gene_pairs"],
                        "median_events_per_pair": np.nan,
                        "maximum_events_per_pair": np.nan,
                        "unique_p_values": np.nan,
                        "minimum_p_value": np.nan,
                        "nominal_0_05": np.nan,
                        "nominal_0_001": np.nan,
                    })
                    metrics.append(reference_metric)
                    details.append(reference_details)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_csv(args.output_dir / "split_reproducibility.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(args.output_dir / "event_diagnostics.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(details, ignore_index=True).to_csv(args.output_dir / "gene_pvalues.tsv.gz", sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
