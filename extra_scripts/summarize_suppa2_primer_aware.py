#!/usr/bin/env python3
"""Split primer-aware SUPPA2 tests and summarize fold reproducibility."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from tealeaf.sc.ds_benchmark import benjamini_hochberg, simes_pvalue


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def split_table(table, fold):
    return table.loc[table.fold.astype(str).eq(str(fold))].copy()


def aggregate(table):
    local = table.loc[table.gene_id.notna() & table.p_value.notna()].copy()
    local["gene_id"] = local.gene_id.astype(str).str.split(".").str[0]
    return (
        local.groupby(["gene_id", "contrast_id"], sort=False)["p_value"]
        .agg(p_value=simes_pvalue, n_events="size")
        .reset_index()
    )


def summarize_split(fold0, fold1):
    first, second = aggregate(fold0), aggregate(fold1)
    keys = sorted(set(zip(first.gene_id, first.contrast_id)) & set(zip(second.gene_id, second.contrast_id)))
    first = first.set_index(["gene_id", "contrast_id"]).loc[keys]
    second = second.set_index(["gene_id", "contrast_id"]).loc[keys]
    p0, p1 = first.p_value.to_numpy(float), second.p_value.to_numpy(float)
    q0, q1 = benjamini_hochberg(p0), benjamini_hochberg(p1)
    conjunction_q = benjamini_hochberg(np.maximum(p0, p1))
    rep0 = q0 <= 0.05
    rep1 = q1 <= 0.05
    held = []
    if rep0.any():
        held.append(float(np.mean(p1[rep0] <= 0.05)))
    if rep1.any():
        held.append(float(np.mean(p0[rep1] <= 0.05)))
    rho = float(spearmanr(-np.log10(np.maximum(p0, 1e-300)), -np.log10(np.maximum(p1, 1e-300))).statistic) if len(keys) > 1 else np.nan
    return pd.DataFrame([{
        "method": "SUPPA2 (primer aware)",
        "shared_gene_contrasts": len(keys),
        "fold0_bh": int(rep0.sum()),
        "fold1_bh": int(rep1.sum()),
        "replicated_bh": int((conjunction_q <= 0.05).sum()),
        "heldout_nominal_replication": float(np.mean(held)) if held else np.nan,
        "spearman_logp": rho,
    }]), pd.DataFrame({
        "gene_id": [key[0] for key in keys],
        "contrast_id": [key[1] for key in keys],
        "fold0_p_value": p0,
        "fold1_p_value": p1,
        "fold0_q_value": q0,
        "fold1_q_value": q1,
        "conjunction_q_value": conjunction_q,
    })


def main():
    args = parse_args()
    table = pd.read_csv(args.tests, sep="\t", compression="infer", low_memory=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    full = split_table(table, "full_data")
    fold0 = split_table(table, 0)
    fold1 = split_table(table, 1)
    for name, value in (("full_data_tests.tsv.gz", full), ("split_data_fold0_tests.tsv.gz", fold0), ("split_data_fold1_tests.tsv.gz", fold1)):
        value.to_csv(args.output_dir / name, sep="\t", index=False, compression="gzip")
    metrics, details = summarize_split(fold0, fold1)
    metrics.to_csv(args.output_dir / "split_data_metrics.tsv", sep="\t", index=False)
    details.to_csv(args.output_dir / "split_data_gene_contrasts.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame([{
        "method": "SUPPA2 (primer aware)",
        "tests": len(full),
        "events": full.feature_id.nunique(),
        "contrasts": full.contrast_id.nunique(),
        "bh_events": int((full.q_value <= 0.05).sum()),
    }]).to_csv(args.output_dir / "full_data_summary.tsv", sep="\t", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
