#!/usr/bin/env python3
"""Compare native SUPPA2 and Tealeaf event directions on the same Tilgner audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-tests", required=True, type=Path)
    parser.add_argument("--hybrid-tests", required=True, type=Path)
    parser.add_argument("--tilgner-audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-rank", type=int, default=200)
    parser.add_argument("--minimum-depth", type=float, default=20.0)
    return parser.parse_args()


def load_tests(path):
    table = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    table = table.loc[
        table.effect.eq("cell_type") & table.p_value.notna(),
        ["contrast_id", "feature_id", "p_value", "effect_size"],
    ].copy()
    table["p_value"] = pd.to_numeric(table.p_value, errors="coerce")
    table["effect_size"] = pd.to_numeric(table.effect_size, errors="coerce")
    table = table.dropna(subset=["p_value", "effect_size"])
    return table.drop_duplicates(["contrast_id", "feature_id"], keep="first")


def main():
    args = parse_args()
    native = load_tests(args.native_tests).rename(columns={
        "p_value": "native_p_value",
        "effect_size": "native_effect",
    })
    hybrid = load_tests(args.hybrid_tests).rename(columns={
        "p_value": "hybrid_p_value",
        "effect_size": "hybrid_effect",
    })
    audit = pd.read_csv(
        args.tilgner_audit, sep="\t", compression="infer", low_memory=False
    )
    audit = audit.loc[
        audit.method.eq("SUPPA2 (full data)")
        & audit.mapping_complete.astype(str).str.lower().eq("true")
        & audit.minimum_pooled_depth.ge(args.minimum_depth)
    ].copy()
    keys = ["contrast_id", "feature_id"]
    joined = audit.merge(native, on=keys, how="inner", validate="one_to_one")
    joined = joined.merge(hybrid, on=keys, how="inner", validate="one_to_one")
    native_sign = np.sign(joined.native_effect.to_numpy(float))
    joined["long_read_pooled_sign"] = native_sign * np.where(
        joined.pooled_replicated.astype(bool), 1.0, -1.0
    )
    joined["long_read_rep1_sign"] = native_sign * np.sign(
        joined.replicate_1_dot_product.to_numpy(float)
    )
    joined["long_read_rep2_sign"] = native_sign * np.sign(
        joined.replicate_2_dot_product.to_numpy(float)
    )
    joined = joined.loc[joined.long_read_pooled_sign.ne(0)].copy()
    joined["strict_eligible"] = (
        joined.minimum_replicate_depth.ge(args.minimum_depth / 2)
        & joined.long_read_rep1_sign.ne(0)
        & joined.long_read_rep2_sign.ne(0)
    )

    ranked = []
    for method, p_column, effect_column in (
        ("SUPPA2 (full data)", "native_p_value", "native_effect"),
        ("Tealeaf EC; SUPPA2 event definitions", "hybrid_p_value", "hybrid_effect"),
    ):
        local = joined.copy()
        local["method"] = method
        local["p_value"] = local[p_column]
        local["effect"] = np.sign(local[effect_column].to_numpy(float))
        local["pooled_agreement"] = (
            local.effect.to_numpy(float) == local.long_read_pooled_sign.to_numpy(float)
        )
        local["both_replicates_agreement"] = (
            (local.effect.to_numpy(float) == local.long_read_rep1_sign.to_numpy(float))
            & (local.effect.to_numpy(float) == local.long_read_rep2_sign.to_numpy(float))
        )
        local = local.sort_values(
            ["p_value", "contrast_id", "feature_id"], kind="stable"
        ).reset_index(drop=True)
        local["rank"] = np.arange(1, len(local) + 1)
        local["pooled_cumulative_agreement"] = (
            local.pooled_agreement.cumsum() / local["rank"]
        )
        strict_values = local.loc[local.strict_eligible, "both_replicates_agreement"]
        strict_counts = np.zeros(len(local), dtype=float)
        strict_denominators = np.zeros(len(local), dtype=float)
        eligible_mask = local.strict_eligible.to_numpy()
        strict_counts[eligible_mask] = (
            local.loc[eligible_mask, "both_replicates_agreement"].astype(int).cumsum()
        )
        strict_denominators[eligible_mask] = np.arange(1, len(strict_values) + 1)
        local["strict_cumulative_agreement"] = np.divide(
            strict_counts,
            strict_denominators,
            out=np.full(len(local), np.nan),
            where=strict_denominators > 0,
        )
        ranked.append(local)
    result = pd.concat(ranked, ignore_index=True)
    summary = []
    for method, group in result.groupby("method", sort=False):
        top = group.head(args.top_rank)
        strict = top.loc[top.strict_eligible]
        summary.append({
            "method": method,
            "n_evaluable": len(group),
            "top_rank": min(args.top_rank, len(top)),
            "pooled_n": len(top),
            "pooled_agree": int(top.pooled_agreement.sum()),
            "pooled_rate": float(top.pooled_agreement.mean()) if len(top) else np.nan,
            "strict_n": len(strict),
            "strict_agree": int(strict.both_replicates_agreement.sum()),
            "strict_rate": float(strict.both_replicates_agreement.mean()) if len(strict) else np.nan,
        })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "rank_agreement_by_rank.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(summary).to_csv(args.output_dir / "rank_agreement_summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
