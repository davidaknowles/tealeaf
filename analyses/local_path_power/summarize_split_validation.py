#!/usr/bin/env python3
"""Summarize exact-event split replication and paired-label null calibration."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold0", required=True, type=Path)
    parser.add_argument("--fold1", required=True, type=Path)
    parser.add_argument("--null", type=Path)
    parser.add_argument("--effect-fold0", type=Path)
    parser.add_argument("--effect-fold1", type=Path)
    parser.add_argument("--candidate-cache", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def read_results(path):
    table = pd.read_csv(path, sep="\t")
    converged = table["null_converged"] & table["alternative_converged"]
    return table.loc[converged].drop_duplicates("test_id", keep="last")


def main():
    args = parse_args()
    fold0 = read_results(args.fold0)
    fold1 = read_results(args.fold1)
    test_ids = None
    if args.candidate_cache is not None:
        with args.candidate_cache.open("rb") as handle:
            candidates = pickle.load(handle)["candidates"]
        test_ids = {candidate[0] for candidate in candidates}
        fold0 = fold0.loc[fold0["test_id"].isin(test_ids)]
        fold1 = fold1.loc[fold1["test_id"].isin(test_ids)]
    common = fold0[["test_id", "lrt_p_value"]].merge(
        fold1[["test_id", "lrt_p_value"]], on="test_id", suffixes=("_0", "_1")
    )
    common["conjunction_p_value"] = common[
        ["lrt_p_value_0", "lrt_p_value_1"]
    ].max(axis=1)
    common["conjunction_fdr"] = benjamini_hochberg(
        common["conjunction_p_value"].to_numpy()
    )
    both_nominal = common["conjunction_p_value"] <= 0.05
    discoveries = common["conjunction_fdr"] <= 0.05
    summary = {
        "fold0_tests": len(fold0),
        "fold1_tests": len(fold1),
        "exact_common_tests": len(common),
        "both_halves_nominal_0.05": int(both_nominal.sum()),
        "conjunction_bh_0.05": int(discoveries.sum()),
        "split_p_value_spearman": float(
            spearmanr(
                -np.log10(common["lrt_p_value_0"].clip(lower=1e-300)),
                -np.log10(common["lrt_p_value_1"].clip(lower=1e-300)),
            ).statistic
        ),
    }
    if args.null is not None:
        null = read_results(args.null)
        if test_ids is not None:
            null = null.loc[null["test_id"].isin(test_ids)]
        summary.update(
            {
                "null_converged_tests": len(null),
                "null_rate_0.05": float((null["lrt_p_value"] <= 0.05).mean()),
                "null_rate_0.01": float((null["lrt_p_value"] <= 0.01).mean()),
                "null_rate_0.001": float((null["lrt_p_value"] <= 0.001).mean()),
                "null_bh_0.05": int(
                    (
                        benjamini_hochberg(null["lrt_p_value"].to_numpy())
                        <= 0.05
                    ).sum()
                ),
            }
        )
    if args.effect_fold0 is not None and args.effect_fold1 is not None:
        effect0 = read_results(args.effect_fold0)[["test_id", "tested_effect"]]
        effect1 = read_results(args.effect_fold1)[["test_id", "tested_effect"]]
        effects = common.merge(effect0, on="test_id").merge(
            effect1, on="test_id", suffixes=("_0", "_1")
        )
        evaluated = effects.loc[effects["conjunction_p_value"] <= 0.05]
        selected = effects.loc[effects["conjunction_fdr"] <= 0.05]
        for label, table in (("nominal", evaluated), ("bh", selected)):
            summary[f"{label}_effects_available"] = len(table)
            summary[f"{label}_same_sign"] = int(
                (np.sign(table["tested_effect_0"]) == np.sign(table["tested_effect_1"])).sum()
            )
            summary[f"{label}_effect_spearman"] = (
                float(spearmanr(table["tested_effect_0"], table["tested_effect_1"]).statistic)
                if len(table) >= 2
                else np.nan
            )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
