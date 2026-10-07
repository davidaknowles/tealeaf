#!/usr/bin/env python3
"""Calibrate and assess EC count nulls without treating label permutations as truth."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.merge_paired_path_test import add_calibration_strata, empirical_null_calibration, moderate_scalar_tests
from tealeaf.sc.ds_benchmark import benjamini_hochberg


def validate_requested_trials(table, settings):
    """Preserve and validate a declared nuisance-stress trial family."""
    if "requested_ids" not in settings:
        return table
    keys = ["test_id", "draw", "strategy"]
    expected = {(test_id, draw, strategy) for test_id in settings["requested_ids"] for draw in range(settings["draws"]) for strategy in settings["expected_strategies"]}
    if table.duplicated(keys).any() or set(map(tuple, table[keys].values)) != expected:
        raise ValueError("missing or duplicate requested nuisance-stress trials")
    table = table.copy()
    failed = ~table.converged.astype(str).str.lower().eq("true")
    table.loc[failed, "p_value"] = 1.
    if not np.isfinite(table.p_value).all() or not table.p_value.between(0, 1).all():
        raise ValueError("invalid successful nuisance-stress p-value")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--moderate-variances", action="store_true")
    args = parser.parse_args()
    paths = sorted(args.cache.glob("shard_*/observed.tsv"))
    if len(paths) != 16:
        raise ValueError(f"expected 16 shards, found {len(paths)}")
    settings = [json.loads((path.parent / "settings.json").read_text()) for path in paths]
    if any(item != settings[0] for item in settings[1:]):
        raise ValueError("count-null settings differ between shards")
    observed = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    null_tables = []
    for path in paths:
        try:
            null_tables.append(pd.read_csv(path.parent / "null.tsv.gz", sep="\t"))
        except pd.errors.EmptyDataError:
            if "requested_ids" not in settings[0]:
                raise
    if not null_tables:
        raise ValueError("no successful null-training fits available")
    null = pd.concat(null_tables, ignore_index=True)
    observed = validate_requested_trials(observed, settings[0])
    observed = observed.loc[np.isfinite(observed.p_value)].copy()
    null = null.loc[np.isfinite(null.p_value)].copy()
    results, summaries = [], []
    for (draw, strategy), local in observed.groupby(["draw", "strategy"]):
        training = null.loc[null.draw.eq(draw) & null.strategy.eq(strategy)].copy()
        if args.moderate_variances:
            local, training = moderate_scalar_tests(local, training)
        table = add_calibration_strata(local.copy().reset_index(drop=True), 100)
        if strategy.startswith("maximum-coordinate"):
            table["calibration_stratum"] += "|paths=" + table.n_paths.astype(str)
        calibrated, _ = empirical_null_calibration(table, training)
        if "requested_ids" in settings[0]:
            failed = ~calibrated.converged.astype(str).str.lower().eq("true")
            calibrated.loc[failed, ["p_value", "raw_p_value"]] = 1.
        calibrated["draw"] = draw
        results.append(calibrated)
        summaries.append({"draw": draw, "strategy": strategy, "n_tests": len(calibrated), "n_converged": int(calibrated.converged.sum()), "raw_reject_0_05": calibrated.raw_p_value.le(.05).mean(), "calibrated_reject_0_05": calibrated.p_value.le(.05).mean(), "calibrated_reject_0_01": calibrated.p_value.le(.01).mean(), "BH_discoveries": int(np.sum(benjamini_hochberg(calibrated.p_value.to_numpy(float)) <= .05))})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tests = pd.concat(results, ignore_index=True)
    tests.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output_dir / "draw_summary.tsv", sep="\t", index=False, na_rep="NA")
    aggregate = tests.groupby("strategy", observed=True).agg(n_tests=("p_value", "size"), n_blocks=("block_id", "nunique"), n_converged=("converged", "sum"), raw_reject_0_05=("raw_p_value", lambda p: p.le(.05).mean()), calibrated_reject_0_05=("p_value", lambda p: p.le(.05).mean()), calibrated_reject_0_01=("p_value", lambda p: p.le(.01).mean())).reset_index()
    aggregate.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    if "converged" in null:
        flagged = null.copy()
        flagged["fit_available"] = flagged.converged.astype(str).str.lower().eq("true")
        flagged.groupby("strategy").agg(n_null_trials=("p_value", "size"), n_converged=("fit_available", "sum"), failed_fraction=("fit_available", lambda values: (~values).mean())).reset_index().to_csv(args.output_dir / "null_fit_summary.tsv", sep="\t", index=False)
    if "alternative_concentration" in tests:
        fitted = tests.loc[tests.converged.astype(str).str.lower().eq("true")].copy()
        fitted["grid_boundary"] = fitted.profile_boundary.astype(str).str.lower().eq("true")
        fitted.groupby("strategy").agg(n_fitted=("p_value", "size"), median_precision=("alternative_concentration", "median"), grid_boundary_fraction=("grid_boundary", "mean"), median_effective_depth=("effective_depth_median", "median")).reset_index().to_csv(args.output_dir / "fit_precision_summary.tsv", sep="\t", index=False)
    # Overall nominal calibration can hide failure in individual dimensions or
    # repeated rejection of the same null hypothesis across count draws.
    diagnostics = {
        "n_tests": ("p_value", "size"),
        "n_converged": ("converged", "sum"),
        "raw_reject_0_05": ("raw_p_value", lambda p: p.le(.05).mean()),
        "calibrated_reject_0_05": ("p_value", lambda p: p.le(.05).mean()),
        "calibrated_reject_0_01": ("p_value", lambda p: p.le(.01).mean()),
    }
    tests.groupby(["strategy", "degrees_of_freedom", "n_subjects"]).agg(**diagnostics).reset_index().to_csv(args.output_dir / "stratum_summary.tsv", sep="\t", index=False, na_rep="NA")
    tests.groupby(["strategy", "test_id", "block_id"]).agg(**diagnostics).reset_index().to_csv(args.output_dir / "hypothesis_summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": settings[0], "null": settings[0]["null"], "assessment": "True EC count nulls, distinct from label-permutation nulls", "calibration": "32 within-subject permutations/sign flips per simulated test and draw; own test excluded from pooled empirical calibration", "variance_moderation": args.moderate_variances, "limitation": "Fixed compatibility, no biological cell-type heteroscedasticity, no EB reselection; block-wise simulations do not preserve a coherent joint-gene null", "interpretation": "Reject fractions diagnose conditional testing calibration; per-draw BH counts do not certify joint real-data FDR"}, indent=2) + "\n")
    print(aggregate.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
