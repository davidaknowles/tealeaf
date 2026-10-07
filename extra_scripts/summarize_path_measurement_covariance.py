#!/usr/bin/env python3
"""Collate fixed-truth covariance diagnostics, separating incomplete resampling."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def block_summaries(observations):
    """Summarize each block before aggregating, not each abundant subject/type."""
    metrics = ["fisher_variance_ratio", "scalar_variance_ratio", "scalar_to_fisher_trace_ratio", "bias_squared_norm", "empirical_trace", "mean_squared_error", "effective_depth_median", "minimum_true_usage"]
    return observations.groupby(["strategy", "test_id", "block_id", "gene_id", "n_paths", "complete_bootstrap"]).agg(n_observations=("subject", "size"), **{metric: (metric, "median") for metric in metrics}).reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    paths = sorted(args.cache.glob("shard_*/settings.json"))
    if len(paths) != 16 or {path.parent.name for path in paths} != {f"shard_{index}" for index in range(16)}:
        raise ValueError("expected 16 completed contiguous shards")
    settings = [json.loads(path.read_text()) for path in paths]
    if any(setting != settings[0] for setting in settings[1:]):
        raise ValueError("inconsistent diagnostic settings")
    status = pd.concat([pd.read_csv(path.parent / "status.tsv", sep="\t") for path in paths], ignore_index=True)
    if status.duplicated(["test_id", "strategy"]).any():
        raise ValueError("duplicate block/strategy diagnostics")
    observations = []
    for path in paths:
        try:
            observations.append(pd.read_csv(path.parent / "observations.tsv.gz", sep="\t"))
        except pd.errors.EmptyDataError:
            pass
    if not observations:
        raise ValueError("no observations with estimable bootstrap covariance")
    observations = pd.concat(observations, ignore_index=True)
    blocks = block_summaries(observations)
    summary = []
    for strategy, denominator in status.groupby("strategy"):
        for scope in ("complete bootstrap", "all available, conditional on success"):
            selected = blocks.loc[blocks.strategy.eq(strategy)].copy()
            if scope == "complete bootstrap":
                selected = selected.loc[selected.complete_bootstrap]
            summary.append({"strategy": strategy, "scope": scope, "n_blocks_requested": len(denominator), "n_blocks_available": len(selected), "count_resamples_requested": denominator.bootstrap_requested.sum(), "count_resamples_converged": denominator.bootstrap_converged.sum(), "n_observations": int(selected.n_observations.sum()), **{metric: selected[metric].median() for metric in ("fisher_variance_ratio", "scalar_variance_ratio", "scalar_to_fisher_trace_ratio", "bias_squared_norm", "empirical_trace", "mean_squared_error", "effective_depth_median", "minimum_true_usage")}})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    observations.to_csv(args.output_dir / "observations.tsv.gz", sep="\t", index=False, na_rep="NA")
    status.to_csv(args.output_dir / "status.tsv", sep="\t", index=False)
    blocks.to_csv(args.output_dir / "blocks.tsv", sep="\t", index=False, na_rep="NA")
    summary = pd.DataFrame(summary)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": settings[0], "interpretation": "Measurement covariance diagnostic, not a biological p-value or real-data endpoint", "primary_scope": "All count resamples converge, complete bootstrap only, block-balanced median", "secondary_scope": "Available observations condition on successful quantification, failures visible in status table", "ratio_definition": "Median across subject/type observations within block, then median across blocks", "comparison": "Predicted Fisher and scalar-matched proportion covariance versus repeated actual EC-count variance, with latent transcript mixture fixed", "baseline": "Pooled and subject baselines refitted on each count resample; oracle subject nuisance is diagnostic only", "production_changes": False}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
