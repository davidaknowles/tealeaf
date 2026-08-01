#!/usr/bin/env python3
"""Merge EC-GLMM simulation shards and compute null-calibrated power."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    table = pd.concat([pd.read_csv(path, sep="\t") for path in args.inputs])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "replicates.tsv", sep="\t", index=False)
    summary = []
    for (family, method), group in table.groupby(["family", "method"]):
        null = group.loc[group["effect"] == 0, "evidence_gain"].dropna()
        threshold = float(np.quantile(null, 0.95)) if len(null) else np.nan
        for effect, records in group.groupby("effect"):
            errors = records["estimate"] - records["effect"]
            ess = records["alternative_importance_ess"].dropna()
            summary.append({
                "family": family,
                "method": method,
                "effect": effect,
                "replicates": len(records),
                "null_95_threshold": threshold,
                "calibrated_rejection": float(
                    np.mean(records["evidence_gain"] > threshold)
                ),
                "mean_bias": float(errors.mean()),
                "rmse": float(np.sqrt(np.mean(np.square(errors)))),
                "convergence": float(np.mean(
                    records["null_converged"] & records["alternative_converged"]
                )),
                "median_seconds": float(records["seconds"].median()),
                "median_alternative_ess": float(
                    ess.median()
                ) if len(ess) else np.nan,
                "median_alternative_observation_noise_sd": float(
                    records.get(
                        "alternative_observation_noise_sd",
                        pd.Series(0.0, index=records.index),
                    ).median()
                ),
            })
    summary = pd.DataFrame(summary)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary.to_dict(orient="records"), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
