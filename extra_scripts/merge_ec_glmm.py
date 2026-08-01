#!/usr/bin/env python3
"""Merge genome-wide EC-GLMM shards and calibrate gains by simulation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def benjamini_hochberg(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    adjusted = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--simulation-replicates", type=Path)
    args = parser.parse_args()
    tables = []
    failures = []
    for directory in args.inputs:
        path = directory / "ec_glmm_fits.tsv"
        if path.is_file() and path.stat().st_size:
            tables.append(pd.read_csv(path, sep="\t"))
        failure_path = directory / "failures.json"
        if failure_path.is_file():
            failures.extend(json.loads(failure_path.read_text()))
    if not tables:
        raise ValueError("no EC GLMM fit tables found")
    table = pd.concat(tables, ignore_index=True)
    table = table.drop_duplicates(["gene", "method", "contrast"], keep="last")
    if args.simulation_replicates is not None:
        simulation = pd.read_csv(args.simulation_replicates, sep="\t")
        for method, positions in table.groupby("method").groups.items():
            positions = np.asarray(list(positions))
            calibration_mask = (
                (table.loc[positions, "contrast"].to_numpy() == "condition")
                & (table.loc[positions, "n_isoforms"].to_numpy() == 2)
            )
            positions = positions[calibration_mask]
            if not len(positions):
                continue
            null = simulation.loc[
                (simulation["method"] == method) & (simulation["effect"] == 0),
                "evidence_gain",
            ].dropna().to_numpy()
            if not len(null):
                continue
            observed = table.loc[positions, "evidence_gain"].to_numpy()
            p_values = np.asarray([
                (1 + np.sum(null >= value)) / (len(null) + 1)
                for value in observed
            ])
            table.loc[positions, "simulation_p_value"] = p_values
            table.loc[positions, "simulation_fdr"] = benjamini_hochberg(p_values)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "ec_glmm_fits.tsv", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    summary = []
    for (contrast, method), records in table.groupby(["contrast", "method"]):
        summary.append({
            "contrast": contrast,
            "method": method,
            "genes": int(len(records)),
            "convergence": float(np.mean(
                records["null_converged"] & records["alternative_converged"]
            )),
            "median_evidence_gain": float(records["evidence_gain"].median()),
            "positive_evidence_gain": float(np.mean(records["evidence_gain"] > 0)),
            "median_tested_coefficient_norm": float(
                records["tested_coefficient_norm"].median()
            ),
            "median_seconds": float(records["seconds"].median()),
            "simulation_nominal_0.05": int(
                (
                    records.get(
                        "simulation_p_value", pd.Series(dtype=float)
                    )
                    <= 0.05
                ).sum()
            ),
            "simulation_fdr_0.05": int(
                (records.get("simulation_fdr", pd.Series(dtype=float)) <= 0.05).sum()
            ),
        })
    summary = pd.DataFrame(summary)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary.to_dict(orient="records"), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
