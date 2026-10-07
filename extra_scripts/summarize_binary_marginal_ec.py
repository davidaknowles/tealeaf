#!/usr/bin/env python3
"""Collate binary direct-EC accuracy/runtime pilots, not endpoint benchmarks."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    files = sorted(args.cache.glob("shard_*/settings.json"))
    if not files:
        raise ValueError("no completed shards")
    settings = [json.loads(path.read_text()) for path in files]
    expected = settings[0].get("shard_count", 4)
    if len(files) != expected or {path.parent.name for path in files} != {f"shard_{i}" for i in range(expected)}:
        raise ValueError("expected complete contiguous pilot shards")
    if any(item != settings[0] for item in settings[1:]):
        raise ValueError("inconsistent pilot settings")
    table = pd.concat([pd.read_csv(path.parent / "tests.tsv.gz", sep="\t") for path in files], ignore_index=True)
    scenarios = settings[0].get("requested_scenarios", [{"scenario": "observed", "draw": -1}, {"scenario": "common composition", "draw": 0}, {"scenario": "biological precision 20", "draw": 0}])
    if "draw" not in table:
        table["draw"] = np.where(table.scenario.eq("observed"), -1, 0)
    keys = {(test_id, item["scenario"], item["draw"]) for test_id in settings[0]["selected_ids"] for item in scenarios}
    if table.duplicated(["test_id", "scenario", "draw"]).any() or set(zip(table.test_id, table.scenario, table.draw)) != keys:
        raise ValueError("missing or duplicate selected pilot cases")
    good = table.converged.astype(str).str.lower().eq("true")
    table["converged"] = good
    table.loc[~good, "p_value"] = 1.
    # Invalid fitted effects cannot enter a future assessment by accident.
    if "standardized_means" in table:
        table.loc[~good, "standardized_means"] = np.nan
    for column in ("runtime_seconds", "quadrature_error"):
        if column not in table:
            table[column] = np.nan
    summary = table.groupby("scenario").agg(n_requested=("p_value", "size"), n_converged=("converged", "sum"), native_reject_0_05=("p_value", lambda values: values.le(.05).sum()), native_reject_rate_0_05=("p_value", lambda values: values.le(.05).mean()), median_runtime_seconds=("runtime_seconds", "median"), max_quadrature_error=("quadrature_error", "max")).reset_index()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": settings[0], "scope": "Independently selected binary runtime/accuracy or count-null panel, not a full-dimensional calibration or split/LR assessment", "selection": "covered production-eligible tests, seed311381, without discovery or LR outcome filtering", "failure_policy": "All requested cases retained, failed p=1 and direction unavailable", "production_changes": False}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
