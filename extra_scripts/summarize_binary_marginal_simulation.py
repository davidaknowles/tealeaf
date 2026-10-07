#!/usr/bin/env python3
"""Collate exact-count random-subject controls without dropping failed trials."""

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
    expected = settings[0]["shard_count"]
    if len(files) != expected or {item["shard_index"] for item in settings} != set(range(expected)):
        raise ValueError("incomplete shard family")
    shared = [{key: value for key, value in item.items() if key not in ("output_dir", "shard_index")} for item in settings]
    if any(item != shared[0] for item in shared[1:]):
        raise ValueError("inconsistent simulation settings")
    data = pd.concat([pd.read_csv(path.parent / "tests.tsv.gz", sep="\t") for path in files], ignore_index=True).sort_values("draw")
    if data.draw.duplicated().any() or set(data.draw) != set(range(settings[0]["draws"])):
        raise ValueError("duplicate or missing simulated trials")
    good = data.converged.astype(str).str.lower().eq("true")
    data.loc[~good, "p_value"] = 1.
    summary = {"n_trials": len(data), "n_converged": int(good.sum()), "reject_0_05": data.p_value.le(.05).mean(), "reject_0_01": data.p_value.le(.01).mean(), "null_precision_median": data.loc[good, "null_concentration"].median(), "alternative_precision_median": data.loc[good, "alternative_concentration"].median(), "null_subject_sd_median": data.loc[good, "null_subject_sd"].median(), "alternative_subject_sd_median": data.loc[good, "alternative_subject_sd"].median(), "n_precision_upper_boundary": int(data.loc[good, "alternative_concentration"].ge(9999).sum()), "runtime_median_seconds": data.runtime_seconds.median(), "max_quadrature_error_converged": data.loc[good, "quadrature_error"].max()}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.DataFrame([summary]).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": shared[0], "failure_policy": "All requested trials retained, nonconverged or inaccurate-quadrature fits p=1", "comparison": "Synthetic control only, no claim about actual-EC calibration or split/LR performance", "production_changes": False}, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
