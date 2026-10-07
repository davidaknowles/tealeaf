#!/usr/bin/env python3
"""Separate fitted-test and pooled-calibration losses on completed split audits."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.summarize_path_reporting_omnibus import split_omnibus_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(args.assessment / f"omnibus_{fold}_tests.tsv.gz", sep="\t") for fold in (0, 1)]
    native = [frame.assign(p_value=frame.raw_p_value) for frame in frames]
    summary = split_omnibus_summary(native, args.output_dir, prefix="native_split")
    diagnostics = []
    for fold, frame in enumerate(frames):
        for strategy, selected in frame.loc[frame.strategy.str.startswith("joint DM")].groupby("strategy"):
            fitted = selected.loc[selected.fit_available]
            diagnostics.append({"fold": fold, "strategy": strategy, "tests": len(selected), "fitted": len(fitted), "native_nominal": int(selected.raw_p_value.le(.05).sum()), "calibrated_nominal": int(selected.p_value.le(.05).sum()), "median_effective_depth": fitted.effective_depth_median.median(), "multinomial_null_fraction": np.isinf(fitted.null_concentration).mean(), "multinomial_alternative_fraction": np.isinf(fitted.alternative_concentration).mean(), "finite_alternative_concentration_median": fitted.loc[np.isfinite(fitted.alternative_concentration), "alternative_concentration"].median()})
    pd.DataFrame(diagnostics).to_csv(args.output_dir / "fit_diagnostics.tsv", sep="\t", index=False)
    print(summary.loc[summary.strategy.str.startswith("joint DM")].to_string(index=False), flush=True)
    print(pd.DataFrame(diagnostics).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
