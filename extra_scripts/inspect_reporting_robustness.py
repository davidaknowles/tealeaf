#!/usr/bin/env python3
"""Separate tiny reporting directions from reproducible effect magnitudes."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    table = pd.read_csv(args.input, sep="\t", low_memory=False)
    table = table.loc[table.event_BH_union & table.direction_agrees.notna()].copy()
    rows = []
    for (comparison, strategy), local in table.groupby(["comparison", "strategy"]):
        for discovery_fold in (0, 1):
            norm = local[f"norm_{discovery_fold}"]
            for cutoff in (0., np.sqrt(2) * .01, np.sqrt(2) * .05, np.sqrt(2) * .10):
                selected = local.loc[norm >= cutoff]
                rows.append({"comparison": comparison, "strategy": strategy, "discovery_fold": discovery_fold, "minimum_discovery_effect_L2": cutoff, "n_selected": len(selected), "agreement": selected.direction_agrees.mean(), "median_discovery_effect": selected[f"norm_{discovery_fold}"].median(), "median_held_effect": selected[f"norm_{1 - discovery_fold}"].median()})
        rows.append({"comparison": comparison, "strategy": strategy, "discovery_fold": "both, diagnostic only", "minimum_discovery_effect_L2": 1e-12, "n_selected": len(local), "agreement": local.direction_agrees.mean(), "near_zero_either_fraction": (local.norm_0.le(1e-12) | local.norm_1.le(1e-12)).mean()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
