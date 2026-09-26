#!/usr/bin/env python3
"""Keep the cell-type contrasts from a primer-aware SUPPA2 result table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.input, sep="\t", compression="infer", low_memory=False)
    table["effect"] = table.contrast_id.astype(str).str.startswith("cell_type__").map({True: "cell_type", False: "condition"})
    table = table.loc[table.effect.eq("cell_type")].copy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(table):,} cell-type tests")


if __name__ == "__main__":
    main()
