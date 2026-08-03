#!/usr/bin/env python3
"""Validate Split-seq STARsolo barcode recovery and junction counts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc.junction_benchmark import read_starsolo_sj


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--cell-metadata", type=Path, required=True)
    parser.add_argument("--min-overlap", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    counts, barcodes, junctions = read_starsolo_sj(args.run_dir)
    metadata = pd.read_csv(args.cell_metadata, sep="\t", dtype=str)
    expected = set(metadata.loc[(metadata.run == args.run) | (metadata.run == "*"), "barcode"])
    observed = set(barcodes)
    report = {
        "run": args.run,
        "observed_barcodes": len(observed),
        "expected_barcodes": len(expected),
        "overlapping_barcodes": len(observed & expected),
        "junctions": len(junctions),
        "junction_umis": int(counts.sum()),
        "cells_with_junction_umis": int((counts.getnnz(axis=1) > 0).sum()),
    }
    if report["overlapping_barcodes"] < args.min_overlap:
        raise ValueError(f"insufficient metadata overlap: {report}")
    if report["junction_umis"] <= 0:
        raise ValueError(f"no junction UMIs: {report}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
