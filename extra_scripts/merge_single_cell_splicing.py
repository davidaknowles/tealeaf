#!/usr/bin/env python3
"""Merge sharded single-cell splicing results and recompute global FDR."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = []
    summaries = []
    estimates_path = args.output_dir / "single_cell_estimates.jsonl.gz"
    with gzip.open(estimates_path, "wt") as output:
        for shard in args.shards:
            with gzip.open(
                shard / "single_cell_estimates.jsonl.gz",
                "rt",
            ) as source:
                for line in source:
                    output.write(line)
            table_path = shard / "differential_splicing.tsv"
            if table_path.stat().st_size:
                tables.append(pd.read_csv(table_path, sep="\t"))
            summaries.append(
                json.loads((shard / "validation_summary.json").read_text())
            )

    table = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    if len(table):
        table["fdr"] = pd.NA
        for method, positions in table.groupby("covariance_method").groups.items():
            table.loc[positions, "fdr"] = benjamini_hochberg(
                table.loc[positions, "p_value"].to_numpy()
            )
        table["fdr"] = table["fdr"].astype(float)
        table = table.sort_values("p_value")
    table.to_csv(
        args.output_dir / "differential_splicing.tsv",
        sep="\t",
        index=False,
    )

    permutation_tests = sum(
        row.get("single_cell_permutation_tests", 0) for row in summaries
    )
    permutation_small = sum(
        row.get("single_cell_permutation_tests", 0)
        * row.get("single_cell_permutation_p_lt_0.05", 0.0)
        for row in summaries
    )
    summary = {
        "seconds_max_shard": max(row["seconds"] for row in summaries),
        "cells": max(row["cells"] for row in summaries),
        "blocks": sum(row["blocks"] for row in summaries),
        "estimates": sum(row["estimates"] for row in summaries),
        "converged": sum(row["converged"] for row in summaries),
        "reliable": sum(row["reliable"] for row in summaries),
        "tests": int(len(table)),
        "fdr_0.05": (
            int((table["fdr"] <= 0.05).sum()) if len(table) else 0
        ),
        "single_cell_permutation_tests": permutation_tests,
        "single_cell_permutation_p_lt_0.05": (
            permutation_small / permutation_tests
            if permutation_tests
            else None
        ),
    }
    summary["converged_fraction"] = (
        summary["converged"] / max(summary["estimates"], 1)
    )
    summary["reliable_fraction"] = (
        summary["reliable"] / max(summary["estimates"], 1)
    )
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
