#!/usr/bin/env python3
"""Normalize cell labels, conditions, and subjects for junction benchmarks."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def read_mapping(path: Path, name: str) -> pd.DataFrame:
    table = pd.read_csv(path, header=None, names=["cell_id", name], dtype=str)
    if table.cell_id.duplicated().any():
        raise ValueError(f"duplicate cell IDs in {path}")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path)
    parser.add_argument("--conditions", type=Path)
    parser.add_argument("--subjects", type=Path)
    parser.add_argument(
        "--combined-groups",
        type=Path,
        help="two-column BARCODE,CELL_TYPE__CONDITION__SUBJECT table",
    )
    parser.add_argument(
        "--default-run",
        default="*",
        help="run for unprefixed IDs; '*' matches every run",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--combined-output",
        type=Path,
        help="optional two-column cell,cell_type__condition__subject table",
    )
    args = parser.parse_args()

    if args.combined_groups is not None:
        if any(value is not None for value in (args.labels, args.conditions, args.subjects)):
            parser.error("--combined-groups cannot be combined with separate metadata tables")
        rows = []
        with open(args.combined_groups, newline="") as handle:
            for cell_id, value in csv.reader(handle):
                cell_type, condition, subject = value.rsplit("__", 2)
                rows.append((cell_id, cell_type, condition, subject))
        table = pd.DataFrame(rows, columns=["cell_id", "cell_type", "condition", "subject"])
    else:
        if any(value is None for value in (args.labels, args.conditions, args.subjects)):
            parser.error("separate metadata requires --labels, --conditions, and --subjects")
        table = read_mapping(args.labels, "cell_type")
        table = table.merge(read_mapping(args.conditions, "condition"), on="cell_id")
        table = table.merge(read_mapping(args.subjects, "subject"), on="cell_id")

    split = table.cell_id.str.partition(":")
    has_run = split[1].eq(":")
    table["run"] = split[0].where(has_run, args.default_run)
    table["barcode"] = split[2].where(has_run, split[0])
    output = table[["run", "barcode", "cell_type", "condition", "subject"]]
    if output[["run", "barcode"]].duplicated().any():
        raise ValueError("normalized run/barcode keys are not unique")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    if args.combined_output is not None:
        args.combined_output.parent.mkdir(parents=True, exist_ok=True)
        cell_id = output["barcode"].where(
            output["run"].eq("*"), output["run"] + ":" + output["barcode"]
        )
        group = output["cell_type"] + "__" + output["condition"] + "__" + output["subject"]
        pd.DataFrame({"cell_id": cell_id, "group": group}).to_csv(
            args.combined_output, header=False, index=False
        )
    print(
        f"cells={len(output)} cell_types={output.cell_type.nunique()} "
        f"conditions={output.condition.nunique()} subjects={output.subject.nunique()}"
    )


if __name__ == "__main__":
    main()
