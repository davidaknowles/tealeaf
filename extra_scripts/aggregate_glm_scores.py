#!/usr/bin/env python3
"""Combine independently computed GLM label-score summaries."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


def read_fit_names(path):
    with open(path, newline="") as handle:
        return [row[0] for row in csv.reader(handle, delimiter="\t") if row]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits-file", required=True, type=Path)
    parser.add_argument("--score-output", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    fieldnames = []
    for name in read_fit_names(args.fits_file):
        summary = args.score_output / name / "label_score_summary.csv"
        with open(summary, newline="") as handle:
            fit_rows = list(csv.DictReader(handle))
        if len(fit_rows) != 1:
            raise ValueError(f"expected one summary row in {summary}")
        row = fit_rows[0]
        if row.get("name") != name:
            raise ValueError(f"summary name mismatch in {summary}")
        rows.append(row)
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)

    args.score_output.mkdir(parents=True, exist_ok=True)
    output = args.score_output / "label_score_summary.csv"
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    with open(temporary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output)
    print(f"combined {len(rows)} fit summaries in {output}")


if __name__ == "__main__":
    main()
