#!/usr/bin/env python3
"""Validate merged alevin matrices, metadata overlap, and weighted sidecar."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from tealeaf.data.alevin import validate_alevin_quantification


def read_first_column(path, delimiter=","):
    with open(path, newline="") as handle:
        return [row[0] for row in csv.reader(handle, delimiter=delimiter) if row]


def read_manifest_runs(path):
    with open(path, newline="") as handle:
        return {
            row["run_accession"]
            for row in csv.DictReader(handle, delimiter="\t")
            if int(row["read"]) == 1
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--reference-labels", type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--min-cell-umis", type=float, default=500)
    parser.add_argument("--min-reference-overlap", type=int, default=1)
    parser.add_argument("--max-total-molecules", type=float)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    report = validate_alevin_quantification(
        args.alevin_dir,
        expected_prefixes=(
            read_manifest_runs(args.manifest) if args.manifest is not None else None
        ),
        reference_ids=(
            read_first_column(args.reference_labels)
            if args.reference_labels is not None else None
        ),
        primer_pair_file=args.primer_pairs,
        min_cell_umis=args.min_cell_umis,
        min_reference_overlap=args.min_reference_overlap,
        max_total_molecules=args.max_total_molecules,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
