#!/usr/bin/env python3
"""Convert a metadata table into tealeaf label and group tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def mapping(value):
    source, separator, target = value.partition("=")
    if not separator or not source or not target:
        raise argparse.ArgumentTypeError("batch maps must have SOURCE=TARGET form")
    return source, target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--cell-barcode-column", default="cell_barcode")
    parser.add_argument("--batch-column", default="Batch")
    parser.add_argument("--label-column", default="annotation")
    parser.add_argument("--group-column", default="CaseNum")
    parser.add_argument("--batch-map", action="append", required=True, type=mapping)
    parser.add_argument("--labels-output", required=True, type=Path)
    parser.add_argument("--groups-output", required=True, type=Path)
    args = parser.parse_args()
    published_to_internal = {
        published: internal for internal, published in args.batch_map
    }
    if len(published_to_internal) != len(args.batch_map):
        raise ValueError("published batch names must be unique")

    labels = {}
    groups = {}
    with open(args.metadata, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            args.cell_barcode_column,
            args.batch_column,
            args.label_column,
            args.group_column,
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"metadata is missing columns: {sorted(missing)}")
        for row in reader:
            published_batch = row[args.batch_column]
            internal_batch = published_to_internal.get(published_batch)
            if internal_batch is None:
                continue
            barcode = row[args.cell_barcode_column]
            label = row[args.label_column]
            group = row[args.group_column]
            if not barcode or not label:
                continue
            cell_id = f"{internal_batch}:{barcode}"
            previous = labels.setdefault(cell_id, label)
            if previous != label:
                raise ValueError(f"conflicting labels for {cell_id}")
            if group:
                previous_group = groups.setdefault(cell_id, group)
                if previous_group != group:
                    raise ValueError(f"conflicting groups for {cell_id}")
    if not labels:
        raise ValueError("no metadata rows matched the batch mapping")

    args.labels_output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.labels_output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(sorted(labels.items()))
    args.groups_output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.groups_output, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(sorted(groups.items()))
    print(f"labels={len(labels)} groups={len(groups)}")


if __name__ == "__main__":
    main()
