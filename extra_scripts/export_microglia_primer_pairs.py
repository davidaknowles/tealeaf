#!/usr/bin/env python3
"""Export paired primer barcodes and reference labels from microglia AnnData."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


def decode_column(obs, name):
    value = obs[name]
    if isinstance(value, h5py.Dataset):
        array = value[:]
    else:
        codes = value["codes"][:]
        categories = value["categories"][:]
        array = np.empty(len(codes), dtype=object)
        valid = codes >= 0
        array[~valid] = None
        array[valid] = categories[codes[valid]]
    return np.asarray([
        item.decode() if isinstance(item, bytes) else item for item in array
    ], dtype=object)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--groups", required=True, type=Path)
    args = parser.parse_args()

    with h5py.File(args.h5ad) as handle:
        obs = handle["obs"]
        polydt = decode_column(obs, "CB_polydT")
        ranhex = decode_column(obs, "CB_ranhex")
        labels = decode_column(obs, "cluster_name")
        sublibraries = decode_column(obs, "sublibrary")
    if len({len(polydt), len(ranhex), len(labels), len(sublibraries)}) != 1:
        raise ValueError("AnnData observation columns have inconsistent lengths")

    unique_pairs = {}
    for index, (poly, hexa, label, sublibrary) in enumerate(
        zip(polydt, ranhex, labels, sublibraries)
    ):
        if not poly or not hexa:
            continue
        record = unique_pairs.setdefault(
            (poly, hexa),
            {
                "first_index": index,
                "labels": set(),
                "sublibraries": set(),
                "source_count": 0,
            },
        )
        record["source_count"] += 1
        if label:
            record["labels"].add(label)
        if sublibrary:
            record["sublibraries"].add(sublibrary)

    args.pairs.parent.mkdir(parents=True, exist_ok=True)
    retained_labels = 0
    with (
        open(args.pairs, "w", newline="") as pair_handle,
        open(args.labels, "w", newline="") as label_handle,
        open(args.groups, "w", newline="") as group_handle,
    ):
        pair_writer = csv.writer(pair_handle, delimiter="\t")
        label_writer = csv.writer(label_handle)
        group_writer = csv.writer(group_handle)
        pair_writer.writerow(("cell_id", "polydt_barcode", "ranhex_barcode"))
        for (poly, hexa), record in unique_pairs.items():
            if record["source_count"] != 1:
                continue
            cell_id = f"microglia_{record['first_index']}"
            pair_writer.writerow((cell_id, poly, hexa))
            if len(record["labels"]) == 1 and len(record["sublibraries"]) == 1:
                label = next(iter(record["labels"]))
                sublibrary = next(iter(record["sublibraries"]))
                label_writer.writerow((cell_id, label))
                group_writer.writerow((cell_id, f"x__x__{sublibrary}"))
                retained_labels += 1
    print(f"source_records={len(polydt)}")
    unambiguous = sum(record["source_count"] == 1 for record in unique_pairs.values())
    print(f"distinct_pair_keys={len(unique_pairs)}")
    print(f"unambiguous_pairs={unambiguous}")
    print(f"excluded_reused_pair_keys={len(unique_pairs) - unambiguous}")
    print(f"labeled_pairs={retained_labels}")


if __name__ == "__main__":
    main()
