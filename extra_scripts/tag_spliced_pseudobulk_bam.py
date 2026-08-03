#!/usr/bin/env python3
"""Retain one spliced alignment per cell/UMI/locus and tag its pseudobulk."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pysam

from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--cell-metadata", type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-deduplicate", action="store_true")
    args = parser.parse_args()

    barcode_to_sample = {}
    with open(args.cell_metadata, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["run"] not in (args.run, "*"):
                continue
            clean_cell_type = "".join(
                character if character.isalnum() or character in "_.-" else "_"
                for character in row["cell_type"]
            )
            barcode_to_sample[row["barcode"]] = row["subject"] + "__" + clean_cell_type
    if not barcode_to_sample:
        raise ValueError(f"no cell metadata for run {args.run}")

    with pysam.AlignmentFile(args.input, "rb") as source:
        header = source.header.to_dict()
        samples = sorted(set(barcode_to_sample.values()))
        header["RG"] = [{"ID": sample, "SM": sample} for sample in samples]
        args.output.parent.mkdir(parents=True, exist_ok=True)
        seen = set()
        retained = 0
        with pysam.AlignmentFile(args.output, "wb", header=header) as target:
            for alignment in source.fetch(until_eof=True):
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    continue
                if alignment.get_tag("NH") != 1 or not any(
                    operation == 3 for operation, _ in (alignment.cigartuples or ())
                ):
                    continue
                if not alignment.has_tag("CB") or not alignment.has_tag("UB"):
                    continue
                barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
                sample = barcode_to_sample.get(barcode)
                if sample is None:
                    continue
                if not args.no_deduplicate:
                    key = (
                        barcode,
                        alignment.get_tag("UB"),
                        alignment.reference_id,
                        alignment.reference_start,
                        alignment.cigarstring,
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                alignment.set_tag("RG", sample, value_type="Z", replace=True)
                target.write(alignment)
                retained += 1
    if not retained:
        raise ValueError("no labeled spliced alignments were retained")
    print(f"retained={retained} pseudobulks={len(samples)}")


if __name__ == "__main__":
    main()
