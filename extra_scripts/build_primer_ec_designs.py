#!/usr/bin/env python3
"""Build separate fixed weighted EC designs for paired primer chemistries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import sc_utils


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    args = parser.parse_args()

    with open(args.alevin_dir / "quants_mat_rows.txt") as handle:
        barcodes = [line.strip() for line in handle]
    barcode_to_row = {barcode: index for index, barcode in enumerate(barcodes)}
    if len(barcode_to_row) != len(barcodes):
        raise ValueError("alevin barcode rows are not unique")
    groups = np.full(len(barcodes), -1, dtype=np.int8)
    with open(args.pairs, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            poly = barcode_to_row.get(row["polydt_barcode"])
            hexa = barcode_to_row.get(row["ranhex_barcode"])
            if poly is not None:
                groups[poly] = 0
            if hexa is not None:
                groups[hexa] = 1
    membership = sp.load_npz(args.alevin_dir / "gene_eqclass.npz")
    caches = [
        args.alevin_dir / "gene_eqclass_fixed_weights_polydt.npz",
        args.alevin_dir / "gene_eqclass_fixed_weights_ranhex.npz",
    ]
    designs = sc_utils.grouped_ec_probability_matrices(
        args.alevin_dir / "gene_eqclass_probs.tsv.gz",
        membership,
        groups,
        2,
        cache_files=caches,
    )
    for name, path, design in zip(("polydT", "ranhex"), caches, designs):
        print(
            f"{name}: matched_rows={np.count_nonzero(groups == (name == 'ranhex'))} "
            f"shape={design.shape} nnz={design.nnz} cache={path}"
        )


if __name__ == "__main__":
    main()
