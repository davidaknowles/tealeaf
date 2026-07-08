"""Aggregate a cell-by-transcript sparse matrix into pseudobulk samples.

This helper is intended for tealeaf single-cell runs where transcript-level
quantification already exists, for example alevin-fry ``quants_mat.npz``.  It
creates the preprocessed files expected by ``tealeaf-sc --preprocessed``:

    {outprefix}pseudo_spliced_count.npz
    {outprefix}pseudo_spliced_TPM.npz
    {outprefix}pseudo_spliced_cols.txt
    {outprefix}pseudo_rows.txt

The barcode-to-group file should contain two comma-separated columns:

    barcode,group
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_lines(path: Path) -> list[str]:
    with path.open() as handle:
        return [line.rstrip("\n") for line in handle]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate a cell-by-transcript sparse matrix by barcode group."
    )
    parser.add_argument("--matrix", required=True, type=Path,
                        help="Input sparse matrix .npz with cells as rows and transcripts as columns.")
    parser.add_argument("--rows", required=True, type=Path,
                        help="Cell barcode names, one per matrix row.")
    parser.add_argument("--cols", required=True, type=Path,
                        help="Transcript names, one per matrix column.")
    parser.add_argument("--barcode-groups", required=True, type=Path,
                        help="CSV without header: barcode,group.")
    parser.add_argument("-o", "--outprefix", required=True, type=Path,
                        help="Output prefix for tealeaf-sc --preprocessed inputs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    barcodes = read_lines(args.rows)
    barcode_to_row = {barcode: idx for idx, barcode in enumerate(barcodes)}

    groups_df = pd.read_csv(
        args.barcode_groups,
        header=None,
        names=["barcode", "group"],
        dtype=str,
    )
    groups_df = groups_df[groups_df["barcode"].isin(barcode_to_row)].copy()
    groups_df["row"] = groups_df["barcode"].map(barcode_to_row)

    # If a barcode appears more than once after upstream processing, keep one
    # barcode/group pair. Conflicting labels should be resolved before this step.
    groups_df = groups_df.drop_duplicates(["barcode", "group"])

    group_names = sorted(groups_df["group"].unique())
    group_to_idx = {group: idx for idx, group in enumerate(group_names)}
    groups_df["group_idx"] = groups_df["group"].map(group_to_idx)

    print(f"Loaded {len(barcodes)} matrix row names")
    print(f"Using {len(groups_df)} barcode/group rows across {len(group_names)} groups")

    matrix = sp.load_npz(args.matrix).tocsr()
    print(f"Loaded matrix shape={matrix.shape} nnz={matrix.nnz}")

    selector = sp.coo_matrix(
        (
            np.ones(len(groups_df), dtype=np.float64),
            (groups_df["group_idx"].to_numpy(), groups_df["row"].to_numpy()),
        ),
        shape=(len(group_names), matrix.shape[0]),
    ).tocsr()

    pseudobulk = selector @ matrix
    pseudobulk = pseudobulk.tocsr()
    print(f"Aggregated matrix shape={pseudobulk.shape} nnz={pseudobulk.nnz}")

    outprefix = args.outprefix
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    sp.save_npz(f"{outprefix}pseudo_spliced_count.npz", pseudobulk)
    sp.save_npz(f"{outprefix}pseudo_spliced_TPM.npz", pseudobulk)

    with Path(f"{outprefix}pseudo_rows.txt").open("w") as output:
        output.write("\n".join(group_names) + "\n")

    cols = read_lines(args.cols)
    with Path(f"{outprefix}pseudo_spliced_cols.txt").open("w") as output:
        output.write("\n".join(cols) + "\n")


if __name__ == "__main__":
    main()
