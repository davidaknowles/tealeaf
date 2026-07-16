#!/usr/bin/env python3
"""Compare column-normalized binary and fixed weighted EC designs."""

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import sc_utils


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--weighted", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    membership = sp.load_npz(args.membership).tocsr()
    weighted = sp.load_npz(args.weighted).tocsr()
    binary = sc_utils.glm_design_matrix(
        membership,
        np.ones(membership.shape[1]),
        parameterization="phi",
        design="binary",
    )
    if weighted.shape != binary.shape:
        raise ValueError(f"shape mismatch: weighted={weighted.shape}, binary={binary.shape}")
    weighted = weighted.tocsr()
    weighted.eliminate_zeros()
    binary.eliminate_zeros()

    overlap_nnz = weighted.sign().multiply(binary.sign()).nnz
    union_nnz = weighted.nnz + binary.nnz - overlap_nnz
    delta = weighted - binary
    delta.eliminate_zeros()
    abs_delta = np.abs(delta.data)
    # Include entries that agree exactly when summarizing differences over the
    # union of both sparse supports.
    union_abs_delta = np.zeros(union_nnz, dtype=float)
    union_abs_delta[:abs_delta.size] = abs_delta
    weighted_norm = np.sqrt(weighted.multiply(weighted).sum())
    binary_norm = np.sqrt(binary.multiply(binary).sum())
    denominator = weighted_norm * binary_norm
    metrics = {
        "shape": list(weighted.shape),
        "binary_nnz": int(binary.nnz),
        "weighted_nnz": int(weighted.nnz),
        "support_overlap_nnz": int(overlap_nnz),
        "support_union_nnz": int(union_nnz),
        "fraction_binary_support_retained": float(overlap_nnz / binary.nnz),
        "fraction_weighted_support_outside_binary": float(
            (weighted.nnz - overlap_nnz) / weighted.nnz
        ),
        "cosine_similarity": float(weighted.multiply(binary).sum() / denominator),
        "relative_frobenius_difference": float(
            np.sqrt(delta.multiply(delta).sum()) / binary_norm
        ),
        "mean_absolute_difference": float(union_abs_delta.mean()),
        "median_absolute_difference": float(np.median(union_abs_delta)),
        "max_absolute_difference": float(union_abs_delta.max()),
        "fraction_abs_difference_gt_1e-4": float(
            np.mean(union_abs_delta > 1e-4)
        ),
        "fraction_abs_difference_gt_1e-3": float(
            np.mean(union_abs_delta > 1e-3)
        ),
        "fraction_abs_difference_gt_1e-2": float(
            np.mean(union_abs_delta > 1e-2)
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
