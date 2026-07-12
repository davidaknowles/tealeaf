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
    if not np.array_equal(weighted.indptr, binary.indptr) or not np.array_equal(
        weighted.indices, binary.indices
    ):
        raise ValueError("weighted and binary designs do not have identical sparse support")

    delta = weighted.data - binary.data
    abs_delta = np.abs(delta)
    weighted_norm = np.linalg.norm(weighted.data)
    binary_norm = np.linalg.norm(binary.data)
    denominator = weighted_norm * binary_norm
    metrics = {
        "shape": list(weighted.shape),
        "nnz": int(weighted.nnz),
        "cosine_similarity": float(np.dot(weighted.data, binary.data) / denominator),
        "relative_frobenius_difference": float(np.linalg.norm(delta) / binary_norm),
        "mean_absolute_difference": float(abs_delta.mean()),
        "median_absolute_difference": float(np.median(abs_delta)),
        "max_absolute_difference": float(abs_delta.max()),
        "fraction_abs_difference_gt_1e-4": float(np.mean(abs_delta > 1e-4)),
        "fraction_abs_difference_gt_1e-3": float(np.mean(abs_delta > 1e-3)),
        "fraction_abs_difference_gt_1e-2": float(np.mean(abs_delta > 1e-2)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
