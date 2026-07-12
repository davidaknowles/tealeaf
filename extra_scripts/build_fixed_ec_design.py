#!/usr/bin/env python3
"""Precompute a shared weighted EC design from an alevin-fry sidecar."""

import argparse
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from tealeaf.sc import sc_utils


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--probabilities", type=Path, default=None)
    parser.add_argument("--cache", type=Path, default=None)
    args = parser.parse_args()

    probability_file = args.probabilities or args.alevin_dir / "gene_eqclass_probs.tsv.gz"
    cache_file = args.cache or args.alevin_dir / "gene_eqclass_fixed_weights.npz"
    membership = sp.load_npz(args.alevin_dir / "gene_eqclass.npz")
    design = sc_utils.averaged_ec_probability_matrix(
        probability_file, membership, cache_file
    )
    column_sums = design.sum(axis=0)
    print(f"wrote {cache_file}")
    print(f"shape={design.shape} nnz={design.nnz}")
    print(f"nonempty_transcripts={np.count_nonzero(np.asarray(column_sums).ravel())}")


if __name__ == "__main__":
    main()
