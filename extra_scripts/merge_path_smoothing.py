#!/usr/bin/env python3
"""Merge path-smoothing evidence and select cross-fitted EB strengths."""

import argparse
from pathlib import Path

from extra_scripts.estimate_path_smoothing import select_cross_fitted_alpha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--minimum-genes", type=int, default=25)
    args = parser.parse_args()
    select_cross_fitted_alpha(
        [path for path in args.inputs if path.is_file() and path.stat().st_size],
        args.output,
        folds=args.folds,
        minimum_genes=args.minimum_genes,
    )


if __name__ == "__main__":
    main()
