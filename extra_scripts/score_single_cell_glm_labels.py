#!/usr/bin/env python3
"""Score one or more tealeaf cell-factor fits against reference labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc import representation_scoring


def named_prefix(value):
    name, separator, prefix = value.partition("=")
    if not separator or not name or not prefix:
        raise argparse.ArgumentTypeError("fits must have the form NAME=OUTPUT_PREFIX")
    return name, Path(prefix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", action="append", required=True, type=named_prefix)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--groups", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--pca-components", type=int, default=30)
    parser.add_argument("--silhouette-sample-size", type=int, default=10_000)
    args = parser.parse_args()

    labels = pd.read_csv(
        args.labels, header=None, names=["cell_id", "label"], dtype=str
    )
    groups = None
    if args.groups is not None:
        group_table = pd.read_csv(
            args.groups, header=None, names=["cell_id", "combined"], dtype=str
        )
        split = group_table["combined"].str.rsplit("__", n=2, expand=True)
        if split.shape[1] != 3:
            raise ValueError("group values must have label__condition__sample format")
        groups = group_table[["cell_id"]].assign(group=split.iloc[:, 2])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for name, prefix in args.fit:
        try:
            cell_ids, factors = representation_scoring.load_factor_representation(prefix)
            aligned = representation_scoring.align_reference_labels(
                cell_ids, factors, labels, groups
            )
            report, folds = representation_scoring.score_representation(
                aligned[0], aligned[1], aligned[2], n_splits=args.folds,
                pca_components=args.pca_components,
                silhouette_sample_size=args.silhouette_sample_size,
            )
            report["status"] = "ok"
        except (OSError, ValueError) as error:
            report = {"status": "invalid", "error": str(error)}
            folds = pd.DataFrame()
        report.update(name=name, fit_prefix=str(prefix))
        representation_scoring.write_score(
            report, folds, args.output_dir / f"{name}_"
        )
        summaries.append(report)
        print(json.dumps(report, indent=2))
    pd.DataFrame(summaries).to_csv(args.output_dir / "label_score_summary.csv", index=False)


if __name__ == "__main__":
    main()
