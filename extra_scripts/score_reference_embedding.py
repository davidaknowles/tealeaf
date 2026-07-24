#!/usr/bin/env python3
"""Score a published reference embedding on mapped, eligible cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc import representation_scoring


def _read_cells(path):
    with open(path) as handle:
        cells = [line.strip() for line in handle if line.strip()]
    if len(cells) != len(set(cells)):
        raise ValueError(f"eligible cell IDs are duplicated in {path}")
    return set(cells)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding", required=True, type=Path)
    parser.add_argument("--cell-map", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--groups", type=Path)
    parser.add_argument("--eligible-cells", type=Path)
    parser.add_argument("--name", default="published_reference")
    parser.add_argument("--representation", default="published_embedding")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--silhouette-sample-size", type=int, default=10_000)
    args = parser.parse_args()

    embedding = pd.read_csv(args.embedding, sep="\t")
    if embedding.shape[1] < 2:
        raise ValueError("reference embedding must contain an ID and coordinates")
    source_column = embedding.columns[0]
    embedding[source_column] = embedding[source_column].astype(str)
    if embedding[source_column].duplicated().any():
        raise ValueError("reference embedding contains duplicate cell IDs")

    cell_map = pd.read_csv(args.cell_map, dtype=str)
    required = {"source_cell_id", "target_cell_id"}
    if not required.issubset(cell_map.columns):
        raise ValueError(f"cell map must contain columns {sorted(required)}")
    if (
        cell_map["source_cell_id"].duplicated().any()
        or cell_map["target_cell_id"].duplicated().any()
    ):
        raise ValueError("cell map must be one-to-one")
    mapped = embedding.merge(
        cell_map,
        left_on=source_column,
        right_on="source_cell_id",
        how="inner",
        validate="one_to_one",
    )
    if len(mapped) != len(embedding):
        raise ValueError(
            f"only {len(mapped)} of {len(embedding)} embedding rows have cell mappings"
        )
    if args.eligible_cells is not None:
        eligible = _read_cells(args.eligible_cells)
        mapped = mapped[mapped["target_cell_id"].isin(eligible)].copy()
    if mapped.empty:
        raise ValueError("no mapped reference cells pass eligibility filtering")

    coordinate_columns = [
        column
        for column in embedding.columns
        if column != source_column
    ]
    coordinates = mapped[coordinate_columns].to_numpy(dtype=np.float32)
    if not np.isfinite(coordinates).all():
        raise ValueError("reference embedding contains nonfinite coordinates")

    labels = pd.read_csv(
        args.labels, header=None, names=["cell_id", "label"], dtype=str
    )
    groups = None
    if args.groups is not None:
        groups = pd.read_csv(
            args.groups, header=None, names=["cell_id", "group"], dtype=str
        )
    positions, aligned_labels, aligned_groups, aligned_ids = (
        representation_scoring.align_reference_metadata(
            mapped["target_cell_id"].to_numpy(), labels, groups
        )
    )
    coordinates = coordinates[positions]
    report, folds = representation_scoring.score_embedding(
        coordinates,
        aligned_labels,
        aligned_groups,
        n_splits=args.folds,
        silhouette_sample_size=args.silhouette_sample_size,
    )
    report.update(
        status="ok",
        name=args.name,
        representation=args.representation,
        standard_analysis_baseline=True,
        n_published_cells=int(len(embedding)),
        n_mapped_eligible_cells=int(len(mapped)),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = args.output_dir / f"{args.name}_"
    np.savez_compressed(
        f"{output_prefix}embedding.npz",
        embedding=coordinates,
    )
    np.savetxt(f"{output_prefix}embedding_cells.txt", aligned_ids, fmt="%s")
    representation_scoring.write_score(report, folds, output_prefix)
    pd.DataFrame([report]).to_csv(
        args.output_dir / "label_score_summary.csv", index=False
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
