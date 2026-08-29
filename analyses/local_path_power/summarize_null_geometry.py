#!/usr/bin/env python3
"""Measure local-path contrast leakage in the global isoform null basis."""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import group_metadata, local_gene_data
from tealeaf.sc import differential, ec_block_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def screened_data(data_cache, settings):
    with data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(handle)
    metadata = group_metadata(groups)
    subject_folds = settings.get("subject_folds")
    if subject_folds is not None:
        folds = pd.read_csv(subject_folds, sep="\t", dtype={"subject": str})
        selected = set(folds.loc[folds.fold.eq(settings["subject_fold"]), "subject"].astype(str))
        keep = metadata["mouse"].astype(str).isin(selected).to_numpy()
        metadata = metadata.loc[keep].reset_index(drop=True)
        counts = tuple(matrix[keep] for matrix in counts)
    if str(settings["test_effect"]).startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        cell_types = represented.index[represented >= int(settings["min_celltype_mice"])]
        keep = metadata["cell_type"].isin(cell_types).to_numpy()
        counts = tuple(matrix[keep] for matrix in counts)
    return counts, gene_ecs, designs


def main():
    args = parse_args()
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    counts, gene_ecs, designs = screened_data(args.data_cache, cached["settings"])
    weight_cache = {}
    records = []
    for candidate in cached["candidates"]:
        test_id, block_id, gene_id, gene, transcripts, path_index, _, rows, _, _ = candidate
        key = (int(gene), tuple(rows), tuple(transcripts))
        if key not in weight_cache:
            local_counts = tuple(matrix[rows] for matrix in counts)
            base, _, _ = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(rows), 1)),
                np.arange(len(rows)),
                drop_zero=False,
            )
            weight_cache[key] = ec_block_glmm.pooled_isoform_weights(base)
        weights = weight_cache[key]
        path_index = np.asarray(path_index, dtype=int)
        tested_basis, nuisance_basis = ec_block_glmm.block_effect_bases(path_index)
        jacobian = differential.path_logratio_jacobian(weights, path_index)[:, :-1]
        tested_derivative = jacobian @ tested_basis
        nuisance_derivative = jacobian @ nuisance_basis
        tested_norm = float(np.linalg.norm(tested_derivative))
        nuisance_norm = float(np.linalg.norm(nuisance_derivative))
        records.append({
            "test_id": test_id,
            "block_id": block_id,
            "gene_id": gene_id,
            "n_isoforms": len(path_index),
            "n_paths": len(np.unique(path_index[path_index >= 0])),
            "n_background_isoforms": int(np.sum(path_index < 0)),
            "n_global_null_directions": nuisance_basis.shape[1],
            "null_path_derivative_norm": nuisance_norm,
            "tested_path_derivative_norm": tested_norm,
            "relative_null_leakage": nuisance_norm / tested_norm if tested_norm > 0 else np.nan,
            "max_abs_null_path_derivative": float(
                np.max(np.abs(nuisance_derivative), initial=0.0)
            ),
        })
    table = pd.DataFrame(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    leakage = table["relative_null_leakage"].dropna()
    print(json.dumps({
        "tests": len(table),
        "tests_with_global_nuisance": int(table["n_global_null_directions"].gt(0).sum()),
        "tests_with_nonzero_leakage": int(table["relative_null_leakage"].gt(1e-10).sum()),
        "median_relative_leakage": float(leakage.median()),
        "p90_relative_leakage": float(leakage.quantile(0.9)),
        "maximum_relative_leakage": float(leakage.max()),
    }, indent=2))


if __name__ == "__main__":
    main()
