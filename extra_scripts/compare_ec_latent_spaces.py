#!/usr/bin/env python3
"""Compare matched isoform-space and pre-collapsed path-space EC GLMM fits."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import benjamini_hochberg
from tealeaf.sc.ds_benchmark import (
    aggregate_gene_pvalues,
    shared_gene_reproducibility,
)

METHODS = ("Isoform-space EC GLMM", "Path-space EC GLMM")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isoform-discovery", required=True, type=Path)
    parser.add_argument("--path-discovery", required=True, type=Path)
    parser.add_argument("--isoform-folds", required=True, nargs=2, type=Path)
    parser.add_argument("--path-folds", required=True, nargs=2, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_results(path):
    table = pd.read_csv(path, sep="\t")
    table["p_value"] = table.get("p_value", table["lrt_p_value"])
    table["joint_converged"] = (
        table["null_converged"].astype(bool)
        & table["alternative_converged"].astype(bool)
        & table["p_value"].notna()
    )
    return table


def discovery_comparison(isoform, path):
    rows = []
    for method, table in zip(METHODS, (isoform, path)):
        eligible = table.loc[table.joint_converged].copy()
        q_value = benjamini_hochberg(eligible.p_value.to_numpy())
        rows.append(
            {
                "universe": "method-specific converged",
                "method": method,
                "tests": len(table),
                "eligible_tests": len(eligible),
                "convergence": float(table.joint_converged.mean()),
                "nominal_0.05": int(np.sum(eligible.p_value <= 0.05)),
                "bh_0.05": int(np.sum(q_value <= 0.05)),
                "bh_rate": (
                    float(np.mean(q_value <= 0.05)) if len(q_value) else np.nan
                ),
                "shared_bh_overlap": np.nan,
                "shared_bh_jaccard": np.nan,
            }
        )

    shared = set(isoform.loc[isoform.joint_converged, "test_id"]) & set(
        path.loc[path.joint_converged, "test_id"]
    )
    calls = []
    for method, table in zip(METHODS, (isoform, path)):
        eligible = table.loc[table.test_id.isin(shared)].copy()
        eligible = (
            eligible.drop_duplicates("test_id")
            .set_index("test_id")
            .loc[sorted(shared)]
        )
        q_value = benjamini_hochberg(eligible.p_value.to_numpy())
        selected = set(eligible.index[q_value <= 0.05])
        calls.append(selected)
        rows.append(
            {
                "universe": "common converged tests",
                "method": method,
                "tests": len(shared),
                "eligible_tests": len(shared),
                "convergence": 1.0,
                "nominal_0.05": int(np.sum(eligible.p_value <= 0.05)),
                "bh_0.05": len(selected),
                "bh_rate": len(selected) / len(shared) if shared else np.nan,
                "shared_bh_overlap": np.nan,
                "shared_bh_jaccard": np.nan,
            }
        )
    overlap = len(calls[0] & calls[1])
    union = len(calls[0] | calls[1])
    for row in rows[-2:]:
        row["shared_bh_overlap"] = overlap
        row["shared_bh_jaccard"] = overlap / union if union else np.nan
    return pd.DataFrame(rows)


def replication_comparison(isoform_folds, path_folds):
    tables = [*isoform_folds, *path_folds]
    shared = set.intersection(
        *(set(table.loc[table.joint_converged, "test_id"]) for table in tables)
    )
    fold_tables = []
    for fold in range(2):
        methods = []
        for method, table in zip(
            METHODS,
            (isoform_folds[fold], path_folds[fold]),
        ):
            local = table.loc[
                table.test_id.isin(shared), ["gene_id", "p_value"]
            ].copy()
            local["method"] = method
            methods.append(
                aggregate_gene_pvalues(local[["method", "gene_id", "p_value"]])
            )
        fold_tables.append(pd.concat(methods, ignore_index=True))
    metrics, topk, genes = shared_gene_reproducibility(
        fold_tables,
        reference_method=METHODS[0],
    )
    metrics["shared_block_tests"] = len(shared)
    topk["shared_block_tests"] = len(shared)
    genes["shared_block_tests"] = len(shared)
    return metrics, topk, genes


def main():
    args = parse_args()
    isoform_discovery = read_results(args.isoform_discovery)
    path_discovery = read_results(args.path_discovery)
    isoform_folds = [read_results(path) for path in args.isoform_folds]
    path_folds = [read_results(path) for path in args.path_folds]
    discovery = discovery_comparison(isoform_discovery, path_discovery)
    replication, topk, genes = replication_comparison(
        isoform_folds, path_folds
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    discovery.to_csv(args.output_dir / "discovery.tsv", sep="\t", index=False)
    replication.to_csv(
        args.output_dir / "replication.tsv", sep="\t", index=False
    )
    topk.to_csv(
        args.output_dir / "replication_topk.tsv", sep="\t", index=False
    )
    genes.to_csv(
        args.output_dir / "replication_genes.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
