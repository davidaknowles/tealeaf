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
    parser.add_argument("--wald-discovery", type=Path)
    parser.add_argument("--wald-folds", nargs=2, type=Path)
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


def discovery_comparison(tables, methods):
    rows = []
    tested_universe = set.union(*(
        set(table["test_id"]) for table in tables
    ))
    for method, table in zip(methods, tables):
        eligible = table.loc[table.joint_converged].copy()
        q_value = benjamini_hochberg(eligible.p_value.to_numpy())
        rows.append(
            {
                "universe": "method-specific converged",
                "method": method,
                "tests": len(tested_universe),
                "eligible_tests": len(eligible),
                "convergence": (
                    len(eligible) / len(tested_universe)
                    if tested_universe
                    else np.nan
                ),
                "nominal_0.05": int(np.sum(eligible.p_value <= 0.05)),
                "bh_0.05": int(np.sum(q_value <= 0.05)),
                "bh_rate": (
                    float(np.mean(q_value <= 0.05)) if len(q_value) else np.nan
                ),
                "shared_bh_overlap": np.nan,
                "shared_bh_jaccard": np.nan,
            }
        )

    shared = set.intersection(*(
        set(table.loc[table.joint_converged, "test_id"])
        for table in tables
    ))
    calls = []
    for method, table in zip(methods, tables):
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
    for row, selected in zip(rows[-len(tables):], calls):
        overlap = len(calls[0] & selected)
        union = len(calls[0] | selected)
        row["shared_bh_overlap"] = overlap
        row["shared_bh_jaccard"] = overlap / union if union else np.nan
    return pd.DataFrame(rows)


def replication_comparison(method_folds, methods):
    tables = [table for folds in method_folds for table in folds]
    shared = set.intersection(
        *(set(table.loc[table.joint_converged, "test_id"]) for table in tables)
    )
    fold_tables = []
    for fold in range(2):
        aggregated = []
        for method, folds in zip(methods, method_folds):
            table = folds[fold]
            local = table.loc[
                table.test_id.isin(shared), ["gene_id", "p_value"]
            ].copy()
            local["method"] = method
            aggregated.append(
                aggregate_gene_pvalues(local[["method", "gene_id", "p_value"]])
            )
        fold_tables.append(pd.concat(aggregated, ignore_index=True))
    metrics, topk, genes = shared_gene_reproducibility(
        fold_tables,
        reference_method=methods[0],
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
    methods = list(METHODS)
    discovery_tables = [isoform_discovery, path_discovery]
    method_folds = [isoform_folds, path_folds]
    if (args.wald_discovery is None) != (args.wald_folds is None):
        raise ValueError("Wald discovery and folds must be supplied together")
    if args.wald_discovery is not None:
        methods.append("Estimate-once path Wald")
        discovery_tables.append(read_results(args.wald_discovery))
        method_folds.append([read_results(path) for path in args.wald_folds])
    discovery = discovery_comparison(discovery_tables, methods)
    replication, topk, genes = replication_comparison(
        method_folds, methods
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
