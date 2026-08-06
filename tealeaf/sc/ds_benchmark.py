"""Common-unit and split-subject benchmarks for differential splicing."""

from __future__ import annotations

from collections import defaultdict
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from tealeaf.sc.junction_benchmark import benjamini_hochberg


def simes_pvalue(pvalues) -> float:
    values = np.sort(np.asarray(pvalues, dtype=float))
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan
    return float(
        min(1.0, np.min(values * len(values) / np.arange(1, len(values) + 1)))
    )


def aggregate_gene_pvalues(table: pd.DataFrame) -> pd.DataFrame:
    """Combine feature p-values into one Simes p-value per method and gene."""
    required = {"method", "gene_id", "p_value"}
    if missing := required - set(table):
        raise ValueError(f"feature table is missing {sorted(missing)}")
    local = table.loc[table.gene_id.notna() & table.p_value.notna()].copy()
    local["gene_id"] = local.gene_id.astype(str).str.split(".").str[0]
    result = (
        local.groupby(["method", "gene_id"], sort=False)["p_value"]
        .agg(p_value=simes_pvalue, n_features="size")
        .reset_index()
    )
    return result


def leafcutter_cluster_gene_map(
    junction_table: str | Path, counts_path: str | Path
) -> pd.DataFrame:
    """Map LeafCutter clusters to a unique gene through bundle junctions."""
    requested: dict[tuple[str, int, int], str] = {}
    opener = gzip.open if str(counts_path).endswith(".gz") else open
    with opener(counts_path, "rt") as handle:
        next(handle, None)
        for line in handle:
            feature = line.split(None, 1)[0]
            chromosome, start, end, cluster = feature.split(":", 3)
            requested[(chromosome, int(start), int(end))] = (
                f"{chromosome}:{cluster}"
            )

    cluster_genes: dict[str, set[str]] = defaultdict(set)
    opener = gzip.open if str(junction_table).endswith(".gz") else open
    with opener(junction_table, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature = requested.get(
                (row["chromosome"], int(row["start"]), int(row["end"]))
            )
            gene = row.get("gene_id")
            if feature is not None and gene:
                cluster_genes[feature].add(gene.split(".", 1)[0])
    rows = [
        (feature, next(iter(genes)))
        for feature, genes in cluster_genes.items()
        if len(genes) == 1
    ]
    return pd.DataFrame(rows, columns=["feature_id", "gene_id"])


def shared_gene_reproducibility(
    fold_tables: list[pd.DataFrame],
    top_k=(25, 50, 100, 200, 500),
    reference_method: str = "Tealeaf EC GLMM",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score reference/comparator pairs on matched cross-fold gene sets."""
    if len(fold_tables) != 2:
        raise ValueError("exactly two fold tables are required")
    methods = sorted(set(fold_tables[0].method) & set(fold_tables[1].method))
    if reference_method not in methods:
        raise ValueError(f"reference method {reference_method!r} is absent")
    paired_rows = []
    metric_rows = []
    top_rows = []
    for comparator in (method for method in methods if method != reference_method):
        comparison = f"{reference_method} vs {comparator}"
        gene_sets = [
            set(table.loc[table.method.eq(method), "gene_id"])
            for table in fold_tables
            for method in (reference_method, comparator)
        ]
        shared = set.intersection(*gene_sets)
        for method in (reference_method, comparator):
            first = fold_tables[0].loc[
                fold_tables[0].method.eq(method)
                & fold_tables[0].gene_id.isin(shared),
                ["gene_id", "p_value"],
            ].set_index("gene_id")
            second = fold_tables[1].loc[
                fold_tables[1].method.eq(method)
                & fold_tables[1].gene_id.isin(shared),
                ["gene_id", "p_value"],
            ].set_index("gene_id")
            genes = sorted(shared)
            p1 = first.loc[genes, "p_value"].to_numpy(dtype=float)
            p2 = second.loc[genes, "p_value"].to_numpy(dtype=float)
            conjunction = np.maximum(p1, p2)
            conjunction_q = benjamini_hochberg(conjunction)
            q1 = benjamini_hochberg(p1)
            q2 = benjamini_hochberg(p2)
            rank_correlation = (
                float(
                    spearmanr(
                        -np.log10(np.maximum(p1, 1e-300)),
                        -np.log10(np.maximum(p2, 1e-300)),
                    ).statistic
                )
                if len(genes) > 1
                else np.nan
            )
            discoveries_1 = q1 <= 0.05
            discoveries_2 = q2 <= 0.05
            replication_rates = []
            if discoveries_1.any():
                replication_rates.append(float(np.mean(p2[discoveries_1] <= 0.05)))
            if discoveries_2.any():
                replication_rates.append(float(np.mean(p1[discoveries_2] <= 0.05)))
            metric_rows.append({
                "comparison": comparison,
                "method": method,
                "shared_genes": len(genes),
                "fold0_bh": int(discoveries_1.sum()),
                "fold1_bh": int(discoveries_2.sum()),
                "replicated_bh": int(np.sum(conjunction_q <= 0.05)),
                "heldout_nominal_replication": (
                    float(np.mean(replication_rates)) if replication_rates else np.nan
                ),
                "spearman_logp": rank_correlation,
            })
            for k in top_k:
                effective = min(int(k), len(genes))
                selected_1 = set(np.asarray(genes)[np.argsort(p1)[:effective]])
                selected_2 = set(np.asarray(genes)[np.argsort(p2)[:effective]])
                top_rows.append({
                    "comparison": comparison,
                    "method": method,
                    "k": effective,
                    "overlap": len(selected_1 & selected_2),
                    "overlap_fraction": (
                        len(selected_1 & selected_2) / effective if effective else np.nan
                    ),
                })
            paired_rows.extend(
                {
                    "comparison": comparison,
                    "method": method,
                    "gene_id": gene,
                    "fold0_p_value": a,
                    "fold1_p_value": b,
                    "conjunction_p_value": c,
                    "conjunction_q_value": q,
                }
                for gene, a, b, c, q in zip(
                    genes, p1, p2, conjunction, conjunction_q
                )
            )
    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(top_rows),
        pd.DataFrame(paired_rows),
    )
