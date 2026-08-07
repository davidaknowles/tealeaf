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


def filter_grouped_subject_records(
    grouped: dict,
    selected_subjects,
    *,
    subject_key: str = "mouse",
) -> dict:
    """Restrict grouped record dictionaries to a subject set."""
    selected = {str(value) for value in selected_subjects}
    result = {}
    for feature, records in grouped.items():
        retained = [
            record
            for record in records
            if str(record[subject_key]) in selected
        ]
        if retained:
            result[feature] = retained
    return result


def compositional_pairwise_table(
    table: pd.DataFrame,
    *,
    method: str,
    min_paired_subjects: int = 0,
) -> pd.DataFrame:
    """Normalize paired compositional block tests for common-unit scoring."""
    required = {
        "contrast",
        "gene_id",
        "p_value",
        "converged",
        "n_mice",
    }
    if missing := required - set(table):
        raise ValueError(f"compositional table is missing {sorted(missing)}")
    eligible = table["converged"].astype(bool) & table["n_mice"].ge(
        int(min_paired_subjects)
    )
    if "inference_eligible" in table:
        eligible &= table["inference_eligible"].astype(bool)
    result = table.loc[
        eligible,
        ["gene_id", "contrast", "p_value"],
    ].copy()
    levels = result["contrast"].str.split("_vs_", n=1, expand=True)
    if levels.shape[1] != 2 or levels.isna().any(axis=None):
        raise ValueError("pairwise contrasts must have LEVEL_A_vs_LEVEL_B form")
    result["level_a"] = levels[0].to_numpy()
    result["level_b"] = levels[1].to_numpy()
    result.insert(0, "method", method)
    return result.drop(columns="contrast")


def calibrate_pvalues_from_null(
    table: pd.DataFrame,
    null_table: pd.DataFrame,
) -> pd.DataFrame:
    """Map method-specific p-values through pooled empirical null CDFs."""
    required = {"method", "p_value"}
    if missing := required - set(table):
        raise ValueError(f"test table is missing {sorted(missing)}")
    if missing := required - set(null_table):
        raise ValueError(f"null table is missing {sorted(missing)}")
    result = table.copy()
    result["raw_p_value"] = result["p_value"]
    for method, indices in result.groupby("method").groups.items():
        pool = np.sort(
            pd.to_numeric(
                null_table.loc[null_table.method.eq(method), "p_value"],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
        )
        if not len(pool):
            continue
        values = pd.to_numeric(
            result.loc[indices, "p_value"], errors="coerce"
        ).to_numpy(dtype=float)
        finite = np.isfinite(values)
        values[finite] = (
            1 + np.searchsorted(pool, values[finite], side="right")
        ) / (len(pool) + 1)
        result.loc[indices, "p_value"] = values
    return result


def normalized_junction_pairwise_table(
    table: pd.DataFrame,
    *,
    method: str,
    min_paired_subjects: int = 0,
) -> pd.DataFrame:
    """Normalize a common-schema junction table for pairwise scoring."""
    required = {
        "effect",
        "gene_id",
        "level_a",
        "level_b",
        "p_value",
    }
    if missing := required - set(table):
        raise ValueError(f"junction table is missing {sorted(missing)}")
    eligible = table.effect.eq("cell_type") & table.gene_id.notna()
    if "converged" in table:
        eligible &= table.converged.astype(bool)
    if "n_subjects" in table:
        eligible &= table.n_subjects.ge(int(min_paired_subjects))
    result = table.loc[
        eligible,
        ["gene_id", "level_a", "level_b", "p_value"],
    ].copy()
    result.insert(0, "method", method)
    return result


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


def aggregate_feature_pvalues(table: pd.DataFrame) -> pd.DataFrame:
    """Combine pairwise p-values into one Simes p-value per feature."""
    required = {"method", "feature_id", "gene_id", "p_value"}
    if missing := required - set(table):
        raise ValueError(f"pairwise table is missing {sorted(missing)}")
    local = table.loc[
        table.feature_id.notna()
        & table.gene_id.notna()
        & table.p_value.notna()
    ].copy()
    local["gene_id"] = local.gene_id.astype(str).str.split(".").str[0]
    return (
        local.groupby(["method", "feature_id", "gene_id"], sort=False)[
            "p_value"
        ]
        .agg(p_value=simes_pvalue, n_pairwise="size")
        .reset_index()
    )


def aggregate_gene_pair_pvalues(table: pd.DataFrame) -> pd.DataFrame:
    """Combine feature tests within each gene and unordered level pair."""
    required = {"method", "gene_id", "level_a", "level_b", "p_value"}
    if missing := required - set(table):
        raise ValueError(f"pairwise table is missing {sorted(missing)}")
    local = table.loc[
        table.gene_id.notna()
        & table.level_a.notna()
        & table.level_b.notna()
        & table.p_value.notna()
    ].copy()
    local["gene_id"] = local.gene_id.astype(str).str.split(".").str[0]
    levels = np.sort(
        local[["level_a", "level_b"]].astype(str).to_numpy(), axis=1
    )
    local["pair_id"] = levels[:, 0] + "||" + levels[:, 1]
    return (
        local.groupby(["method", "gene_id", "pair_id"], sort=False)[
            "p_value"
        ]
        .agg(p_value=simes_pvalue, n_features="size")
        .reset_index()
    )


def shared_pair_gene_reproducibility(
    fold_tables: list[pd.DataFrame],
    top_k=(25, 50, 100, 200, 500),
    reference_method: str = "Tealeaf paired EC GLMM",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score methods after matching gene-by-level-pair hypotheses."""
    if len(fold_tables) != 2:
        raise ValueError("exactly two fold tables are required")
    paired = [aggregate_gene_pair_pvalues(table) for table in fold_tables]
    methods = sorted(set(paired[0].method) & set(paired[1].method))
    if reference_method not in methods:
        raise ValueError(f"reference method {reference_method!r} is absent")
    metrics = []
    topk = []
    genes = []
    for comparator in (method for method in methods if method != reference_method):
        key_sets = []
        for table in paired:
            for method in (reference_method, comparator):
                local = table.loc[table.method.eq(method), ["gene_id", "pair_id"]]
                key_sets.append(set(map(tuple, local.to_numpy())))
        shared = set.intersection(*key_sets)
        fold_genes = []
        for table in paired:
            method_tables = []
            for method in (reference_method, comparator):
                local = table.loc[table.method.eq(method)].copy()
                keys = list(zip(local.gene_id, local.pair_id))
                local = local.loc[[key in shared for key in keys]]
                method_tables.append(
                    aggregate_gene_pvalues(local[["method", "gene_id", "p_value"]])
                )
            fold_genes.append(pd.concat(method_tables, ignore_index=True))
        local_metrics, local_topk, local_genes = shared_gene_reproducibility(
            fold_genes,
            top_k=top_k,
            reference_method=reference_method,
        )
        local_metrics["shared_gene_pairs"] = len(shared)
        local_topk["shared_gene_pairs"] = len(shared)
        local_genes["shared_gene_pairs"] = len(shared)
        metrics.append(local_metrics)
        topk.append(local_topk)
        genes.append(local_genes)
    return (
        pd.concat(metrics, ignore_index=True),
        pd.concat(topk, ignore_index=True),
        pd.concat(genes, ignore_index=True),
    )


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
