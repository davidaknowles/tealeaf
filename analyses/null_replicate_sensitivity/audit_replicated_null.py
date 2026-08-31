"""Propagate paired-path null families through split reproducibility scoring."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import benjamini_hochberg
from extra_scripts.score_ds_reproducibility import load_pairwise_fold
from tealeaf.sc.ds_benchmark import aggregate_gene_pair_pvalues, simes_pvalue


TEALEAF_LABEL = "Tealeaf direct local path"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold-dir", action="append", required=True, type=Path)
    parser.add_argument("--leafcutter-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def add_pair_id(table: pd.DataFrame) -> pd.DataFrame:
    result = table.copy()
    result["gene_id"] = result.gene_id.astype(str).str.split(".").str[0]
    levels = np.sort(
        result[["level_a", "level_b"]].astype(str).to_numpy(), axis=1
    )
    result["pair_id"] = levels[:, 0] + "||" + levels[:, 1]
    return result


def shared_pair_sets(folds: list[pd.DataFrame]) -> dict[str, set[tuple[str, str]]]:
    paired = [aggregate_gene_pair_pvalues(table) for table in folds]
    methods = sorted(set(paired[0].method) & set(paired[1].method))
    result = {}
    for comparator in (method for method in methods if method != TEALEAF_LABEL):
        key_sets = []
        for table in paired:
            for method in (TEALEAF_LABEL, comparator):
                local = table.loc[
                    table.method.eq(method), ["gene_id", "pair_id"]
                ]
                key_sets.append(set(map(tuple, local.to_numpy())))
        result[comparator] = set.intersection(*key_sets)
    return result


def load_null_pairs(path: Path) -> pd.DataFrame:
    observed = pd.read_csv(
        path / "tealeaf" / "paired_path.tsv",
        sep="\t",
        usecols=["test_id", "gene_id", "level_a", "level_b"],
    )
    observed = add_pair_id(observed).drop(columns=["level_a", "level_b"])
    null = pd.read_csv(
        path / "tealeaf" / "paired_path_null.tsv.gz",
        sep="\t",
        usecols=["test_id", "replicate", "p_value"],
    ).merge(observed, on="test_id", validate="many_to_one")
    return (
        null.groupby(["replicate", "gene_id", "pair_id"], sort=False)[
            "p_value"
        ]
        .agg(simes_pvalue)
        .reset_index()
    )


def aggregate_shared_genes(
    table: pd.DataFrame, shared: set[tuple[str, str]]
) -> pd.Series:
    keys = pd.MultiIndex.from_frame(table[["gene_id", "pair_id"]])
    local = table.loc[keys.isin(shared)]
    return local.groupby(["replicate", "gene_id"], sort=False)["p_value"].agg(
        simes_pvalue
    )


def audit_family(folds: list[pd.Series], replicate: int) -> int:
    genes = [table.xs(replicate, level="replicate") for table in folds]
    conjunction = pd.concat(genes, axis=1, join="inner").max(axis=1)
    return int(np.sum(benjamini_hochberg(conjunction.to_numpy()) <= 0.05))


def main() -> None:
    args = parse_args()
    if len(args.fold_dir) != 2:
        raise ValueError("provide exactly two --fold-dir values")
    leaf_map = pd.read_csv(args.leafcutter_map, sep="\t")
    observed = [
        load_pairwise_fold(
            path,
            leaf_map,
            "tealeaf",
            TEALEAF_LABEL,
            "laplace_multinomial",
            0,
            4,
            "paired_path",
            "comparison_majiq_min3_cov3",
        )
        for path in args.fold_dir
    ]
    shared = shared_pair_sets(observed)
    null = [load_null_pairs(path) for path in args.fold_dir]
    records = []
    for comparator, pairs in shared.items():
        null_genes = [aggregate_shared_genes(table, pairs) for table in null]
        calls = np.array(
            [audit_family(null_genes, replicate) for replicate in range(128)]
        )
        for null_replicates in (32, 64, 128):
            local_calls = calls[:null_replicates]
            records.append(
                {
                    "null_replicates": null_replicates,
                    "comparator": comparator,
                    "shared_gene_pairs": len(pairs),
                    "families_with_bh": int(np.sum(local_calls > 0)),
                    "mean_bh": float(np.mean(local_calls)),
                    "max_bh": int(np.max(local_calls)),
                }
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
