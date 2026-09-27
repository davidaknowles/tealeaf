#!/usr/bin/env python3
"""Score primer-aware SUPPA2 on the Tealeaf matched split universe."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tealeaf.sc.ds_benchmark import shared_pair_gene_reproducibility


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tealeaf-shards", action="append", required=True, type=Path, help="Fold directory containing shard_*/paired_path.tsv; repeat twice.")
    parser.add_argument("--primer-tests", action="append", required=True, type=Path, help="Primer-aware fold test table; repeat twice.")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def load_tealeaf(path):
    files = sorted(path.glob("shard_*/paired_path.tsv"))
    if not files:
        raise FileNotFoundError(f"no paired_path.tsv shards under {path}")
    table = pd.concat((pd.read_csv(file, sep="\t") for file in files), ignore_index=True)
    table = table.loc[table.method.eq("local_path") & table.converged.astype(str).str.lower().eq("true")].copy()
    table["method"] = "Tealeaf direct local path"
    return table[["method", "gene_id", "level_a", "level_b", "p_value"]]


def load_primer(path):
    table = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    table = table.loc[table.effect.eq("cell_type") & table.p_value.notna()].copy()
    table["method"] = "SUPPA2 (primer aware)"
    return table[["method", "gene_id", "level_a", "level_b", "p_value"]]


def main():
    args = parse_args()
    if len(args.tealeaf_shards) != 2 or len(args.primer_tests) != 2:
        raise ValueError("provide exactly two Tealeaf shard directories and two primer test tables")
    folds = [pd.concat([load_tealeaf(tealeaf), load_primer(primer)], ignore_index=True) for tealeaf, primer in zip(args.tealeaf_shards, args.primer_tests)]
    metrics, topk, genes = shared_pair_gene_reproducibility(folds, reference_method="Tealeaf direct local path")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.output_dir / "split_data_tealeaf_comparison_metrics.tsv", sep="\t", index=False)
    topk.to_csv(args.output_dir / "split_data_tealeaf_comparison_topk.tsv", sep="\t", index=False)
    genes.to_csv(args.output_dir / "split_data_tealeaf_comparison_genes.tsv.gz", sep="\t", index=False, compression="gzip")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
