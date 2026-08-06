#!/usr/bin/env python3
"""Score split-subject DS results on a shared gene universe."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import pandas as pd

from tealeaf.sc.ds_benchmark import (
    aggregate_feature_pvalues,
    aggregate_gene_pvalues,
    leafcutter_cluster_gene_map,
    shared_pair_gene_reproducibility,
    shared_gene_reproducibility,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold-dir", action="append", type=Path, required=True)
    parser.add_argument("--junction-table", type=Path, required=True)
    parser.add_argument("--leafcutter-counts", type=Path, required=True)
    parser.add_argument("--leafcutter-map", type=Path, required=True)
    parser.add_argument("--tealeaf-subdir", default="tealeaf")
    parser.add_argument("--tealeaf-label", default="Tealeaf EC GLMM")
    parser.add_argument("--tealeaf-fit-method", default="laplace_multinomial")
    parser.add_argument("--min-median-gene-umis", type=float, default=0.0)
    parser.add_argument("--match-pairs", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_fold(
    path: Path,
    leaf_map: pd.DataFrame,
    tealeaf_subdir: str,
    tealeaf_label: str,
    tealeaf_fit_method: str,
    min_median_gene_umis: float,
) -> pd.DataFrame:
    external = pd.read_csv(path / "comparison/cell_type_simes.tsv.gz", sep="\t")
    external.loc[external.method.eq("LeafCutter"), "gene_id"] = external.loc[
        external.method.eq("LeafCutter"), "feature_id"
    ].map(leaf_map.set_index("feature_id").gene_id)
    tealeaf = pd.read_csv(
        path / tealeaf_subdir / "ec_block_glmm.tsv", sep="\t"
    )
    tealeaf = tealeaf.loc[
        tealeaf.method.eq(tealeaf_fit_method)
        & tealeaf.median_gene_umis.ge(min_median_gene_umis)
        & tealeaf.null_converged
        & tealeaf.alternative_converged,
        ["block_id", "gene_id", "p_value"],
    ].rename(columns={"block_id": "feature_id"})
    tealeaf.insert(0, "method", tealeaf_label)
    tealeaf = aggregate_feature_pvalues(tealeaf)
    features = pd.concat(
        [external[["method", "feature_id", "gene_id", "p_value"]], tealeaf],
        ignore_index=True,
    )
    return aggregate_gene_pvalues(features)


def load_pairwise_fold(
    path: Path,
    leaf_map: pd.DataFrame,
    tealeaf_subdir: str,
    tealeaf_label: str,
    tealeaf_fit_method: str,
    min_median_gene_umis: float,
) -> pd.DataFrame:
    external = pd.read_csv(
        path / "comparison/all_tests.tsv.gz", sep="\t", low_memory=False
    )
    external = external.loc[external.effect.eq("cell_type")].copy()
    external.loc[external.method.eq("LeafCutter"), "gene_id"] = external.loc[
        external.method.eq("LeafCutter"), "feature_id"
    ].map(leaf_map.set_index("feature_id").gene_id)
    tealeaf = pd.read_csv(
        path / tealeaf_subdir / "ec_block_glmm.tsv", sep="\t"
    )
    tealeaf = tealeaf.loc[
        tealeaf.method.eq(tealeaf_fit_method)
        & tealeaf.median_gene_umis.ge(min_median_gene_umis)
        & tealeaf.null_converged
        & tealeaf.alternative_converged,
        ["gene_id", "level_a", "level_b", "p_value"],
    ].copy()
    tealeaf.insert(0, "method", tealeaf_label)
    return pd.concat(
        [
            external[["method", "gene_id", "level_a", "level_b", "p_value"]],
            tealeaf,
        ],
        ignore_index=True,
    )


def main():
    args = parse_args()
    if len(args.fold_dir) != 2:
        raise ValueError("provide exactly two --fold-dir values")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.leafcutter_map.exists():
        leaf_map = pd.read_csv(args.leafcutter_map, sep="\t")
    else:
        leaf_map = leafcutter_cluster_gene_map(
            args.junction_table, args.leafcutter_counts
        )
        args.leafcutter_map.parent.mkdir(parents=True, exist_ok=True)
        leaf_map.to_csv(args.leafcutter_map, sep="\t", index=False)
    if args.match_pairs:
        folds = [
            load_pairwise_fold(
                path,
                leaf_map,
                args.tealeaf_subdir,
                args.tealeaf_label,
                args.tealeaf_fit_method,
                args.min_median_gene_umis,
            )
            for path in args.fold_dir
        ]
        metrics, topk, genes = shared_pair_gene_reproducibility(
            folds, reference_method=args.tealeaf_label
        )
    else:
        folds = [
            load_fold(
                path,
                leaf_map,
                args.tealeaf_subdir,
                args.tealeaf_label,
                args.tealeaf_fit_method,
                args.min_median_gene_umis,
            )
            for path in args.fold_dir
        ]
        metrics, topk, genes = shared_gene_reproducibility(
            folds, reference_method=args.tealeaf_label
        )
    metrics.to_csv(args.output_dir / "gene_metrics.tsv", sep="\t", index=False)
    topk.to_csv(args.output_dir / "topk_overlap.tsv", sep="\t", index=False)
    genes.to_csv(
        args.output_dir / "gene_reproducibility.tsv.gz", sep="\t", index=False
    )
    try:
        from plotnine import (
            aes,
            facet_wrap,
            geom_line,
            geom_point,
            ggplot,
            labs,
            theme_bw,
        )
    except ImportError:
        warnings.warn("plotnine is unavailable; skipping top-K overlap plot")
    else:
        plot = (
            ggplot(topk, aes("k", "overlap_fraction", color="method"))
            + geom_line()
            + geom_point()
            + facet_wrap("comparison", scales="free_x")
            + theme_bw()
            + labs(x="Top genes per fold", y="Cross-fold overlap fraction")
        )
        plot.save(args.output_dir / "topk_overlap.pdf", width=6.5, height=4.2)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
