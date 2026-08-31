"""Summarize calibrated-result stability across Monte Carlo null budgets."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats


def parse_result(value: str) -> tuple[str, int, Path]:
    analysis, replicates, path = value.split("=", 2)
    return analysis, int(replicates), Path(path)


def parse_reproducibility(value: str) -> tuple[int, Path]:
    replicates, path = value.split("=", 1)
    return int(replicates), Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result",
        action="append",
        required=True,
        type=parse_result,
        metavar="ANALYSIS=R=DIR",
    )
    parser.add_argument(
        "--reproducibility",
        action="append",
        type=parse_reproducibility,
        metavar="R=DIR",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def summarize_result(analysis: str, replicates: int, path: Path) -> dict:
    summary = pd.read_csv(path / "summary.tsv", sep="\t").iloc[0].to_dict()
    families = pd.read_csv(path / "null_families.tsv", sep="\t")
    families_with_bh = int(families["bh_0.05"].gt(0).sum())
    return {
        "analysis": analysis,
        "null_replicates": replicates,
        **summary,
        "null_families_with_bh": families_with_bh,
        "null_family_bh_rate": families_with_bh / len(families),
    }


def compare_tables(analysis: str, replicates: int, baseline: Path, path: Path) -> dict:
    columns = ["test_id", "p_value", "fdr", "converged"]
    left = pd.read_csv(baseline / "paired_path.tsv", sep="\t", usecols=columns)
    right = pd.read_csv(path / "paired_path.tsv", sep="\t", usecols=columns)
    merged = left.merge(
        right,
        on="test_id",
        suffixes=("_32", "_candidate"),
        validate="one_to_one",
    )
    finite = np.isfinite(merged.p_value_32) & np.isfinite(merged.p_value_candidate)
    log_left = -np.log10(np.maximum(merged.loc[finite, "p_value_32"], np.finfo(float).tiny))
    log_right = -np.log10(np.maximum(merged.loc[finite, "p_value_candidate"], np.finfo(float).tiny))
    bh_left = set(merged.loc[merged.fdr_32.le(0.05), "test_id"])
    bh_right = set(merged.loc[merged.fdr_candidate.le(0.05), "test_id"])
    union = bh_left | bh_right
    return {
        "analysis": analysis,
        "null_replicates": replicates,
        "shared_tests": len(merged),
        "spearman_log10_p": scipy.stats.spearmanr(log_left, log_right).statistic,
        "median_absolute_p_change": float(
            np.median(
                np.abs(
                    merged.loc[finite, "p_value_32"]
                    - merged.loc[finite, "p_value_candidate"]
                )
            )
        ),
        "bh_32": len(bh_left),
        "bh_candidate": len(bh_right),
        "bh_shared": len(bh_left & bh_right),
        "bh_gained": len(bh_right - bh_left),
        "bh_lost": len(bh_left - bh_right),
        "bh_jaccard": len(bh_left & bh_right) / len(union) if union else 1.0,
    }


def main() -> None:
    args = parse_args()
    grouped: dict[str, dict[int, Path]] = {}
    summaries = []
    for analysis, replicates, path in args.result:
        grouped.setdefault(analysis, {})[replicates] = path
        summaries.append(summarize_result(analysis, replicates, path))
    comparisons = []
    for analysis, results in grouped.items():
        if 32 not in results:
            raise ValueError(f"{analysis} has no R=32 baseline")
        for replicates, path in sorted(results.items()):
            if replicates != 32:
                comparisons.append(compare_tables(analysis, replicates, results[32], path))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).sort_values(
        ["analysis", "null_replicates"]
    ).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(comparisons).sort_values(
        ["analysis", "null_replicates"]
    ).to_csv(args.output_dir / "stability.tsv", sep="\t", index=False)
    if args.reproducibility:
        reproducibility = []
        for replicates, path in args.reproducibility:
            table = pd.read_csv(path / "gene_metrics.tsv", sep="\t")
            table.insert(0, "null_replicates", replicates)
            reproducibility.append(table)
        pd.concat(reproducibility, ignore_index=True).sort_values(
            ["comparison", "method", "null_replicates"]
        ).to_csv(
            args.output_dir / "reproducibility.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
