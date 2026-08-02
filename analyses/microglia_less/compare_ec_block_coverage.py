#!/usr/bin/env python3
"""Compare calibrated EC-block tests across gene-coverage rules."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from extra_scripts.run_differential_splicing import benjamini_hochberg


METHODS = {
    "multinomial": "multinomial_full",
    "logistic_normal": "multinomial_noise_full",
    "dirichlet_multinomial": "dirichlet_multinomial_full",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--differential-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def load_table(root, coverage, method):
    path = (
        root
        / f"ec_block_direct{coverage}_{method}_bootstrap5_merged"
        / "ec_block_glmm.tsv"
    )
    return pd.read_csv(path, sep="\t").set_index("block_id")


def eligible(table):
    return (
        table["null_converged"]
        & table["alternative_converged"]
        & table["p_value"].notna()
    )


def compare_method(strict, permissive, name):
    common = strict.index.intersection(permissive.index)
    common_set = set(common)
    strict_significant = set(strict.index[strict["fdr"] <= 0.05])
    permissive_significant = set(
        permissive.index[permissive["fdr"] <= 0.05]
    )
    strict_common = strict_significant & common_set
    permissive_common = permissive_significant & common_set
    paired = (
        strict.loc[common, "p_value"].notna()
        & permissive.loc[common, "p_value"].notna()
    )
    paired_blocks = common[paired]
    permissive_shared = permissive.loc[common]
    permissive_shared_eligible = eligible(permissive_shared)
    shared_fdr = benjamini_hochberg(
        permissive_shared.loc[permissive_shared_eligible, "p_value"].to_numpy()
    )
    return {
        "method": name,
        "strict_tests": len(strict),
        "permissive_tests": len(permissive),
        "shared_tests": len(common),
        "strict_converged": int(eligible(strict).sum()),
        "permissive_converged": int(eligible(permissive).sum()),
        "strict_nominal": int(
            (strict.loc[eligible(strict), "p_value"] <= 0.05).sum()
        ),
        "permissive_nominal": int(
            (permissive.loc[eligible(permissive), "p_value"] <= 0.05).sum()
        ),
        "strict_fdr": len(strict_significant),
        "permissive_fdr": len(permissive_significant),
        "permissive_shared_family_fdr": int((shared_fdr <= 0.05).sum()),
        "significant_both": len(strict_common & permissive_common),
        "strict_only_shared": len(strict_common - permissive_common),
        "permissive_only_shared": len(permissive_common - strict_common),
        "permissive_only_universe": len(
            permissive_significant - common_set
        ),
        "p_value_spearman": spearmanr(
            strict.loc[paired_blocks, "p_value"],
            permissive.loc[paired_blocks, "p_value"],
        ).statistic,
        "statistic_spearman": spearmanr(
            strict.loc[paired_blocks, "statistic"],
            permissive.loc[paired_blocks, "statistic"],
        ).statistic,
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    discoveries = {"strict": {}, "permissive": {}}
    for name, method in METHODS.items():
        strict = load_table(args.differential_root, "25_min30", method)
        permissive = load_table(args.differential_root, "10", method)
        rows.append(compare_method(strict, permissive, name))
        discoveries["strict"][name] = set(
            strict.index[strict["fdr"] <= 0.05]
        )
        discoveries["permissive"][name] = set(
            permissive.index[permissive["fdr"] <= 0.05]
        )
    pd.DataFrame(rows).to_csv(
        args.output_dir / "coverage_comparison.tsv", sep="\t", index=False
    )
    overlap_rows = []
    for coverage, sets in discoveries.items():
        mn = sets["multinomial"]
        ln = sets["logistic_normal"]
        dm = sets["dirichlet_multinomial"]
        overlap_rows.append({
            "coverage": coverage,
            "multinomial_logistic_normal": len(mn & ln),
            "multinomial_only": len(mn - ln),
            "logistic_normal_only": len(ln - mn),
            "all_three": len(mn & ln & dm),
        })
    pd.DataFrame(overlap_rows).to_csv(
        args.output_dir / "method_overlap.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
