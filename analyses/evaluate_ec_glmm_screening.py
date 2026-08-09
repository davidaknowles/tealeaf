#!/usr/bin/env python3
"""Evaluate conservative path-ranking screens for exact EC-GLMM tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ec-glmm", required=True, type=Path)
    parser.add_argument("--path-screen", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    exact = pd.read_csv(args.ec_glmm, sep="\t").drop_duplicates("test_id")
    screen = pd.read_csv(args.path_screen, sep="\t").drop_duplicates("test_id")
    exact_p = exact.lrt_p_value.where(
        exact.null_converged.astype(bool) & exact.alternative_converged.astype(bool),
        1.0,
    ).fillna(1.0)
    exact = exact.assign(exact_p_value=exact_p)
    full_calls = set(
        exact.loc[benjamini_hochberg(exact_p.to_numpy()) <= 0.05, "test_id"]
    )
    screen_ids = set(screen.test_id)
    mandatory = set(exact.test_id) - screen_ids
    scores = {
        "path_lrt": ("p_value", True),
        "score_mixture_bf": ("mixture_log_bayes_factor", False),
        "bic_bf": ("bic_log_bayes_factor", False),
    }
    rows = []
    for method, (column, ascending) in scores.items():
        ranked = screen.sort_values(column, ascending=ascending).test_id.to_numpy()
        for fraction in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
            retained_count = int(np.ceil(fraction * len(ranked)))
            retained = mandatory | set(ranked[:retained_count])
            screened_p = exact.exact_p_value.where(exact.test_id.isin(retained), 1.0)
            calls = set(
                exact.loc[benjamini_hochberg(screened_p.to_numpy()) <= 0.05, "test_id"]
            )
            rows.append(
                {
                    "screen": method,
                    "screen_fraction": fraction,
                    "screen_eligible_tests": len(ranked),
                    "mandatory_exact_tests": len(mandatory),
                    "exact_tests_run": len(retained),
                    "exact_test_fraction": len(retained) / len(exact),
                    "full_ec_glmm_discoveries": len(full_calls),
                    "screened_ec_glmm_discoveries": len(calls),
                    "discovery_recall": (
                        len(calls & full_calls) / len(full_calls)
                        if full_calls
                        else np.nan
                    ),
                    "screened_calls_not_in_full_analysis": len(calls - full_calls),
                }
            )
    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
