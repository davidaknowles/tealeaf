#!/usr/bin/env python3
"""Select plain or overdispersed EC GLMMs by the null-model BIC."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.junction_benchmark import benjamini_hochberg


def select_hybrid(plain: pd.DataFrame, noise: pd.DataFrame) -> pd.DataFrame:
    """Return test rows from the null-BIC-selected observation model."""
    if plain.test_id.duplicated().any() or noise.test_id.duplicated().any():
        raise ValueError("each input must have one row per test_id")
    common = plain.set_index("test_id").index.intersection(
        noise.set_index("test_id").index, sort=False
    )
    plain = plain.set_index("test_id").loc[common].copy()
    noise = noise.set_index("test_id").loc[common].copy()
    independent_samples = np.maximum(noise.n_samples.to_numpy(dtype=float) / 2, 2)
    bic_gain = (
        2
        * (
            plain.null_objective.to_numpy(dtype=float)
            - noise.null_objective.to_numpy(dtype=float)
        )
        - np.log(independent_samples)
    )
    use_noise = bic_gain > 0
    result = plain.copy()
    shared_columns = result.columns.intersection(noise.columns)
    result.loc[use_noise, shared_columns] = noise.loc[use_noise, shared_columns]
    result["method"] = "laplace_multinomial_hybrid_bic"
    result["selected_observation_model"] = np.where(
        use_noise, "logistic_normal_noise", "multinomial"
    )
    result["null_noise_bic_gain"] = bic_gain
    joint_convergence = (
        plain.null_converged
        & plain.alternative_converged
        & noise.null_converged
        & noise.alternative_converged
    )
    result["null_converged"] = joint_convergence
    result["alternative_converged"] = joint_convergence
    result["fdr"] = np.nan
    eligible = joint_convergence & result.p_value.notna()
    result.loc[eligible, "fdr"] = benjamini_hochberg(
        result.loc[eligible, "p_value"].to_numpy(dtype=float)
    )
    return result.reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plain", required=True, type=Path)
    parser.add_argument("--noise", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    result = select_hybrid(
        pd.read_csv(args.plain, sep="\t"),
        pd.read_csv(args.noise, sep="\t"),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False)
    result.loc[result.fdr <= 0.05].to_csv(
        args.output_dir / "significant_ec_block_glmm.tsv", sep="\t", index=False
    )
    pd.DataFrame([{
        "method": "laplace_multinomial_hybrid_bic",
        "events": len(result),
        "convergence": float(
            (result.null_converged & result.alternative_converged).mean()
        ),
        "noise_selected": int(
            result.selected_observation_model.eq("logistic_normal_noise").sum()
        ),
        "fdr_0.05": int((result.fdr <= 0.05).sum()),
    }]).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
