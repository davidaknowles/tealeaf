#!/usr/bin/env python3
"""Select one path-smoothing strength with a calibrated two-groups model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.special


def beta_uniform_mixture(table: pd.DataFrame) -> dict[str, float]:
    """Fit a gene-equal beta-uniform mixture to calibrated p-values."""
    required = {"gene_id", "p_value"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"p-value table is missing columns: {sorted(missing)}")
    retained = table.loc[
        table.p_value.notna() & table.p_value.between(0, 1),
        ["gene_id", "p_value"],
    ]
    if retained.empty:
        raise ValueError("no finite p-values were supplied")
    p_values = np.clip(retained.p_value.to_numpy(dtype=float), 1e-12, 1 - 1e-12)
    gene_codes, genes = pd.factorize(retained.gene_id.astype(str))
    gene_counts = np.bincount(gene_codes)

    def objective(parameters):
        null_weight = scipy.special.expit(parameters[0])
        alternative_shape = scipy.special.expit(parameters[1])
        density = (
            null_weight
            + (1 - null_weight)
            * alternative_shape
            * p_values ** (alternative_shape - 1)
        )
        log_density = np.log(np.maximum(density, 1e-300))
        per_gene = np.bincount(gene_codes, weights=log_density) / gene_counts
        return -float(per_gene.mean())

    candidates = []
    for initial in ((0.0, -1.0), (2.0, -2.0), (-2.0, -0.5)):
        fit = scipy.optimize.minimize(
            objective,
            initial,
            method="Nelder-Mead",
            options={"maxiter": 2000, "xatol": 1e-9, "fatol": 1e-10},
        )
        candidates.append(fit)
    fit = min(candidates, key=lambda value: value.fun)
    return {
        "mean_log_evidence": -float(fit.fun),
        "null_weight": float(scipy.special.expit(fit.x[0])),
        "alternative_shape": float(scipy.special.expit(fit.x[1])),
        "genes": int(len(genes)),
        "tests": int(len(retained)),
    }


def bh_count(p_values, level=0.05):
    """Return the number of Benjamini-Hochberg rejections."""
    ordered = np.sort(np.asarray(p_values, dtype=float))
    if not len(ordered):
        return 0
    accepted = ordered <= float(level) * np.arange(1, len(ordered) + 1) / len(ordered)
    return int(np.flatnonzero(accepted)[-1] + 1) if accepted.any() else 0


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result",
        nargs="+",
        action="append",
        metavar=("STRENGTH", "DIRECTORY"),
        required=True,
        help="A strength followed by one or more independent assessment directories.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    observed_by_strength = {}
    null_by_strength = {}
    metadata = []
    for specification in args.result:
        if len(specification) < 2:
            raise ValueError("each --result needs a strength and at least one directory")
        strength = float(specification[0])
        observed_parts = []
        null_parts = []
        for assessment, value in enumerate(specification[1:]):
            directory = Path(value)
            observed = pd.read_csv(directory / "paired_path.tsv", sep="\t")
            observed = observed.loc[observed.converged.astype(bool)].copy()
            observed_parts.append(observed)
            null = pd.read_csv(directory / "paired_path_null.tsv.gz", sep="\t")
            null = null.merge(
                observed[["test_id", "gene_id"]], on="test_id", how="inner"
            )
            null["assessment"] = assessment
            null_parts.append(null)
            metadata.append({
                "strength": strength,
                "assessment": assessment,
                "path_prior_center": str(observed.path_prior_center.iloc[0]),
                "path_prior_family": str(
                    observed.path_prior_family.iloc[0]
                    if "path_prior_family" in observed
                    else "dirichlet"
                ),
                "path_pseudocount_scaling": str(
                    observed.path_pseudocount_scaling.iloc[0]
                ),
            })
        observed_by_strength[strength] = pd.concat(
            observed_parts, ignore_index=True
        )
        null_by_strength[strength] = pd.concat(null_parts, ignore_index=True)
    metadata = pd.DataFrame(metadata)
    for column in (
        "path_prior_center",
        "path_prior_family",
        "path_pseudocount_scaling",
    ):
        values = metadata[column].unique()
        if len(values) != 1:
            raise ValueError(f"all inputs must use one {column}")
    if metadata.path_prior_center.iloc[0] != "uniform":
        raise ValueError("testing-level EB requires uniform-centered inputs")

    curve = []
    for strength, table in observed_by_strength.items():
        curve.append({"strength": strength, **beta_uniform_mixture(table)})
    curve = pd.DataFrame(curve).sort_values("strength").reset_index(drop=True)
    selected = curve.sort_values(
        ["mean_log_evidence", "strength"], ascending=[False, True]
    ).iloc[0]

    replicate_sets = [
        set(table.replicate.unique()) for table in null_by_strength.values()
    ]
    common_replicates = sorted(set.intersection(*replicate_sets))
    null_rows = []
    assessments = sorted(metadata.assessment.unique())
    for replicate in common_replicates:
        candidate_scores = {}
        for strength, table in null_by_strength.items():
            local = table.loc[table.replicate.eq(replicate)]
            candidate_scores[strength] = beta_uniform_mixture(local)[
                "mean_log_evidence"
            ]
        strength = max(
            candidate_scores,
            key=lambda value: (candidate_scores[value], -value),
        )
        selected_null = null_by_strength[strength]
        row = {"replicate": int(replicate), "selected_strength": strength}
        for assessment in assessments:
            p_values = selected_null.loc[
                selected_null.replicate.eq(replicate)
                & selected_null.assessment.eq(assessment),
                "p_value",
            ]
            row[f"assessment_{assessment}_rejection_0.05"] = float(
                p_values.lt(0.05).mean()
            )
            row[f"assessment_{assessment}_bh_0.05"] = bh_count(p_values)
        null_rows.append(row)
    null_audit = pd.DataFrame(null_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    curve.to_csv(args.output_dir / "evidence_curve.tsv", sep="\t", index=False)
    null_audit.to_csv(args.output_dir / "null_selection.tsv", sep="\t", index=False)
    result = {
        "selection_method": "gene_equal_beta_uniform_mixture",
        "selection_scope": "global",
        "strength": float(selected.strength),
        "mean_log_evidence": float(selected.mean_log_evidence),
        "null_weight": float(selected.null_weight),
        "alternative_shape": float(selected.alternative_shape),
        "path_prior_center": str(metadata.path_prior_center.iloc[0]),
        "path_prior_family": str(metadata.path_prior_family.iloc[0]),
        "path_pseudocount_scaling": str(
            metadata.path_pseudocount_scaling.iloc[0]
        ),
        "null_replicates": int(len(null_audit)),
        "null_families_with_bh": int(
            null_audit.filter(like="_bh_0.05").gt(0).any(axis=1).sum()
        ),
    }
    (args.output_dir / "selection.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
