#!/usr/bin/env python3
"""Compare cell-resolved and pseudobulk splice-path inference."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats

from tealeaf.sc import differential


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk-estimates", required=True, type=Path)
    parser.add_argument("--single-cell-estimates", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-condition-replicates", type=int, default=3)
    parser.add_argument("--effect-size", type=float, default=0.5)
    return parser.parse_args()


def read_estimates(path, method, block_filter=None):
    records = {}
    seen_blocks = set()
    with gzip.open(path, "rt") as source:
        for line in source:
            row = json.loads(line)
            if block_filter is not None and row["block_id"] not in block_filter:
                continue
            seen_blocks.add(row["block_id"])
            if method == "pseudobulk":
                reliable = row["fit_converged"] and row["conditional_reliable"]
                covariance = row["conditional_covariance"]
            else:
                reliable = row["fit_converged"] and row["reliable"]
                covariance = row["covariance"]
            if not reliable:
                continue
            records[(row["block_id"], row["group"])] = {
                "value": np.asarray(row["path_logratios"], dtype=float),
                "covariance": np.asarray(covariance, dtype=float),
                "cluster": row["cluster"],
                "condition": row["condition"],
            }
    return records, seen_blocks


def fit_matched_tests(pseudobulk, single_cell, min_replicates, effect_size):
    grouped = defaultdict(list)
    for key in sorted(set(pseudobulk) & set(single_cell)):
        block_id, group = key
        aggregate = pseudobulk[key]
        resolved = single_cell[key]
        if (
            aggregate["cluster"] != resolved["cluster"]
            or aggregate["condition"] != resolved["condition"]
        ):
            raise ValueError(f"metadata disagree for {block_id}, {group}")
        grouped[(block_id, aggregate["cluster"])].append(
            (group, aggregate, resolved)
        )

    rows = []
    for (block_id, cluster), records in grouped.items():
        conditions = sorted({
            record[1]["condition"] for record in records
        })
        counts = {
            condition: sum(
                record[1]["condition"] == condition for record in records
            )
            for condition in conditions
        }
        conditions = [
            condition
            for condition in conditions
            if counts[condition] >= int(min_replicates)
        ]
        records = [
            record
            for record in records
            if record[1]["condition"] in conditions
        ]
        if len(conditions) < 2 or len(records) <= len(conditions):
            continue
        condition_index = {
            condition: index for index, condition in enumerate(conditions)
        }
        design = np.zeros((len(records), len(conditions)), dtype=float)
        design[:, 0] = 1.0
        for row, record in enumerate(records):
            index = condition_index[record[1]["condition"]]
            if index:
                design[row, index] = 1.0
        values = {}
        covariances = {}
        fits = {}
        for method, position in (("pseudobulk", 1), ("single_cell", 2)):
            values[method] = np.asarray([
                record[position]["value"] for record in records
            ])
            covariances[method] = np.asarray([
                record[position]["covariance"] for record in records
            ])
            fits[method] = differential.multivariate_gls_test(
                values[method],
                covariances[method],
                design,
                tested_columns=range(1, len(conditions)),
            )
        row = {
            "block_id": block_id,
            "cluster": cluster,
            "conditions": ",".join(conditions),
            "n_samples": len(records),
            "logratio_dimension": values["pseudobulk"].shape[1],
        }
        for method, fit in fits.items():
            row[f"{method}_p_value"] = fit["p_value"]
            row[f"{method}_biological_variance"] = fit[
                "biological_variance"
            ]
        if len(conditions) == 2 and row["logratio_dimension"] == 1:
            cutoff = scipy.stats.chi2.ppf(0.95, 1)
            for method, fit in fits.items():
                variance = fit["coefficient_covariance"][1, 1]
                noncentrality = float(effect_size) ** 2 / variance
                row[f"{method}_coefficient_se"] = np.sqrt(variance)
                row[f"{method}_power"] = scipy.stats.ncx2.sf(
                    cutoff, 1, noncentrality
                )
        rows.append(row)
    return pd.DataFrame(rows)


def finite_summary(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "q10": float(np.quantile(values, 0.1)),
        "q90": float(np.quantile(values, 0.9)),
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    single_cell, target_blocks = read_estimates(
        args.single_cell_estimates, "single_cell"
    )
    pseudobulk, _ = read_estimates(
        args.pseudobulk_estimates,
        "pseudobulk",
        block_filter=target_blocks,
    )
    common = sorted(set(pseudobulk) & set(single_cell))
    covariance_ratios = [
        np.trace(single_cell[key]["covariance"])
        / np.trace(pseudobulk[key]["covariance"])
        for key in common
    ]
    response_differences = [
        np.linalg.norm(
            pseudobulk[key]["value"] - single_cell[key]["value"]
        )
        / np.sqrt(len(pseudobulk[key]["value"]))
        for key in common
    ]
    tests = fit_matched_tests(
        pseudobulk,
        single_cell,
        args.min_condition_replicates,
        args.effect_size,
    )
    tests.to_csv(
        args.output_dir / "matched_condition_tests.tsv",
        sep="\t",
        index=False,
    )
    two_condition = tests.dropna(
        subset=["pseudobulk_power", "single_cell_power"]
    )
    summary = {
        "pseudobulk_reliable_estimates": len(pseudobulk),
        "single_cell_reliable_estimates": len(single_cell),
        "common_reliable_estimates": len(common),
        "single_cell_over_pseudobulk_covariance_trace": finite_summary(
            covariance_ratios
        ),
        "response_rms_difference": finite_summary(response_differences),
        "matched_condition_tests": len(tests),
        "pseudobulk_nominal_p_lt_0.05": int(
            (tests["pseudobulk_p_value"] < 0.05).sum()
        ),
        "single_cell_nominal_p_lt_0.05": int(
            (tests["single_cell_p_value"] < 0.05).sum()
        ),
        "p_value_spearman": float(scipy.stats.spearmanr(
            tests["pseudobulk_p_value"],
            tests["single_cell_p_value"],
        ).statistic),
        "two_condition_two_path_tests": len(two_condition),
        "effect_size": args.effect_size,
        "pseudobulk_mean_power": float(
            two_condition["pseudobulk_power"].mean()
        ),
        "single_cell_mean_power": float(
            two_condition["single_cell_power"].mean()
        ),
        "median_single_cell_over_pseudobulk_coefficient_se": float(
            np.median(
                two_condition["single_cell_coefficient_se"]
                / two_condition["pseudobulk_coefficient_se"]
            )
        ),
    }
    (args.output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
