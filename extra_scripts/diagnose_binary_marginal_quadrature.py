"""Describe integration failures and compare matched numerical backends."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--prior-pilot", type=Path, required=True)
    parser.add_argument("--shared-pilot", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    table = pd.read_csv(args.diagnostics / "tests.tsv.gz", sep="\t")
    table["converged"] = table.converged.astype(str).str.lower().eq("true")
    tolerance = 1e-3
    table["subject_rule_failed"] = table.subject_only_quadrature_error.gt(tolerance)
    table["path_rule_failed"] = table.path_only_quadrature_error.gt(tolerance)
    table["minimum_path_support_stratum"] = pd.cut(table.true_median_subject_minimum_path_mean, [0, .01, .05, .2, .5], include_lowest=True).astype(str)
    def aggregate(keys):
        return table.groupby(keys, observed=True).agg(n_requested=("p_value", "size"), n_converged=("converged", "sum"), native_reject_0_05=("p_value", lambda values: values.le(.05).sum()), n_subject_rule_failed=("subject_rule_failed", "sum"), n_path_rule_failed=("path_rule_failed", "sum"), median_true_minimum_path_mean=("true_median_subject_minimum_path_mean", "median"), median_proposal_effective_depth=("proposal_effective_depth_median", "median"), median_gene_depth=("gene_ec_depth_median", "median"), median_true_block_mass=("true_block_mass_median", "median"), max_subject_quadrature_error=("subject_only_quadrature_error", "max"), max_path_quadrature_error=("path_only_quadrature_error", "max")).reset_index()
    summary = aggregate(["scenario"])
    support = aggregate(["scenario", "minimum_path_support_stratum"])
    prior = pd.read_csv(args.prior_pilot / "tests.tsv.gz", sep="\t")
    shared = pd.read_csv(args.shared_pilot / "tests.tsv.gz", sep="\t")
    keys = ["test_id", "scenario", "draw"]
    if prior.duplicated(keys).any() or shared.duplicated(keys).any() or set(map(tuple, prior[keys].values)) != set(map(tuple, shared[keys].values)):
        raise ValueError("pilot families differ or contain duplicates")
    matched = prior.merge(shared, on=keys, suffixes=("_prior", "_shared"), validate="one_to_one")
    matched["runtime_ratio"] = matched.runtime_seconds_prior / matched.runtime_seconds_shared
    for column in ("p_value", "statistic"):
        matched[f"absolute_{column}_difference"] = np.abs(matched[f"{column}_prior"] - matched[f"{column}_shared"])
    matched["max_absolute_usage_difference"] = [np.max(np.abs(np.asarray(json.loads(left)) - np.asarray(json.loads(right)))) for left, right in zip(matched.standardized_means_prior, matched.standardized_means_shared)]
    timing = matched.groupby("scenario").agg(n_matched=("test_id", "size"), median_runtime_ratio=("runtime_ratio", "median"), max_absolute_p_difference=("absolute_p_value_difference", "max"), max_absolute_statistic_difference=("absolute_statistic_difference", "max"), max_absolute_usage_difference=("max_absolute_usage_difference", "max")).reset_index()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, result in (("quadrature_summary", summary), ("support_strata", support), ("matched_numerical_pilots", matched), ("timing_summary", timing)):
        result.to_csv(args.output_dir / f"{name}.tsv", sep="\t", index=False, na_rep="NA")
    print(summary.to_string(index=False))
    print(support.to_string(index=False))
    print(timing.to_string(index=False))


if __name__ == "__main__":
    main()
