#!/usr/bin/env python3
"""Independent Gaussian biological-null audit of omnibus sensitivities."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import softmax

from extra_scripts.audit_path_reporting_omnibus import omnibus_statistics
from extra_scripts.summarize_path_reporting_omnibus import calibrate_omnibus
from extra_scripts.audit_wild_omnibus import wild_tests
from tealeaf.sc.differential import helmert_basis, paired_mean_test
from tealeaf.sc.ds_benchmark import simes_pvalue, cauchy_pvalue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-tests", type=int, default=600)
    args = parser.parse_args()
    rng = np.random.default_rng(8721)
    labels = np.tile(np.arange(4), 12)
    subjects = np.repeat(np.arange(12), 4).astype(str)
    observed, null, combination_rows = [], [], []
    for number in range(args.n_tests):
        dimension = 1 if number % 2 else 3
        heteroscedastic = number % 4 >= 2
        is_null = number % 6 != 0
        covariance = .7 * np.ones((dimension, dimension)) + .3 * np.eye(dimension)
        noise = rng.multivariate_normal(np.zeros(dimension), covariance, len(labels))
        if heteroscedastic:
            noise *= np.array([.5, 1., 2., 4.])[labels, None]
        values = rng.normal(size=(12, dimension))[subjects.astype(int)] + noise
        if not is_null:
            values += labels[:, None] * 1.2
        proportions = softmax(values @ helmert_basis(dimension + 1).T, axis=1)
        header = {"test_id": f"sim_{number}", "n_subjects": 12, "converged": True, "is_null": is_null, "heteroscedastic": heteroscedastic, "dimension": dimension}
        by_subject = values.reshape(12, 4, dimension)
        pair_pvalues = np.array([paired_mean_test(by_subject[:, b] - by_subject[:, a])["p_value"] for a in range(4) for b in range(a + 1, 4)])
        for strategy, pvalue in (("Bonferroni paired omnibus", min(1., len(pair_pvalues) * pair_pvalues.min())), ("Simes paired omnibus", simes_pvalue(pair_pvalues)), ("ACAT paired omnibus", cauchy_pvalue(pair_pvalues))):
            combination_rows.append({**header, "strategy": strategy, "p_value": pvalue, "raw_p_value": pvalue})
        stats = omnibus_statistics({32.: values, 1.: values}, proportions, labels, subjects)
        observed.extend({**header, "strategy": name, **{key: result[key] for key in ("p_value", "statistic", "degrees_of_freedom")}} for name, result in stats.items())
        wild_observed, wild_null = wild_tests(values, labels, subjects, rng)
        observed.append({**header, "strategy": "CR2 maximum-coordinate wild ILR A32", **wild_observed})
        null.extend({**header, "strategy": "CR2 maximum-coordinate wild ILR A32", **record} for record in wild_null)
        for replicate in range(64):
            permuted = labels.reshape(12, 4).copy()
            for row in permuted:
                rng.shuffle(row)
            stats = omnibus_statistics({32.: values, 1.: values}, proportions, permuted.ravel(), subjects)
            null.extend({**header, "strategy": name, "replicate": replicate, **{key: result[key] for key in ("p_value", "statistic", "degrees_of_freedom")}} for name, result in stats.items())
    table = pd.DataFrame(observed)
    calibrated, held, held_summary = calibrate_omnibus(table, pd.DataFrame(null))
    calibrated = pd.concat([calibrated, pd.DataFrame(combination_rows)], ignore_index=True)
    summaries = []
    for (strategy, hetero, dimension), local in calibrated.groupby(["strategy", "heteroscedastic", "dimension"]):
        for truth, subset in local.groupby("is_null"):
            summaries.append({"strategy": strategy, "heteroscedastic": hetero, "dimension": dimension, "is_null": truth, "n_tests": len(subset), "analytic_reject_0_05": subset.raw_p_value.le(.05).mean(), "calibrated_reject_0_05": subset.p_value.le(.05).mean(), "analytic_reject_0_01": subset.raw_p_value.le(.01).mean(), "calibrated_reject_0_01": subset.p_value.le(.01).mean()})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(args.output_dir / "gaussian_null_summary.tsv", sep="\t", index=False)
    held_summary.to_csv(args.output_dir / "gaussian_held_permutations.tsv", sep="\t", index=False)
    calibrated.to_csv(args.output_dir / "gaussian_tests.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
