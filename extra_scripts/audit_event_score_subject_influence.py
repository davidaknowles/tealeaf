#!/usr/bin/env python3
"""Audit complete binary score archives, never replace production tests."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.reassess_event_score_archive import read_table
from tealeaf.sc.path_score_mixed import MODEL_VERSION, binary_score_components_from_records, binary_score_components_to_proportions, binary_score_subject_influence, mixed_score_test


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", type=Path, required=True)
    parser.add_argument("--merged-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=64)
    parser.add_argument("--threshold", type=float, default=1e-5)
    parser.add_argument("--compare-proportions", action="store_true", help="Additional conditional score-coordinate sensitivity, not regenerated sign-null inference.")
    args = parser.parse_args()
    if args.shard_count < 1 or not 0 < args.threshold < .1:
        raise ValueError("positive shard count and raw threshold below .1 required")
    observed_path = args.merged_dir / "paired_path.tsv"
    observed = read_table(observed_path)
    if observed.test_id.duplicated().any():
        raise ValueError("unique complete-cohort hypothesis identities required")
    complete = observed.complete_subject_fits.astype(str).str.lower().eq("true")
    strong = observed.loc[complete & observed.raw_p_value.le(args.threshold)].copy()
    strong["panel"] = "small observed raw p-value"
    weak = observed.loc[complete & observed.raw_p_value.between(.1, .9)].copy()
    rng, controls = np.random.default_rng(47852), []
    for size, group in strong.groupby("n_subjects", sort=True):
        pool = weak.loc[weak.n_subjects.eq(size)]
        if len(pool) < len(group):
            raise ValueError("not enough distinct weak controls with matching subject count")
        controls.append(pool.sample(n=len(group), random_state=int(rng.integers(2**31))))
    if not controls:
        raise ValueError("no strong-tail parents to assess")
    control = pd.concat(controls)
    control["panel"] = "subject-count-matched weak control"
    selected = pd.concat([strong, control], ignore_index=True).set_index("test_id", drop=False)
    requested, seen, trials, subjects, declared = set(selected.index), set(), [], [], 0
    for shard_index in range(args.shard_count):
        shard = args.cohort_root / f"shard_{shard_index}"
        summary = json.loads((shard / "summary.json").read_text())
        settings = json.loads((shard / "settings.json").read_text())
        if summary["completed"] + summary["failures"] != summary["tests_in_shard"] or settings["model_version"] != MODEL_VERSION or settings["arguments"]["information_metric"] != "reference":
            raise ValueError("incomplete or incompatible source cohort")
        declared += summary["tests_in_shard"]
        source = read_table(shard / "paired_path.tsv")
        local_ids = requested & set(source.test_id)
        if not local_ids:
            continue
        if seen & local_ids:
            raise ValueError("duplicate requested identities across shards")
        seen |= local_ids
        source = source.set_index("test_id")
        archives = read_table(shard / "subject_scores.tsv.gz")
        null = read_table(shard / "paired_path_null.tsv.gz")
        null = null.loc[null.test_id.isin(local_ids)]
        null_groups = dict(tuple(null.groupby("test_id", sort=False)))
        for test_id, rows in archives.loc[archives.test_id.isin(local_ids)].groupby("test_id", sort=False):
            original, chosen = source.loc[test_id], selected.loc[test_id]
            header = dict(test_id=test_id, panel=chosen.panel, gene_id=chosen.gene_id, source_raw_p_value=float(original.p_value), merged_raw_p_value=float(chosen.raw_p_value), calibrated_p_value=float(chosen.p_value), median_gene_umis=float(original.median_gene_umis), n_isoforms=int(original.n_isoforms), n_ecs=int(original.n_ecs))
            try:
                components = binary_score_components_from_records(rows.to_dict("records"))
                records, influence = binary_score_subject_influence(components)
                draws = null_groups[test_id]
                if len(draws) != 32 or set(draws.replicate) != set(range(32)):
                    raise ValueError("exactly 32 unique original sign draws required")
                values = draws.p_value.to_numpy(float)
                replay_match = np.isclose(influence["fitted_p_value"], original.p_value, rtol=3e-6, atol=1e-300)
                status = "ok" if replay_match else "replay_mismatch"
                diagnostic = {**header, **influence, "status": status, "n_original_sign_draws": len(values), "n_sign_draws_below_threshold": int((values <= args.threshold).sum()), "n_sign_draws_below_observed": int((values <= original.p_value).sum()), "minimum_sign_p_value": float(values.min()), "maximum_sign_p_value": float(values.max()), "identical_sign_p_values": bool((values == values[0]).all()), "leave_dominant_out_available": False}
                if args.compare_proportions:
                    try:
                        _, comparison = binary_score_subject_influence(binary_score_components_to_proportions(components))
                        diagnostic.update({"proportion_" + name: value for name, value in comparison.items()})
                        diagnostic["proportion_status"] = "ok"
                    except (ValueError, np.linalg.LinAlgError) as error:
                        diagnostic.update(proportion_status="diagnostic_error", proportion_error=repr(error))
                if influence["n_informative_subjects"] >= 5:
                    keep = components.subject_ids != influence["dominant_subject"]
                    try:
                        leave = mixed_score_test(components.scores[keep], components.information[keep], components.biological_shapes[keep], reference_information=components.reference_information[keep], scalar_fast=True)
                        diagnostic.update(leave_dominant_out_available=True, leave_dominant_out_p_value=float(leave["p_value"]), leave_dominant_out_mean=float(leave["mean_difference"][0]), leave_dominant_out_subjects=int(leave["n_subjects"]))
                    except (ValueError, np.linalg.LinAlgError) as error:
                        diagnostic["leave_dominant_out_error"] = repr(error)
                trials.append(diagnostic)
                subjects.extend({**header, **record} for record in records)
            except (ValueError, KeyError, np.linalg.LinAlgError) as error:
                trials.append({**header, "status": "diagnostic_error", "error": repr(error)})
        recorded = {row["test_id"] for row in trials}
        for test_id in sorted(local_ids - recorded):
            trials.append(dict(test_id=test_id, panel=selected.at[test_id, "panel"], status="missing_archive"))
        print(f"shard {shard_index}, {len(trials)}/{len(selected)} diagnostic requests", flush=True)
    if seen != requested or declared != len(observed):
        raise ValueError("selected archive coverage or complete cohort size differs from assessment")
    table = pd.DataFrame(trials)
    if table.test_id.duplicated().any() or set(table.test_id) != requested:
        raise ValueError("retain every diagnostic request exactly once")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output_dir / "selected_cases.tsv.gz", sep="\t", index=False)
    table.to_csv(args.output_dir / "diagnostics.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(subjects).to_csv(args.output_dir / "subject_influence.tsv.gz", sep="\t", index=False, na_rep="NA")
    summaries = []
    for panel, local in table.groupby("panel"):
        valid = local.loc[local.status.eq("ok")]
        summaries.append(dict(panel=panel, requested=len(local), valid=len(valid), failures=len(local) - len(valid), median_maximum_precision_share=valid.maximum_precision_share.median() if len(valid) else np.nan, median_effective_weighted_subjects=valid.effective_weighted_subjects.median() if len(valid) else np.nan, precision_share_gt_half=int(valid.maximum_precision_share.gt(.5).sum()) if len(valid) else 0, precision_share_gt_90_percent=int(valid.maximum_precision_share.gt(.9).sum()) if len(valid) else 0, all_32_sign_draws_below_threshold=int(valid.n_sign_draws_below_threshold.eq(32).sum()) if len(valid) else 0, identical_sign_p_values=int(valid.identical_sign_p_values.sum()) if len(valid) else 0))
        if args.compare_proportions:
            comparison = valid.loc[valid.proportion_status.eq("ok")] if len(valid) else valid
            summaries[-1].update(proportion_valid=len(comparison), proportion_failures=len(local) - len(comparison), median_proportion_precision_share=comparison.proportion_maximum_precision_share.median() if len(comparison) else np.nan, median_proportion_effective_subjects=comparison.proportion_effective_weighted_subjects.median() if len(comparison) else np.nan, proportion_precision_share_gt_90_percent=int(comparison.proportion_maximum_precision_share.gt(.9).sum()) if len(comparison) else 0, proportion_diagnostic_F_le_05=int(comparison.proportion_fitted_p_value.le(.05).sum()) if len(comparison) else 0, proportion_diagnostic_F_le_threshold=int(comparison.proportion_fitted_p_value.le(args.threshold).sum()) if len(comparison) else 0, changed_mean_sign=int((comparison.proportion_fitted_mean * comparison.fitted_mean < 0).sum()) if len(comparison) else 0, proportion_mean_outside_simplex_difference=int(comparison.proportion_fitted_mean.abs().gt(np.sqrt(2.)).sum()) if len(comparison) else 0)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    receipt = dict(seed=47852, threshold=args.threshold, declared_cohort_tests=declared, requested=len(selected), observed_sha256=hashlib.sha256(observed_path.read_bytes()).hexdigest(), source=str(args.cohort_root), selection="all complete observed raw-tail tests; distinct weak controls matched on informative subject count, not LR outcomes", scope="conditional precision and dominant-subject omission diagnostic; omission refits Gaussian score aggregation, not EC counts; not alternate production inference, calibrated FDR or complete-method LR comparison", failure_policy="every selected case retained, with errors or replay mismatches explicit", coordinate_comparison="proportion versus ILR common-effect sensitivity; original sign draws are NOT regenerated for proportion coordinates" if args.compare_proportions else "none", changes="none to source counts, fitted scores, null draws, calibration or testing families")
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
