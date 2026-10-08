#!/usr/bin/env python3
"""Frozen gene-balanced global variance pilot, not production inference."""

import argparse
import hashlib
import json
from pathlib import Path
import zlib

import numpy as np
import pandas as pd

from extra_scripts.reassess_event_score_archive import read_table
from tealeaf.sc.path_score_mixed import MODEL_VERSION
from tealeaf.sc.score_variance_pooling import scalar_score_panel, fit_shared_biological_variance


def freeze_training_cases(observed):
    """One declared hypothesis per gene, selected before fit availability."""
    if observed.test_id.duplicated().any() or observed.gene_id.isna().any():
        raise ValueError("unique declared identities and gene labels required")
    frozen = observed.copy()
    frozen["training_hash"] = frozen.test_id.map(lambda key: hashlib.sha256(("shared-variance-47852|" + str(key)).encode()).hexdigest())
    return frozen.sort_values(["gene_id", "training_hash", "test_id"], kind="stable").drop_duplicates("gene_id").copy()


def prepare_panel(records, ids):
    """Ordered complete archives; preserve original subject order for signs."""
    frames = []
    for index, test_id in enumerate(ids):
        subjects = records[test_id].copy()
        if subjects.subject.duplicated().any():
            raise ValueError("duplicate subject identity in archived score test")
        subjects["group"] = index
        frames.append(subjects)
    table = pd.concat(frames, ignore_index=True)
    return scalar_score_panel(table.score, table.information, table.biological_shape, table.reference_information, table.group), table


def signed_values(panel, table, ids, seed, replicate):
    """Same original per-test sign SeedSequence, sampled before rank masking."""
    counts = table.groupby("group", sort=True).size().to_numpy()
    signs = np.concatenate([np.random.default_rng(np.random.SeedSequence((seed, zlib.crc32(str(test_id).encode()), replicate))).choice((-1., 1.), size=int(count)) for test_id, count in zip(ids, counts)])
    return panel.values * signs[panel.source_positions]


def real_pilot(args):
    if args.output_dir.exists():
        raise ValueError("use a new output directory, preserve earlier pilots")
    observed_path = args.merged_dir / "paired_path.tsv"
    observed = read_table(observed_path)
    complete = observed.complete_subject_fits.astype(str).str.lower().eq("true")
    training = freeze_training_cases(observed)
    training["available"] = training.complete_subject_fits.astype(str).str.lower().eq("true")
    training_ids = sorted(training.loc[training.available, "test_id"])
    selected = read_table(args.diagnostic_dir / "selected_cases.tsv.gz")
    if selected.test_id.duplicated().any() or not set(selected.test_id) <= set(observed.loc[complete, "test_id"]):
        raise ValueError("original complete diagnostic requests required")
    wanted = set(training_ids) | set(selected.test_id)
    records, identities, seed, declared = {}, [], None, 0
    for index in range(args.shard_count):
        shard = args.cohort_root / f"shard_{index}"
        summary = json.loads((shard / "summary.json").read_text())
        settings = json.loads((shard / "settings.json").read_text())
        source = read_table(shard / "paired_path.tsv")
        failures = json.loads((shard / "failures.json").read_text())
        if summary["completed"] != len(source) or summary["failures"] != len(failures) or len(source) + len(failures) != summary["tests_in_shard"] or settings["model_version"] != MODEL_VERSION or settings["arguments"]["information_metric"] != "reference":
            raise ValueError("incomplete or incompatible source cohort")
        if seed is not None and seed != settings["arguments"]["seed"]:
            raise ValueError("null seeds differ within source cohort")
        seed = settings["arguments"]["seed"]
        declared += summary["tests_in_shard"]
        identities.extend(source.test_id)
        identities.extend(row["test_id"] for row in failures)
        contexts = read_table(shard / "score_contexts.tsv.gz").set_index("test_id", verify_integrity=True)
        subjects = read_table(shard / "subject_scores.tsv.gz")
        for test_id, local in subjects.loc[subjects.test_id.isin(wanted)].groupby("test_id", sort=False):
            if test_id in records or test_id not in contexts.index or len(local) != contexts.at[test_id, "n_expected_subjects"] or contexts.at[test_id, "score_coordinate"] != "ilr":
                raise ValueError("duplicated, incomplete or non-ILR score archive")
            records[test_id] = local
        print(f"shard {index}, restored {len(records)}/{len(wanted)} requests", flush=True)
    if declared != len(observed) or len(set(identities)) != len(identities) or set(identities) != set(observed.test_id) or set(records) != wanted:
        raise ValueError("whole source family and frozen archive requests must be complete")
    if len(training_ids) < 100:
        raise ValueError("at least 100 available prespecified training genes required")
    train_panel, train_table = prepare_panel(records, training_ids)
    diagnostic_ids = sorted(selected.test_id)
    diagnostic_panel, diagnostic_table = prepare_panel(records, diagnostic_ids)
    variance = fit_shared_biological_variance(train_panel)
    results = diagnostic_panel.evaluate(variance["biological_variance"])
    output = selected.set_index("test_id").loc[diagnostic_ids].reset_index()
    for field, values in results.items():
        output["pooled_" + field] = values
    original = read_table(args.diagnostic_dir / "diagnostics.tsv.gz").set_index("test_id")
    if not original.loc[diagnostic_ids, "status"].eq("ok").all():
        raise ValueError("original diagnostic replay must cover all selected cases")
    for field in ("maximum_precision_share", "effective_weighted_subjects", "fitted_p_value", "fitted_mean", "n_sign_draws_below_threshold"):
        output["original_" + field] = original.loc[diagnostic_ids, field].to_numpy()
    null, hyperparameters = [], []
    for replicate in range(32):
        train_values = signed_values(train_panel, train_table, training_ids, seed, replicate)
        null_variance = fit_shared_biological_variance(train_panel, values=train_values)
        hyperparameters.append(dict(replicate=replicate, **null_variance))
        values = signed_values(diagnostic_panel, diagnostic_table, diagnostic_ids, seed, replicate)
        signed = diagnostic_panel.evaluate(null_variance["biological_variance"], values=values)
        null.extend(dict(test_id=test_id, replicate=replicate, p_value=float(probability), biological_variance=null_variance["biological_variance"]) for test_id, probability in zip(diagnostic_ids, signed["p_value"]))
    null_table = pd.DataFrame(null)
    output["pooled_sign_draws_le_1e5"] = null_table.assign(tail=null_table.p_value.le(1e-5)).groupby("test_id")["tail"].sum().reindex(diagnostic_ids).to_numpy()
    output["pooled_mean_sign_changed"] = output.pooled_mean_difference * output.original_fitted_mean < 0
    summaries = []
    for name, local in output.groupby("panel"):
        summaries.append(dict(panel=name, requested=len(local), original_native_le_05=int(local.original_fitted_p_value.le(.05).sum()), pooled_native_le_05=int(local.pooled_p_value.le(.05).sum()), pooled_native_le_1e5=int(local.pooled_p_value.le(1e-5).sum()), original_dominant_gt_90=int(local.original_maximum_precision_share.gt(.9).sum()), pooled_dominant_gt_90=int(local.pooled_maximum_precision_share.gt(.9).sum()), original_all_32_sign_draws_le_1e5=int(local.original_n_sign_draws_below_threshold.eq(32).sum()), pooled_all_32_sign_draws_le_1e5=int(local.pooled_sign_draws_le_1e5.eq(32).sum()), pooled_mean_sign_changed=int(local.pooled_mean_sign_changed.sum()), median_pooled_effective_subjects=local.pooled_effective_weighted_subjects.median()))
    args.output_dir.mkdir(parents=True)
    training.to_csv(args.output_dir / "training_cases.tsv.gz", sep="\t", index=False, na_rep="NA")
    output.to_csv(args.output_dir / "diagnostics.tsv.gz", sep="\t", index=False, na_rep="NA")
    null_table.to_csv(args.output_dir / "sign_diagnostics.tsv.gz", sep="\t", index=False)
    pd.DataFrame(hyperparameters).to_csv(args.output_dir / "null_variances.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(fitted_variance=variance, declared_tests=declared, requested_training_genes=len(training), available_training_genes=len(training_ids), unavailable_selected_training_genes=int((~training.available).sum()), source=str(args.cohort_root), observed_sha256=hashlib.sha256(observed_path.read_bytes()).hexdigest(), training="one declared hypothesis per gene by fixed SHA256 key, including failures; never replace an unavailable selected case", diagnostic_selection="same previous complete strong-tail and subject-count-matched weak-control panels, no LR selection", null="original 32 per-test subject-sign seeds, with global variance re-estimated in each whole training null family", scope="selected-panel variance/influence pilot, not empirical-tail calibration, complete split/LR power, bounded PSI reporting or joint count-null certification", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(variance, indent=2), flush=True)
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", type=Path, required=True)
    parser.add_argument("--merged-dir", type=Path, required=True)
    parser.add_argument("--diagnostic-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=64)
    args = parser.parse_args()
    if args.shard_count < 1:
        parser.error("positive shard count required")
    real_pilot(args)


if __name__ == "__main__":
    main()
