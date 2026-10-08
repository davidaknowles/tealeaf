#!/usr/bin/env python3
"""Reassess validated actual-count trials with gene-balanced variance pooling."""

import argparse
import hashlib
import json
from pathlib import Path
import zlib

import numpy as np
import pandas as pd

from extra_scripts.audit_shared_score_variance import freeze_training_cases, prepare_panel
from extra_scripts.merge_paired_path_test import add_calibration_strata, empirical_null_calibration
from extra_scripts.reassess_event_score_archive import read_table
from tealeaf.sc.score_variance_pooling import fit_shared_biological_variance
from tealeaf.sc.path_score_mixed import MODEL_VERSION


def count_signed_values(panel, table, ids, draw, replicate):
    """Original count-null sequential sign stream, sampled before rank masking."""
    counts = table.groupby("group", sort=True).size().to_numpy()
    signs = []
    for test_id, count in zip(ids, counts):
        rng = np.random.default_rng(381924 + zlib.crc32(str(test_id).encode()) + 1721 * draw)
        for _ in range(replicate + 1):
            value = rng.choice((-1., 1.), size=int(count))
        signs.append(value)
    return panel.values * np.concatenate(signs)[panel.source_positions]


def assess(args):
    if args.output_dir.exists():
        raise ValueError("use a new output directory")
    manifest_path = args.source_assessment / "manifest.json"
    source = json.loads(manifest_path.read_text())
    trials_path = args.source_assessment / "trials.tsv.gz"
    trials = read_table(trials_path)
    recipe = source["source_recipe"]
    expected = {(test_id, draw) for test_id in recipe["requested_ids"] for draw in range(recipe["draws"])}
    if len(expected) != len(trials) or set(zip(trials.test_id, trials.draw)) != expected or trials.duplicated(["test_id", "draw"]).any() or len(recipe["expected_strategies"]) != 1 or recipe["score_coordinate"] != "ilr" or recipe["information_metric"] != "reference" or recipe["mixed_score_version"] != MODEL_VERSION:
        raise ValueError("complete compatible validated null family required")
    subjects = []
    for receipt in source["shards"]:
        shard = Path(source["source"]) / f"shard_{receipt['shard']}"
        diagnostics = shard / "subject_null_diagnostics.tsv.gz"
        observed = shard / "observed.tsv"
        if hashlib.sha256(diagnostics.read_bytes()).hexdigest() != receipt["diagnostics_sha256"] or hashlib.sha256(observed.read_bytes()).hexdigest() != receipt["observed_sha256"]:
            raise ValueError("actual count score archive changed after source assessment")
        subjects.append(read_table(diagnostics))
    subjects = pd.concat(subjects, ignore_index=True)
    if subjects.duplicated(["test_id", "draw", "subject"]).any():
        raise ValueError("unique actual count subject archives required")
    outcomes, nulls, variances, frozen = [], [], [], []
    for draw, local in trials.groupby("draw", sort=True):
        local = local.copy()
        complete = local.complete_subject_fits.astype(str).str.lower().eq("true")
        training = freeze_training_cases(local)
        training["available"] = training.complete_subject_fits.astype(str).str.lower().eq("true")
        frozen.append(training.assign(count_draw=draw))
        ids = sorted(local.loc[complete, "test_id"])
        training_ids = sorted(training.loc[training.available, "test_id"])
        if len(training_ids) < 4:
            raise ValueError("insufficient prespecified training genes for count-null pooling")
        records = dict(tuple(subjects.loc[subjects.draw.eq(draw) & subjects.test_id.isin(ids)].groupby("test_id", sort=False)))
        if set(records) != set(ids):
            raise ValueError("missing archive for a successful original trial")
        for row in local.loc[complete].itertuples(index=False):
            if len(records[row.test_id]) != row.n_expected_subjects or row.n_fitted_subjects != row.n_expected_subjects:
                raise ValueError("incomplete original successful count fit")
        train_panel, train_table = prepare_panel(records, training_ids)
        panel, table = prepare_panel(records, ids)
        expected_subjects = local.set_index("test_id").loc[ids, "n_subjects"].to_numpy()
        if not np.array_equal(expected_subjects, panel.n_subjects):
            raise ValueError("count-null information rank changed")
        fitted = fit_shared_biological_variance(train_panel)
        variances.append(dict(draw=int(draw), replicate="observed", requested_training_genes=len(training), available_training_genes=len(training_ids), **fitted))
        result = panel.evaluate(fitted["biological_variance"])
        local["original_native_p_value"] = local.p_value
        local["p_value"] = 1.
        positions = local.set_index("test_id").loc[ids].index
        for field, values in result.items():
            lookup = pd.Series(values, index=positions)
            local["pooled_" + field] = local.test_id.map(lookup)
        local.loc[complete, "p_value"] = local.loc[complete, "pooled_p_value"]
        local["degrees_of_freedom"] = 1
        local["n_subjects"] = local.n_subjects.where(complete, 0)
        draws = []
        for replicate in range(32):
            training_values = count_signed_values(train_panel, train_table, training_ids, int(draw), replicate)
            fitted_null = fit_shared_biological_variance(train_panel, values=training_values)
            variances.append(dict(draw=int(draw), replicate=replicate, requested_training_genes=len(training), available_training_genes=len(training_ids), **fitted_null))
            signed = count_signed_values(panel, table, ids, int(draw), replicate)
            probabilities = panel.evaluate(fitted_null["biological_variance"], values=signed)["p_value"]
            draws.extend(dict(test_id=test_id, draw=int(draw), replicate=replicate, p_value=float(probability)) for test_id, probability in zip(ids, probabilities))
        calibrated, null = empirical_null_calibration(add_calibration_strata(local, 100), pd.DataFrame(draws))
        if calibrated.p_value.isna().any():
            raise ValueError("missing pooled calibration probability")
        outcomes.append(calibrated)
        nulls.append(null)
        print(f"draw {draw}, available training genes={len(training_ids)}, pooled variance={fitted['biological_variance']:.6g}", flush=True)
    output = pd.concat(outcomes, ignore_index=True)
    summaries = []
    for label, field in (("source event-specific variance native F", "original_native_p_value"), ("pooled variance native F", "raw_p_value"), ("pooled variance sign-calibrated", "p_value")):
        summaries.append(dict(strategy=label, requested=len(output), usable=int(output.complete_subject_fits.astype(str).str.lower().eq("true").sum()), rejected_05=int(output[field].le(.05).sum()), rejected_01=int(output[field].le(.01).sum()), rejected_001=int(output[field].le(.001).sum())))
    args.output_dir.mkdir(parents=True)
    output.to_csv(args.output_dir / "trials.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.concat(nulls, ignore_index=True).to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False)
    pd.concat(frozen, ignore_index=True).to_csv(args.output_dir / "training_cases.tsv.gz", sep="\t", index=False)
    pd.DataFrame(variances).to_csv(args.output_dir / "variances.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    receipt = dict(source_assessment=str(args.source_assessment), source_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(), source_trials_sha256=hashlib.sha256(trials_path.read_bytes()).hexdigest(), requested=len(output), source_recipe=recipe, original_fit_failures_retained_at_p1=True, training="one declared hypothesis per gene by the same fixed hash; unavailable choices not replaced; separately within each count draw", null="original count-null sign streams, global variance re-estimated for all32 null families of each draw", scope="complete conditional actual-count null reassessment; count draws are not independent biological events, training panel is smaller than real data, and no joint-gene FDR or split/LR power is certified", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-assessment", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    assess(parser.parse_args())


if __name__ == "__main__":
    main()
