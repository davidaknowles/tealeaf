#!/usr/bin/env python3
"""Experimental geometry-stratified calibration, preserving whole test families."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_event_input_controls import fold_table
from extra_scripts.assess_paired_inference_audit import split_assessment
from extra_scripts.merge_paired_path_test import empirical_null_calibration
from extra_scripts.reassess_event_score_archive import read_table
from extra_scripts.run_differential_splicing import benjamini_hochberg
from extra_scripts.plot_tilgner_method_replication import _rank_table
from tealeaf.sc.path_score_mixed import MODEL_VERSION, binary_information_geometry
from tealeaf.sc.empirical_null import leave_parent_out_cdf
from tealeaf.sc.replication_audit import ranked_direction_summary, ranked_category_summary, reexpress_event_directions


def calibrate_tables(observed, null, geometry):
    """Split existing strata by frozen score/sign-invariant geometry classes.

    No class is dropped or combined with another geometry class. Small pools
    retain their limited CDF resolution, rather than falling back to a mixed
    geometry pool. Original failed tests stay p1 in the original family.
    """
    if observed.test_id.duplicated().any() or geometry.test_id.duplicated().any():
        raise ValueError("unique observed and geometry identities required")
    complete = observed.complete_subject_fits.astype(str).str.lower().eq("true").to_numpy()
    if set(geometry.test_id) != set(observed.loc[complete, "test_id"]):
        raise ValueError("geometry must cover every and only complete inferential test")
    if not geometry.geometry_class.isin(["balanced", "intermediate", "dominated"]).all():
        raise ValueError("unrecognized frozen geometry class")
    if null.duplicated(["test_id", "replicate"]).any() or set(null.test_id) != set(geometry.test_id) or not null.groupby("test_id").replicate.apply(lambda values: set(values) == set(range(32))).all():
        raise ValueError("exact 32 unique original sign draws required for every complete test")
    if not np.isfinite(observed.raw_p_value).all() or not observed.raw_p_value.between(0., 1.).all() or not np.isfinite(null.raw_p_value).all() or not null.raw_p_value.between(0., 1.).all():
        raise ValueError("finite original analytic probabilities required")
    table = observed.merge(geometry, on="test_id", how="left", validate="one_to_one", sort=False)
    if not np.array_equal(table.loc[complete, "n_subjects"], table.loc[complete, "n_informative_subjects"]):
        raise ValueError("geometry information rank differs from original test")
    table["legacy_calibrated_p_value"] = table.p_value
    table["legacy_calibration_stratum"] = table.calibration_stratum
    table["geometry_class"] = table.geometry_class.fillna("failed")
    table["calibration_stratum"] = table.calibration_stratum + "|geometry=" + table.geometry_class
    table["p_value"] = table.raw_p_value
    draws = null.copy()
    draws["p_value"] = draws.raw_p_value
    calibrated, calibrated_null = empirical_null_calibration(table, draws)
    calibrated.loc[~complete, "p_value"] = 1.
    if calibrated.loc[complete, "p_value"].isna().any():
        raise ValueError("missing calibrated probability on a complete test")
    calibrated["fdr"] = benjamini_hochberg(calibrated.p_value)
    return calibrated, calibrated_null


def cohort(args):
    if args.output_dir.exists():
        raise ValueError("use a new output directory, preserve earlier assessments")
    table_path = args.merged_dir / "paired_path.tsv"
    observed = read_table(table_path)
    complete = observed.complete_subject_fits.astype(str).str.lower().eq("true")
    wanted = set(observed.loc[complete, "test_id"])
    lookup = observed.set_index("test_id")
    rows, all_ids, declared = [], [], 0
    for index in range(args.shard_count):
        root = args.cohort_root / f"shard_{index}"
        summary = json.loads((root / "summary.json").read_text())
        settings = json.loads((root / "settings.json").read_text())
        source = read_table(root / "paired_path.tsv")
        failures = json.loads((root / "failures.json").read_text())
        if len(source) != summary["completed"] or len(failures) != summary["failures"] or len(source) + len(failures) != summary["tests_in_shard"]:
            raise ValueError("incomplete source shard")
        if settings["model_version"] != MODEL_VERSION or settings["arguments"]["information_metric"] != "reference" or settings["arguments"].get("score_coordinate", "ilr") != "ilr":
            raise ValueError("incompatible source fitting recipe")
        declared += summary["tests_in_shard"]
        all_ids.extend(source.test_id)
        all_ids.extend(row["test_id"] for row in failures)
        local = source.set_index("test_id")
        contexts = read_table(root / "score_contexts.tsv.gz").set_index("test_id", verify_integrity=True)
        archive = read_table(root / "subject_scores.tsv.gz")
        for test_id, subjects in archive.loc[archive.test_id.isin(wanted)].groupby("test_id", sort=False):
            if test_id not in local.index or not np.isclose(local.at[test_id, "p_value"], lookup.at[test_id, "raw_p_value"], rtol=3e-6, atol=1e-300):
                raise ValueError("merged analytic tail differs from source test")
            if subjects.subject.duplicated().any():
                raise ValueError("duplicated subject archive identity")
            if test_id not in contexts.index or contexts.at[test_id, "score_coordinate"] != "ilr" or len(subjects) != contexts.at[test_id, "n_expected_subjects"]:
                raise ValueError("incomplete or non-ILR archived count fit")
            result = binary_information_geometry(subjects.information, subjects.biological_shape, subjects.reference_information)
            rows.append(dict(test_id=test_id, **result))
        print(f"shard {index}, {len(rows)}/{len(wanted)} complete-test geometries", flush=True)
    if declared != len(observed) or len(set(all_ids)) != len(all_ids) or set(all_ids) != set(observed.test_id):
        raise ValueError("complete source and assessment families differ")
    geometry = pd.DataFrame(rows)
    null_path = args.merged_dir / "paired_path_null.tsv.gz"
    null = read_table(null_path)
    calibrated, draws = calibrate_tables(observed, null, geometry)
    args.output_dir.mkdir(parents=True)
    geometry.to_csv(args.output_dir / "geometry.tsv.gz", sep="\t", index=False)
    calibrated.to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False, na_rep="NA")
    draws.to_csv(args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False)
    summaries = []
    for stratum, local in calibrated.groupby("calibration_stratum"):
        training = draws.loc[draws.calibration_stratum.eq(stratum)]
        summaries.append(dict(stratum=stratum, geometry_class=local.geometry_class.iloc[0], requested=len(local), original_raw_le_1e5=int(local.raw_p_value.le(1e-5).sum()), legacy_calibrated_le_1e5=int(local.legacy_calibrated_p_value.le(1e-5).sum()), geometry_calibrated_le_1e5=int(local.p_value.le(1e-5).sum()), training_draws=len(training), training_raw_le_1e5=int(training.raw_p_value.le(1e-5).sum()), geometry_global_BH_calls=int(local.fdr.le(.05).sum())))
    pd.DataFrame(summaries).to_csv(args.output_dir / "stratum_summary.tsv", sep="\t", index=False)
    receipt = dict(source=str(args.cohort_root), merged_source=str(args.merged_dir), declared_tests=declared, complete_tests=len(geometry), source_shards=args.shard_count, observed_sha256=hashlib.sha256(table_path.read_bytes()).hexdigest(), null_sha256=hashlib.sha256(null_path.read_bytes()).hexdigest(), grid=dict(points=81, relative_log_variance_min=-24., relative_log_variance_max=16., endpoints=[0., "infinity"], scale="median log(1/(I*B))"), classes="balanced max share <=.5; intermediate <=.9; dominated >.9", training="original 32 sign draws, original dimension/subject strata crossed with frozen geometry; exact leave-own-test-out plus-one CDF", small_pools="no fallback across geometry classes; limited empirical tail resolution remains", family="complete original family; every original failed test remains p1", changes="calibration strata only; count fits, scores, analytic tails, reporting, original strata, null signs and family unchanged", selection="no p-value, fitted heterogeneity, effect direction or LR outcomes used for geometry assignment", scope="experimental, requires actual-count-null and both replication endpoints before adoption", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")


def splits(args):
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=False)
    method = "Reference-sequence EC score, geometry-stratified sign calibration"
    paths = [root / "paired_path.tsv" for root in (args.fold0, args.fold1)]
    folds = [fold_table(path, method) for path in paths]
    split_assessment(folds, repo, args.output_dir, method)
    alternate = args.output_dir / "testing_directions"
    alternate.mkdir()
    split_assessment([fold_table(path, method, testing_direction=True) for path in paths], repo, alternate, method + ", efficient-score ILR direction")
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(cohorts=[json.loads((root / "manifest.json").read_text()) for root in (args.fold0, args.fold1)], scope="complete split assessment, frozen published gene/pair families; no LR claim", production_changes=False), indent=2) + "\n")


def count_null(args):
    """Same actual counts, comparing both calibrations from one frozen real pool."""
    if args.output_dir.exists():
        raise ValueError("use a new output directory")
    training = read_table(args.training_dir / "paired_path.tsv")
    signs = read_table(args.training_dir / "paired_path_null.tsv.gz")
    training = training.loc[training.complete_subject_fits.astype(str).str.lower().eq("true")]
    by_size = training.groupby("n_subjects").legacy_calibration_stratum.unique()
    if any(len(value) != 1 for value in by_size):
        raise ValueError("training subject count has inconsistent original strata")
    base_lookup = {int(size): values[0] for size, values in by_size.items()}
    trials, receipts, expected, reference = [], [], set(), None
    for index in range(args.shard_count):
        shard = args.null_root / f"shard_{index}"
        settings = json.loads((shard / "settings.json").read_text())
        if settings.get("mixed_score_version") != MODEL_VERSION or settings.get("score_coordinate") != "ilr" or settings.get("information_metric") != "reference" or len(settings["expected_strategies"]) != 1:
            raise ValueError("incompatible actual-count-null recipe")
        if reference is not None and settings != reference:
            raise ValueError("count-null recipe differs across shards")
        reference = settings
        local = read_table(shard / "observed.tsv")
        diagnostics = read_table(shard / "subject_null_diagnostics.tsv.gz")
        declared = {(parent, draw) for parent in settings["requested_ids"][index::args.shard_count] for draw in range(settings["draws"])}
        actual = set(zip(local.test_id, local.draw))
        if actual != declared or local.duplicated(["test_id", "draw"]).any() or len(set(settings["requested_ids"])) != len(settings["requested_ids"]) or expected & declared:
            raise ValueError("missing, duplicated or foreign requested count-null trials")
        expected |= declared
        groups = dict(tuple(diagnostics.groupby(["test_id", "draw"], sort=False))) if len(diagnostics) else {}
        for row in local.to_dict("records"):
            complete = str(row["converged"]).lower() == "true"
            if complete:
                subjects = groups.get((row["test_id"], row["draw"]))
                if subjects is None or subjects.subject.duplicated().any() or len(subjects) != row["n_expected_subjects"] or row["n_fitted_subjects"] != row["n_expected_subjects"]:
                    raise ValueError("incomplete count-null subject archive")
                geometry = binary_information_geometry(subjects.information, subjects.biological_shape, subjects.reference_information)
                if geometry["n_informative_subjects"] != row["n_subjects"]:
                    raise ValueError("count-null rank differs from original fit")
                row.update(geometry)
            else:
                row.update(geometry_class="failed")
                row["p_value"] = 1.
            base = base_lookup.get(int(row["n_subjects"]), "unavailable") if complete else "failed"
            row.update(complete_subject_fits=complete, baseline_stratum=base, geometry_stratum=base + "|geometry=" + row["geometry_class"])
            trials.append(row)
        receipts.append(dict(shard=index, trials=len(local), observed_sha256=hashlib.sha256((shard / "observed.tsv").read_bytes()).hexdigest(), diagnostics_sha256=hashlib.sha256((shard / "subject_null_diagnostics.tsv.gz").read_bytes()).hexdigest()))
    if len(expected) != len(reference["requested_ids"]) * reference["draws"]:
        raise ValueError("complete null trial family required")
    table = pd.DataFrame(trials)
    original_strata = signs.calibration_stratum.str.split("|geometry=", n=1, regex=False).str[0]
    for label, target, train in (("baseline", "baseline_stratum", original_strata), ("geometry", "geometry_stratum", signs.calibration_stratum)):
        probabilities, counts = leave_parent_out_cdf(table.p_value, table.test_id, table[target], signs.raw_p_value, signs.test_id, train)
        table[label + "_calibrated_p_value"] = probabilities
        table[label + "_training_draws"] = counts
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "trials.tsv.gz", sep="\t", index=False, na_rep="NA")
    summaries = []
    for label, field in (("native", "p_value"), ("baseline real-pool calibration", "baseline_calibrated_p_value"), ("geometry real-pool calibration", "geometry_calibrated_p_value")):
        summaries.append(dict(calibration=label, requested=len(table), usable=int(table.complete_subject_fits.sum()), failures=int((~table.complete_subject_fits).sum()), rejected_05=int(table[field].le(.05).sum()), rejected_01=int(table[field].le(.01).sum()), rejected_001=int(table[field].le(.001).sum())))
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    receipt = dict(source=str(args.null_root), source_recipe=reference, shards=receipts, training=str(args.training_dir), training_manifest=json.loads((args.training_dir / "manifest.json").read_text()), requested=len(table), complete_subject_fits=int(table.complete_subject_fits.sum()), unavailable_baseline_training=int((table.complete_subject_fits & table.baseline_training_draws.eq(0)).sum()), unavailable_geometry_training=int((table.complete_subject_fits & table.geometry_training_draws.eq(0)).sum()), calibration="frozen complete real fold0 sign pool; baseline versus geometry strata on identical analytic trial tails; own real parent excluded, all its sign draws", scope="actual-count conditional-null rejection rates, not independent biological replicates or joint-gene FDR certification", selection="unchanged frozen random hypotheses and every failed draw at p1", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def replace_fixed_mapping_probabilities(mapped, tests):
    """Replace calibration only, rejecting any changed analytic tails or IDs."""
    keys = ["feature_id", "contrast_id"]
    if mapped.duplicated(keys).any() or tests.duplicated(keys).any():
        raise ValueError("unique unchanged association identities required")
    replacement = tests[keys + ["p_value", "raw_p_value", "legacy_calibrated_p_value", "fdr"]].rename(columns={field: "geometry_" + field for field in ("p_value", "raw_p_value", "legacy_calibrated_p_value", "fdr")})
    result = mapped.merge(replacement, on=keys, how="left", validate="one_to_one", indicator=True, sort=False)
    if not result._merge.eq("both").all() or not np.allclose(result.raw_p_value, result.geometry_raw_p_value, rtol=3e-6, atol=1e-300) or not np.allclose(result.p_value, result.geometry_legacy_calibrated_p_value, rtol=3e-6, atol=1e-300):
        raise ValueError("source mapping is missing tests or differs from original inference")
    result["p_value"], result["fdr"] = result.geometry_p_value, result.geometry_fdr
    return result.drop(columns=["_merge", *replacement.columns.difference(keys)])


def long_read(args):
    args.output_dir.mkdir(parents=True, exist_ok=False)
    tests = read_table(args.cohort_dir / "paired_path.tsv")
    source_path = args.source_assessment / "lr_mapping.tsv.gz"
    mapped = replace_fixed_mapping_probabilities(read_table(source_path), tests)
    method = "Reference-sequence EC score, geometry-stratified sign calibration"
    mapped["method"] = method
    mapped.to_csv(args.output_dir / "lr_mapping.tsv.gz", sep="\t", index=False, na_rep="NA")
    valid = mapped.mapping_complete.astype(str).str.lower().eq("true") & mapped.minimum_pooled_depth.ge(20) & mapped.pooled_replicated.notna()
    eligible = mapped.loc[valid].copy()
    eligible["pooled_replicated"] = eligible.pooled_replicated.astype(str).str.lower().eq("true")
    score = reexpress_event_directions(eligible, tests, "test_ilr_effect_size", method + ", efficient-score ILR direction")
    for label, local in (("usage", eligible), ("score", score)):
        root = args.output_dir / label
        root.mkdir()
        ranked = _rank_table(local, len(local))
        pd.DataFrame(ranked_direction_summary(ranked)).to_csv(root / "lr_rank_summary.tsv", sep="\t", index=False)
        pd.DataFrame(ranked_category_summary(ranked)).to_csv(root / "lr_event_type_composition.tsv", sep="\t", index=False)
        ranked.head(200).to_csv(root / "lr_rank.tsv.gz", sep="\t", index=False)
    receipt = dict(cohort=json.loads((args.cohort_dir / "manifest.json").read_text()), source_assessment=str(args.source_assessment), mapping_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(), mapping="same complete all-tested source association mapping and unchanged count fits, independent A1 reports and LR counts; calibration reranked without an FDR or overlap screen", eligible_associations=len(eligible), score_direction_available=int(score.direction_available.sum()), scope="candidate-own calibrated ranking, separate fixed usage and efficient-score directions; no mixing different variants' endpoints", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    collect = sub.add_parser("cohort")
    collect.add_argument("--cohort-root", type=Path, required=True)
    collect.add_argument("--merged-dir", type=Path, required=True)
    collect.add_argument("--output-dir", type=Path, required=True)
    collect.add_argument("--shard-count", type=int, default=64)
    assess = sub.add_parser("splits")
    assess.add_argument("--fold0", type=Path, required=True)
    assess.add_argument("--fold1", type=Path, required=True)
    assess.add_argument("--output-dir", type=Path, required=True)
    validate = sub.add_parser("count-null")
    validate.add_argument("--training-dir", type=Path, required=True)
    validate.add_argument("--null-root", type=Path, required=True)
    validate.add_argument("--output-dir", type=Path, required=True)
    validate.add_argument("--shard-count", type=int, default=16)
    lr = sub.add_parser("long-read")
    lr.add_argument("--cohort-dir", type=Path, required=True)
    lr.add_argument("--source-assessment", type=Path, required=True)
    lr.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "cohort":
        if args.shard_count < 1:
            parser.error("positive shard count required")
        cohort(args)
    elif args.mode == "splits":
        splits(args)
    elif args.mode == "count-null":
        if args.shard_count < 1:
            parser.error("positive shard count required")
        count_null(args)
    else:
        long_read(args)


if __name__ == "__main__":
    main()
