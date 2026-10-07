"""Complete-family split and fresh all-tested LR assessments of input controls."""

import argparse
import json
from pathlib import Path
import pickle
import subprocess
import sys

import numpy as np
import pandas as pd

from extra_scripts.assess_paired_inference_audit import split_assessment
from extra_scripts.plot_tilgner_method_replication import _rank_table
from extra_scripts.run_ec_block_glmm import group_metadata
from extra_scripts.evaluate_suppa2_statistics import normalize_pairs
from tealeaf.sc.replication_audit import complete_paired_fits, coverage_correlation, ranked_direction_summary


def guard_completed_shards(root, output, concentration, shard_count=32):
    """Never mistake a partial array or partial subject fit for a full family."""
    summaries = []
    all_ids = []
    for index in range(shard_count):
        shard = root / f"shard_{index}"
        summary = json.loads((shard / "summary.json").read_text())
        table = pd.read_csv(shard / "paired_path.tsv", sep="\t")
        failures = json.loads((shard / "failures.json").read_text())
        if len(table) != summary["completed"] or len(failures) != summary["failures"] or len(table) + len(failures) != summary["tests_in_shard"]:
            raise ValueError("shard completion does not match its declared family")
        if not table.path_pseudocount.eq(concentration).all() or not table.profile_event_mass.astype(str).str.lower().eq("true").all() or not table.report_pseudocount.eq(1).all():
            raise ValueError("input-control fitting recipe changed within the array")
        complete, reporting = complete_paired_fits(table)
        table["complete_subject_fits"] = complete
        table["complete_reporting_fits"] = reporting
        table.loc[~complete, ["p_value", "statistic"]] = [1., 0.]
        table.loc[~complete, "converged"] = False
        table.loc[~reporting, ["effect_size", "report_psi_effect", "report_ilr_effect"]] = np.nan
        try:
            null = pd.read_csv(shard / "paired_path_null.tsv.gz", sep="\t")
        except pd.errors.EmptyDataError:
            null = pd.DataFrame(columns=["test_id", "replicate", "p_value"])
        null = null.loc[null.test_id.isin(table.loc[complete, "test_id"])].copy()
        if len(null) and (null.duplicated(["test_id", "replicate"]).any() or not null.groupby("test_id").replicate.apply(lambda values: set(values) == set(range(32))).all()):
            raise ValueError("incomplete or duplicated sign-null realizations")
        if set(null.test_id) != set(table.loc[complete, "test_id"]):
            raise ValueError("complete inferential tests must have their declared null draws")
        for failure in failures:
            event, factor, first, second = failure["test_id"].rsplit("|", 3)
            if factor != "cell_type" or not event.startswith("SUPPA2:"):
                raise ValueError("unrecognized failed event identifier")
            table = pd.concat([table, pd.DataFrame([{"test_id": failure["test_id"], "block_id": event.removeprefix("SUPPA2:"), "feature_id": event, "gene_id": event.removeprefix("SUPPA2:").split(";", 1)[0], "level_a": first, "level_b": second, "contrast_id": f"cell_type__{first}__{second}", "effect": "cell_type", "converged": False, "n_samples": 0, "n_subjects": 0, "n_paths": 2, "degrees_of_freedom": 1, "p_value": 1., "statistic": 0., "mean_difference_norm": np.nan, "path_pseudocount": concentration, "path_pseudocount_scaling": "total", "complete_subject_fits": False, "complete_reporting_fits": False}])], ignore_index=True)
        all_ids.extend(table.test_id)
        target = output / f"shard_{index}"
        target.mkdir(parents=True, exist_ok=True)
        table.to_csv(target / "paired_path.tsv", sep="\t", index=False)
        null.to_csv(target / "paired_path_null.tsv.gz", sep="\t", index=False)
        (target / "failures.json").write_text(json.dumps(failures) + "\n")
        summaries.append({**summary, "shard_index": index, "complete_subject_fits": int(complete.sum()), "complete_reporting_fits": int(reporting.sum())})
    if len(set(all_ids)) != len(all_ids):
        raise ValueError("duplicate requested tests across shards")
    return summaries


def fold_table(path, method, testing_direction=False):
    table = pd.read_csv(path, sep="\t")
    table = normalize_pairs(table)
    table["coverage"] = table.median_gene_umis
    field = "test_ilr_effect_size" if testing_direction else "effect_size"
    effect = table[field].where(table.complete_subject_fits, np.nan)
    table["effect_vector"] = effect.map(lambda value: [float(value)])
    table["effect_features"] = [["inclusion"] for _ in range(len(table))]
    table["method"] = method
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--split-only", action="store_true", help="Assess both complete subject halves without waiting for full-data fits.")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    control = args.cache / f"{args.source}_paired"
    with (control / "prepared.pkl").open("rb") as handle:
        groups = pickle.load(handle)[0]
    metadata = group_metadata(groups)
    if metadata.duplicated(["mouse", "cell_type"]).any():
        raise ValueError("n_samples/2 guard requires one pseudobulk per subject/type")
    method = f"Hybrid, {args.source} input control"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cohorts = []
    for fold in ((0, 1) if args.split_only else (0, 1, "full")):
        root = args.cache / (f"full/{args.source}" if fold == "full" else f"split/{args.source}/fold{fold}")
        staged = args.cache / f"assess/{args.source}/{fold}/guarded"
        checks = guard_completed_shards(root, staged, 64 if fold == "full" else 32)
        merged = staged.parent / "merged"
        command = [sys.executable, str(repo / "extra_scripts/merge_paired_path_test.py"), "--shards", *map(str, [staged / f"shard_{index}" for index in range(32)]), "--output-dir", str(merged), "--calibration", "empirical", "--retain-failed-family"]
        if fold == "full":
            command.append("--moderate-variances")
        subprocess.run(command, check=True)
        cohorts.append({"fold": fold, "root": str(root), "merged": str(merged), "shards": checks})
    folds = [fold_table(Path(cohorts[fold]["merged"]) / "paired_path.tsv", method) for fold in (0, 1)]
    split_assessment(folds, repo, args.output_dir, method)
    test_output = args.output_dir / "testing_directions"
    test_output.mkdir(exist_ok=True)
    split_assessment([fold_table(Path(cohorts[fold]["merged"]) / "paired_path.tsv", method, True) for fold in (0, 1)], repo, test_output, method + ", testing ILR directions")
    correlations = [{"fold": fold, "p_scale": column, **coverage_correlation(table[column], table.coverage, table[["n_subjects"]])} for fold, table in enumerate(folds) for column in ("p_value", "raw_p_value")]
    pd.DataFrame(correlations).to_csv(args.output_dir / "coverage_correlations.tsv", sep="\t", index=False)
    if args.split_only:
        (args.output_dir / "manifest.json").write_text(json.dumps({"control": json.loads((control / "manifest.json").read_text()), "cohorts": cohorts, "family": "complete split arrays, failed subject fits retained at p1; fixed published matched gene-pair split universes", "scope": "split-only assessment; no full-data or LR claim", "null_limitation": "independent event sign flips are training calibration, not an actual biological count-null validation", "production_changes": False}, indent=2) + "\n")
        return
    full = pd.read_csv(Path(cohorts[2]["merged"]) / "paired_path.tsv", sep="\t")
    full = full.loc[full.complete_subject_fits & full.complete_reporting_fits].copy()
    full["method"] = method
    tests = args.output_dir / "full_data_tests.tsv.gz"
    full.to_csv(tests, sep="\t", index=False)
    mapping = args.output_dir / "lr_mapping.tsv.gz"
    subprocess.run([sys.executable, str(repo / "extra_scripts/assess_event_tilgner_replication.py"), "--tests", str(tests), "--events", str(args.event_catalog), "--tilgner-matrix", str(args.matrix_dir), "--gtf", str(args.gtf), "--output", str(mapping), "--summary", str(args.output_dir / "lr_depth_summary.tsv"), "--top-per-contrast", "0"], check=True)
    mapped = pd.read_csv(mapping, sep="\t")
    valid = mapped.mapping_complete.astype(str).str.lower().eq("true") & mapped.minimum_pooled_depth.ge(20) & mapped.pooled_replicated.notna()
    eligible = mapped.loc[valid].copy()
    eligible["pooled_replicated"] = eligible.pooled_replicated.astype(str).str.lower().eq("true")
    ranked = _rank_table(eligible, len(eligible))
    pd.DataFrame(ranked_direction_summary(ranked)).to_csv(args.output_dir / "lr_rank_summary.tsv", sep="\t", index=False)
    ranked.head(200).to_csv(args.output_dir / "lr_rank.tsv.gz", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({"control": json.loads((control / "manifest.json").read_text()), "cohorts": cohorts, "family": "all screened events, failed subject fits retained at p1; fixed published matched gene-pair split universes", "LR": "fresh all-complete-tested mapping with unchanged source depth and zero-effect policy; no FDR filter", "null_limitation": "independent event sign flips are training calibration, not an actual biological count-null validation", "production_changes": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
