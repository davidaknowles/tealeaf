#!/usr/bin/env python3
"""Joint split/LR assessment of prespecified subject-paired estimators."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_split_coverage_direction import eligible_reference
from extra_scripts.summarize_path_reporting_omnibus import prepare_effects, reporting_summary, rank_split_reports, lr_pairwise_summary, load_shards
from extra_scripts.assess_tilgner_long_read_replication import normalized_difference, vector_agreement
from tealeaf.sc.replication_audit import ranked_direction_summary


def read_reports(root, expected):
    paths = sorted(root.glob("shard_*/observed.tsv"))
    if len(paths) != expected:
        raise ValueError(f"expected {expected} completed shards in {root}, found {len(paths)}")
    table = pd.concat([pd.read_csv(path, sep="\t", low_memory=False) for path in paths], ignore_index=True)
    if table.duplicated(["test_id", "strategy"]).any():
        raise ValueError("duplicate reporting estimates")
    table["report_fallback"] = table.report_fallback.fillna(False).astype(bool)
    return table


def cohort_settings(root):
    settings = [json.loads(path.read_text())["candidate_settings"] for path in sorted(root.glob("shard_*/settings.json"))]
    if not settings or any(item != settings[0] for item in settings[1:]):
        raise ValueError(f"incomplete or inconsistent candidate metadata in {root}")
    return settings[0]


def check_reporting_control(new, control, tolerance=5e-5):
    """Reject cohort/model drift, allowing small optimizer roundoff only."""
    first = new.loc[new.strategy.eq("subject arithmetic mean A1"), ["test_id", "effect", "n_subjects"]]
    second = control.loc[control.strategy.eq("subject-mean A1"), ["test_id", "effect", "n_subjects"]]
    same = first.merge(second, on="test_id", suffixes=("_new", "_control"), validate="one_to_one", how="outer", indicator=True)
    if not same["_merge"].eq("both").all() or not same.n_subjects_new.eq(same.n_subjects_control).all():
        raise ValueError("arithmetic control has different hypotheses or subject counts")
    errors = np.array([np.max(np.abs(np.asarray(json.loads(a)) - np.asarray(json.loads(b)))) for a, b in zip(same.effect_new, same.effect_control)])
    if not np.isfinite(errors).all() or np.max(errors, initial=0) > tolerance:
        raise ValueError("arithmetic control effect differences exceed optimizer tolerance")
    return {"n_tests": len(same), "maximum_absolute_component_difference": float(np.max(errors, initial=0)), "median_absolute_component_difference": float(np.median(errors)), "tolerance": tolerance}


def lr_endpoint_audit(effects, source, tests, output):
    """Keep the pooled endpoint primary; inspect replicate weighting separately."""
    source = source.loc[source.mapping_complete.astype(str).str.lower().eq("true") & source.minimum_pooled_depth.ge(20) & source.pooled_replicated.notna()].copy()
    source = source.merge(tests[["test_id", "p_value", "raw_p_value", "statistic"]], on="test_id", validate="one_to_one")
    external, endpoint_records = {}, []
    for row in source.itertuples(index=False):
        a = np.asarray([json.loads(getattr(row, f"counts_a_rep{i}")) for i in (1, 2)])
        b = np.asarray([json.loads(getattr(row, f"counts_b_rep{i}")) for i in (1, 2)])
        individual = np.asarray([normalized_difference(first, second) for first, second in zip(a, b)])
        pooled = normalized_difference(a.sum(axis=0), b.sum(axis=0))
        equal = individual.mean(axis=0)
        dot, cosine = vector_agreement(pooled, equal)
        between, _ = vector_agreement(individual[0], individual[1])
        external[row.test_id] = {"pooled": pooled, "equal biological replicates": equal, "replicate 1": individual[0], "replicate 2": individual[1]}
        endpoint_records.append({"test_id": row.test_id, "gene_id": row.gene_id, "block_id": row.block_id, "minimum_replicate_depth": row.minimum_replicate_depth, "pooled_vs_equal_agree": dot > 0 if np.isfinite(dot) else np.nan, "pooled_vs_equal_cosine": cosine, "between_replicates_agree": between > 0 if np.isfinite(between) else np.nan, "celltype_difference_in_rep1_read_weight": a[0].sum() / a.sum() - b[0].sum() / b.sum()})
    endpoint_table = pd.DataFrame(endpoint_records)
    endpoint_table.to_csv(output / "lr_endpoint_consistency.tsv", sep="\t", index=False, na_rep="NA")
    records, ranks = [], []
    common = set.intersection(*[set(local.loc[local.converged, "test_id"]) for _, local in effects.groupby("strategy")]) & set(source.test_id)
    for strategy, local in effects.groupby("strategy"):
        observed = source.merge(local[["test_id", "effect", "converged", "report_fallback"]], on="test_id", how="left", validate="one_to_one")
        for endpoint in ("pooled", "equal biological replicates", "replicate 1", "replicate 2"):
            observed = observed.copy()
            observed["pooled_replicated"] = [vector_agreement(np.asarray(json.loads(row.effect)), external[row.test_id][endpoint])[0] > 0 if isinstance(row.effect, str) and row.converged and np.isfinite(vector_agreement(np.asarray(json.loads(row.effect)), external[row.test_id][endpoint])[0]) else np.nan for row in observed.itertuples(index=False)]
            observed["effect_norm"] = observed.effect.map(lambda value: np.linalg.norm(json.loads(value)) if isinstance(value, str) else np.nan)
            for scope, subset in (("own evaluable", observed), ("common estimators", observed.loc[observed.test_id.isin(common)])):
                evaluable = subset.loc[subset.pooled_replicated.notna()].copy()
                for ordering, columns, ascending in (("calibrated then raw", ["p_value", "raw_p_value", "statistic", "test_id"], [True, True, False, True]), ("continuous raw", ["raw_p_value", "p_value", "statistic", "test_id"], [True, True, False, True]), ("effect magnitude, not significance", ["effect_norm", "raw_p_value", "test_id"], [False, True, True])):
                    ranked = evaluable.sort_values(columns, ascending=ascending, kind="stable").copy()
                    ranked["rank"] = np.arange(1, len(ranked) + 1)
                    ranked["method"] = f"{strategy}; {endpoint}; {scope}; {ordering}"
                    for summary in ranked_direction_summary(ranked):
                        records.append({**summary, "strategy": strategy, "endpoint": endpoint, "scope": scope, "ranking": ordering, "n_eligible": len(evaluable), "n_selected": len(subset), "n_fallback": int(evaluable.report_fallback.fillna(False).sum())})
                    columns = ["test_id", "gene_id", "block_id", "rank", "method", "pooled_replicated", "report_fallback", "effect_norm"]
                    ranks.append(ranked.head(200)[columns].assign(endpoint=endpoint, scope=scope, ranking=ordering))
    pd.DataFrame(records).to_csv(output / "lr_endpoint_rank_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(ranks, ignore_index=True).to_csv(output / "lr_endpoint_rank.tsv.gz", sep="\t", index=False, na_rep="NA")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--previous-cache", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--free-cache", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cohort_checks = []
    for fold in (0, 1, "full"):
        folder = "pairwise_full" if fold == "full" else f"pairwise_fold{fold}"
        roots = [args.cache / folder, args.previous_cache / folder]
        if args.free_cache:
            roots.append(args.free_cache / folder)
        cohorts = [cohort_settings(root) for root in roots]
        if any(cohort != cohorts[0] for cohort in cohorts[1:]):
            raise ValueError(f"candidate settings differ across reporting controls in {folder}")
        cohort_checks.append({"fold": fold, "candidate_settings": cohorts[0]})
    references = [eligible_reference(args.run_root / f"junction_benchmark/reproducibility/fold{fold}/tealeaf_paired_path_total_a32_production/paired_path.tsv") for fold in (0, 1)]
    controls = [load_shards(args.previous_cache / f"pairwise_fold{fold}") for fold in (0, 1)]
    control_checks = []
    for fold in (0, 1):
        new = read_reports(args.cache / f"pairwise_fold{fold}", 32)
        control_checks.append({"fold": fold, **check_reporting_control(new, controls[fold])})
    effects = [prepare_effects(pd.concat([read_reports(args.cache / f"pairwise_fold{fold}", 32), controls[fold]], ignore_index=True)) for fold in (0, 1)]
    if args.free_cache:
        for fold in (0, 1):
            free = read_reports(args.free_cache / f"pairwise_fold{fold}", 32)
            free["strategy"] = "Free isoforms; " + free.strategy
            effects[fold] = pd.concat([effects[fold], prepare_effects(free)], ignore_index=True)
    masks = pd.read_csv(repo / "analyses/split_coverage_direction/event_direction_agreement.tsv.gz", sep="\t")
    masks = masks.loc[masks.method.eq("Tealeaf")].drop(columns=["direction_agrees", "cosine", "nonzero_components", "agreeing_components"])
    joined = effects[0].merge(effects[1], on=["gene_id", "pair_id", "block_id", "strategy"], suffixes=("_0", "_1"), validate="one_to_one")
    summary = reporting_summary(joined, masks, args.output_dir, compact_details=True)
    rank_split_reports(effects, references, args.output_dir)
    full = pd.concat([read_reports(args.cache / "pairwise_full", 16), load_shards(args.previous_cache / "pairwise_full")], ignore_index=True)
    control_checks.append({"fold": "full", **check_reporting_control(full, full)})
    (args.output_dir / "control_checks.json").write_text(json.dumps(control_checks, indent=2) + "\n")
    if args.free_cache:
        free = read_reports(args.free_cache / "pairwise_full", 16)
        free["strategy"] = "Free isoforms; " + free.strategy
        full = pd.concat([full, free], ignore_index=True)
    strength = float(json.loads((args.run_root / "differential/local_path_full_testing_eb_pairwise/selection/selection.json").read_text())["strength"])
    tests = pd.read_csv(args.run_root / f"differential/local_path_full_testing_eb_pairwise/a{strength:g}/paired_path.tsv", sep="\t")
    source = pd.read_csv(repo / "analyses/tilgner_long_read_replication/pairwise_replication.tsv", sep="\t")
    lr_pairwise_summary(full, source, tests, args.output_dir)
    lr_endpoint_audit(full, source, tests, args.output_dir)
    split = summary.loc[summary.comparison.eq("SUPPA2") & summary.selection.eq("event BH union")]
    lr = pd.read_csv(args.output_dir / "lr_reporting_summary.tsv", sep="\t")
    lr = lr.loc[lr.cutoff.eq(100) & lr.ranking.eq("continuous raw")]
    joint = split.merge(lr[["strategy", "agreement", "normalized_auc", "n_available"]].rename(columns={"agreement": "lr_top100_agreement", "normalized_auc": "lr_A100"}), on="strategy", validate="one_to_one")
    joint.to_csv(args.output_dir / "joint_endpoint_summary.tsv", sep="\t", index=False, na_rep="NA")
    print(joint[["strategy", "n_selected", "n_evaluable", "agreement", "rho_effect_components", "lr_A100", "lr_top100_agreement"]].to_string(index=False), flush=True)
    (args.output_dir / "manifest.json").write_text(json.dumps({"assessment": "Paired reporting estimates against independent subject halves and external long reads", "cohort_checks": cohort_checks, "production_tests_changed": False, "concentrations": [1, 4, 16, 32, 64], "free_isoform_sensitivity": bool(args.free_cache), "free_isoform_concentrations": [1, 32, 64] if args.free_cache else [], "selection": "frozen production matched event-BH union, with selected and evaluable denominators", "estimator_selection_uses_long_reads": False, "reporting_estimators": ["arithmetic", "geometric", "paired median CLR", "paired harmonic depth", "paired isotropic REML"], "external_primary": "original pooled long-read direction", "external_secondary": ["equal biological replicates", "replicate 1", "replicate 2"], "interpretation": "Exploratory grid, not a newly independently validated production strategy; tests and effect reporting must both be assessed"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
