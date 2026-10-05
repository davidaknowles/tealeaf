#!/usr/bin/env python3
"""Audit the Table 1 matched universes without changing their tests."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.evaluate_suppa2_statistics import normalize_pairs, grouped_pvalues, metric_row
from tealeaf.sc.ds_benchmark import benjamini_hochberg
from tealeaf.sc.differential import helmert_basis
from tealeaf.sc.replication_audit import coverage_correlation, aligned_direction
from tealeaf.sc.junction_benchmark import JunctionBundle
from extra_scripts.assess_tilgner_junction_replication import (
    leafcutter_groups, scquint_groups, mean_composition_batch, majiq_feature_ids,
)


def eligible_reference(path):
    table = pd.read_csv(path, sep="\t")
    table = table.loc[table.converged & table.n_subjects.ge(4) & table.method.eq("local_path")].copy()
    table["feature_id"] = table.block_id
    table["method"] = "Tealeaf"
    return normalize_pairs(table)


def load_comparators(repro, repo, fold):
    table = pd.read_csv(repro / f"fold{fold}/comparison_majiq_min3_cov3/all_tests.tsv.gz", sep="\t", low_memory=False)
    table = table.loc[table.effect.eq("cell_type")].copy()
    mapping = pd.read_csv(repro / "leafcutter_cluster_gene.tsv.gz", sep="\t").set_index("feature_id").gene_id
    leaf = table.method.eq("LeafCutter")
    table.loc[leaf, "gene_id"] = table.loc[leaf, "feature_id"].map(mapping)
    print(f"external methods fold={fold}: {table.method.unique().tolist()}", flush=True)
    result = {method: normalize_pairs(local) for method, local in table.groupby("method") if method in ("LeafCutter", "MAJIQ Heterogen", "scQuint", "rMATS")}
    base = repo / "analyses/comparator_suppa_rmats"
    for label, path in (("SUPPA2", base / f"suppa2/split_data_matched_exact_fold{fold}_tests.tsv.gz"),
                        ("SUPPA2 primer-aware", base / f"suppa2_primer_aware/split_data_matched_exact_fold{fold}_tests.tsv.gz"),
                        ("SUPPA2/Tealeaf hybrid", base / f"suppa2_tealeaf_hybrid/split_data/fold{fold}_tests.tsv.gz")):
        result[label] = normalize_pairs(pd.read_csv(path, sep="\t", low_memory=False))
    return result


def subset(table, shared):
    return table.loc[[(gene, pair) in shared for gene, pair in zip(table.gene_id, table.pair_id)]].copy()


def add_coverage(table, reference):
    depth = reference.groupby(["gene_id", "pair_id"], as_index=False).agg(coverage=("median_gene_umis", "median"), reference_n_subjects=("n_subjects", "median"))
    return table.merge(depth, on=["gene_id", "pair_id"], how="left", validate="many_to_one")


def gene_table(table):
    pairs = grouped_pvalues(table, ["gene_id", "pair_id"], "simes")
    genes = grouped_pvalues(pairs, ["gene_id"], "simes")
    depth = table.groupby(["gene_id", "pair_id"]).coverage.median().groupby("gene_id").median().rename("coverage")
    counts = table.groupby("gene_id").agg(n_features=("p_value", "size"), n_subjects=("reference_n_subjects", "median"))
    return genes.merge(depth, on="gene_id").merge(counts, on="gene_id")


def signed_junctions(table, bundle, benchmark, fold, method):
    """Recover split effects, retaining the original event tests and identities.

    LeafCutter/scQuint use equal-pseudobulk junction proportions within their
    tested groups. MAJIQ uses each tested edge's native median PSI difference.
    Missing effects are retained as unevaluable, not removed before event BH.
    """
    local = table.copy()
    local["effect_vector"] = [[] for _ in range(len(local))]
    local["effect_features"] = [[] for _ in range(len(local))]
    if method == "MAJIQ Heterogen":
        for contrast, selected in local.groupby("contrast_id"):
            raw = pd.read_csv(benchmark / f"reproducibility/fold{fold}/majiq_min3_cov3/tests/raw/{contrast}.tsv", sep="\t", comment="#")
            features = majiq_feature_ids(raw)
            a, b = selected[["level_a", "level_b"]].iloc[0]
            delta = raw[f"{b}-raw_psi_quantile_0.500"] - raw[f"{a}-raw_psi_quantile_0.500"]
            effects = dict(zip(features, delta))
            for index, row in selected.iterrows():
                local.at[index, "effect_vector"] = [float(effects.get(row.feature_id, np.nan))]
                local.at[index, "effect_features"] = [str(row.feature_id)]
        return local
    groups = (leafcutter_groups(local, bundle, benchmark / "leafcutter/clustering/benchmark_perind_numers.counts.gz")
              if method == "LeafCutter" else scquint_groups(local, bundle))
    sample_lookup = dict(zip(bundle.samples.sample_id, range(len(bundle.samples))))
    contrasts = {item["contrast_id"]: item for item in json.loads((benchmark / f"reproducibility/fold{fold}/contrasts.json").read_text())}
    effects = {}
    batches = {}
    for key, group in groups.items():
        batches.setdefault((group["contrast_id"], len(group["indices"])), []).append((key, group))
    for (contrast, _), batch in batches.items():
        manifest = contrasts[contrast]
        indices = [group["indices"] for _, group in batch]
        first = mean_composition_batch(bundle.counts, [sample_lookup[sample] for sample in manifest["samples_a"]], indices)
        second = mean_composition_batch(bundle.counts, [sample_lookup[sample] for sample in manifest["samples_b"]], indices)
        for (key, group), delta in zip(batch, second - first):
            effects[key] = (delta.tolist(), group["indices"].astype(str).tolist())
    for index, row in local.iterrows():
        vector, features = effects.get((method, row.contrast_id, row.feature_id), ([], []))
        local.at[index, "effect_vector"] = vector
        local.at[index, "effect_features"] = features
    return local


def signed_reference(table, paths, diagnostic_path=None):
    fits = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    if fits.test_id.duplicated().any():
        raise ValueError("duplicate refitted tests")
    local = table.merge(fits[["test_id", "mean_difference_norm", "mean_difference", "path_signatures"]], on="test_id", how="left", suffixes=("", "_refit"), validate="one_to_one")
    missing = local.mean_difference.isna()
    error = np.abs(local.mean_difference_norm - local.mean_difference_norm_refit)
    if diagnostic_path is not None:
        check = local[["test_id", "gene_id", "median_gene_umis", "mean_difference_norm", "mean_difference_norm_refit"]].copy()
        check["absolute_error"] = error
        check["relative_error"] = error / np.maximum(local.mean_difference_norm, 1e-8)
        check.to_csv(diagnostic_path, sep="\t", index=False)
    if missing.any():
        raise ValueError(f"direction refit mismatch, missing={missing.sum()}, maximum_norm_error={error.max()}, error_quantiles={error.quantile([.5, .9, .99]).to_dict()}, exceeding_1e6={(error > 1e-6).sum()}, exceeding_1e4={(error > 1e-4).sum()}, exceeding_1e2={(error > .01).sum()}")
    local["refit_verified"] = error <= 1e-6
    local["effect_vector"] = [((helmert_basis(len(json.loads(value)) + 1) @ np.array(json.loads(value))).tolist() if verified else []) for value, verified in zip(local.mean_difference, local.refit_verified)]
    local["effect_features"] = [([json.dumps(item, sort_keys=True) for item in json.loads(value)] if verified else []) for value, verified in zip(local.path_signatures, local.refit_verified)]
    return local


def direction_rows(tables, genes, comparison, method):
    # Recompute event BH on each complete matched fold before intersection.
    keyed = []
    for fold, table in enumerate(tables):
        local = table.copy()
        if "effect_vector" not in local:
            field = "test_ilr_effect_size" if method == "SUPPA2/Tealeaf hybrid" else "effect_size"
            if field not in local:
                return [], {"comparison": comparison, "method": method, "status": "signed effects not available"}
            if not pd.to_numeric(local[field], errors="coerce").notna().any():
                return [], {"comparison": comparison, "method": method, "status": "signed effects not available"}
            local["effect_vector"] = local[field].map(lambda value: [float(value)])
        # All source tables use b-minus-a; normalize reversed level order.
        flipped = local.level_a.astype(str) > local.level_b.astype(str)
        for index in local.index[flipped]:
            local.at[index, "effect_vector"] = [-value for value in local.at[index, "effect_vector"]]
        local["event_q"] = benjamini_hochberg(local.p_value.to_numpy(float))
        local["published_q"] = local.fdr if "fdr" in local else np.nan
        keep = ["gene_id", "pair_id", "feature_id", "p_value", "event_q", "coverage", "effect_vector"]
        keep += ["published_q"]
        if "n_paths" in local:
            keep.append("n_paths")
        if "effect_features" in local:
            keep.append("effect_features")
        if local.duplicated(["gene_id", "pair_id", "feature_id"]).any():
            raise ValueError(f"duplicate event tests, {method}")
        keyed.append(local[keep])
    joined = keyed[0].merge(keyed[1], on=["gene_id", "pair_id", "feature_id"], suffixes=("_0", "_1"), validate="one_to_one")
    significant_genes = set()
    for table in genes:
        significant_genes.update(table.loc[benjamini_hochberg(table.p_value.to_numpy(float)) <= .05, "gene_id"])
    records = []
    incompatible = 0
    for row in joined.itertuples(index=False):
        try:
            features0 = getattr(row, "effect_features_0", None)
            features1 = getattr(row, "effect_features_1", None)
            direction = aligned_direction(row.effect_vector_0, row.effect_vector_1, features0, features1)
        except ValueError:
            incompatible += 1
            direction = {"direction_agrees": np.nan, "cosine": np.nan, "nonzero_components": 0, "agreeing_components": 0}
        records.append({"comparison": comparison, "method": method, "gene_id": row.gene_id, "pair_id": row.pair_id, "feature_id": row.feature_id,
                        "p_value_0": row.p_value_0, "p_value_1": row.p_value_1, "event_q_0": row.event_q_0, "event_q_1": row.event_q_1,
                        "coverage": min(row.coverage_0, row.coverage_1), "gene_BH_union": row.gene_id in significant_genes,
                        "event_BH_union": min(row.event_q_0, row.event_q_1) <= .05, "event_BH_intersection": max(row.event_q_0, row.event_q_1) <= .05,
                        "event_BH_union_two_path": min(row.event_q_0, row.event_q_1) <= .05 and getattr(row, "n_paths_0", 2) == 2,
                        "event_BH_union_multi_path": min(row.event_q_0, row.event_q_1) <= .05 and getattr(row, "n_paths_0", 2) > 2,
                        "published_event_BH_union": row.published_q_0 <= .05 or row.published_q_1 <= .05,
                        "fold0_event_BH": row.event_q_0 <= .05, "fold1_event_BH": row.event_q_1 <= .05, **direction})
    return records, {"comparison": comparison, "method": method, "status": "computed", "common_events": len(joined), "incompatible_or_nonfinite": incompatible}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--directions-dir", type=Path)
    parser.add_argument("--junction-directions", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    repro = args.run_root / "junction_benchmark/reproducibility"
    ref = [eligible_reference(repro / f"fold{k}/tealeaf_paired_path_total_a32_production/paired_path.tsv") for k in (0, 1)]
    refit_checks = []
    if args.directions_dir:
        ref = [signed_reference(table, sorted((args.directions_dir / f"fold{k}").glob("shard_*/paired_path.tsv")), args.output_dir / f"fold{k}_refit_norm_diagnostics.tsv") for k, table in enumerate(ref)]
        refit_checks = [{"fold": k, "eligible_tests": len(table), "verified_tests": int(table.refit_verified.sum()), "excluded_tests": int((~table.refit_verified).sum()), "maximum_effect_norm_error": np.abs(table.mean_difference_norm - table.mean_difference_norm_refit).max(), "maximum_allowed_error": 1e-6, "baseline_revision": "d31da47"} for k, table in enumerate(ref)]
        pd.DataFrame(refit_checks).to_csv(args.output_dir / "direction_refit_checks.tsv", sep="\t", index=False)
    others = [load_comparators(repro, args.repo_root, k) for k in (0, 1)]
    bundle = JunctionBundle.load(args.run_root / "junction_benchmark/pseudobulk_junctions") if args.junction_directions else None
    correlations, details, directions, statuses, gene_details, null_correlations, metrics = [], [], [], [], [], [], []
    nulls = []
    for k in (0, 1):
        local = pd.read_csv(repro / f"fold{k}/tealeaf_paired_path_total_a32_production/paired_path_null.tsv.gz", sep="\t")
        local = local.merge(ref[k][["test_id", "gene_id", "pair_id", "median_gene_umis", "n_subjects", "n_paths"]], on="test_id", validate="many_to_one")
        nulls.append(local)
    for comparison in sorted(set(others[0]) & set(others[1])):
        shared = set.intersection(*(set(zip(t.gene_id, t.pair_id)) for t in [*ref, others[0][comparison], others[1][comparison]]))
        for method in ("Tealeaf", comparison):
            tables = [add_coverage(subset((ref[k] if method == "Tealeaf" else others[k][comparison]), shared), ref[k]) for k in (0, 1)]
            if bundle is not None and method in ("LeafCutter", "scQuint", "MAJIQ Heterogen"):
                tables = [signed_junctions(table, bundle, args.run_root / "junction_benchmark", k, method) for k, table in enumerate(tables)]
            gene_folds = []
            for k, table in enumerate(tables):
                gene = gene_table(table)
                gene_folds.append(gene)
                gene_details.append(gene.assign(comparison=comparison, method=method, fold=k))
                for unit, data, controls in (("event", table, table[["reference_n_subjects"]]), ("gene", gene, gene[["n_subjects", "n_features"]])):
                    correlations.append({"comparison": comparison, "method": method, "fold": k, "unit": unit, **coverage_correlation(data.p_value, data.coverage, controls)})
                if method == "Tealeaf":
                    for replicate, null in subset(nulls[k], shared).groupby("replicate"):
                        null = null.rename(columns={"median_gene_umis": "coverage", "n_subjects": "reference_n_subjects"})
                        for unit, data in (("event", null), ("gene", gene_table(null))):
                            null_correlations.append({"comparison": comparison, "fold": k, "replicate": replicate, "unit": unit, **coverage_correlation(data.p_value, data.coverage)})
            if method != "Tealeaf" or args.directions_dir:
                records, status = direction_rows(tables, gene_folds, comparison, method)
                directions.extend(records)
                statuses.append(status)
            metric, _ = metric_row(comparison, method, "simes", "simes", [grouped_pvalues(table, ["gene_id", "pair_id"], "simes") for table in tables])
            metrics.append(metric)
            details.extend([table[["gene_id", "pair_id", "feature_id", "p_value", "coverage"]].assign(comparison=comparison, method=method, fold=k) for k, table in enumerate(tables)])
    pd.DataFrame(correlations).to_csv(args.output_dir / "coverage_correlations.tsv", sep="\t", index=False)
    pd.DataFrame(null_correlations).to_csv(args.output_dir / "null_coverage_correlations.tsv", sep="\t", index=False)
    pd.concat(gene_details).to_csv(args.output_dir / "gene_coverage_pvalues.tsv.gz", sep="\t", index=False)
    pd.concat(details).to_csv(args.output_dir / "event_coverage_pvalues.tsv.gz", sep="\t", index=False)
    statuses.extend({"comparison": "rMATS (historical)", "method": method, "status": "historical split event tables unavailable; full-data effects not substituted"} for method in ("Tealeaf", "rMATS"))
    pd.DataFrame(statuses).to_csv(args.output_dir / "direction_status.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(metrics).to_csv(args.output_dir / "table1_gene_metrics.tsv", sep="\t", index=False)
    # Coverage quintiles refer to the common event set, not a selected tail.
    direction = pd.DataFrame(directions)
    summaries = []
    if len(direction):
        direction["coverage_bin"] = direction.groupby(["comparison", "method"]).coverage.transform(lambda depth: pd.qcut(depth, 5, labels=False, duplicates="drop"))
        for (comparison, method), group in direction.groupby(["comparison", "method"]):
            for selection in ("event_BH_union", "event_BH_union_two_path", "event_BH_union_multi_path", "published_event_BH_union", "event_BH_intersection", "gene_BH_union", "fold0_event_BH", "fold1_event_BH"):
                selected = group.loc[group[selection]]
                for bin_id, local in [("all", selected), *[(str(int(bin_id)), local) for bin_id, local in selected.groupby("coverage_bin")]]:
                    finite = local.direction_agrees.notna()
                    summaries.append({"comparison": comparison, "method": method, "selection": selection, "coverage_bin": bin_id, "n_events": len(local), "n_genes": local.gene_id.nunique(), "n_direction_evaluable": int(finite.sum()), "agree": int(local.loc[finite, "direction_agrees"].sum()), "agreement": local.loc[finite, "direction_agrees"].mean(), "gene_mean_agreement": local.loc[finite].groupby("gene_id").direction_agrees.mean().mean(), "median_cosine": local.cosine.median(), "median_coverage": local.coverage.median()})
        direction.to_csv(args.output_dir / "event_direction_agreement.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "direction_summary.tsv", sep="\t", index=False, na_rep="NA")
    # Fixed-concentration null coverage check, distinct from selector-aware FDR.
    bins = []
    for k, local in enumerate(nulls):
        depth = ref[k].set_index("test_id").median_gene_umis
        categories = pd.qcut(depth, 5, labels=False, duplicates="drop")
        local["coverage_bin"] = local.test_id.map(categories)
        for scope, scoped in (("all", local), ("two_path", local.loc[local.n_paths.eq(2)]), ("multi_path", local.loc[local.n_paths.gt(2)])):
            for bin_id, values in scoped.groupby("coverage_bin"):
                bins.append({"fold": k, "scope": scope, "coverage_bin": int(bin_id), "n_events": values.test_id.nunique(), "null_tests": len(values), "median_coverage": values.median_gene_umis.median(), "reject_0_05": (values.p_value <= .05).mean(), "reject_0_01": (values.p_value <= .01).mean(), "reject_0_001": (values.p_value <= .001).mean()})
    pd.DataFrame(bins).to_csv(args.output_dir / "null_coverage_bins.tsv", sep="\t", index=False)
    print(pd.DataFrame(correlations).to_string(index=False), flush=True)
    print("NULL CORRELATIONS, fixed concentration 32", flush=True)
    print(pd.DataFrame(null_correlations).groupby(["comparison", "fold", "unit"]).rho_p_coverage.agg(["mean", "min", "max"]).to_string(), flush=True)
    print(pd.DataFrame(summaries).query("coverage_bin == 'all'").to_string(index=False) if summaries else "No direction estimates", flush=True)


if __name__ == "__main__":
    main()
