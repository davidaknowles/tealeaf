#!/usr/bin/env python3
"""Assess independently fitted reporting effects and omnibus alternatives."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from extra_scripts.audit_split_coverage_direction import eligible_reference
from extra_scripts.merge_paired_path_test import add_calibration_strata, empirical_null_calibration
from extra_scripts.assess_tilgner_long_read_replication import normalized_difference, vector_agreement
from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, load_blocks, block_feature_rows
from tealeaf.sc.ds_benchmark import benjamini_hochberg
from tealeaf.sc.ds_benchmark import simes_pvalue, cauchy_pvalue
from tealeaf.sc.replication_audit import aligned_direction, ranked_direction_summary


KEYS = ["gene_id", "pair_id", "block_id"]


def load_shards(root, name="observed.tsv", expected=None):
    paths = sorted(root.glob(f"shard_*/{name}"))
    if expected is None:
        expected = 32 if root.name in ("pairwise_fold0", "pairwise_fold1") else 16
    if len(paths) != expected:
        raise ValueError(f"expected {expected} completed shards, found {len(paths)} in {root}")
    table = pd.concat([pd.read_csv(path, sep="\t", low_memory=False) for path in paths], ignore_index=True)
    if root.name.startswith("pairwise") and name == "observed.tsv":
        table = add_dirichlet_fallback(table)
    return table


def add_dirichlet_fallback(table):
    """Evaluate operational pooling on the full denominator, including failures."""
    baseline = table.loc[table.strategy.eq("subject-mean A1")].copy()
    pooled = table.loc[table.strategy.eq("effective-count Dirichlet pooling"), ["test_id", "effect", "converged"]].rename(columns={"effect": "pooled_effect", "converged": "pooled_converged"})
    fallback = baseline.merge(pooled, on="test_id", how="left", validate="one_to_one")
    valid = fallback.pooled_converged.astype(str).str.lower().eq("true") & fallback.pooled_effect.notna()
    for index in fallback.index[valid]:
        effect = np.asarray(json.loads(fallback.at[index, "pooled_effect"]))
        valid.loc[index] = np.isfinite(effect).all() and np.linalg.norm(effect) > 1e-12
    fallback.loc[valid, "effect"] = fallback.loc[valid, "pooled_effect"]
    fallback["report_fallback"] = ~valid
    fallback["strategy"] = "Dirichlet pooling with subject-mean fallback"
    fallback = fallback.drop(columns=["pooled_effect", "pooled_converged"])
    table = table.copy()
    table["report_fallback"] = False
    return pd.concat([table, fallback], ignore_index=True)


def prepare_effects(table):
    table = table.copy()
    table["gene_id"] = table.gene_id.str.split(".").str[0]
    levels = np.sort(table[["level_a", "level_b"]].astype(str).to_numpy(), axis=1)
    table["pair_id"] = levels[:, 0] + "||" + levels[:, 1]
    table["effect_vector"] = table.effect.map(lambda value: np.asarray(json.loads(value), dtype=float))
    for index in table.index[table.level_a.astype(str) > table.level_b.astype(str)]:
        table.at[index, "effect_vector"] *= -1
    table["features"] = table.path_signatures.map(lambda value: [json.dumps(item, sort_keys=True) for item in json.loads(value)])
    return table


def reporting_summary(joined, masks, output, compact_details=False):
    summaries, details = [], []
    for strategy, local in joined.groupby("strategy"):
        records = []
        for row in local.itertuples(index=False):
            record = {key: getattr(row, key) for key in KEYS}
            valid = bool(row.converged_0 and row.converged_1)
            try:
                if not valid:
                    raise ValueError("failed reporting fit")
                metrics = aligned_direction(row.effect_vector_0, row.effect_vector_1, row.features_0, row.features_1)
                order = [row.features_1.index(feature) for feature in row.features_0]
                first, second = row.effect_vector_0, row.effect_vector_1[order]
                metrics.update({"norm_0": np.linalg.norm(first), "norm_1": np.linalg.norm(second), "difference_norm": np.linalg.norm(first - second), "effect_0": json.dumps(first.tolist()), "effect_1": json.dumps(second.tolist())})
            except ValueError:
                metrics = {"direction_agrees": np.nan, "cosine": np.nan, "effect_0": "[]", "effect_1": "[]"}
            records.append({**record, "strategy": strategy, "report_fallback_either": row.report_fallback_0 or row.report_fallback_1, **metrics})
        directions = pd.DataFrame(records)
        for comparison, mask in masks.groupby("comparison"):
            selected = mask.rename(columns={"feature_id": "block_id"}).merge(directions, on=KEYS, how="left", validate="one_to_one")
            details.append(selected.assign(strategy=strategy))
            scopes = [("all common", np.ones(len(selected), dtype=bool)), ("event BH union", selected.event_BH_union), ("event BH intersection", selected.event_BH_intersection), ("event BH union, two paths", selected.event_BH_union_two_path), ("event BH union, multi paths", selected.event_BH_union_multi_path)]
            for scope, condition in scopes:
                subset = selected.loc[condition].copy()
                finite = subset.direction_agrees.notna()
                observed = subset.loc[finite]
                norms = observed.dropna(subset=["norm_0", "norm_1"])
                norms_rho = float(spearmanr(norms.norm_0, norms.norm_1).statistic) if len(norms) > 2 else np.nan
                arrays = [(np.asarray(json.loads(a)), np.asarray(json.loads(b))) for a, b in zip(observed.effect_0, observed.effect_1)]
                effects_rho = float(spearmanr(np.concatenate([a for a, b in arrays]), np.concatenate([b for a, b in arrays])).statistic) if arrays else np.nan
                summaries.append({"comparison": comparison, "strategy": strategy, "selection": scope, "n_selected": len(subset), "n_evaluable": int(finite.sum()), "n_agree": int(observed.direction_agrees.sum()), "agreement": observed.direction_agrees.mean(), "median_cosine": observed.cosine.median(), "rho_effect_components": effects_rho, "rho_effect_norms": norms_rho, "median_absolute_difference": norms.difference_norm.median(), "n_fallback_either": int(subset.report_fallback_either.fillna(False).sum())})
    pd.DataFrame(summaries).to_csv(output / "split_reporting_summary.tsv", sep="\t", index=False, na_rep="NA")
    detail_table = pd.concat(details, ignore_index=True)
    if compact_details:
        columns = [*KEYS, "strategy", "report_fallback_either", "direction_agrees", "cosine", "nonzero_components", "agreeing_components", "norm_0", "norm_1", "difference_norm", "effect_0", "effect_1"]
        detail_table = detail_table[[column for column in columns if column in detail_table]].drop_duplicates([*KEYS, "strategy"])
    detail_table.to_csv(output / "split_reporting_directions.tsv.gz", sep="\t", index=False, na_rep="NA")
    return pd.DataFrame(summaries)


def rank_split_reports(effects, references, output):
    summaries = []
    merged = effects[0].merge(effects[1], on=[*KEYS, "strategy"], suffixes=("_0", "_1"), validate="one_to_one")
    for fold in (0, 1):
        tests = references[fold].rename(columns={"p_value": "discovery_p", "raw_p_value": "discovery_raw", "fdr": "discovery_q"})
        held = references[1 - fold].rename(columns={"p_value": "held_p", "fdr": "held_q"})
        for strategy, local in merged.groupby("strategy"):
            local = local.merge(tests[[*KEYS, "discovery_p", "discovery_raw", "discovery_q"]], on=KEYS, validate="one_to_one").merge(held[[*KEYS, "held_p", "held_q"]], on=KEYS, validate="one_to_one")
            for ranking, columns in (("calibrated then raw", ["discovery_p", "discovery_raw"]), ("continuous raw", ["discovery_raw", "discovery_p"])):
                ranked = local.sort_values([*columns, *KEYS], kind="stable")
                for cutoff in (100, 200):
                    selected = ranked.head(cutoff)
                    agreements, cosines = [], []
                    for row in selected.itertuples(index=False):
                        try:
                            if not (row.converged_0 and row.converged_1):
                                raise ValueError("failed fit")
                            metric = aligned_direction(row.effect_vector_0, row.effect_vector_1, row.features_0, row.features_1)
                            agreements.append(metric["direction_agrees"])
                            cosines.append(metric["cosine"])
                        except ValueError:
                            agreements.append(np.nan)
                    summaries.append({"discovery_fold": fold, "strategy": strategy, "ranking": ranking, "cutoff": cutoff, "n_selected": len(selected), "n_evaluable": int(np.isfinite(np.asarray(agreements, float)).sum()), "agreement": np.nanmean(agreements), "median_cosine": np.nanmedian(cosines), "held_nominal_replication": selected.held_p.le(.05).mean(), "held_BH_replication": selected.held_q.le(.05).mean()})
    pd.DataFrame(summaries).to_csv(output / "split_top_rank_summary.tsv", sep="\t", index=False, na_rep="NA")


def calibrate_omnibus(observed, null, training_families=32):
    calibrated, held_frames, audits = [], [], []
    for strategy, local in observed.groupby("strategy"):
        local = local.copy().reset_index(drop=True)
        pool = null.loc[null.strategy.eq(strategy)].copy().reset_index(drop=True)
        training = pool.loc[pool.replicate.lt(training_families)].copy()
        held = pool.loc[pool.replicate.ge(training_families)].copy()
        local = add_calibration_strata(local, 100)
        if strategy.startswith(("maximum-coordinate", "CR2")):
            dimension = local.n_paths if "n_paths" in local else local.dimension + 1
            local["calibration_stratum"] += "|paths=" + dimension.astype(str)
        fitted, _ = empirical_null_calibration(local, training)
        fitted["fdr"] = benjamini_hochberg(fitted.p_value.to_numpy(float))
        calibrated.append(fitted)
        strata = fitted.set_index("test_id").calibration_stratum
        held["calibration_stratum"] = held.test_id.map(strata)
        training["calibration_stratum"] = training.test_id.map(strata)
        held["raw_p_value"] = held.p_value
        for stratum, positions in held.groupby("calibration_stratum").groups.items():
            train = training.loc[training.calibration_stratum.eq(stratum)]
            all_values = np.sort(train.p_value.to_numpy(float))
            own = {test_id: np.sort(frame.p_value.to_numpy(float)) for test_id, frame in train.groupby("test_id")}
            for test_id, own_positions in held.loc[positions].groupby("test_id").groups.items():
                raw = held.loc[own_positions, "raw_p_value"].to_numpy(float)
                local_values = own.get(test_id, np.empty(0))
                held.loc[own_positions, "p_value"] = (1 + np.searchsorted(all_values, raw, side="right") - np.searchsorted(local_values, raw, side="right")) / (1 + len(all_values) - len(local_values))
        held_frames.append(held)
        for replicate, frame in held.groupby("replicate"):
            audits.append({"strategy": strategy, "replicate": replicate, "n_tests": len(frame), "reject_0_05": frame.p_value.le(.05).mean(), "reject_0_01": frame.p_value.le(.01).mean(), "reject_0_001": frame.p_value.le(.001).mean(), "BH_discoveries": int(np.sum(benjamini_hochberg(frame.p_value.to_numpy(float)) <= .05))})
    return pd.concat(calibrated, ignore_index=True), pd.concat(held_frames, ignore_index=True), pd.DataFrame(audits)


def paired_combination_omnibus(paired, template):
    """Combine complete paired-testing families, without a discovery screen.

    Bonferroni is valid under any dependence if marginal p-values are valid.
    Simes needs suitable positive dependence; ACAT has an asymptotic tail
    justification. Neither is certified here by independently seeded pairwise
    null families, which do not preserve shared-subject contrast dependence.
    """
    paired = paired.loc[paired.converged.astype(str).str.lower().eq("true") & paired.n_subjects.ge(4) & np.isfinite(paired.p_value)].copy()
    controls = template.loc[template.strategy.str.startswith("null-variance Wald")]
    prototypes = (controls if len(controls) else template).drop_duplicates("block_id").set_index("block_id")
    rows = []
    combinations = {"Bonferroni paired omnibus": lambda p: min(1., len(p) * np.min(p)), "Simes paired omnibus": simes_pvalue, "ACAT paired omnibus": cauchy_pvalue}
    for block_id, family in paired.groupby("block_id"):
        if block_id not in prototypes.index:
            continue
        for strategy, combination in combinations.items():
            p_value = combination(family.p_value.to_numpy(float))
            raw = combination(family.raw_p_value.to_numpy(float))
            record = prototypes.loc[block_id].to_dict()
            record.update({"block_id": block_id, "strategy": strategy, "p_value": p_value, "raw_p_value": raw, "statistic": -np.log10(max(raw, 1e-300)), "degrees_of_freedom": np.nan, "n_paired_contrasts": len(family), "calibration_stratum": "combination of archived marginal tests, not refitted"})
            rows.append(record)
    result = pd.DataFrame(rows)
    result["fdr"] = result.groupby("strategy").p_value.transform(lambda p: benjamini_hochberg(p.to_numpy(float)))
    return result


def split_omnibus_summary(folds, output, prefix="split_omnibus"):
    summaries, directions = [], []
    for strategy in sorted(set(folds[0].strategy) & set(folds[1].strategy)):
        first, second = [frame.loc[frame.strategy.eq(strategy)] for frame in folds]
        common = first.merge(second, on="block_id", suffixes=("_0", "_1"), validate="one_to_one")
        p0, p1 = common.p_value_0.to_numpy(float), common.p_value_1.to_numpy(float)
        q0, q1 = benjamini_hochberg(p0), benjamini_hochberg(p1)
        conjunction = benjamini_hochberg(np.maximum(p0, p1))
        records = []
        for row in common.itertuples(index=False):
            try:
                features0 = [json.dumps(item, sort_keys=True) for item in json.loads(row.path_signatures_0)]
                features1 = [json.dumps(item, sort_keys=True) for item in json.loads(row.path_signatures_1)]
                if set(features0) != set(features1):
                    raise ValueError("incompatible path definitions")
                levels0, levels1 = json.loads(row.levels_0), json.loads(row.levels_1)
                common_levels = sorted(set(levels0) & set(levels1))
                if len(common_levels) < 2:
                    raise ValueError("fewer than two common cell types")
                a = np.asarray(json.loads(row.adjusted_effects_0))[[levels0.index(level) for level in common_levels]]
                b = np.asarray(json.loads(row.adjusted_effects_1))[[levels1.index(level) for level in common_levels]][:, [features1.index(feature) for feature in features0]]
                a -= a.mean(axis=0)
                b -= b.mean(axis=0)
                agreement = aligned_direction(a.ravel(), b.ravel())
            except ValueError:
                agreement = {"direction_agrees": np.nan, "cosine": np.nan}
            records.append({"block_id": row.block_id, "strategy": strategy, **agreement})
        direction = pd.DataFrame(records)
        direction["BH_union"] = (q0 <= .05) | (q1 <= .05)
        directions.append(direction)
        finite = direction.BH_union & direction.direction_agrees.notna()
        held = [p1[q0 <= .05], p0[q1 <= .05]]
        summaries.append({"strategy": strategy, "common_blocks": len(common), "fold0_BH": int(np.sum(q0 <= .05)), "fold1_BH": int(np.sum(q1 <= .05)), "replicated_BH": int(np.sum(conjunction <= .05)), "held_nominal_replication": np.mean([np.mean(values <= .05) for values in held if len(values)]) if any(len(values) for values in held) else np.nan, "rho_logp": spearmanr(-np.log10(p0), -np.log10(p1)).statistic, "union_selected": int(direction.BH_union.sum()), "direction_evaluable": int(finite.sum()), "direction_agreement": direction.loc[finite, "direction_agrees"].mean(), "median_cosine": direction.loc[finite, "cosine"].median()})
    pd.DataFrame(summaries).to_csv(output / f"{prefix}_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(directions, ignore_index=True).to_csv(output / f"{prefix}_directions.tsv.gz", sep="\t", index=False, na_rep="NA")
    return pd.DataFrame(summaries)


def lr_pairwise_summary(effects, source, tests, output):
    source = source.loc[source.mapping_complete.astype(str).str.lower().eq("true") & source.minimum_pooled_depth.ge(20) & source.pooled_replicated.notna()].copy()
    source = source.merge(tests[["test_id", "p_value", "raw_p_value", "statistic"]], on="test_id", validate="one_to_one")
    source["lr_effect"] = [normalized_difference(sum(np.asarray(json.loads(row[f"counts_a_rep{i}"])) for i in (1, 2)), sum(np.asarray(json.loads(row[f"counts_b_rep{i}"])) for i in (1, 2))) for _, row in source.iterrows()]
    summaries, ranked_all = [], []
    for strategy, frame in effects.groupby("strategy"):
        local = source.merge(frame[["test_id", "effect", "converged", "report_fallback"]].rename(columns={"report_fallback": "estimator_fallback"}), on="test_id", how="left", validate="one_to_one")
        valid = local.converged.astype(str).str.lower().eq("true") & local.effect.notna()
        local["report_fallback"] = ~valid | local.estimator_fallback.fillna(False)
        for index in local.index[valid]:
            effect = np.asarray(json.loads(local.at[index, "effect"]))
            dot, _ = vector_agreement(effect, local.at[index, "lr_effect"])
            if np.isfinite(dot) and np.linalg.norm(effect) > 0:
                local.at[index, "pooled_replicated"] = dot > 0
            else:
                local.at[index, "report_fallback"] = True
        for ranking, columns in (("calibrated then raw", ["p_value", "raw_p_value", "statistic", "test_id"]), ("continuous raw", ["raw_p_value", "p_value", "statistic", "test_id"])):
            ranked = local.sort_values(columns, ascending=[True, True, False, True], kind="stable").copy()
            ranked["rank"] = np.arange(1, len(ranked) + 1)
            ranked["method"] = f"{strategy}; {ranking}"
            rows = ranked_direction_summary(ranked)
            for row in rows:
                row.update({"strategy": strategy, "ranking": ranking, "n_fallback": int(ranked.head(row["cutoff"]).report_fallback.sum())})
            summaries.extend(rows)
            ranked_all.append(ranked.head(200).drop(columns="lr_effect"))
    pd.DataFrame(summaries).to_csv(output / "lr_reporting_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(ranked_all, ignore_index=True).to_csv(output / "lr_reporting_rank.tsv.gz", sep="\t", index=False, na_rep="NA")


def lr_omnibus_summary(tests, reports, production_tests, blocks_path, matrix_path, gtf_path, output, paired_tests=None):
    matrix, features, columns = read_tilgner_matrix(matrix_path, gtf_path)
    blocks = load_blocks(blocks_path)
    groups = {level: local.column.to_numpy(int) for level, local in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby("tealeaf_cell_type")}
    reporting = tests.loc[tests.strategy.str.startswith("null-variance Wald")].drop_duplicates("block_id")[["block_id", "levels", "path_signatures", "adjusted_effects"]].assign(reporting="subject-blocked A1")
    pooled = reports.loc[reports.converged.astype(str).str.lower().eq("true")].rename(columns={"strategy": "reporting"})
    reporting = pd.concat([reporting, pooled[["block_id", "levels", "path_signatures", "adjusted_effects", "reporting"]]], ignore_index=True)
    supported_pairs = {}
    if paired_tests is not None:
        paired = paired_tests.loc[paired_tests.converged.astype(str).str.lower().eq("true") & paired_tests.n_subjects.ge(4) & np.isfinite(paired_tests.raw_p_value)]
        paired = paired.loc[paired.level_a.isin(groups) & paired.level_b.isin(groups)].sort_values(["raw_p_value", "p_value", "test_id"], kind="stable")
        supported_pairs = {row.block_id: (row.level_a, row.level_b) for row in paired.drop_duplicates("block_id").itertuples(index=False)}
        supported = reporting.loc[reporting.reporting.eq("subject-blocked A1") & reporting.block_id.isin(supported_pairs)].copy()
        supported["reporting"] = "subject-blocked A1, paired-evidence contrast"
        reporting = pd.concat([reporting, supported], ignore_index=True)
    mapped, count_cache = [], {}
    for row in reporting.itertuples(index=False):
        signatures = json.loads(row.path_signatures)
        levels = json.loads(row.levels)
        effects = np.asarray(json.loads(row.adjusted_effects))
        represented = sorted(set(levels) & set(groups))
        if len(represented) < 2 or row.block_id not in blocks:
            continue
        candidates = [(np.linalg.norm(effects[levels.index(b)] - effects[levels.index(a)]), a, b) for i, a in enumerate(represented) for b in represented[i + 1:]]
        # Select from short-read effects BEFORE mapping/depth/direction checks.
        norm, a, b = sorted(candidates, key=lambda item: (-item[0], item[1], item[2]))[0]
        if row.reporting == "subject-blocked A1, paired-evidence contrast":
            a, b = supported_pairs[row.block_id]
            if a not in levels or b not in levels:
                continue
            norm = np.linalg.norm(effects[levels.index(b)] - effects[levels.index(a)])
        key = (row.block_id, row.path_signatures)
        if key not in count_cache:
            local_features = block_feature_rows(blocks[row.block_id], signatures, features)
            path_rows = {path: group.row.to_numpy(int) for path, group in local_features.groupby("path_number")}
            complete = set(path_rows) == set(range(1, len(signatures) + 1))
            counts = {level: np.array([matrix[path_rows.get(path, np.array([], dtype=int))][:, indices].sum() for path in range(1, len(signatures) + 1)], dtype=float) for level, indices in groups.items()}
            count_cache[key] = complete, counts
        complete, counts = count_cache[key]
        external = normalized_difference(counts[a], counts[b])
        effect = effects[levels.index(b)] - effects[levels.index(a)]
        dot, cosine = vector_agreement(effect, external)
        minimum_depth = min(counts[a].sum(), counts[b].sum())
        eligible = complete and minimum_depth >= 20 and np.isfinite(dot) and norm > 0 and np.linalg.norm(external) > 0
        mapped.append({"block_id": row.block_id, "reporting": row.reporting, "level_a": a, "level_b": b, "effect_norm": norm, "mapping_complete": complete, "minimum_depth": minimum_depth, "eligible": eligible, "pooled_replicated": dot > 0 if eligible else np.nan, "cosine": cosine})
    mapped = pd.DataFrame(mapped)
    control = production_tests[["block_id", "p_value", "raw_p_value", "statistic", "fdr"]].copy().assign(strategy="archived production omnibus")
    all_tests = pd.concat([tests[["block_id", "p_value", "raw_p_value", "statistic", "fdr", "strategy"]], control], ignore_index=True)
    all_tests = all_tests.loc[np.isfinite(all_tests.p_value) & np.isfinite(all_tests.raw_p_value)]
    common_tested = set.intersection(*[set(local.block_id) for _, local in all_tests.groupby("strategy")])
    common_evaluable = set.intersection(*[set(local.loc[local.eligible, "block_id"]) for _, local in mapped.groupby("reporting")])
    common_blocks = common_tested & common_evaluable
    summaries, ranks = [], []
    for (strategy, report), local in all_tests.merge(mapped, on="block_id", validate="many_to_many").groupby(["strategy", "reporting"]):
        for scope, scoped in (("all tested", local), ("BH discoveries", local.loc[local.fdr.lt(.05)]), ("common tested and LR eligible", local.loc[local.block_id.isin(common_blocks)])):
            eligible = scoped.loc[scoped.eligible].copy()
            for ranking, columns in (("calibrated then raw", ["p_value", "raw_p_value", "statistic", "block_id"]), ("continuous raw", ["raw_p_value", "p_value", "statistic", "block_id"])):
                ranked = eligible.sort_values(columns, ascending=[True, True, False, True], kind="stable").copy()
                ranked["rank"] = np.arange(1, len(ranked) + 1)
                ranked["method"] = f"{strategy}; {report}; {ranking}; {scope}"
                for summary in ranked_direction_summary(ranked):
                    summaries.append({**summary, "strategy": strategy, "reporting": report, "ranking": ranking, "selection": scope})
                ranks.append(ranked.head(200))
    pd.DataFrame(summaries).to_csv(output / "lr_omnibus_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(ranks, ignore_index=True).to_csv(output / "lr_omnibus_rank.tsv.gz", sep="\t", index=False, na_rep="NA")
    mapped.to_csv(output / "lr_omnibus_mapping.tsv.gz", sep="\t", index=False, na_rep="NA")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--reporting-only", action="store_true")
    parser.add_argument("--include-wild", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    references = [eligible_reference(args.run_root / f"junction_benchmark/reproducibility/fold{fold}/tealeaf_paired_path_total_a32_production/paired_path.tsv") for fold in (0, 1)]
    effects = [prepare_effects(load_shards(args.cache / f"pairwise_fold{fold}")) for fold in (0, 1)]
    masks = pd.read_csv(repo / "analyses/split_coverage_direction/event_direction_agreement.tsv.gz", sep="\t")
    masks = masks.loc[masks.method.eq("Tealeaf")].drop(columns=["direction_agrees", "cosine", "nonzero_components", "agreeing_components"])
    joined = effects[0].merge(effects[1], on=[*KEYS, "strategy"], suffixes=("_0", "_1"), validate="one_to_one")
    summary = reporting_summary(joined, masks, args.output_dir)
    rank_split_reports(effects, references, args.output_dir)
    full_effects = load_shards(args.cache / "pairwise_full")
    strength = float(json.loads((args.run_root / "differential/local_path_full_testing_eb_pairwise/selection/selection.json").read_text())["strength"])
    tests = pd.read_csv(args.run_root / f"differential/local_path_full_testing_eb_pairwise/a{strength:g}/paired_path.tsv", sep="\t")
    lr_pairwise_summary(full_effects, pd.read_csv(repo / "analyses/tilgner_long_read_replication/pairwise_replication.tsv", sep="\t"), tests, args.output_dir)
    print(summary.query("selection == 'event BH union'").to_string(index=False), flush=True)
    if args.reporting_only:
        return
    calibrated_folds, audits, full_tests = [], [], None
    for fold in (0, 1, "full"):
        root = args.cache / (f"omnibus_fold{fold}" if fold != "full" else "omnibus_full")
        observed, null = load_shards(root), load_shards(root, "null.tsv.gz")
        if args.include_wild:
            wild_root = args.cache / (f"wild_omnibus_fold{fold}" if fold != "full" else "wild_omnibus_full")
            observed = pd.concat([observed, load_shards(wild_root)], ignore_index=True)
            null = pd.concat([null, load_shards(wild_root, "null.tsv.gz")], ignore_index=True)
        calibrated, held, audit = calibrate_omnibus(observed, null)
        # Every statistic uses one explicit A1 reporting control. Alphabetical
        # method ordering must not silently substitute the wild test's A32
        # descriptive arrays for another method's A1 arrays.
        report_control = calibrated.loc[calibrated.strategy.str.startswith("null-variance Wald")].set_index("block_id")
        for field in ("levels", "path_signatures", "adjusted_effects"):
            value = calibrated.block_id.map(report_control[field])
            available = value.notna()
            calibrated.loc[available, field] = value.loc[available]
        if fold == "full":
            paired_table = tests
        else:
            paired_table = references[fold]
        calibrated = pd.concat([calibrated, paired_combination_omnibus(paired_table, calibrated)], ignore_index=True)
        calibrated.to_csv(args.output_dir / f"omnibus_{fold}_tests.tsv.gz", sep="\t", index=False, na_rep="NA")
        audit["fold"] = fold
        audits.append(audit)
        if fold != "full":
            calibrated_folds.append(calibrated)
        else:
            full_tests = calibrated
    pd.concat(audits, ignore_index=True).to_csv(args.output_dir / "held_null_summary.tsv", sep="\t", index=False, na_rep="NA")
    omnibus_summary = split_omnibus_summary(calibrated_folds, args.output_dir)
    common = set.intersection(*[set(local.block_id) for fold in calibrated_folds for _, local in fold.groupby("strategy")])
    split_omnibus_summary([frame.loc[frame.block_id.isin(common)] for frame in calibrated_folds], args.output_dir, prefix="split_omnibus_matched")
    omnibus_strength = float(json.loads((args.run_root / "differential/local_path_full_testing_eb_omnibus_production/selection/selection.json").read_text())["strength"])
    production_omnibus = pd.read_csv(args.run_root / f"differential/local_path_full_testing_eb_omnibus_production/a{omnibus_strength:g}/paired_path.tsv", sep="\t")
    lr_omnibus_summary(full_tests, load_shards(args.cache / "omnibus_pooled"), production_omnibus, args.run_root / "differential/gencode_vM32_splice_blocks.json.gz", args.matrix_dir, args.gtf, args.output_dir, paired_tests=tests)
    print(omnibus_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
