#!/usr/bin/env python3
"""Frozen matched split benchmarks and complete LR remapping for prototypes."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.evaluate_suppa2_statistics import normalize_pairs, grouped_pvalues, metric_row
from extra_scripts.audit_split_coverage_direction import direction_rows
from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, load_blocks, block_feature_rows, normalized_difference, vector_agreement, load_path_usage
from tealeaf.sc.differential import helmert_basis
from tealeaf.sc.replication_audit import ranked_direction_summary


def fit_table(path):
    table = pd.read_csv(path, sep="\t", low_memory=False)
    # Preserve the original eligible family when a prototype cannot fit a
    # hypothesis: failed prototype tests contribute p=1, not a smaller family.
    table.loc[~table.converged.astype(bool) | table.n_subjects.lt(4), "p_value"] = 1.
    table["feature_id"] = table.block_id
    table = normalize_pairs(table)
    table["coverage"] = table.median_gene_umis
    table["effect_vector"] = [((helmert_basis(len(json.loads(value)) + 1) @ np.asarray(json.loads(value))).tolist() if value != "[]" else []) for value in table.mean_difference]
    table["effect_features"] = table.path_signatures.map(lambda value: [json.dumps(item, sort_keys=True) for item in json.loads(value)])
    return table


def split_assessment(folds, repo, output, model):
    archived = pd.read_csv(repo / "analyses/split_coverage_direction/event_coverage_pvalues.tsv.gz", sep="\t")
    archived = archived.loc[archived.method.eq("Tealeaf")]
    metrics, directions, genes = [], [], []
    for comparison, local in archived.groupby("comparison"):
        shared = sorted(set(zip(local.gene_id, local.pair_id)))
        universe = pd.DataFrame(shared, columns=["gene_id", "pair_id"])
        tables, pairs, gene_tables = [], [], []
        for fold in folds:
            table = fold.merge(universe, on=["gene_id", "pair_id"], validate="many_to_one")
            paired = grouped_pvalues(table, ["gene_id", "pair_id"], "simes")
            paired = universe.merge(paired, on=["gene_id", "pair_id"], how="left", validate="one_to_one").fillna({"p_value": 1.})
            tables.append(table)
            pairs.append(paired)
            gene_tables.append(grouped_pvalues(paired, ["gene_id"], "simes"))
        metric, detail = metric_row(comparison, model, "simes", "simes", pairs)
        records, status = direction_rows(tables, gene_tables, comparison, model)
        record_table = pd.DataFrame(records)
        selected = record_table.loc[record_table.event_BH_union] if len(record_table) else pd.DataFrame()
        evaluable = selected.loc[selected.direction_agrees.notna()] if len(selected) else pd.DataFrame()
        metric.update({"direction_selected": len(selected), "direction_evaluable": len(evaluable), "direction_agreement": evaluable.direction_agrees.mean() if len(evaluable) else np.nan, "common_events": status.get("common_events", 0), "failed_tests_fold0": int((tables[0].p_value.eq(1.) & ~tables[0].converged).sum()), "failed_tests_fold1": int((tables[1].p_value.eq(1.) & ~tables[1].converged).sum())})
        metrics.append(metric)
        directions.extend(records)
        genes.append(detail)
    pd.DataFrame(metrics).to_csv(output / "split_metrics.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(directions).to_csv(output / "split_directions.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.concat(genes, ignore_index=True).to_csv(output / "split_gene_tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    print(pd.DataFrame(metrics).to_string(index=False), flush=True)


def long_read_assessment(tests, root, matrix_dir, gtf, block_path, output, model, usage=None, effect_vectors=None):
    matrix, features, columns = read_tilgner_matrix(matrix_dir, gtf)
    blocks = load_blocks(block_path)
    groups = {(level, int(rep)): frame.column.to_numpy(int) for (level, rep), frame in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"])}
    represented = set(level for level, rep in groups)
    tests = tests.loc[tests.converged & tests.n_subjects.ge(4) & tests.level_a.isin(represented) & tests.level_b.isin(represented)].copy()
    if usage is None and effect_vectors is None:
        complete_reporting = "inference_backend" in tests and tests.inference_backend.eq("null-corrected").all()
        usage = load_path_usage(root, set(tests.test_id), require_complete=complete_reporting).set_index(["test_id", "cell_type", "path_number"]).proportion
    # Retain original matrix row indices while avoiding a whole-annotation
    # string scan for every tested contrast of the same block.
    feature_groups = {gene: frame for gene, frame in features.groupby("stable_gene_id")}
    count_cache, records = {}, []
    for row in tests.itertuples(index=False):
        if row.block_id not in blocks:
            continue
        signatures = json.loads(row.path_signatures)
        key = (row.block_id, row.path_signatures)
        if key not in count_cache:
            gene = row.gene_id.split(".")[0]
            local = block_feature_rows(blocks[row.block_id], signatures, feature_groups.get(gene, features.iloc[:0]))
            path_rows = {path: frame.row.to_numpy(int) for path, frame in local.groupby("path_number")}
            complete = set(path_rows) == set(range(1, len(signatures) + 1))
            counts = {group: np.asarray([matrix[path_rows.get(path, np.empty(0, int))][:, indices].sum() for path in range(1, len(signatures) + 1)], dtype=float) for group, indices in groups.items()}
            count_cache[key] = complete, counts
        complete, counts = count_cache[key]
        complete = complete and all((level, rep) in groups for level in (row.level_a, row.level_b) for rep in (1, 2))
        a = np.asarray([counts.get((row.level_a, rep), np.zeros(len(signatures))) for rep in (1, 2)])
        b = np.asarray([counts.get((row.level_b, rep), np.zeros(len(signatures))) for rep in (1, 2)])
        if effect_vectors is None:
            short_a = np.asarray([usage.get((row.test_id, row.level_a, path), np.nan) for path in range(1, len(signatures) + 1)])
            short_b = np.asarray([usage.get((row.test_id, row.level_b, path), np.nan) for path in range(1, len(signatures) + 1)])
            effect = short_b - short_a
        else:
            effect = np.asarray(effect_vectors.get(row.test_id, np.full(len(signatures), np.nan)), dtype=float)
        external = normalized_difference(a.sum(axis=0), b.sum(axis=0))
        dot, cosine = vector_agreement(effect, external)
        minimum = min(a.sum(), b.sum())
        eligible = complete and minimum >= 20 and np.isfinite(effect).all() and np.isfinite(dot)
        records.append({"test_id": row.test_id, "block_id": row.block_id, "gene_id": row.gene_id, "level_a": row.level_a, "level_b": row.level_b, "n_paths": len(signatures), "p_value": row.p_value, "raw_p_value": row.raw_p_value, "statistic": row.statistic, "fdr": row.fdr, "mapping_complete": complete, "minimum_pooled_depth": minimum, "eligible": eligible, "pooled_replicated": dot > 0 if eligible else np.nan, "cosine": cosine, "effect": json.dumps(effect.tolist()), **{f"counts_{level}_rep{rep + 1}": json.dumps(value.tolist()) for level, values in (("a", a), ("b", b)) for rep, value in enumerate(values)}})
    mapped = pd.DataFrame(records)
    mapped.to_csv(output / "lr_mapping.tsv.gz", sep="\t", index=False, na_rep="NA")
    summaries, ranks = [], []
    for scope, selected in (("all tested", mapped), ("BH discoveries", mapped.loc[mapped.fdr.le(.05)])):
        eligible = selected.loc[selected.eligible]
        for ranking, columns, ascending in (("calibrated then raw", ["p_value", "raw_p_value", "statistic", "test_id"], [True, True, False, True]), ("continuous raw", ["raw_p_value", "p_value", "statistic", "test_id"], [True, True, False, True])):
            ranked = eligible.sort_values(columns, ascending=ascending, kind="stable").copy()
            ranked["rank"] = np.arange(1, len(ranked) + 1)
            ranked["method"] = model
            for summary in ranked_direction_summary(ranked):
                summaries.append({**summary, "scope": scope, "ranking": ranking})
            ranks.append(ranked.head(200).assign(scope=scope, ranking=ranking))
    pd.DataFrame(summaries).to_csv(output / "lr_rank_summary.tsv", sep="\t", index=False, na_rep="NA")
    pd.concat(ranks, ignore_index=True).to_csv(output / "lr_rank.tsv.gz", sep="\t", index=False, na_rep="NA")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--minimum-gene-umis", type=float, default=25.)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cohort_checks = []
    for fold in (0, 1, "full"):
        summaries = list((args.cache / f"{args.model}_{fold}").glob("shard_*/summary.json"))
        if len(summaries) != 32:
            raise ValueError("expected 32 completed inference shards per cohort")
        cohort = None
        for summary in summaries:
            settings = json.loads(summary.read_text()).get("candidate_settings")
            if settings is None or settings.get("min_gene_umis") != args.minimum_gene_umis:
                raise ValueError(f"candidate cohort is missing or differs from the required gene-count threshold in {summary}")
            if settings.get("subject_fold") != (None if fold == "full" else fold):
                raise ValueError("inference subject fold does not match its assessment cohort")
            if cohort is not None and settings != cohort:
                raise ValueError("candidate settings differ within an inference cohort")
            cohort = settings
        cohort_checks.append({"fold": fold, "candidate_settings": cohort})
    folds = [fit_table(args.cache / f"{args.model}_{fold}/merged/paired_path.tsv") for fold in (0, 1)]
    split_assessment(folds, repo, args.output_dir, args.model)
    full = pd.read_csv(args.cache / f"{args.model}_full/merged/paired_path.tsv", sep="\t")
    corrected = "inference_backend" in full and full.inference_backend.eq("null-corrected").all()
    if corrected:
        reporting_folds = []
        for fold, original in enumerate(folds):
            local = original.copy()
            usage = load_path_usage(args.cache / f"{args.model}_{fold}", set(local.test_id), require_complete=True).set_index(["test_id", "cell_type", "path_number"]).proportion
            effects = []
            for row in local.itertuples(index=False):
                size = len(json.loads(row.path_signatures))
                first = np.asarray([usage.get((row.test_id, row.level_a, path), np.nan) for path in range(1, size + 1)])
                second = np.asarray([usage.get((row.test_id, row.level_b, path), np.nan) for path in range(1, size + 1)])
                vector = second - first
                effects.append(-vector if str(row.level_a) > str(row.level_b) else vector)
            local["effect_vector"] = effects
            reporting_folds.append(local)
        report_output = args.output_dir / "reporting_A1"
        report_output.mkdir(exist_ok=True)
        split_assessment(reporting_folds, repo, report_output, f"{args.model}, free A1 reported direction")
    long_read_assessment(full, args.cache / f"{args.model}_full", args.matrix_dir, args.gtf, args.block_cache, args.output_dir, args.model)
    (args.output_dir / "manifest.json").write_text(json.dumps({"model": args.model, "cohort_checks": cohort_checks, "production_changes": False, "split_universe": "Frozen published comparator-matched gene/pair families, unavailable prototype families assigned p=1", "split_testing": "Same gene/pair Simes, gene Simes, conjunction and BH as Table 1", "variance_moderation": True, "calibration_families": 32, "long_read": "Complete finite converged tested pairwise family, freshly remapped without discovery or historical 704-event restriction", "long_read_reporting": "Independent free-transcript A1 subject arithmetic mean; any failed reporting aggregate invalidates the entire effect" if corrected else "Testing-strength subject-mean fitted usage; score backend uses one-step softmax proxies", "split_direction": "Null-corrected linear-proportion testing response" if corrected else "Testing response", "primary_rank": "calibrated then raw", "interpretation": "Exploratory, requires matched endpoint checks and biological-null validation before promotion"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
