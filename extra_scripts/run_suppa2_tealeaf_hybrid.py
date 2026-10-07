#!/usr/bin/env python3
"""Test SUPPA2 event definitions with Tealeaf's local-path EC model.

Each SUPPA2 event's included and excluded transcript sets define two tested
paths with fixed pooled within-class mixtures. With ``--profile-event-mass``,
their combined mass is estimated per subject and cell type; otherwise it is
fixed at the pooled value for a legacy sensitivity analysis. Other
EC-supported isoforms form one fixed-mixture nuisance component. Primer-specific EC
counts are fitted with Tealeaf's paired path estimator and paired test; signed
subject-label nulls are drawn independently for each event/contrast and emitted
in the production path format for pooled empirical calibration. The reported
testing effect is the subject-mean included-minus-excluded ILR difference,
level b minus level a, equal to the mean inclusion-logit difference divided
by sqrt(2). Optional descriptive refits also export mean PSI differences
under a separate reporting pseudocount, without changing test statistics.

The opt-in experimental mixed-score backend instead retains every supported
transcript, fits an unpenalized shared-path subject null with type-specific
within-path nuisance, propagates efficient-score information into REML and an
approximate small-sample F test, and reports independent weak-prior PSI fits.
It never applies the default backend's target smoothing or variance moderation.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    covered_celltype_pairwise_designs,
    local_test_design,
    modeled_gene_umis,
)
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_paired_path_test import filtered_inputs, signed_null_p_value
from tealeaf.sc import ec_block_glmm
from tealeaf.sc.event_paths import collapse_event_nuisance
from tealeaf.sc.path_score_mixed import shared_path_score_components, aggregate_path_scores, paired_score_reporting, signed_path_score_p_value, binary_subject_score_records, MODEL_VERSION
from tealeaf.sc.replication_audit import complete_cluster_fit

MIXED_SCORE_EMPTY_COLUMNS = ("test_id", "path_pseudocount", "profile_event_mass", "report_pseudocount", "inference_backend", "model_version", "converged", "n_subjects", "n_expected_subjects", "n_fitted_subjects", "n_samples", "complete_reporting_fits", "report_n_subjects", "effect_size", "report_psi_effect", "p_value", "statistic")


def canonical(value):
    return str(value).split(".", 1)[0]


def event_path_index(transcripts, features, included, excluded):
    """Map event membership to a binary path vector; return None if incomplete."""
    feature_ids = [canonical(features[index]) for index in transcripts]
    included_ids = {canonical(value) for value in included.split(",") if value}
    excluded_ids = {canonical(value) for value in excluded.split(",") if value}
    if not included_ids or not excluded_ids or included_ids & excluded_ids:
        return None
    positions = {}
    for index, transcript in enumerate(feature_ids):
        positions.setdefault(transcript, []).append(index)
    if not (included_ids | excluded_ids) <= positions.keys():
        return None
    path_index = np.full(len(transcripts), -1, dtype=int)
    for transcript in included_ids:
        path_index[positions[transcript]] = 0
    for transcript in excluded_ids:
        path_index[positions[transcript]] = 1
    return path_index


def supported_gene_transcripts(gene, gene_transcripts, gene_ecs, designs):
    """Retain gene isoforms represented in at least one primer EC map."""
    transcripts = np.asarray(gene_transcripts[gene], dtype=int)
    ecs = np.asarray(gene_ecs[gene], dtype=int)
    supported = np.zeros(len(transcripts), dtype=bool)
    for design in designs:
        local = design[ecs][:, transcripts]
        supported |= np.asarray(local.sum(axis=0)).ravel() > 0
    return transcripts[supported]


def partition_event_tests(tests, shard_count):
    """Shard by gene/context so one event family stays together."""
    groups = {}
    for test in tests:
        candidate = test[0]
        key = (int(candidate[3]), tuple(np.asarray(candidate[7], dtype=int)))
        groups.setdefault(key, []).append(test)
    shards = [[] for _ in range(int(shard_count))]
    loads = [0] * int(shard_count)
    for group in sorted(groups.values(), key=len, reverse=True):
        index = min(range(len(shards)), key=loads.__getitem__)
        shards[index].extend(group)
        loads[index] += len(group)
    return shards


def mixed_event_record(base, path_index, labels, clusters, baseline, event, gene_id, tested_levels, totals, n_ecs, args, *, score_archive=None, components=None):
    """Experimental full-transcript event score, not fixed-share paired t."""
    if components is None:
        components = shared_path_score_components(base, path_index, labels, clusters, baseline=baseline, max_iter=args.max_iter, null_multistart=args.null_multistart, reporting_concentration=args.report_pseudocount)
    expected = len(np.unique(clusters))
    event_id = str(event.event_id)
    test_id = f"SUPPA2:{event_id}|cell_type|{'|'.join(tested_levels)}"
    if score_archive is not None:
        contexts, subject_scores = score_archive
        contexts.append(dict(test_id=test_id, block_id=event_id, feature_id=event.feature_id, gene_id=gene_id, event_type=event.event_type, level_a=tested_levels[0], level_b=tested_levels[1], n_expected_subjects=expected, n_samples=len(labels), n_isoforms=base.n_isoforms, n_ecs=n_ecs, median_gene_umis=float(np.median(totals)), baseline_event_mass=float(baseline[path_index >= 0].sum()), report_pseudocount=args.report_pseudocount, score_coordinate=components.score_coordinate, model_version=MODEL_VERSION))
        subject_scores.extend(binary_subject_score_records(components, test_id))
    information_metric = getattr(args, "information_metric", "absolute")
    result = aggregate_path_scores(components, information_metric=information_metric)
    complete = complete_cluster_fit(result, expected)
    report = paired_score_reporting(components)
    row = dict(test_id=test_id, block_id=event_id, gene_id=gene_id, contrast="cell_type_pairwise", contrast_id=f"cell_type__{tested_levels[0]}__{tested_levels[1]}", effect="cell_type", level_a=tested_levels[0], level_b=tested_levels[1], method="Tealeaf EC mixed score; SUPPA2 event definitions", inference_backend="mixed-score", model_version=MODEL_VERSION, path_pseudocount=0., path_prior_center="none", path_pseudocount_scaling="total", retain_uncertainty=True, uncertainty_scale=1., n_paths=2, n_isoforms=base.n_isoforms, n_source_isoforms=base.n_isoforms, n_ecs=n_ecs, n_samples=len(labels), n_expected_subjects=expected, n_subjects=result["n_subjects"], n_fitted_subjects=result["n_fitted_subjects"], degrees_of_freedom=result["degrees_of_freedom"], denominator_degrees_of_freedom=result["denominator_degrees_of_freedom"], median_gene_umis=float(np.median(totals)), statistic=result["statistic"] if complete else 0., p_value=result["p_value"] if complete else 1., chi_square_p_value=result["chi_square_p_value"], biological_variance=result["biological_variance"], restricted_objective=result["restricted_objective"], residual_inflation=result["residual_inflation"], converged=complete, complete_subject_fits=complete, complete_reporting_fits=report["complete"], mean_difference_norm=float(np.linalg.norm(result["mean_difference"])), effect_size=float(report["effect"][0]), effect_coordinate="psi", test_ilr_effect_size=float(result["mean_difference"][0]), report_psi_effect=float(report["effect"][0]), report_n_subjects=report["n_reported_subjects"], report_pseudocount=args.report_pseudocount, profile_event_mass=True, baseline_event_mass=float(baseline[path_index >= 0].sum()), event_type=event.event_type, event_id=event_id, feature_id=event.feature_id)
    nulls = []
    row["information_metric"] = information_metric
    if complete:
        test_hash = zlib.crc32(test_id.encode("utf-8"))
        for replicate in range(args.null_replicates):
            rng = np.random.default_rng(np.random.SeedSequence((args.seed, test_hash, replicate)))
            nulls.append(dict(test_id=test_id, block_id=event_id, replicate=replicate, p_value=signed_path_score_p_value(components, rng, information_metric=information_metric)))
    usage = []
    if args.export_path_usage:
        for subject, reports in zip(components.subject_ids, components.reporting_proportions):
            for level, proportions in reports:
                usage.append(dict(test_id=test_id, feature_id=event.feature_id, contrast_id=row["contrast_id"], subject=subject, cell_type=tested_levels[int(level)], inclusion=float(proportions[0])))
    return row, nulls, usage


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--event-catalog", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260927)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--path-pseudocount", type=float, default=32.0)
    parser.add_argument("--path-prior-center", choices=("uniform", "baseline"), default="uniform")
    parser.add_argument("--path-pseudocount-scaling", choices=("per_path", "total"), default="total")
    parser.add_argument("--max-tests", type=int, help="limit tests for a smoke run")
    parser.add_argument("--scan-all-events", action="store_true", help="screen the full supported gene transcript set rather than the block candidate gene set")
    parser.add_argument("--profile-event-mass", action="store_true", help="Estimate included-plus-excluded mass separately for each subject and cell type.")
    parser.add_argument("--report-pseudocount", type=float, help="Optional separate smoothing strength for descriptive effect estimates; testing is unchanged.")
    parser.add_argument("--export-path-usage", action="store_true")
    parser.add_argument("--inference", choices=("paired", "mixed-score"), default="paired", help="Experimental mixed score profiles all within-path transcript shares; the default paired procedure is unchanged.")
    parser.add_argument("--null-multistart", action="store_true", help="Experimental mixed-score null fit with pooled and interior starts.")
    parser.add_argument("--export-score-components", action=argparse.BooleanOptionalAction, default=True, help="Archive fitted binary subject scores for inexpensive rank-rule reassessment; no inference change.")
    parser.add_argument("--information-metric", choices=("absolute", "reference"), default="absolute", help="Experimental mixed-score numerical rank rule, not a change to effect units or priors.")
    return parser.parse_args()


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    if args.inference == "mixed-score" and (args.report_pseudocount is None or not args.profile_event_mass):
        raise ValueError("mixed score requires explicit independent reporting strength and profiled event mass")
    if args.null_multistart and args.inference != "mixed-score":
        raise ValueError("null multistart is only available for mixed score")
    if args.information_metric != "absolute" and args.inference != "mixed-score":
        raise ValueError("target-normalized information rank is only available for mixed score")
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    if settings.get("test_effect") != "cell_type_pairwise":
        raise ValueError("the hybrid currently requires paired cell-type candidates")
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(
        args.data_cache, settings
    )
    with args.data_cache.open("rb") as handle:
        _, _, genes, gene_transcripts, _, _ = pickle.load(handle)
    features = args.features.read_text().splitlines()
    if len(features) != max(int(indices.max()) for indices in gene_transcripts if len(indices)) + 1:
        raise ValueError("feature list does not align with prepared transcript indices")
    catalog = pd.read_csv(
        args.event_catalog, sep="\t", compression="infer", dtype=str
    ).fillna("")
    events_by_gene = {}
    for event in catalog.itertuples(index=False):
        events_by_gene.setdefault(canonical(event.gene_id), []).append(event)

    # A gene/subject-fold/contrast context is shared by all standard events for
    # that gene. Collapse the many annotation blocks to avoid duplicate tests.
    tests = []
    represented_event_ids = set()
    contexts = {}
    if args.scan_all_events:
        gene_lookup = {canonical(gene): index for index, gene in enumerate(genes)}
        screening_counts = tuple(matrix.tocsc() for matrix in counts)
        candidate_genes = sorted({
            gene_lookup[gene_id]
            for gene_id in events_by_gene
            if gene_id in gene_lookup
        })
        for gene in candidate_genes:
            ecs = np.asarray(gene_ecs[gene], dtype=int)
            if not len(ecs) or len(ecs) > int(settings.get("max_ecs", 128)):
                continue
            transcripts = supported_gene_transcripts(
                gene, gene_transcripts, gene_ecs, designs
            )
            if len(transcripts) < 2:
                continue
            gene_umis = modeled_gene_umis(
                screening_counts, designs, ecs, transcripts
            )
            coverage_specs = covered_celltype_pairwise_designs(
                metadata,
                gene_umis,
                min_gene_umis=float(settings.get("min_gene_umis", 10.0)),
                min_samples=int(settings.get("min_gene_samples", 0)),
                min_celltype_mice=int(settings.get("min_celltype_mice", 3)),
            )
            gene_id = str(genes[gene])
            for coverage, _ in coverage_specs:
                rows, _, _, _, _, levels = coverage
                candidate = (
                    "", "", gene_id, gene, transcripts, None, [], rows, None, tuple(levels)
                )
                context_key = (gene, tuple(np.asarray(rows, dtype=int)), tuple(levels))
                contexts.setdefault(context_key, candidate)
                for event in events_by_gene.get(canonical(gene_id), ()):
                    path_index = event_path_index(
                        transcripts, features, event.included, event.excluded
                    )
                    if path_index is not None:
                        tests.append((candidate, event, path_index))
                        represented_event_ids.add(str(event.event_id))
    else:
        for candidate in cached["candidates"]:
            _, _, _, gene, _, _, _, rows, _, levels = candidate
            key = (int(gene), tuple(np.asarray(rows, dtype=int)), tuple(levels))
            contexts.setdefault(key, candidate)
        supported_transcript_cache = {}
        for (gene, rows, levels), candidate in contexts.items():
            if gene not in supported_transcript_cache:
                supported_transcript_cache[gene] = supported_gene_transcripts(
                    gene, gene_transcripts, gene_ecs, designs
                )
            transcripts = candidate[4]
            gene_id = canonical(genes[gene])
            for event in events_by_gene.get(gene_id, ()):
                path_index = event_path_index(
                    transcripts, features, event.included, event.excluded
                )
                if path_index is not None:
                    tests.append((candidate, event, path_index))
                    represented_event_ids.add(str(event.event_id))
    skipped_events = len(catalog) - len(represented_event_ids)
    tests = partition_event_tests(tests, args.shard_count)[args.shard_index]
    if args.max_tests is not None:
        tests = tests[: args.max_tests]

    row_lookup = {}
    observed, null, failures, usage = [], [], [], []
    score_archive = ([], []) if args.inference == "mixed-score" and args.export_score_components else None
    baseline_cache = {}
    active_score_context, score_component_cache = None, {}
    started = time.perf_counter()
    for test_number, (candidate, event, path_index) in enumerate(tests):
        if args.inference == "mixed-score" and test_number % 100 == 0:
            print(f"mixed-score shard {args.shard_index}, {test_number}/{len(tests)} requested tests, {len(observed)} completed records, {len(failures)} exceptions, {time.perf_counter() - started:.1f}s", flush=True)
        test_id, _, gene_id, gene, _, _, _, rows, _, tested_levels = candidate
        event_id = str(event.event_id)
        event_test_id = f"SUPPA2:{event_id}|cell_type|{'|'.join(tested_levels)}"
        try:
            local_metadata, _, labels = local_test_design(
                metadata, rows, tested_levels, "cell_type_pairwise"
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata.mouse.astype(str).to_numpy()
            transcripts = candidate[4]
            base, _, totals = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(local_metadata), 1)),
                clusters,
                drop_zero=False,
            )
            cache_key = (gene, tuple(rows))
            if cache_key not in baseline_cache:
                baseline_cache[cache_key] = ec_block_glmm.pooled_isoform_weights(base)
            if args.inference == "mixed-score":
                if cache_key != active_score_context:
                    active_score_context, score_component_cache = cache_key, {}
                partition = tuple(path_index)
                if partition not in score_component_cache:
                    score_component_cache[partition] = shared_path_score_components(base, path_index, labels, clusters, baseline=baseline_cache[cache_key], max_iter=args.max_iter, null_multistart=args.null_multistart, reporting_concentration=args.report_pseudocount)
                record, generated_null, generated_usage = mixed_event_record(base, path_index, labels, clusters, baseline_cache[cache_key], event, gene_id, tested_levels, totals, len(gene_ecs[gene]), args, score_archive=score_archive, components=score_component_cache[partition])
                observed.append(record)
                null.extend(generated_null)
                usage.extend(generated_usage)
                continue
            fit_base, fit_path_index, fit_baseline = collapse_event_nuisance(
                base, path_index, baseline_cache[cache_key]
            )
            result = ec_block_glmm.paired_path_test(
                fit_base,
                fit_path_index,
                labels,
                clusters,
                baseline=fit_baseline,
                max_iter=args.max_iter,
                path_pseudocount=args.path_pseudocount,
                path_prior_center=args.path_prior_center,
                path_pseudocount_scaling=args.path_pseudocount_scaling,
                profile_event_mass=args.profile_event_mass,
            )
            values = result["differences"]
            covariances = result["difference_covariances"]
            mean_event_ilr_difference = float(values.mean(axis=0)[0]) if len(values) else np.nan
            report = result
            if args.report_pseudocount is not None:
                report = ec_block_glmm.paired_path_test(fit_base, fit_path_index, labels, clusters, baseline=fit_baseline, max_iter=args.max_iter, path_pseudocount=args.report_pseudocount, path_prior_center=args.path_prior_center, path_pseudocount_scaling=args.path_pseudocount_scaling, profile_event_mass=args.profile_event_mass)
            report_values = report["differences"]
            report_ilr = float(report_values.mean(axis=0)[0]) if len(report_values) else np.nan
            psi_differences = np.array([pair[1].path_proportions[0] - pair[0].path_proportions[0] for pair in report["path_fits"]])
            report_psi = float(psi_differences.mean()) if len(psi_differences) else np.nan
            if args.export_path_usage:
                for subject, fits in zip(report["subject_ids"], report["path_fits"]):
                    for level, fit in zip(report["levels"], fits):
                        usage.append({"test_id": event_test_id, "feature_id": event.feature_id, "contrast_id": f"cell_type__{tested_levels[0]}__{tested_levels[1]}", "subject": subject, "cell_type": tested_levels[int(level)], "inclusion": float(fit.path_proportions[0]), "event_mass": float(fit.theta[:2].sum()), "ilr": float(fit.path_logratios[0])})
            observed.append({
                "test_id": event_test_id,
                "block_id": event_id,
                "gene_id": gene_id,
                "contrast": "cell_type_pairwise",
                "contrast_id": f"cell_type__{tested_levels[0]}__{tested_levels[1]}",
                "effect": "cell_type",
                "level_a": tested_levels[0],
                "level_b": tested_levels[1],
                "method": "Tealeaf EC; SUPPA2 event definitions",
                "path_pseudocount": args.path_pseudocount,
                "path_prior_center": args.path_prior_center,
                "path_pseudocount_scaling": args.path_pseudocount_scaling,
                "retain_uncertainty": False,
                "uncertainty_scale": 0.0,
                "n_paths": 2,
                "n_isoforms": fit_base.n_isoforms,
                "n_source_isoforms": base.n_isoforms,
                "n_ecs": len(gene_ecs[gene]),
                "n_samples": result.get("n_observations", len(local_metadata)),
                "n_subjects": result["n_subjects"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "median_gene_umis": float(np.median(totals)),
                "statistic": result["statistic"],
                "p_value": result["p_value"],
                "biological_variance": result.get("biological_variance", np.nan),
                "restricted_objective": result.get("restricted_objective", np.nan),
                "converged": result["converged"],
                "mean_difference_norm": float(np.linalg.norm(values.mean(axis=0))) if len(values) else 0.0,
                "effect_size": report_psi if args.report_pseudocount is not None else mean_event_ilr_difference,
                "effect_coordinate": "psi" if args.report_pseudocount is not None else "ilr",
                "test_ilr_effect_size": mean_event_ilr_difference,
                "report_ilr_effect": report_ilr,
                "report_psi_effect": report_psi,
                "report_n_subjects": report["n_subjects"],
                "report_pseudocount": args.report_pseudocount if args.report_pseudocount is not None else args.path_pseudocount,
                "profile_event_mass": args.profile_event_mass,
                "baseline_event_mass": float(fit_baseline[:2].sum()),
                "event_type": event.event_type,
                "event_id": event_id,
                "feature_id": event.feature_id,
            })
            if result["converged"]:
                test_hash = zlib.crc32(event_test_id.encode("utf-8"))
                for replicate in range(args.null_replicates):
                    rng = np.random.default_rng(
                        np.random.SeedSequence((args.seed, test_hash, replicate))
                    )
                    null.append({
                        "test_id": event_test_id,
                        "block_id": event_id,
                        "replicate": replicate,
                        "p_value": signed_null_p_value(
                            values, covariances, rng, False, 0.0
                        ),
                    })
        except Exception as error:  # retain failed IDs for coverage accounting
            failures.append({"test_id": event_test_id, "error": repr(error)})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(observed)
    if args.inference == "mixed-score" and table.empty:
        table = pd.DataFrame(columns=MIXED_SCORE_EMPTY_COLUMNS)
    table.to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(null).to_csv(args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False)
    if args.export_path_usage:
        pd.DataFrame(usage).to_csv(args.output_dir / "path_usage.tsv.gz", sep="\t", index=False)
    if score_archive is not None:
        pd.DataFrame(score_archive[0]).to_csv(args.output_dir / "score_contexts.tsv.gz", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(score_archive[1]).to_csv(args.output_dir / "subject_scores.tsv.gz", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidate_contexts": len(contexts),
        "catalog_events": len(catalog),
        "represented_events": len(represented_event_ids),
        "skipped_events": skipped_events,
        "tests_in_shard": len(tests),
        "completed": len(observed),
        "failures": len(failures),
        "elapsed_seconds": time.perf_counter() - started,
    }, indent=2) + "\n")
    if args.inference == "mixed-score":
        (args.output_dir / "settings.json").write_text(json.dumps({"arguments": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}, "candidate_settings": settings, "model_version": MODEL_VERSION, "target_prior": "none", "nuisance_pseudocount": 1e-4, "production_changes": False}, indent=2) + "\n")
    print(f"wrote {len(observed):,} tests and {len(null):,} nulls; skipped {skipped_events:,} unsupported catalogue events")


if __name__ == "__main__":
    main()
