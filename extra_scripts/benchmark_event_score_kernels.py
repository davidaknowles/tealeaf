"""Measure equivalent event-score kernels, without ranking or selecting hits."""

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, local_test_design, modeled_gene_umis
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts, event_path_index
from tealeaf.sc.ec_glmm import subset_gene_data
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_bias import SharedPathNullProblem
from tealeaf.sc.path_score_mixed import mixed_score_test, shared_path_score_components


def timed(function, repetitions):
    function()
    started = time.perf_counter()
    for _ in range(repetitions):
        result = function()
    return result, (time.perf_counter() - started) / repetitions


def compare_scalar(scores, information, shapes, reference, context, rng, null_draws):
    """Identical scores and sign draws, include unavailable cases rather than selecting them."""
    rows = []
    signs = np.vstack([np.ones(len(scores)), rng.choice((-1., 1.), size=(null_draws, len(scores)))])
    for metric in ("absolute", "reference"):
        target = reference if metric == "reference" else None
        for draw, sign in enumerate(signs):
            arguments = (scores * sign[:, None], information, shapes)
            results, elapsed, errors = [], [], []
            for fast in (False, True):
                start = time.perf_counter()
                try:
                    results.append(mixed_score_test(*arguments, reference_information=target, scalar_fast=fast))
                    errors.append(None)
                except (ValueError, np.linalg.LinAlgError) as error:
                    results.append(None)
                    errors.append(str(error))
                elapsed.append(time.perf_counter() - start)
            row = dict(context=context, metric=metric, draw=draw, n_subjects=len(scores), generic_seconds=elapsed[0], scalar_seconds=elapsed[1], same_availability=(results[0] is None) == (results[1] is None), generic_error=errors[0], scalar_error=errors[1])
            if results[0] is not None and results[1] is not None:
                for key in ("p_value", "statistic", "mean_difference", "mean_covariance", "restricted_objective", "biological_variance"):
                    first, second = np.asarray(results[0][key]), np.asarray(results[1][key])
                    row[key + "_max_absolute_difference"] = float(np.max(np.abs(first - second)))
                    row[key + "_matches"] = bool(np.allclose(first, second, rtol=3e-6, atol=1e-8))
                row["p_generic"], row["p_scalar"] = results[0]["p_value"], results[1]["p_value"]
            rows.append(row)
    return rows


def summarize(output):
    inputs = pd.read_csv(output / "input_selection.tsv", sep="\t")
    objective = pd.read_csv(output / "null_objective.tsv", sep="\t")
    scalar = pd.read_csv(output / "scalar_reml.tsv.gz", sep="\t")
    complete = scalar.p_value_matches.notna()
    row = dict(contexts=len(inputs), input_median_speedup=float(inputs.speedup.median()), input_total_speedup=float(inputs.row_first_seconds.sum() / inputs.column_first_seconds.sum()), objective_median_speedup=float(objective.speedup.median()), dense_fit_converged=int(objective.dense_fit_converged.sum()), vector_fit_converged=int(objective.vector_fit_converged.sum()), fit_max_absolute_objective_difference=float((objective.dense_fit_objective - objective.vector_fit_objective).abs().max()), scalar_trials=len(scalar), scalar_complete=int(complete.sum()), same_availability=bool(scalar.same_availability.all()), scalar_total_speedup=float(scalar.generic_seconds.sum() / scalar.scalar_seconds.sum()), scalar_median_speedup=float((scalar.generic_seconds / scalar.scalar_seconds).median()))
    for key in ("p_value", "statistic", "mean_difference", "mean_covariance", "restricted_objective", "biological_variance"):
        row[key + "_all_match"] = bool(scalar.loc[complete, key + "_matches"].all())
        row[key + "_max_absolute_difference"] = float(scalar[key + "_max_absolute_difference"].max())
    (output / "summary.json").write_text(json.dumps(row, indent=2) + "\n")
    print(json.dumps(row, indent=2), flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--alias-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--contexts", type=int, default=12)
    parser.add_argument("--null-draws", type=int, default=32)
    args = parser.parse_args()
    if args.contexts < 1 or args.null_draws < 1:
        raise ValueError("positive context and null draw counts required")
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(args.cache / "original_binary_paired/prepared.pkl", settings)
    features = (args.cache / "original_binary_paired/features.txt").read_text().splitlines()
    screening = tuple(value.tocsc() for value in counts)
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    audit = pd.read_csv(args.alias_audit, sep="\t")
    audit = audit.loc[audit.source.eq("original_binary")].drop_duplicates("gene_id")
    # Fixed-random genes plus dimension extremes, never significance or LR support.
    selected = pd.concat([audit.sample(min(args.contexts, len(audit)), random_state=817), audit.nsmallest(1, "transcripts"), audit.nlargest(1, "transcripts")]).drop_duplicates("gene_id")
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    input_rows, objective_rows, scalar_rows, failures = [], [], [], []
    rng = np.random.default_rng(751)
    for context, record in enumerate(selected.itertuples()):
        gene = lookup[canonical(record.gene_id)]
        transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
        umis = modeled_gene_umis(screening, designs, gene_ecs[gene], transcripts)
        specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
        coverage = specs[context % len(specs)][0]
        rows, levels = coverage[0], coverage[-1]
        local, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
        subjects = local.mouse.astype(str).to_numpy()
        fixed = np.ones((len(rows), 1))
        common = (designs, transcripts, gene_ecs[gene], fixed, subjects)
        old, old_time = timed(lambda: subset_gene_data(tuple(value[rows] for value in counts), *common, drop_zero=False), 10)
        new, new_time = timed(lambda: subset_gene_data(screening, *common, rows=rows, drop_zero=False), 10)
        for first, second in zip(old[0].counts + old[0].compatibility + (old[1], old[2]), new[0].counts + new[0].compatibility + (new[1], new[2])):
            np.testing.assert_array_equal(first, second)
        input_rows.append(dict(gene_id=record.gene_id, feature_id=record.feature_id, transcripts=len(transcripts), n_ecs=len(gene_ecs[gene]), n_rows=len(rows), row_first_seconds=old_time, column_first_seconds=new_time, speedup=old_time / new_time, bit_identical=True))
        event = catalog.loc[record.feature_id]
        paths = event_path_index(transcripts, features, event.included, event.excluded)
        base = new[0]
        baseline = pooled_isoform_weights(base)
        subject = np.unique(subjects)[0]
        local_counts = tuple(np.array([matrix[(subjects == subject) & (labels == level)].sum(axis=0) for level in np.unique(labels)]) for matrix in base.counts)
        problem = SharedPathNullProblem(local_counts, base.compatibility, baseline, paths)
        parameters = problem.initial + rng.normal(scale=.15, size=problem.dimension)
        dense, dense_time = timed(lambda: problem.objective(parameters), 200)
        vector, vector_time = timed(lambda: problem.objective_vector(parameters), 200)
        np.testing.assert_allclose(vector[0], dense[0], rtol=1e-12, atol=1e-8)
        np.testing.assert_allclose(vector[1], dense[1], rtol=1e-8, atol=1e-8)
        fit_rows = dict(gene_id=record.gene_id, feature_id=record.feature_id, transcripts=len(transcripts), dimension=problem.dimension, dense_seconds=dense_time, vector_seconds=vector_time, speedup=dense_time / vector_time, objective_absolute_difference=abs(dense[0] - vector[0]), gradient_max_absolute_difference=float(np.max(np.abs(dense[1] - vector[1]))))
        for method in ("dense", "vector"):
            start = time.perf_counter()
            fit = problem.fit(max_iter=2000, multistart=True, objective_method=method)
            fit_rows[method + "_fit_seconds"] = time.perf_counter() - start
            fit_rows[method + "_fit_objective"] = fit.objective
            fit_rows[method + "_fit_iterations"] = fit.iterations
            fit_rows[method + "_fit_converged"] = fit.converged
        objective_rows.append(fit_rows)
        try:
            components = shared_path_score_components(base, paths, labels, subjects, baseline=baseline, max_iter=2000, null_multistart=True)
            scalar_rows.extend(compare_scalar(components.scores, components.information, components.biological_shapes, components.reference_information, record.feature_id, rng, args.null_draws))
        except (ValueError, np.linalg.LinAlgError) as error:
            failures.append(dict(gene_id=record.gene_id, feature_id=record.feature_id, error=repr(error)))
        print(f"context {context + 1}/{len(selected)}, T={len(transcripts)}, input speedup={old_time / new_time:.2f}, objective speedup={dense_time / vector_time:.2f}", flush=True)
    for subjects in (4, 8, 16, 32):
        for spread in (1., 6.):
            info = np.exp(rng.normal(0., spread, subjects))[:, None, None]
            shapes = np.exp(rng.normal(0., 2., subjects))[:, None, None]
            scores = rng.normal(.3, 1., (subjects, 1)) * info[:, 0]
            reference = info * rng.uniform(1., 100., (subjects, 1, 1))
            scalar_rows.extend(compare_scalar(scores, info, shapes, reference, f"synthetic_M{subjects}_spread{spread}", rng, args.null_draws))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(input_rows).to_csv(args.output_dir / "input_selection.tsv", sep="\t", index=False)
    pd.DataFrame(objective_rows).to_csv(args.output_dir / "null_objective.tsv", sep="\t", index=False)
    pd.DataFrame(scalar_rows).to_csv(args.output_dir / "scalar_reml.tsv.gz", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(selection="fixed-random whole-catalog supported genes and minimum/maximum transcript counts, no significance or LR selection", random_gene_seed=817, numerical_seed=751, null_draws=args.null_draws, candidate_settings=settings, input_assertion="exact count, map, mask and total array equality", purpose="runtime and numerical equivalence, not endpoint performance", production_statistical_changes=False), indent=2) + "\n")
    summarize(args.output_dir)


if __name__ == "__main__":
    main()
