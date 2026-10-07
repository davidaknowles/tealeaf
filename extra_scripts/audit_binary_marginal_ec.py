#!/usr/bin/env python3
"""Runtime/accuracy pilot of binary direct-EC random-subject inference.

Selection is from covered production-eligible binary tests, not discoveries or
LR agreement. This pilot is not a calibrated full-family endpoint comparison.
"""

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_marginal import MODEL_VERSION, binary_marginal_test, prepare_binary_ec_likelihood
from tealeaf.sc.path_marginal import random_subject_objective
from tealeaf.sc.differential import path_proportions
from tealeaf.sc.path_simulation import simulate_counts


def main(*, likelihood_transform=None, quadrature_backend=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=4)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=4)
    parser.add_argument("--subject-nodes", type=int, default=48)
    parser.add_argument("--path-nodes", type=int, default=32)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--draws", type=int, default=1)
    parser.add_argument("--only-null", action="store_true", help="Repeat count-level nulls without observed-data fitting.")
    parser.add_argument("--diagnostics", action="store_true", help="Retain generating support and separate path/subject quadrature errors.")
    args = parser.parse_args()
    if args.blocks < 1 or args.draws < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("invalid panel size or shard")
    scenarios = ([] if args.only_null else [("observed", np.inf, -1)]) + [(name, precision, draw) for draw in range(args.draws) for name, precision in (("common composition", np.inf), ("biological precision 20", 20.))]
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    reference = pd.read_csv(args.reference, sep="\t", low_memory=False)
    eligible = set(reference.loc[reference.converged.astype(str).str.lower().eq("true") & reference.n_subjects.ge(4), "test_id"])
    candidates = [candidate for candidate in cached["candidates"] if candidate[0] in eligible and len(candidate[6]) == 2]
    rng = np.random.default_rng(311381)
    chosen = sorted(rng.choice(len(candidates), min(len(candidates), args.blocks), replace=False))
    selected_ids = [candidates[index][0] for index in chosen]
    candidates = partition_candidates([candidates[index] for index in chosen], args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    output = []
    for candidate in candidates:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id}
        try:
            local_metadata, _, labels = local_test_design(metadata, rows, tested_levels, cached["settings"]["test_effect"])
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            baseline, converged = pooled_isoform_weights(base, max_iter=250, return_status=True)
            if not converged:
                raise ValueError("initial pooled baseline did not converge")
            if args.diagnostics:
                header.update(baseline_min_path_usage=float(path_proportions(baseline, path_index).min()), baseline_block_mass=float(baseline[np.asarray(path_index) >= 0].sum()))
        except (ValueError, np.linalg.LinAlgError) as exception:
            for scenario, _, draw in scenarios:
                output.append({**header, "scenario": scenario, "draw": draw, "p_value": 1., "converged": False, "error": str(exception)})
            continue
        for scenario, concentration, draw in scenarios:
            start = time.monotonic()
            try:
                if scenario == "observed":
                    generated, nuisance = base, baseline
                else:
                    suffix = f" draw{draw}" if draw > 0 else ""
                    rng = np.random.default_rng(311381 + zlib.crc32((test_id + scenario + suffix).encode()))
                    generated = simulate_counts(base, baseline, subjects, rng, .5, labels=labels, path_index=path_index, residual_concentration=None if np.isinf(concentration) else concentration, return_details=args.diagnostics)
                    if args.diagnostics:
                        generated, generating = generated
                    nuisance, converged = pooled_isoform_weights(generated, max_iter=250, return_status=True)
                    if not converged:
                        raise ValueError("count-draw pooled baseline did not converge")
                likelihood = prepare_binary_ec_likelihood(generated, path_index, labels, subjects, baseline=nuisance)
                if likelihood_transform is not None:
                    likelihood = likelihood_transform(likelihood)
                fitted = binary_marginal_test(likelihood, max_iter=args.max_iter, subject_nodes=args.subject_nodes, path_nodes=args.path_nodes)
                result = {key: value for key, value in fitted.items() if key not in ("null_fit", "alternative_fit", "levels", "standardized_means")}
                result.update(levels=json.dumps(fitted["levels"].tolist()), standardized_means=json.dumps(fitted["standardized_means"].tolist()), null_optimizer_success=bool(fitted["null_fit"].success), alternative_optimizer_success=bool(fitted["alternative_fit"].success), exact_integral_rows=int(np.isfinite(likelihood.binomial_terms).all(axis=1).sum()), error="")
                if args.diagnostics:
                    result.update(null_parameters=json.dumps(fitted["null_fit"].x.tolist()), alternative_parameters=json.dumps(fitted["alternative_fit"].x.tolist()), proposal_effective_depth_median=float(np.median(likelihood.proposal_counts.sum(axis=1))), gene_ec_depth_median=float(np.median(sum(values.sum(axis=1) for values in likelihood.counts))))
                    if scenario != "observed":
                        means = np.asarray([path_proportions(weight, path_index) for weight in generating["subject_weights"]])
                        result.update(true_minimum_subject_path_mean=float(means.min()), true_median_subject_minimum_path_mean=float(np.median(means.min(axis=1))), true_block_mass_median=float(np.median(generating["subject_weights"][:, np.asarray(path_index) >= 0].sum(axis=1))))
                    design = np.column_stack([np.ones(len(likelihood.labels)), *[likelihood.labels == level for level in fitted["levels"][1:]]]).astype(float)
                    for name, subject_order, path_order in (("subject_only", 2 * args.subject_nodes, args.path_nodes), ("path_only", args.subject_nodes, 2 * args.path_nodes)):
                        null_value = random_subject_objective(fitted["null_fit"].x, likelihood, design[:, :1], subject_nodes=subject_order, path_nodes=path_order, gradient=False)
                        alternative_value = random_subject_objective(fitted["alternative_fit"].x, likelihood, design, subject_nodes=subject_order, path_nodes=path_order, gradient=False)
                        result[name + "_quadrature_error"] = max(abs(null_value - fitted["null_fit"].fun), abs(alternative_value - fitted["alternative_fit"].fun), 2 * abs((null_value - alternative_value) - (fitted["null_fit"].fun - fitted["alternative_fit"].fun)))
            except (ValueError, np.linalg.LinAlgError) as exception:
                result = {"p_value": 1., "converged": False, "error": str(exception)}
            output.append({**header, "scenario": scenario, "draw": draw, "runtime_seconds": time.monotonic() - start, **result})
            print(f"{test_id}, {scenario}, draw{draw}, converged={result['converged']}, elapsed={output[-1]['runtime_seconds']:.1f}s", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    settings = {"candidate_settings": cached["settings"], "model_version": MODEL_VERSION, "selected_ids": selected_ids, "shard_count": args.shard_count, "requested_scenarios": [{"scenario": scenario, "draw": draw} for scenario, _, draw in scenarios], "subject_nodes": args.subject_nodes, "path_nodes": args.path_nodes, "max_iter": args.max_iter, "failure_policy": "p=1, retained cases", "scope": "binary runtime/accuracy or count-null panel only, no calibrated or full-family endpoint claim", "production_changes": False}
    if quadrature_backend is not None:
        settings["quadrature_backend"] = quadrature_backend
    if args.diagnostics:
        settings["diagnostics"] = "generating support, raw fit coefficients, and subject-only/path-only doubled-order errors; proposals are not true observed path counts"
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")


if __name__ == "__main__":
    main()
