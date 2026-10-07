#!/usr/bin/env python3
"""Fixed-truth EC resampling, covariance matching and subject-nuisance audit."""

import argparse
import json
from pathlib import Path
import pickle
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc.ec_block_glmm import pooled_isoform_weights
from tealeaf.sc.path_pooling import quantify_effective_paths, subject_isoform_baselines
from tealeaf.sc.path_simulation import simulate_counts, resample_counts
from tealeaf.sc.differential import path_proportions


def summarize_bootstrap(records, truth, requested):
    """One record per subject/type, retaining incomplete-bootstrap status."""
    if len(records) < 2:
        return []
    keys = [(str(subject), str(label)) for subject, label in zip(records[0]["subjects"], records[0]["labels"])]
    if any(keys != [(str(subject), str(label)) for subject, label in zip(record["subjects"], record["labels"])] for record in records):
        raise ValueError("conditional resampling changed the observation family")
    values = np.asarray([record["proportions"] for record in records])
    original = np.asarray([record["proportion_covariances"] for record in records])
    scalar = np.asarray([record["scalar_proportion_covariances"] for record in records])
    centered = values - values.mean(axis=0)
    covariance = np.einsum("rni,rnj->nij", centered, centered) / (len(records) - 1)
    rows = []
    for index, key in enumerate(keys):
        expected = truth[key]
        observed = covariance[index]
        variance = np.trace(observed)
        bias = values[:, index].mean(axis=0) - expected
        fisher = original[:, index].mean(axis=0)
        matched = scalar[:, index].mean(axis=0)
        rows.append({"subject": key[0], "label": key[1], "bootstrap_requested": requested, "bootstrap_converged": len(records), "complete_bootstrap": len(records) == requested, "minimum_true_usage": expected.min(), "empirical_trace": variance, "fisher_trace": np.trace(fisher), "scalar_trace": np.trace(matched), "fisher_variance_ratio": np.trace(fisher) / variance if variance > 0 else np.nan, "scalar_variance_ratio": np.trace(matched) / variance if variance > 0 else np.nan, "scalar_to_fisher_trace_ratio": np.trace(matched) / np.trace(fisher) if np.trace(fisher) > 0 else np.nan, "bias_squared_norm": float(bias @ bias), "mean_squared_error": variance * (len(records) - 1) / len(records) + float(bias @ bias), "effective_depth_median": np.median([record["effective_depths"][index] for record in records]), "ec_depth": records[0]["depths"][index], "truth": json.dumps(expected.tolist()), "bootstrap_mean": json.dumps(values[:, index].mean(axis=0).tolist())})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=32)
    parser.add_argument("--bootstrap", type=int, default=32)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--concentration", type=float, default=.25)
    parser.add_argument("--biological-concentration", type=float, default=20.)
    parser.add_argument("--subject-scale", type=float, default=.5)
    args = parser.parse_args()
    if args.bootstrap < 2 or args.blocks < 1:
        parser.error("at least two count resamples and one block required")
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    reference = pd.read_csv(args.reference, sep="\t", low_memory=False)
    eligible = set(reference.loc[reference.converged & reference.n_subjects.ge(4), "test_id"])
    candidates = [candidate for candidate in cached["candidates"] if candidate[0] in eligible]
    selected = np.random.default_rng(309151).choice(len(candidates), min(args.blocks, len(candidates)), replace=False)
    selected_ids = [candidates[index][0] for index in sorted(selected)]
    candidates = partition_candidates([candidates[index] for index in sorted(selected)], args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    outputs, failures, status = [], [], []
    for candidate in candidates:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "n_paths": len(signatures)}
        try:
            local_metadata, _, labels = local_test_design(metadata, rows, tested_levels, cached["settings"]["test_effect"])
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            baseline, converged = pooled_isoform_weights(base, max_iter=250, return_status=True)
            if not converged:
                raise ValueError("starting pooled baseline did not converge")
            rng = np.random.default_rng(309151 + zlib.crc32(test_id.encode()))
            _, details = simulate_counts(base, baseline, subjects, rng, args.subject_scale, labels=labels, path_index=path_index, residual_concentration=args.biological_concentration, return_details=True)
            oracle = dict(zip(details["subject_levels"], details["subject_weights"]))
            truth = {(str(subject), str(label)): path_proportions(weight, path_index) for subject, label, weight in zip(subjects, labels, details["observation_weights"])}
            records = {name: [] for name in ("pooled nuisance", "subject nuisance", "oracle nuisance")}
            for replicate in range(args.bootstrap):
                generated = resample_counts(base, details["observation_weights"], rng)
                for name in records:
                    try:
                        if name == "pooled nuisance":
                            fitted_baseline, converged = pooled_isoform_weights(generated, max_iter=250, return_status=True)
                            if not converged:
                                raise ValueError("pooled baseline did not converge")
                            options = {"baseline": fitted_baseline}
                        else:
                            local = oracle if name == "oracle nuisance" else subject_isoform_baselines(generated, subjects)
                            options = {"baseline": baseline, "subject_baselines": local}
                        records[name].append(quantify_effective_paths(generated, path_index, labels, subjects, concentration=args.concentration, **options))
                    except (ValueError, np.linalg.LinAlgError) as exception:
                        failures.append({**header, "strategy": name, "replicate": replicate, "error": str(exception)})
            for name, values in records.items():
                status.append({**header, "strategy": name, "bootstrap_requested": args.bootstrap, "bootstrap_converged": len(values), "complete_bootstrap": len(values) == args.bootstrap})
                outputs.extend({**header, "strategy": name, **row} for row in summarize_bootstrap(values, truth, args.bootstrap))
        except (ValueError, np.linalg.LinAlgError) as exception:
            failures.append({**header, "error": str(exception)})
            for name in ("pooled nuisance", "subject nuisance", "oracle nuisance"):
                status.append({**header, "strategy": name, "bootstrap_requested": args.bootstrap, "bootstrap_converged": 0, "complete_bootstrap": False})
        print(f"{test_id}, {len(outputs)} observation diagnostics", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output_dir / "observations.tsv.gz", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(status).to_csv(args.output_dir / "status.tsv", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "settings.json").write_text(json.dumps({"candidate_settings": cached["settings"], "selected_test_ids": selected_ids, "blocks_requested": args.blocks, "bootstrap": args.bootstrap, "concentration": args.concentration, "biological_concentration": args.biological_concentration, "subject_scale": args.subject_scale, "seed": 309151, "fixed_truth": True, "limitations": "Conditional diagnostic, not a new test; incomplete bootstrap covariance estimates are conditional on successful quantification"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
