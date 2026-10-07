#!/usr/bin/env python3
"""Assess paired reporting estimators, independently of external outcomes."""

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_paired_path_test import filtered_inputs
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_ec_block_glmm import local_test_design, partition_candidates
from tealeaf.sc import ec_block_glmm, differential
from tealeaf.sc.path_reporting import paired_reporting, proportion_covariance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=32)
    parser.add_argument("--concentrations", default="1;4;16;32;64")
    parser.add_argument("--free-isoforms", action="store_true")
    args = parser.parse_args()
    if args.free_isoforms:
        differential.fit_path_perturbation = differential.fit_free_isoform_paths
    concentrations = [float(value) for value in args.concentrations.replace(",", ";").split(";")]
    if any(value <= 0 or not np.isfinite(value) for value in concentrations):
        raise ValueError("positive reporting concentrations required")
    reference = pd.read_csv(args.reference, sep="\t", low_memory=False)
    if "converged" in reference:
        reference = reference.loc[reference.converged.astype(str).str.lower().eq("true") & reference.n_subjects.ge(4)]
    requested = set(reference.test_id)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    candidates = partition_candidates([candidate for candidate in cached["candidates"] if candidate[0] in requested], args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    outputs, observations, failures, baseline_cache = [], [], [], {}
    started = time.monotonic()
    for number, candidate in enumerate(candidates):
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "path_signatures": json.dumps(signatures), "n_paths": len(signatures), "level_a": levels[0], "level_b": levels[1]}
        try:
            local_metadata, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            key = (gene, tuple(rows), tuple(transcripts))
            if key not in baseline_cache:
                baseline_cache[key] = ec_block_glmm.pooled_isoform_weights(base)
            for concentration in concentrations:
                result = ec_block_glmm.paired_path_test(base, path_index, labels, subjects, baseline=baseline_cache[key], path_pseudocount=concentration, path_pseudocount_scaling="total")
                fits = result["path_fits"]
                proportions = np.asarray([[fit.path_proportions for fit in pair] for pair in fits])
                if len(fits) < 2:
                    failures.append({**header, "concentration": concentration, "error": "fewer than two fitted subject pairs"})
                    continue
                covariances = np.asarray([[proportion_covariance(fit) for fit in pair] for pair in fits])
                for index, pair in enumerate(fits):
                    for level, fit in enumerate(pair):
                        if not fit.covariance.identifiable:
                            covariances[index, level] = np.nan
                depths = np.asarray([[sum(float(values[(subjects == subject) & (labels == level)].sum()) for values in base.counts) for level in result["levels"]] for subject in result["subject_ids"]])
                observations.append({**header, "concentration": concentration, "subjects": json.dumps(result["subject_ids"].tolist()), "proportions": json.dumps(proportions.tolist()), "covariances": json.dumps(covariances.tolist()), "depths": json.dumps(depths.tolist())})
                for name, report in paired_reporting(proportions, covariances, depths).items():
                    weights = report["weights"]
                    effect = report["effect"]
                    outputs.append({**header, "strategy": f"{name} A{concentration:g}", "concentration": concentration, "effect": json.dumps(effect.tolist()), "converged": bool(np.isfinite(effect).all()), "report_fallback": report.get("fallback", False), "n_subjects": len(fits), "between_subject_variance": report.get("between_subject_variance", np.nan), "max_subject_weight": float(weights.max()), "effective_subjects": float(1 / np.square(weights).sum()), "error": report.get("error", "")})
        except (ValueError, np.linalg.LinAlgError) as exception:
            failures.append({**header, "error": str(exception)})
        if number % 25 == 0:
            print(f"{number + 1}/{len(candidates)} candidates, {time.monotonic() - started:.1f}s", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output_dir / "observed.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(observations).to_csv(args.output_dir / "quantified.tsv.gz", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "settings.json").write_text(json.dumps({"candidate_settings": cached["settings"], "n_candidates": len(candidates), "n_failures": len(failures), "concentrations": concentrations, "free_isoforms": args.free_isoforms, "statistical_tests": "frozen production reference; reporting only", "baseline": "independently refitted within this fold", "external_outcomes_used_by_estimator": False}, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
