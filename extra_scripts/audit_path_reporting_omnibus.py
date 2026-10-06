#!/usr/bin/env python3
"""Subject-split reporting and design-matched omnibus sensitivities."""

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
from tealeaf.sc import differential, ec_block_glmm, ec_glmm
from tealeaf.sc.omnibus import regression_omnibus
from tealeaf.sc.path_reporting import dirichlet_pooling, proportion_covariance


def paired_reports(base, baseline, path_index, labels, subjects):
    result = ec_block_glmm.paired_path_test(base, path_index, labels, subjects, baseline=baseline, path_pseudocount=1., path_pseudocount_scaling="total")
    if len(result["subject_ids"]) < 2:
        raise ValueError("fewer than two fitted subject pairs")
    fits = result["path_fits"]
    proportions = np.asarray([[fit.path_proportions for fit in pair] for pair in fits])
    outputs = [{"strategy": "subject-mean A1", "effect": proportions[:, 1].mean(axis=0) - proportions[:, 0].mean(axis=0), "converged": True, "n_subjects": len(fits)}]
    retained = np.isin(subjects, result["subject_ids"])
    paired_data = ec_glmm.ECGLMMData(tuple(values[retained] for values in base.counts), base.compatibility, base.design[retained], base.clusters[retained])
    for name, local, selected_labels, options in (("pooled local", base, labels, {}), ("primer-balanced pooled local", base, labels, {"balance_primers": True}), ("paired-subject pooled local", paired_data, labels[retained], {}), ("paired-subject balanced pooled local", paired_data, labels[retained], {"balance_primers": True})):
        fitted = ec_block_glmm.pooled_path_effect(local, path_index, selected_labels, baseline=baseline, **options)
        outputs.append({"strategy": name, "effect": fitted["difference"], "converged": fitted["converged"], "n_subjects": len(fits)})
    flat_fits = [fit for pair in fits for fit in pair]
    depths = [sum(float(values[(subjects == subject) & (labels == level)].sum()) for values in base.counts) for subject in result["subject_ids"] for level in result["levels"]]
    try:
        if not all(fit.covariance.identifiable for fit in flat_fits):
            raise ValueError("unidentifiable subject measurement covariance")
        pooling = dirichlet_pooling(proportions.reshape(-1, proportions.shape[-1]), np.asarray([proportion_covariance(fit) for fit in flat_fits]), np.tile([0, 1], len(fits)), depths)
        outputs.append({"strategy": "effective-count Dirichlet pooling", "effect": pooling["means"][1] - pooling["means"][0], "converged": pooling["converged"], "n_subjects": len(fits), "concentration": pooling["concentration"], "median_effective_depth": np.median(pooling["effective_depth"]), "max_subject_weight_fraction": max(pooling["subject_precision_weights"][index::2].max() / pooling["subject_precision_weights"][index::2].sum() for index in (0, 1)), "error": ""})
    except (ValueError, np.linalg.LinAlgError) as exception:
        outputs.append({"strategy": "effective-count Dirichlet pooling", "effect": np.full(proportions.shape[-1], np.nan), "converged": False, "n_subjects": len(fits), "error": str(exception)})
    return outputs


def omnibus_statistics(values, proportions, labels, subjects):
    design, tested, _, _ = ec_block_glmm.blocked_multilevel_design(labels, subjects)
    concentration = max(values)
    reference_values = values[concentration]
    stats = {f"null-variance Wald A{int(concentration)}": differential.multivariate_gls_test(reference_values, np.zeros((len(labels), reference_values.shape[1], reference_values.shape[1])), design, tested)}
    for concentration, array in values.items():
        for name, result in regression_omnibus(array, design, tested).items():
            stats[f"{name} ILR A{int(concentration)}"] = result
    array = proportions @ differential.helmert_basis(proportions.shape[1])
    for name, result in regression_omnibus(array, design, tested).items():
        stats[f"{name} proportions A1"] = result
    return stats


def omnibus_reports(base, baseline, path_index, labels, subjects, replicates, seed, test_concentration=32.):
    results = {concentration: ec_block_glmm.blocked_multilevel_path_test(base, path_index, labels, subjects, baseline=baseline, path_pseudocount=concentration, path_pseudocount_scaling="total") for concentration in (test_concentration, 1.)}
    # Quantification failures must not silently create a different design.
    first, second = results[test_concentration], results[1.]
    keys = [(str(subject), int(level)) for subject, level in zip(first["observation_subjects"], first["observation_labels"])]
    second_keys = [(str(subject), int(level)) for subject, level in zip(second["observation_subjects"], second["observation_labels"])]
    common = sorted(set(keys) & set(second_keys))
    if not common:
        raise ValueError("no common quantified subject-level observations")
    indexed = {}
    for concentration, result, local_keys in ((test_concentration, first, keys), (1., second, second_keys)):
        positions = [local_keys.index(key) for key in common]
        indexed[concentration] = result["values"][positions]
        if concentration == 1.:
            proportions = np.asarray([result["path_fits"][position].path_proportions for position in positions])
    observation_subjects = np.array([key[0] for key in common])
    observation_labels = np.array([key[1] for key in common])
    statistics = omnibus_statistics(indexed, proportions, observation_labels, observation_subjects)
    null = []
    rng = np.random.default_rng(seed)
    positions_by_subject = [np.flatnonzero(observation_subjects == subject) for subject in np.unique(observation_subjects)]
    for replicate in range(replicates):
        permuted = observation_labels.copy()
        for positions in positions_by_subject:
            permuted[positions] = rng.permutation(permuted[positions])
        permuted_stats = omnibus_statistics(indexed, proportions, permuted, observation_subjects)
        null.extend({"strategy": name, "replicate": replicate, **{key: result[key] for key in ("p_value", "statistic", "degrees_of_freedom")}} for name, result in permuted_stats.items())
    # Subject-blocked adjusted usage differences share one reference design.
    design, tested, levels, _ = ec_block_glmm.blocked_multilevel_design(observation_labels, observation_subjects)
    coefficients = np.linalg.lstsq(design, proportions, rcond=None)[0][tested]
    adjusted_effects = np.vstack([np.zeros(proportions.shape[1]), coefficients])
    return statistics, null, {"n_subjects": len(np.unique(observation_subjects)), "n_observations": len(common), "levels": levels, "adjusted_effects": adjusted_effects, "n_fit_only_testing": len(set(keys) - set(common)), "n_fit_only_a1": len(set(second_keys) - set(common)), "quantified": {"labels": observation_labels.tolist(), "subjects": observation_subjects.tolist(), "proportions": proportions.tolist(), "values": {str(key): value.tolist() for key, value in indexed.items()}}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--mode", choices=("pairwise", "omnibus"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=32)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--test-concentration", type=float, default=32.)
    parser.add_argument("--profile-mass", action="store_true")
    args = parser.parse_args()
    if args.profile_mass:
        # Audit-process-only override, never the default production fitting path.
        differential.fit_path_perturbation = differential.fit_profiled_path_perturbation
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    candidates = cached["candidates"]
    if args.reference:
        reference = pd.read_csv(args.reference, sep="\t")
        if "converged" in reference:
            reference = reference.loc[reference.converged.astype(str).str.lower().eq("true") & reference.n_subjects.ge(4)]
        requested = set(reference.test_id)
        candidates = [candidate for candidate in candidates if candidate[0] in requested]
    candidates = partition_candidates(candidates, args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    outputs, nulls, failures, baseline_cache, quantified = [], [], [], {}, []
    started = time.monotonic()
    for number, candidate in enumerate(candidates):
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "path_signatures": json.dumps(signatures), "n_paths": len(signatures)}
        try:
            local_metadata, _, labels = local_test_design(metadata, rows, tested_levels, "cell_type_pairwise" if args.mode == "pairwise" else "cell_type")
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            key = (gene, tuple(rows), tuple(transcripts))
            if key not in baseline_cache:
                baseline_cache[key] = ec_block_glmm.pooled_isoform_weights(base)
            baseline = baseline_cache[key]
            if args.mode == "pairwise":
                for report in paired_reports(base, baseline, path_index, labels, subjects):
                    report["effect"] = json.dumps(report["effect"].tolist())
                    outputs.append({**header, "level_a": tested_levels[0], "level_b": tested_levels[1], **report})
            else:
                statistics, null, details = omnibus_reports(base, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()), args.test_concentration)
                quantified.append({**header, "quantified": json.dumps(details.pop("quantified"))})
                adjusted = details.pop("adjusted_effects")
                level_names = [tested_levels[int(level)] for level in details.pop("levels")]
                for name, result in statistics.items():
                    outputs.append({**header, **details, "strategy": name, "p_value": result["p_value"], "statistic": result["statistic"], "degrees_of_freedom": result["degrees_of_freedom"], "converged": True, "levels": json.dumps(level_names), "adjusted_effects": json.dumps(adjusted.tolist())})
                nulls.extend({**header, **details, **row} for row in null)
        except (ValueError, np.linalg.LinAlgError) as exception:
            failures.append({**header, "error": str(exception)})
        if number % 25 == 0:
            print(f"{number + 1}/{len(candidates)} candidates, {time.monotonic() - started:.1f}s", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output_dir / "observed.tsv", sep="\t", index=False, na_rep="NA")
    if nulls:
        pd.DataFrame(nulls).to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False, na_rep="NA")
    if quantified:
        pd.DataFrame(quantified).to_csv(args.output_dir / "quantified.tsv.gz", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    (args.output_dir / "settings.json").write_text(json.dumps({"mode": args.mode, "candidate_settings": cached["settings"], "baseline": "current analytic EC optimizer, refitted within this subject fold", "profile_mass": args.profile_mass, "test_concentration": args.test_concentration, "null_replicates": args.null_replicates, "n_candidates": len(candidates), "n_failures": len(failures)}, default=str, indent=2) + "\n")


if __name__ == "__main__":
    main()
