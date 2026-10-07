#!/usr/bin/env python3
"""Biological-null EC simulations preserving observed depth and primer totals."""

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
from extra_scripts.audit_path_reporting_omnibus import omnibus_reports, joint_dm_reports
from tealeaf.sc import ec_block_glmm, ec_glmm, differential
from tealeaf.sc.path_score import paired_path_score_test, paired_subject_centered_test
from tealeaf.sc.path_simulation import simulate_counts
from tealeaf.sc.path_bias import paired_null_corrected_path_test, null_corrected_path_responses
from tealeaf.sc.omnibus import regression_omnibus
from tealeaf.sc.path_score_mixed import mixed_path_score_test, mixed_score_test, MODEL_VERSION as MIXED_SCORE_VERSION


def mixed_score_statistics(base, baseline, path_index, labels, subjects, replicates, seed):
    """Fixed three reference diagnostics, sharing one fitted EC score model."""
    result = mixed_path_score_test(base, path_index, labels, subjects, baseline=baseline)
    components = result["components"]
    references = {"mixed EC score cluster F": "p_value", "mixed EC score residual F": "residual_F_p_value", "mixed EC score chi-square": "chi_square_p_value"}
    observed = {name: {**result, "p_value": result[key]} for name, key in references.items()}
    null = []
    rng = np.random.default_rng(seed)
    for replicate in range(replicates):
        flipped = components.scores * rng.choice([-1., 1.], size=(len(components.scores), 1))
        tested = mixed_score_test(flipped, components.information, components.biological_shapes)
        for name, key in references.items():
            null.append({"strategy": name, "replicate": replicate, "n_subjects": result["n_subjects"], "p_value": tested[key], "statistic": tested["statistic"], "degrees_of_freedom": tested["degrees_of_freedom"]})
    return observed, null, {"n_subjects": result["n_subjects"], "n_observations": len(labels)}


def paired_statistics(base, baseline, path_index, labels, subjects, replicates, seed, prior_center, include_score=False, null_corrected=False):
    rng = np.random.default_rng(seed)
    observed, null = {}, []
    for concentration in (32., 1.):
        if null_corrected:
            result = paired_null_corrected_path_test(base, path_index, labels, subjects, baseline=baseline, concentration=concentration)
        else:
            result = ec_block_glmm.paired_path_test(base, path_index, labels, subjects, baseline=baseline, path_pseudocount=concentration, path_pseudocount_scaling="total", path_prior_center=prior_center)
        if result["n_subjects"] < 4:
            raise ValueError("fewer than four fitted subject pairs")
        name = f"paired null-corrected proportions A{concentration:g}" if null_corrected else f"paired ILR A{concentration:g}"
        observed[name] = result
        differences = result["differences"]
        for replicate in range(replicates):
            flipped = differences * rng.choice([-1., 1.], size=(len(differences), 1))
            tested = differential.paired_mean_test(flipped)
            null.append({"strategy": name, "replicate": replicate, "n_subjects": result["n_subjects"], **{key: tested[key] for key in ("p_value", "statistic", "degrees_of_freedom")}})
    if include_score:
        configurations = [(1., 1., False), (32., 1., False), (32., 32., False), (32., 1., True)]
        for denominator, anchor, free in configurations:
            name = f"paired EC score D{denominator:g} N{anchor:g}" + (" free null" if free else "")
            result = paired_path_score_test(base, path_index, labels, subjects, baseline=baseline, denominator_concentration=denominator, null_concentration=anchor, free_null=free)
            observed[name] = result
            differences = result["differences"]
            for replicate in range(replicates):
                flipped = differences * rng.choice([-1., 1.], size=(len(differences), 1))
                tested = differential.paired_mean_test(flipped)
                null.append({"strategy": name, "replicate": replicate, "n_subjects": result["n_subjects"], **{key: tested[key] for key in ("p_value", "statistic", "degrees_of_freedom")}})
        for concentration in (1., 32.):
            name = f"paired subject-centered A{concentration:g}"
            result = paired_subject_centered_test(base, path_index, labels, subjects, baseline=baseline, concentration=concentration)
            observed[name] = result
            differences = result["differences"]
            for replicate in range(replicates):
                flipped = differences * rng.choice([-1., 1.], size=(len(differences), 1))
                tested = differential.paired_mean_test(flipped)
                null.append({"strategy": name, "replicate": replicate, "n_subjects": result["n_subjects"], **{key: tested[key] for key in ("p_value", "statistic", "degrees_of_freedom")}})
    return observed, null


def corrected_omnibus_statistics(base, baseline, path_index, labels, subjects, replicates, seed):
    """Reuse blocked-regression diagnostics after null correction, not production."""
    observed, null = {}, []
    for concentration in (32., 1.):
        quantified = null_corrected_path_responses(base, path_index, labels, subjects, baseline=baseline, concentration=concentration)
        values, levels, clusters = quantified["values"], quantified["encoded_labels"], quantified["subjects"]
        if len(np.unique(clusters)) < 4:
            raise ValueError("fewer than four quantified omnibus subject clusters")
        design, tested, _, _ = ec_block_glmm.blocked_multilevel_design(levels, clusters)
        for name, result in regression_omnibus(values, design, tested).items():
            observed[f"null-corrected {name} A{concentration:g}"] = {**result, "converged": True, "n_subjects": len(np.unique(clusters)), "n_observations": len(values)}
        rng = np.random.default_rng(seed)
        positions = [np.flatnonzero(clusters == subject) for subject in np.unique(clusters)]
        for replicate in range(replicates):
            permuted = levels.copy()
            for local in positions:
                permuted[local] = rng.permutation(permuted[local])
            null_design, null_tested, _, _ = ec_block_glmm.blocked_multilevel_design(permuted, clusters)
            for name, result in regression_omnibus(values, null_design, null_tested).items():
                null.append({"strategy": f"null-corrected {name} A{concentration:g}", "replicate": replicate, "n_subjects": len(np.unique(clusters)), **result})
    first = next(iter(observed.values()))
    return observed, null, {"n_subjects": first["n_subjects"], "n_observations": first["n_observations"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=128)
    parser.add_argument("--draws", type=int, default=16)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--subject-scale", type=float, default=.5)
    parser.add_argument("--residual-concentration", type=float, help="Independent subject/type path Dirichlet variation, with the same conditional mean across types.")
    parser.add_argument("--within-path-type-scale", type=float, help="Label-specific within-path transcript shifts preserving the target path null; paired ILR sensitivity only.")
    parser.add_argument("--free-isoforms", action="store_true")
    parser.add_argument("--null-corrected", action="store_true", help="Experimental subtraction of smoothed expected-null path proportions, with type-specific transcript nuisance.")
    parser.add_argument("--mixed-score-only", action="store_true", help="Experimental likelihood-score REML with count uncertainty and type-specific transcript nuisance.")
    parser.add_argument("--mode", choices=("omnibus", "pairwise"), default="omnibus")
    parser.add_argument("--prior-center", choices=("uniform", "baseline"), default="uniform")
    parser.add_argument("--include-score", action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--joint-dm-only", action="store_true")
    group.add_argument("--cox-reid-only", action="store_true")
    parser.add_argument("--known-concentration", type=float, help="Known latent biological concentration control, only for Cox-Reid count-null diagnostics.")
    args = parser.parse_args()
    if args.mixed_score_only and (args.null_corrected or args.include_score or args.joint_dm_only or args.cox_reid_only or args.free_isoforms or args.prior_center != "uniform"):
        parser.error("mixed-score audit cannot be combined with alternative fit overrides")
    if args.null_corrected and (args.include_score or args.joint_dm_only or args.cox_reid_only or args.free_isoforms or args.prior_center != "uniform"):
        parser.error("null correction requires uniform-prior sensitivity, without fitter overrides")
    expected_strategies = [f"paired null-corrected proportions A{concentration}" if args.null_corrected else f"paired ILR A{concentration}" for concentration in (32, 1)]
    if args.null_corrected and args.mode == "omnibus":
        expected_strategies = [f"null-corrected {test} A{concentration}" for concentration in (32, 1) for test in ("trace F", "Pillai", "maximum-coordinate F")]
    if args.mixed_score_only:
        expected_strategies = ["mixed EC score cluster F", "mixed EC score residual F", "mixed EC score chi-square"]
    preserve_failures = args.within_path_type_scale is not None or args.null_corrected or args.mixed_score_only
    if args.within_path_type_scale is not None and ((args.mode != "pairwise" and not args.null_corrected and not args.mixed_score_only) or args.include_score or args.joint_dm_only or args.cox_reid_only):
        parser.error("within-path nuisance stress currently requires the paired ILR-only sensitivity")
    if args.within_path_type_scale is not None and (not np.isfinite(args.within_path_type_scale) or args.within_path_type_scale < 0):
        parser.error("within-path type scale must be finite and nonnegative")
    if args.known_concentration is not None and (not args.cox_reid_only or args.residual_concentration != args.known_concentration):
        parser.error("known latent precision must match the simulated finite concentration and requires --cox-reid-only")
    if args.free_isoforms:
        differential.fit_path_perturbation = differential.fit_free_isoform_paths
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    reference = pd.read_csv(args.reference, sep="\t", low_memory=False)
    eligible = set(reference.loc[reference.converged & reference.n_subjects.ge(4), "test_id"])
    candidates = [candidate for candidate in cached["candidates"] if candidate[0] in eligible]
    # Fixed random subset of the whole tested universe, not the discoveries.
    rng = np.random.default_rng(381924)
    selected = rng.choice(len(candidates), size=min(args.blocks, len(candidates)), replace=False)
    requested_ids = [candidates[index][0] for index in sorted(selected)]
    candidates = partition_candidates([candidates[index] for index in sorted(selected)], args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    outputs, nulls, failures, baseline_cache = [], [], [], {}
    for candidate in candidates:
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        header = {"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "n_paths": len(signatures)}
        if args.within_path_type_scale is not None:
            multiplicities = np.bincount(np.asarray(path_index)[np.asarray(path_index) >= 0], minlength=len(signatures))
            header.update(n_transcripts=len(transcripts), n_multiplet_paths=int((multiplicities > 1).sum()), n_outside_transcripts=int((np.asarray(path_index) < 0).sum()))
        try:
            local_metadata, _, labels = local_test_design(metadata, rows, tested_levels, "cell_type" if args.mode == "omnibus" else "cell_type_pairwise")
            subjects = local_metadata.mouse.astype(str).to_numpy()
            base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), subjects, drop_zero=False)
            key = (gene, tuple(rows), tuple(transcripts))
            if key not in baseline_cache:
                baseline_cache[key] = ec_block_glmm.pooled_isoform_weights(base)
            rng = np.random.default_rng(zlib.crc32(test_id.encode()) + 381924)
            for draw in range(args.draws):
                simulated = simulate_counts(base, baseline_cache[key], subjects, rng, args.subject_scale, labels=labels, path_index=path_index, residual_concentration=args.residual_concentration, within_path_type_scale=args.within_path_type_scale)
                try:
                    # A failed fit must not prevent later, independent count draws.
                    baseline = ec_block_glmm.pooled_isoform_weights(simulated)
                    if args.mixed_score_only:
                        statistics, null, details = mixed_score_statistics(simulated, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()) + draw * 1721)
                    elif args.null_corrected and args.mode == "omnibus":
                        statistics, null, details = corrected_omnibus_statistics(simulated, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()) + draw * 1721)
                    elif args.joint_dm_only or args.cox_reid_only:
                        statistics, null, details = joint_dm_reports(simulated, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()) + draw * 1721, dispersion_method="cox_reid" if args.cox_reid_only else "ml", known_concentration=args.known_concentration)
                    elif args.mode == "omnibus":
                        statistics, null, details = omnibus_reports(simulated, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()) + draw * 1721, prior_center=args.prior_center, include_score=args.include_score)
                    else:
                        statistics, null = paired_statistics(simulated, baseline, path_index, labels, subjects, args.null_replicates, zlib.crc32(test_id.encode()) + draw * 1721, args.prior_center, args.include_score, args.null_corrected)
                        details = {"n_subjects": statistics[expected_strategies[0]]["n_subjects"], "n_observations": 2 * statistics[expected_strategies[0]]["n_subjects"]}
                except (ValueError, np.linalg.LinAlgError) as exception:
                    failures.append({**header, "draw": draw, "error": str(exception)})
                    if preserve_failures:
                        outputs.extend({**header, "draw": draw, "strategy": strategy, "p_value": 1., "statistic": 0., "degrees_of_freedom": len(signatures) - 1, "n_subjects": 0, "n_observations": len(labels), "converged": False, "mean_difference_norm": np.nan} for strategy in expected_strategies)
                    continue
                for name, result in statistics.items():
                    mean_norm = np.linalg.norm(result["differences"].mean(axis=0)) if "differences" in result and len(result["differences"]) else np.nan
                    outputs.append({**header, "draw": draw, "strategy": name, "p_value": result["p_value"], "statistic": result["statistic"], "degrees_of_freedom": result["degrees_of_freedom"], "n_subjects": result.get("n_subjects", details["n_subjects"]), "n_observations": details["n_observations"], "converged": result.get("converged", True), "mean_difference_norm": mean_norm})
                    if args.cox_reid_only:
                        outputs[-1].update({key: result.get(key, np.nan) for key in ("alternative_concentration", "profile_index", "profile_boundary", "residual_degrees_of_freedom", "effective_depth_median")})
                        outputs[-1]["report_error"] = details["report_error"]
                nulls.extend({**header, "draw": draw, "n_subjects": details["n_subjects"], **row} for row in null)
        except (ValueError, np.linalg.LinAlgError) as exception:
            failures.append({**header, "error": str(exception)})
            if preserve_failures:
                recorded = {(row["test_id"], row["draw"], row["strategy"]) for row in outputs}
                outputs.extend({**header, "draw": draw, "strategy": strategy, "p_value": 1., "statistic": 0., "degrees_of_freedom": len(signatures) - 1, "n_subjects": 0, "n_observations": 0, "converged": False, "mean_difference_norm": np.nan} for draw in range(args.draws) for strategy in expected_strategies if (test_id, draw, strategy) not in recorded)
        print(f"{test_id}, cumulative observed={len(outputs)}, failed blocks={len(failures)}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output_dir / "observed.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(nulls).to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    null_description = "common transcript mixture for all cell types within each subject" if args.residual_concentration is None else "independent subject/type Dirichlet path compositions with identical conditional mean across types"
    settings = {"candidate_settings": cached["settings"], "selected_blocks": args.blocks, "draws": args.draws, "null_replicates": args.null_replicates, "subject_scale": args.subject_scale, "residual_concentration": args.residual_concentration, "free_isoforms": args.free_isoforms, "mode": args.mode, "prior_center": args.prior_center, "include_score": args.include_score, "joint_dm_only": args.joint_dm_only, "cox_reid_only": args.cox_reid_only, "known_latent_concentration": args.known_concentration, "null": null_description + ", observed primer/depth totals preserved", "baseline": "refitted for each simulated count draw", "seed": 381924}
    if preserve_failures:
        settings.update(requested_ids=requested_ids, expected_strategies=expected_strategies, null_corrected=args.null_corrected, mixed_score_only=args.mixed_score_only, failure_policy="all requested tests retained with p=1 and unavailable effects")
    if args.mixed_score_only:
        settings.update(mixed_score_version=MIXED_SCORE_VERSION, primary_reference="Modified-Knapp-Hartung-style cluster F; chi-square and residual F diagnostic only", variance_shape="Dirichlet ILR contrast shape, estimated shared REML multiplier per hypothesis")
    if args.within_path_type_scale is not None:
        settings.update(within_path_type_scale=args.within_path_type_scale, null="same target path usage within subject; within-path transcript composition changes across cell types; observed primer/depth totals preserved")
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
