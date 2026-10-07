"""Actual-count null checks for the complete binary-input event controls."""

import argparse
import json
from pathlib import Path
import pickle
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import covered_celltype_pairwise_designs, local_test_design, modeled_gene_umis
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_paired_path_test import filtered_inputs, signed_null_p_value
from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, collapse_event_nuisance, event_path_index, supported_gene_transcripts
from tealeaf.sc import ec_block_glmm
from tealeaf.sc.path_simulation import simulate_counts
from tealeaf.sc.path_score_mixed import shared_path_score_components, aggregate_path_scores, signed_path_score_p_value, MODEL_VERSION as MIXED_SCORE_VERSION
from tealeaf.sc.replication_audit import complete_cluster_fit
from tealeaf.sc.conditional_path_score import binary_fragment_opportunity_kernels
from tealeaf.sc.ec_glmm import ECGLMMData


def trial_header(record):
    """Keep event identity in the shared count-null assessment schema."""
    return {"test_id": record["test_id"], "block_id": record["feature_id"].removeprefix("SUPPA2:"), "gene_id": record["gene_id"], "n_paths": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--events", type=int, default=32)
    parser.add_argument("--draws", type=int, default=2)
    parser.add_argument("--subject-scale", type=float, default=.5)
    parser.add_argument("--concentrations", type=float, nargs="+", default=[32, 64])
    parser.add_argument("--inference", choices=("profiled", "mixed-score"), default="profiled")
    parser.add_argument("--max-iter", type=int, default=300, help="Shared-null iteration limit for the mixed-score diagnostic.")
    parser.add_argument("--null-multistart", action="store_true", help="Compare pooled and interior initializations of the same mixed-score null.")
    parser.add_argument("--score-coordinate", choices=("ilr", "proportion"), default="ilr")
    parser.add_argument("--information-metric", choices=("absolute", "reference"), default="absolute")
    parser.add_argument("--residual-concentration", type=float)
    parser.add_argument("--within-path-type-scale", type=float)
    parser.add_argument("--count-likelihood", choices=("multinomial", "conditional"), default="multinomial")
    parser.add_argument("--ec-opportunity-scale", type=float, default=0.)
    parser.add_argument("--kernel-units", choices=("prepared", "fragment"), default="prepared")
    parser.add_argument("--simulation-kernel-units", choices=("analysis", "prepared", "fragment"), default="analysis")
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid null shard")
    if args.null_multistart and args.inference != "mixed-score":
        raise ValueError("null multistart is only a mixed-score fitting diagnostic")
    if args.score_coordinate != "ilr" and args.inference != "mixed-score":
        raise ValueError("proportion coordinates require mixed-score inference")
    if args.information_metric != "absolute" and args.inference != "mixed-score":
        raise ValueError("reference rank metric requires mixed-score inference")
    if args.count_likelihood != "multinomial" and args.inference != "mixed-score":
        raise ValueError("conditional count likelihood requires mixed-score inference")
    if not np.isfinite(args.subject_scale) or args.subject_scale < 0 or not args.concentrations or len(set(args.concentrations)) != len(args.concentrations) or any(not np.isfinite(value) or value <= 0 for value in args.concentrations):
        raise ValueError("nonnegative subject scale and unique positive concentrations required")
    source = args.cache / f"{args.source}_paired"
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] != 0 or settings["min_gene_umis"] != 25:
        raise ValueError("count-null control requires the frozen production fold0 coverage recipe")
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(source / "prepared.pkl", settings)
    original_designs = designs
    fragment_designs = binary_fragment_opportunity_kernels(designs) if "fragment" in (args.kernel_units, args.simulation_kernel_units) else None
    if args.kernel_units == "fragment":
        designs = fragment_designs
    simulation_designs = designs if args.simulation_kernel_units == "analysis" else (fragment_designs if args.simulation_kernel_units == "fragment" else original_designs)
    if metadata.duplicated(["mouse", "cell_type"]).any():
        raise ValueError("expected one pseudobulk per subject/type")
    features = (source / "features.txt").read_text().splitlines()
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    family = []
    for index in range(32):
        shard = args.cache / f"split/{args.source}/fold0/shard_{index}"
        summary = json.loads((shard / "summary.json").read_text())
        table = pd.read_csv(shard / "paired_path.tsv", sep="\t")
        failures = json.loads((shard / "failures.json").read_text())
        if len(table) + len(failures) != summary["tests_in_shard"]:
            raise ValueError("count-null sampling requires the whole screened family")
        family.extend(table[["test_id", "feature_id", "gene_id", "level_a", "level_b"]].to_dict("records"))
        for row in failures:
            event, _, first, second = row["test_id"].rsplit("|", 3)
            family.append({"test_id": row["test_id"], "feature_id": event, "gene_id": event.removeprefix("SUPPA2:").split(";", 1)[0], "level_a": first, "level_b": second})
    family = sorted(family, key=lambda row: row["test_id"])
    if len({row["test_id"] for row in family}) != len(family):
        raise ValueError("duplicate event hypotheses in the full sampling family")
    selection = np.random.default_rng(381924).choice(len(family), min(args.events, len(family)), replace=False)
    requested = [family[index] for index in sorted(selection)]
    gene_lookup = {canonical(value): index for index, value in enumerate(genes)}
    screening_counts = tuple(value.tocsc() for value in counts)
    expected_strategies = ["hybrid free-transcript EC mixed score" + (", proportion contrast" if args.score_coordinate == "proportion" else "")] if args.inference == "mixed-score" else [f"hybrid profiled ILR A{value:g}" for value in args.concentrations]
    if args.information_metric == "reference":
        expected_strategies = [value + ", target-normalized rank" for value in expected_strategies]
    if args.count_likelihood == "conditional":
        expected_strategies = [value + ", conditional EC opportunities" for value in expected_strategies]
    if args.kernel_units == "fragment":
        expected_strategies = [value + ", fragment kernel" for value in expected_strategies]
    if args.simulation_kernel_units != "analysis":
        expected_strategies = [value + f", {args.simulation_kernel_units}-kernel null" for value in expected_strategies]
    observed, null, failures, contexts, diagnostics = [], [], [], {}, []
    for record in requested[args.shard_index::args.shard_count]:
        test_id = record["test_id"]
        levels = (record["level_a"], record["level_b"])
        header = trial_header(record)
        try:
            gene = gene_lookup[canonical(record["gene_id"])]
            if gene not in contexts:
                transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
                umis = modeled_gene_umis(screening_counts, designs, gene_ecs[gene], transcripts)
                specs = covered_celltype_pairwise_designs(metadata, umis, min_gene_umis=settings["min_gene_umis"], min_samples=settings["min_gene_samples"], min_celltype_mice=settings["min_celltype_mice"])
                contexts[gene] = transcripts, {tuple(spec[0][-1]): spec[0][0] for spec in specs}
            transcripts, rows_by_level = contexts[gene]
            rows = rows_by_level[levels]
            local, _, labels = local_test_design(metadata, rows, levels, "cell_type_pairwise")
            subjects = local.mouse.astype(str).to_numpy()
            n_expected = len(set(subjects[labels == 0]) & set(subjects[labels == 1]))
            event = catalog.loc[record["feature_id"]]
            path_index = event_path_index(transcripts, features, event.included, event.excluded)
            if path_index is None:
                raise ValueError("sampled event no longer has complete transcript support")
            base, _, _ = local_gene_data(tuple(value[rows] for value in counts), designs, transcripts, gene_ecs[gene], np.ones((len(rows), 1)), subjects, drop_zero=False)
            baseline = ec_block_glmm.pooled_isoform_weights(base)
            if args.simulation_kernel_units == "analysis":
                simulation_base, simulation_baseline = base, baseline
            else:
                simulation_base, _, _ = local_gene_data(tuple(value[rows] for value in counts), simulation_designs, transcripts, gene_ecs[gene], np.ones((len(rows), 1)), subjects, drop_zero=False)
                simulation_baseline = ec_block_glmm.pooled_isoform_weights(simulation_base)
            rng = np.random.default_rng(381924 + zlib.crc32(test_id.encode()))
            for draw in range(args.draws):
                generated = simulate_counts(simulation_base, simulation_baseline, subjects, rng, args.subject_scale, labels=labels, path_index=path_index, residual_concentration=args.residual_concentration, within_path_type_scale=args.within_path_type_scale, ec_opportunity_scale=args.ec_opportunity_scale)
                generated = ECGLMMData(generated.counts, base.compatibility, base.design, base.clusters)
                fitted_baseline = ec_block_glmm.pooled_isoform_weights(generated)
                if args.inference == "profiled":
                    collapsed, collapsed_paths, collapsed_baseline = collapse_event_nuisance(generated, path_index, fitted_baseline)
                for concentration, strategy in zip(([0.] if args.inference == "mixed-score" else args.concentrations), expected_strategies):
                    try:
                        if args.inference == "mixed-score":
                            components = shared_path_score_components(generated, path_index, labels, subjects, baseline=fitted_baseline, max_iter=args.max_iter, null_multistart=args.null_multistart, score_coordinate=args.score_coordinate, count_likelihood=args.count_likelihood)
                            for index, subject in enumerate(components.subject_ids):
                                fit = components.null_fits[index]
                                diagnostics.append({**header, "draw": draw, "subject": subject, "null_inclusion": fit.path_proportions[0], "null_event_mass_min": fit.theta[:, path_index >= 0].sum(axis=1).min(), "null_event_mass_max": fit.theta[:, path_index >= 0].sum(axis=1).max(), "information": components.information[index, 0, 0], "score": components.scores[index, 0], "starts": fit.starts, "selected_start": fit.selected_start, "pooled_refit_inclusion": fitted_baseline[path_index == 0].sum() / fitted_baseline[path_index >= 0].sum(), "pooled_refit_event_mass": fitted_baseline[path_index >= 0].sum()})
                            result = aggregate_path_scores(components, information_metric=args.information_metric)
                        else:
                            result = ec_block_glmm.paired_path_test(collapsed, collapsed_paths, labels, subjects, baseline=collapsed_baseline, path_pseudocount=concentration, path_pseudocount_scaling="total", profile_event_mass=True)
                        fitted_subjects = result.get("n_fitted_subjects", result["n_subjects"])
                        complete = complete_cluster_fit(result, n_expected)
                        norm = float(np.linalg.norm(result["mean_difference"] if args.inference == "mixed-score" else result["differences"].mean(axis=0))) if complete else np.nan
                        trial_null = []
                        if complete:
                            signs = np.random.default_rng(381924 + zlib.crc32(test_id.encode()) + 1721 * draw)
                            for replicate in range(32):
                                if args.inference == "mixed-score":
                                    components = result["components"]
                                    null_p = signed_path_score_p_value(components, signs, information_metric=args.information_metric)
                                else:
                                    null_p = signed_null_p_value(result["differences"], result["difference_covariances"], signs, False, 0.)
                                trial_null.append({**header, "draw": draw, "strategy": strategy, "replicate": replicate, "n_subjects": n_expected, "p_value": null_p})
                        observed.append({**header, "draw": draw, "strategy": strategy, "p_value": result["p_value"] if complete else 1., "statistic": result["statistic"] if complete else 0., "degrees_of_freedom": 1, "n_subjects": result["n_subjects"], "n_fitted_subjects": fitted_subjects, "n_expected_subjects": n_expected, "n_observations": len(rows), "mean_difference_norm": norm, "converged": complete})
                        null.extend(trial_null)
                    except (ValueError, np.linalg.LinAlgError) as error:
                        failures.append({**header, "draw": draw, "strategy": strategy, "error": repr(error)})
        except (ValueError, KeyError, np.linalg.LinAlgError) as error:
            failures.append({**header, "error": repr(error)})
        recorded = {(row["draw"], row["strategy"]) for row in observed if row["test_id"] == test_id}
        observed.extend({**header, "draw": draw, "strategy": strategy, "p_value": 1., "statistic": 0., "degrees_of_freedom": 1, "n_subjects": 0, "n_observations": 0, "mean_difference_norm": np.nan, "converged": False} for draw in range(args.draws) for strategy in expected_strategies if (draw, strategy) not in recorded)
        print(f"{test_id}, completed trials={len(observed)}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(observed).to_csv(args.output_dir / "observed.tsv", sep="\t", index=False)
    pd.DataFrame(null, columns=["test_id", "block_id", "gene_id", "n_paths", "draw", "strategy", "replicate", "n_subjects", "p_value"]).to_csv(args.output_dir / "null.tsv.gz", sep="\t", index=False)
    (args.output_dir / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    if args.inference == "mixed-score":
        pd.DataFrame(diagnostics).to_csv(args.output_dir / "subject_null_diagnostics.tsv.gz", sep="\t", index=False)
    manifest = {"source": args.source, "candidate_settings": settings, "requested_ids": [row["test_id"] for row in requested], "expected_strategies": expected_strategies, "draws": args.draws, "null_replicates": 32, "subject_scale": args.subject_scale, "concentrations": [0.] if args.inference == "mixed-score" else args.concentrations, "inference": args.inference, "mixed_score_version": MIXED_SCORE_VERSION if args.inference == "mixed-score" else None, "residual_concentration": args.residual_concentration, "within_path_type_scale": args.within_path_type_scale, "seed": 381924, "selection": "fixed random sample of the whole screened event family, including fit failures; no significance or LR selection", "baseline": "full transcript pooled fit refitted on each generated count draw; profiled backend collapses event classes, mixed-score retains every supported transcript and frees type-specific nuisance", "null": "zero conditional-mean target-path contrast; observed primer totals, compatibility, coverage and subject missingness retained", "failure_policy": "all requested trials at p1 on exceptions or incomplete subject fits"}
    manifest["max_iter"] = args.max_iter if args.inference == "mixed-score" else None
    manifest["null_multistart"] = args.null_multistart
    manifest["score_coordinate"] = args.score_coordinate if args.inference == "mixed-score" else None
    manifest["information_metric"] = args.information_metric if args.inference == "mixed-score" else None
    manifest["completeness"] = "every eligible subject null must fit; mixed-score degrees of freedom count informative clusters separately from fitted clusters"
    manifest.update(count_likelihood=args.count_likelihood, ec_opportunity_scale=args.ec_opportunity_scale, kernel_units=args.kernel_units, simulation_kernel_units=args.simulation_kernel_units)
    (args.output_dir / "settings.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
