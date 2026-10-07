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
    parser.add_argument("--residual-concentration", type=float)
    parser.add_argument("--within-path-type-scale", type=float)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid null shard")
    source = args.cache / f"{args.source}_paired"
    with args.candidate_cache.open("rb") as handle:
        settings = pickle.load(handle)["settings"]
    if settings["subject_fold"] != 0 or settings["min_gene_umis"] != 25:
        raise ValueError("count-null control requires the frozen production fold0 coverage recipe")
    metadata, counts, genes, gene_tx, gene_ecs, designs = filtered_inputs(source / "prepared.pkl", settings)
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
    expected_strategies = [f"hybrid profiled ILR A{value}" for value in (32, 64)]
    observed, null, failures, contexts = [], [], [], {}
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
            rng = np.random.default_rng(381924 + zlib.crc32(test_id.encode()))
            for draw in range(args.draws):
                generated = simulate_counts(base, baseline, subjects, rng, .5, labels=labels, path_index=path_index, residual_concentration=args.residual_concentration, within_path_type_scale=args.within_path_type_scale)
                fitted_baseline = ec_block_glmm.pooled_isoform_weights(generated)
                collapsed, collapsed_paths, collapsed_baseline = collapse_event_nuisance(generated, path_index, fitted_baseline)
                for concentration, strategy in zip((32, 64), expected_strategies):
                    try:
                        result = ec_block_glmm.paired_path_test(collapsed, collapsed_paths, labels, subjects, baseline=collapsed_baseline, path_pseudocount=concentration, path_pseudocount_scaling="total", profile_event_mass=True)
                        complete = result["converged"] and result["n_subjects"] == n_expected and n_expected >= 4
                        norm = float(np.linalg.norm(result["differences"].mean(axis=0))) if complete else np.nan
                        observed.append({**header, "draw": draw, "strategy": strategy, "p_value": result["p_value"] if complete else 1., "statistic": result["statistic"] if complete else 0., "degrees_of_freedom": 1, "n_subjects": result["n_subjects"], "n_observations": len(rows), "mean_difference_norm": norm, "converged": complete})
                        if complete:
                            signs = np.random.default_rng(381924 + zlib.crc32(test_id.encode()) + 1721 * draw)
                            null.extend({**header, "draw": draw, "strategy": strategy, "replicate": replicate, "n_subjects": n_expected, "p_value": signed_null_p_value(result["differences"], result["difference_covariances"], signs, False, 0.)} for replicate in range(32))
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
    manifest = {"source": args.source, "candidate_settings": settings, "requested_ids": [row["test_id"] for row in requested], "expected_strategies": expected_strategies, "draws": args.draws, "null_replicates": 32, "residual_concentration": args.residual_concentration, "within_path_type_scale": args.within_path_type_scale, "seed": 381924, "selection": "fixed random sample of the whole screened event family, including fit failures; no significance or LR selection", "baseline": "full transcript pooled fit refitted on each generated count draw before actual event-class collapse", "null": "zero conditional-mean target-path contrast; observed primer totals, compatibility, coverage and subject missingness retained", "failure_policy": "all requested trials at p1 on exceptions or incomplete subject fits"}
    (args.output_dir / "settings.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
