#!/usr/bin/env python3
"""Test SUPPA2 event definitions with Tealeaf's local-path EC model.

Each SUPPA2 event's included and excluded transcript sets define two tested
paths. Other annotated isoforms are nuisance components. Primer-specific EC
counts are fitted with Tealeaf's paired path estimator and paired test; signed
subject-label nulls are emitted in the same format as the production path
pipeline for matched empirical calibration.
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

from extra_scripts.run_ec_block_glmm import local_test_design
from extra_scripts.run_ec_glmm import local_gene_data
from extra_scripts.run_paired_path_test import filtered_inputs, signed_null_p_value
from tealeaf.sc import ec_block_glmm


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
    return parser.parse_args()


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
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
    contexts = {}
    for candidate in cached["candidates"]:
        _, _, _, gene, _, _, _, rows, _, levels = candidate
        key = (int(gene), tuple(np.asarray(rows, dtype=int)), tuple(levels))
        contexts.setdefault(key, candidate)
    tests = []
    represented_event_ids = set()
    supported_transcript_cache = {}
    for (gene, rows, levels), candidate in contexts.items():
        # Fit the whole supported gene transcript set so a SUPPA2 event is not
        # restricted to whichever local block happened to qualify as a screen.
        if gene not in supported_transcript_cache:
            supported_transcript_cache[gene] = supported_gene_transcripts(
                gene, gene_transcripts, gene_ecs, designs
            )
        transcripts = supported_transcript_cache[gene]
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
    observed, null, failures = [], [], []
    baseline_cache = {}
    started = time.perf_counter()
    for candidate, event, path_index in tests:
        test_id, _, gene_id, gene, _, _, _, rows, _, tested_levels = candidate
        event_id = str(event.event_id)
        event_test_id = f"SUPPA2:{event_id}|cell_type|{'|'.join(tested_levels)}"
        try:
            local_metadata, _, labels = local_test_design(
                metadata, rows, tested_levels, "cell_type_pairwise"
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata.mouse.astype(str).to_numpy()
            transcripts = supported_transcript_cache[gene]
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
            result = ec_block_glmm.paired_path_test(
                base,
                path_index,
                labels,
                clusters,
                baseline=baseline_cache[cache_key],
                max_iter=args.max_iter,
                path_pseudocount=args.path_pseudocount,
                path_prior_center=args.path_prior_center,
                path_pseudocount_scaling=args.path_pseudocount_scaling,
            )
            values = result["differences"]
            covariances = result["difference_covariances"]
            observed.append({
                "test_id": event_test_id,
                "block_id": event_id,
                "gene_id": gene_id,
                "contrast": "cell_type_pairwise",
                "level_a": tested_levels[0],
                "level_b": tested_levels[1],
                "method": "Tealeaf EC; SUPPA2 event definitions",
                "path_pseudocount": args.path_pseudocount,
                "path_prior_center": args.path_prior_center,
                "path_pseudocount_scaling": args.path_pseudocount_scaling,
                "retain_uncertainty": False,
                "uncertainty_scale": 0.0,
                "n_paths": 2,
                "n_isoforms": base.n_isoforms,
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
    pd.DataFrame(observed).to_csv(args.output_dir / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(null).to_csv(args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False)
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
    print(f"wrote {len(observed):,} tests and {len(null):,} nulls; skipped {skipped_events:,} unsupported catalogue events")


if __name__ == "__main__":
    main()
