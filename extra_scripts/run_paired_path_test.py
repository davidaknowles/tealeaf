#!/usr/bin/env python3
"""Fit subject-paired tests to EC-derived local splice-path compositions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time
import zlib

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    partition_candidates,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import differential, ec_block_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--path-pseudocount", type=float, default=0.0)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--max-candidates", type=int)
    return parser.parse_args()


def filtered_inputs(data_cache, settings):
    with data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(
            handle
        )
    metadata = group_metadata(groups)
    subject_folds = settings.get("subject_folds")
    if subject_folds is not None:
        folds = pd.read_csv(subject_folds, sep="\t", dtype={"subject": str})
        selected = set(
            folds.loc[
                folds.fold.eq(settings["subject_fold"]), "subject"
            ].astype(str)
        )
        retained = metadata.mouse.astype(str).isin(selected).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    represented = metadata.groupby("cell_type")["mouse"].nunique()
    levels = represented.index[
        represented >= int(settings["min_celltype_mice"])
    ]
    retained = metadata.cell_type.isin(levels).to_numpy()
    metadata = metadata.loc[retained].reset_index(drop=True)
    counts = tuple(matrix[retained] for matrix in counts)
    return metadata, counts, genes, gene_transcripts, gene_ecs, designs


def signed_null_p_value(differences, rng):
    signs = rng.choice((-1.0, 1.0), size=(len(differences), 1))
    return differential.paired_mean_test(differences * signs)["p_value"]


def main():
    args = parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid shard index")
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached["settings"]
    if settings.get("test_effect") != "cell_type_pairwise":
        raise ValueError("paired path testing requires a pairwise candidate cache")
    candidates = cached["candidates"]
    if args.max_candidates is not None:
        candidates = candidates[: int(args.max_candidates)]
    candidates = partition_candidates(candidates, args.shard_count)[
        args.shard_index
    ]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(
        args.data_cache, settings
    )
    observed_rows = []
    null_rows = []
    failures = []
    pooled_cache = {}
    started = time.perf_counter()
    for candidate in candidates:
        (
            test_id,
            block_id,
            gene_id,
            gene,
            transcripts,
            path_index,
            signatures,
            rows,
            _,
            tested_levels,
        ) = candidate
        try:
            local_metadata, _, labels = local_test_design(
                metadata, rows, tested_levels, "cell_type_pairwise"
            )
            local_counts = tuple(matrix[rows] for matrix in counts)
            clusters = local_metadata.mouse.astype(str).to_numpy()
            base, _, totals = local_gene_data(
                local_counts,
                designs,
                transcripts,
                gene_ecs[gene],
                np.ones((len(local_metadata), 1)),
                clusters,
                drop_zero=False,
            )
            cache_key = (gene, tuple(rows), tuple(transcripts))
            if cache_key not in pooled_cache:
                pooled_cache[cache_key] = ec_block_glmm.pooled_isoform_weights(
                    base
                )
            baseline = pooled_cache[cache_key]
            result = ec_block_glmm.paired_path_test(
                base,
                path_index,
                labels,
                clusters,
                baseline=baseline,
                max_iter=args.max_iter,
                path_pseudocount=args.path_pseudocount,
            )
            differences = result.pop("differences")
            result.pop("path_fits")
            result.pop("subject_ids")
            result.pop("levels")
            observed_rows.append({
                "test_id": test_id,
                "block_id": block_id,
                "gene_id": gene_id,
                "contrast": "cell_type_pairwise",
                "level_a": tested_levels[0],
                "level_b": tested_levels[1],
                "method": "paired_path",
                "path_pseudocount": args.path_pseudocount,
                "n_paths": len(signatures),
                "n_isoforms": base.n_isoforms,
                "n_ecs": len(gene_ecs[gene]),
                "n_samples": len(local_metadata),
                "n_subjects": result["n_subjects"],
                "degrees_of_freedom": result["degrees_of_freedom"],
                "median_gene_umis": float(np.median(totals)),
                "statistic": result["statistic"],
                "p_value": result["p_value"],
                "converged": result["converged"],
                "mean_difference_norm": float(
                    np.linalg.norm(differences.mean(axis=0))
                ) if len(differences) else 0.0,
            })
            if result["converged"]:
                test_hash = zlib.crc32(test_id.encode("utf-8"))
                for replicate in range(args.null_replicates):
                    rng = np.random.default_rng(
                        np.random.SeedSequence((args.seed, test_hash, replicate))
                    )
                    null_rows.append({
                        "test_id": test_id,
                        "block_id": block_id,
                        "replicate": replicate,
                        "p_value": signed_null_p_value(differences, rng),
                    })
        except Exception as error:
            failures.append({"test_id": test_id, "error": repr(error)})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(observed_rows).to_csv(
        args.output_dir / "paired_path.tsv", sep="\t", index=False
    )
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "paired_path_null.tsv.gz", sep="\t", index=False
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    (args.output_dir / "summary.json").write_text(json.dumps({
        "candidates": len(candidates),
        "completed": len(observed_rows),
        "failures": len(failures),
        "elapsed_seconds": time.perf_counter() - started,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
