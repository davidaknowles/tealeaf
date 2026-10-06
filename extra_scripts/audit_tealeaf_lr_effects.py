#!/usr/bin/env python3
"""Audit short-read-only reporting effects while preserving production test p-values."""

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
from extra_scripts.assess_tilgner_long_read_replication import normalized_difference, vector_agreement
from tealeaf.sc import ec_block_glmm, ec_glmm, differential


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--replication", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=8)
    args = parser.parse_args()
    source = pd.read_csv(args.replication, sep="\t", low_memory=False)
    source = source.loc[source.mapping_complete.astype(str).str.lower().eq("true") & source.minimum_pooled_depth.ge(20)].drop_duplicates("test_id").set_index("test_id", verify_integrity=True)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    candidates = [item for item in cached["candidates"] if item[0] in source.index]
    candidates = partition_candidates(candidates, args.shard_count)[args.shard_index]
    metadata, counts, _, _, gene_ecs, designs = filtered_inputs(args.data_cache, cached["settings"])
    rows_out, baseline_cache = [], {}
    started = time.monotonic()
    for number, candidate in enumerate(candidates):
        test_id, block_id, gene_id, gene, transcripts, path_index, signatures, rows, _, tested_levels = candidate
        local_metadata, _, labels = local_test_design(metadata, rows, tested_levels, "cell_type_pairwise")
        base, _, _ = local_gene_data(tuple(matrix[rows] for matrix in counts), designs, transcripts, gene_ecs[gene], np.ones((len(local_metadata), 1)), local_metadata.mouse.astype(str).to_numpy(), drop_zero=False)
        key = (gene, tuple(rows), tuple(transcripts))
        if key not in baseline_cache:
            baseline_cache[key] = ec_block_glmm.pooled_isoform_weights(base)
        baseline = baseline_cache[key]
        record = source.loc[test_id]
        first = [np.asarray(json.loads(record[f"counts_a_rep{i}"]), dtype=float) for i in (1, 2)]
        second = [np.asarray(json.loads(record[f"counts_b_rep{i}"]), dtype=float) for i in (1, 2)]
        external = normalized_difference(sum(first), sum(second))
        external_replicates = [normalized_difference(a, b) for a, b in zip(first, second)]
        for strategy, options in (("pooled local", {}), ("primer-balanced pooled local", {"balance_primers": True}), ("primer 0 pooled local", {"primer": 0}), ("primer 1 pooled local", {"primer": 1}), ("pooled free isoforms", None)):
            converged, error, effect = False, "", np.full(len(signatures), np.nan)
            try:
                if options is None:
                    proportions, statuses = [], []
                    for level in np.unique(labels):
                        selected = labels == level
                        subset = ec_glmm.ECGLMMData(tuple(values[selected] for values in base.counts), base.compatibility, base.design[selected], base.clusters[selected])
                        weights, status = ec_block_glmm.pooled_isoform_weights(subset, return_status=True)
                        proportions.append(differential.path_proportions(weights, path_index))
                        statuses.append(status)
                    effect, converged = proportions[1] - proportions[0], all(statuses)
                else:
                    result = ec_block_glmm.pooled_path_effect(base, path_index, labels, baseline=baseline, **options)
                    effect, converged = result["difference"], result["converged"]
                dot, cosine = vector_agreement(effect, external)
                replicate_dots = [vector_agreement(effect, delta)[0] for delta in external_replicates]
            except (ValueError, np.linalg.LinAlgError) as exception:
                error = str(exception)
                dot, cosine, replicate_dots = np.nan, np.nan, [np.nan, np.nan]
            rows_out.append({"test_id": test_id, "block_id": block_id, "gene_id": gene_id, "level_a": tested_levels[0], "level_b": tested_levels[1], "strategy": strategy, "n_paths": len(signatures), "n_isoforms": len(transcripts), "n_nuisance_isoforms": int(np.sum(np.asarray(path_index) < 0)), "converged": converged, "effect": json.dumps(effect.tolist()), "effect_norm": float(np.linalg.norm(effect)), "pooled_dot_product": dot, "pooled_cosine": cosine, "pooled_replicated": dot > 0 if np.isfinite(dot) and converged else np.nan, "both_replicates_replicated": all(value > 0 for value in replicate_dots) if np.isfinite(replicate_dots).all() and converged else np.nan, "minimum_replicate_depth": record.minimum_replicate_depth, "error": error})
        if number % 20 == 0:
            print(f"{number + 1}/{len(candidates)} tests, elapsed {time.monotonic() - started:.1f}s", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows_out).to_csv(args.output, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
