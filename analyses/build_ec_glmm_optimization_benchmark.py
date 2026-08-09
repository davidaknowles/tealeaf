#!/usr/bin/env python3
"""Build matched EC-GLMM optimization subsets and batching diagnostics."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    treatment_design,
)
from tealeaf.sc import ec_block_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--groups", type=int, default=5)
    parser.add_argument("--blocks-per-group", type=int, default=4)
    return parser.parse_args()


def filtered_metadata(groups, settings):
    metadata = group_metadata(groups)
    test_effect = settings.get("test_effect", "cell_type")
    if str(test_effect).startswith("cell_type"):
        minimum = int(settings.get("min_celltype_mice", 3))
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        retained_levels = represented.index[represented >= minimum]
        metadata = metadata.loc[metadata.cell_type.isin(retained_levels)].reset_index(
            drop=True
        )
    return metadata


def candidate_dimensions(candidate, metadata, gene_ecs, designs, settings):
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
    transcripts = np.asarray(transcripts, dtype=int)
    rows = np.asarray(rows, dtype=int)
    local, nuisance, labels = local_test_design(
        metadata,
        rows,
        tested_levels,
        settings.get("test_effect", "cell_type"),
    )
    supported_ecs = []
    for mapping in designs:
        local_mapping = mapping[gene_ecs[gene]][:, transcripts]
        supported_ecs.append(
            int(np.sum(np.asarray(local_mapping.sum(axis=1)).ravel() > 0))
        )
    tested = treatment_design(labels, len(tested_levels))
    null_tensor, _, tested_count = ec_block_glmm.block_fixed_effect_tensors(
        nuisance, tested, path_index
    )
    alternative_tensor = ec_block_glmm.full_fixed_effect_tensor(
        nuisance, tested, len(transcripts) - 1
    )
    null_coefficients = null_tensor.shape[2]
    alternative_coefficients = alternative_tensor.shape[2]
    clusters = pd.factorize(local.mouse, sort=True)[0]
    cluster_sizes = tuple(sorted(Counter(clusters).values()))
    batch_key = (
        len(rows),
        len(transcripts),
        *supported_ecs,
        len(cluster_sizes),
        cluster_sizes,
        null_coefficients,
    )
    cost = (
        len(rows)
        * len(transcripts)
        * max(sum(supported_ecs), 1)
        * max(null_coefficients, 1)
    )
    return {
        "test_id": test_id,
        "block_id": block_id,
        "gene_id": gene_id,
        "gene_index": gene,
        "samples": len(rows),
        "clusters": len(cluster_sizes),
        "cluster_sizes": json.dumps(cluster_sizes),
        "isoforms": len(transcripts),
        "primer_1_ecs": supported_ecs[0],
        "primer_2_ecs": supported_ecs[1],
        "paths": len(signatures),
        "tested_coefficients": tested_count,
        "null_coefficients": null_coefficients,
        "alternative_coefficients": alternative_coefficients,
        "batch_key": repr(batch_key),
        "ragged_dot_key": repr((len(transcripts), *supported_ecs, null_coefficients)),
        "ragged_optimizer_key": repr((len(transcripts), null_coefficients)),
        "cost_proxy": int(cost),
    }


def quantile_group_selection(group_rows, count):
    ordered = sorted(
        group_rows,
        key=lambda item: np.median([row["cost_proxy"] for row in item[1]]),
    )
    if len(ordered) <= count:
        return ordered
    positions = np.unique(np.rint(np.linspace(0, len(ordered) - 1, count)).astype(int))
    return [ordered[position] for position in positions]


def padding_summary(table):
    """Estimate dense-padding work after cost-sorted shape bucketing."""
    rows = []
    for batch_size in (8, 16, 32, 64):
        actual = 0
        padded = 0
        batches = 0
        for _, bucket in table.groupby("ragged_optimizer_key"):
            bucket = bucket.sort_values("cost_proxy")
            for start in range(0, len(bucket), batch_size):
                batch = bucket.iloc[start : start + batch_size]
                local_work = (
                    batch.samples
                    * (batch.primer_1_ecs + batch.primer_2_ecs)
                    * batch.isoforms
                )
                actual += int(local_work.sum())
                padded += int(
                    len(batch)
                    * batch.samples.max()
                    * (batch.primer_1_ecs.max() + batch.primer_2_ecs.max())
                    * batch.isoforms.iloc[0]
                )
                batches += 1
        rows.append(
            {
                "batch_size": batch_size,
                "batches": batches,
                "actual_likelihood_work_proxy": actual,
                "padded_likelihood_work_proxy": padded,
                "padded_to_actual_ratio": padded / actual,
            }
        )
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.data_cache.open("rb") as handle:
        groups, _, _, _, gene_ecs, designs = pickle.load(handle)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    settings = cached.get("settings", {})
    metadata = filtered_metadata(groups, settings)
    dimensions = [
        candidate_dimensions(candidate, metadata, gene_ecs, designs, settings)
        for candidate in cached["candidates"]
    ]
    table = pd.DataFrame(dimensions)
    bucket_sizes = table.groupby("batch_key").size()
    table["exact_batch_size"] = table.batch_key.map(bucket_sizes)
    ragged_sizes = table.groupby("ragged_optimizer_key").size()
    table["ragged_optimizer_batch_size"] = table.ragged_optimizer_key.map(ragged_sizes)
    ragged_dot_sizes = table.groupby("ragged_dot_key").size()
    table["ragged_dot_batch_size"] = table.ragged_dot_key.map(ragged_dot_sizes)
    table.to_csv(
        args.output_dir / "candidate_batchability.tsv.gz", sep="\t", index=False
    )
    padding_summary(table).to_csv(
        args.output_dir / "padding_efficiency.tsv", sep="\t", index=False
    )
    grouped_candidates = defaultdict(list)
    grouped_rows = defaultdict(list)
    for candidate, row in zip(cached["candidates"], dimensions):
        key = (candidate[3], tuple(np.asarray(candidate[7], dtype=int)))
        grouped_candidates[key].append(candidate)
        grouped_rows[key].append(row)
    reusable = [(key, rows) for key, rows in grouped_rows.items() if len(rows) >= 2]
    selected_groups = quantile_group_selection(reusable, int(args.groups))
    selected_candidates = []
    selected_rows = []
    for key, rows in selected_groups:
        order = np.argsort([row["cost_proxy"] for row in rows])
        if len(order) > int(args.blocks_per_group):
            positions = np.unique(
                np.rint(
                    np.linspace(0, len(order) - 1, int(args.blocks_per_group))
                ).astype(int)
            )
            order = order[positions]
        candidates = grouped_candidates[key]
        selected_candidates.extend(candidates[index] for index in order)
        selected_rows.extend(rows[index] for index in order)
    with (args.output_dir / "candidates.pkl").open("wb") as handle:
        pickle.dump(
            {"settings": settings, "candidates": selected_candidates},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    pd.DataFrame(selected_rows).to_csv(
        args.output_dir / "benchmark_candidates.tsv", sep="\t", index=False
    )
    summary = {
        "candidates": len(table),
        "gene_row_groups": len(grouped_rows),
        "groups_with_multiple_blocks": len(reusable),
        "candidates_in_exact_batches_of_at_least_2": int(
            np.sum(table.exact_batch_size >= 2)
        ),
        "candidates_in_exact_batches_of_at_least_4": int(
            np.sum(table.exact_batch_size >= 4)
        ),
        "maximum_exact_batch_size": int(table.exact_batch_size.max()),
        "candidates_in_ragged_dot_batches_of_at_least_2": int(
            np.sum(table.ragged_dot_batch_size >= 2)
        ),
        "candidates_in_ragged_dot_batches_of_at_least_8": int(
            np.sum(table.ragged_dot_batch_size >= 8)
        ),
        "maximum_ragged_dot_batch_size": int(table.ragged_dot_batch_size.max()),
        "candidates_in_ragged_optimizer_batches_of_at_least_2": int(
            np.sum(table.ragged_optimizer_batch_size >= 2)
        ),
        "candidates_in_ragged_optimizer_batches_of_at_least_16": int(
            np.sum(table.ragged_optimizer_batch_size >= 16)
        ),
        "maximum_ragged_optimizer_batch_size": int(
            table.ragged_optimizer_batch_size.max()
        ),
        "benchmark_groups": len(selected_groups),
        "benchmark_candidates": len(selected_candidates),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
