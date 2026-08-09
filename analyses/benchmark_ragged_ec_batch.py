#!/usr/bin/env python3
"""Benchmark padded and ragged-dot batches for the EC likelihood kernel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time

import numpy as np
import pandas as pd

from extra_scripts.run_ec_block_glmm import (
    group_metadata,
    local_test_design,
    tensor_data,
)
from extra_scripts.run_ec_glmm import local_gene_data
from tealeaf.sc import ec_block_glmm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-cache", required=True, type=Path)
    parser.add_argument("--candidate-cache", required=True, type=Path)
    parser.add_argument("--batchability", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument(
        "--selection", choices=("largest", "max_padding"), default="largest"
    )
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args()


def filter_metadata_and_counts(groups, counts, settings):
    metadata = group_metadata(groups)
    if str(settings.get("test_effect", "cell_type")).startswith("cell_type"):
        represented = metadata.groupby("cell_type")["mouse"].nunique()
        levels = represented.index[
            represented >= int(settings.get("min_celltype_mice", 3))
        ]
        retained = metadata.cell_type.isin(levels).to_numpy()
        metadata = metadata.loc[retained].reset_index(drop=True)
        counts = tuple(matrix[retained] for matrix in counts)
    return metadata, counts


def prepare_candidates(args):
    with args.data_cache.open("rb") as handle:
        groups, counts, _, _, gene_ecs, designs = pickle.load(handle)
    with args.candidate_cache.open("rb") as handle:
        cached = pickle.load(handle)
    metadata, counts = filter_metadata_and_counts(
        groups, counts, cached.get("settings", {})
    )
    table = pd.read_csv(args.batchability, sep="\t")
    if args.selection == "largest":
        selected_key = table.groupby("ragged_dot_key").size().idxmax()
        eligible = table.loc[table.ragged_dot_key.eq(selected_key)].sort_values(
            "cost_proxy"
        )
        if len(eligible) > int(args.batch_size):
            middle = len(eligible) // 2
            half = int(args.batch_size) // 2
            eligible = eligible.iloc[
                middle - half : middle - half + int(args.batch_size)
            ]
    else:
        choices = []
        for key, bucket in table.groupby("ragged_dot_key"):
            bucket = bucket.sort_values("cost_proxy")
            for start in range(0, len(bucket), int(args.batch_size)):
                batch = bucket.iloc[start : start + int(args.batch_size)]
                if len(batch) < 2:
                    continue
                padding = len(batch) * batch.samples.max() / batch.samples.sum()
                choices.append((padding, key, batch))
        _, selected_key, eligible = max(choices, key=lambda value: value[0])
    by_test = {candidate[0]: candidate for candidate in cached["candidates"]}
    prepared = []
    for test_id in eligible.test_id:
        candidate = by_test[test_id]
        _, _, _, gene, transcripts, path_index, _, rows, _, tested_levels = candidate
        local_metadata, nuisance, labels = local_test_design(
            metadata,
            rows,
            tested_levels,
            cached["settings"].get("test_effect", "cell_type"),
        )
        local_counts = tuple(matrix[rows] for matrix in counts)
        clusters = local_metadata.mouse.to_numpy()
        base, _, _ = local_gene_data(
            local_counts,
            designs,
            transcripts,
            gene_ecs[gene],
            np.ones((len(local_metadata), 1)),
            clusters,
            drop_zero=False,
        )
        tested = np.eye(len(tested_levels), dtype=float)[labels, 1:]
        null_tensor, _, _ = ec_block_glmm.block_fixed_effect_tensors(
            nuisance, tested, path_index
        )
        prepared.append(tensor_data(base, null_tensor))
    return selected_key, prepared


def benchmark(prepared, repetitions, seed):
    import jax
    import jax.numpy as jnp
    import jax.scipy as jsp

    jax.config.update("jax_enable_x64", True)
    batch_size = len(prepared)
    isoforms = prepared[0].n_isoforms
    dimension = isoforms - 1
    coefficient_count = prepared[0].fixed_effect_tensor.shape[2]
    sample_sizes = np.asarray([len(data.clusters) for data in prepared], dtype=np.int32)
    maximum_samples = int(sample_sizes.max())
    sample_blocks = np.repeat(np.arange(batch_size), sample_sizes)
    rng = np.random.default_rng(seed)
    initial = rng.normal(scale=0.15, size=(batch_size, coefficient_count))
    padded_tensors = np.zeros(
        (batch_size, maximum_samples, dimension, coefficient_count)
    )
    padded_masks = np.zeros((batch_size, maximum_samples))
    ragged_tensors = []
    padded_counts = []
    ragged_counts = []
    mappings = []
    for primer in range(len(prepared[0].counts)):
        ec_count = prepared[0].counts[primer].shape[1]
        local_padded_counts = np.zeros((batch_size, maximum_samples, ec_count))
        local_ragged_counts = []
        local_mappings = []
        for block, data in enumerate(prepared):
            samples = len(data.clusters)
            if primer == 0:
                padded_tensors[block, :samples] = data.fixed_effect_tensor
                padded_masks[block, :samples] = 1.0
            local_padded_counts[block, :samples] = data.counts[primer]
            local_ragged_counts.append(data.counts[primer])
            local_mappings.append(data.compatibility[primer])
        padded_counts.append(jnp.asarray(local_padded_counts))
        ragged_counts.append(jnp.asarray(np.concatenate(local_ragged_counts)))
        mappings.append(jnp.asarray(np.stack(local_mappings)))
    ragged_tensors = jnp.asarray(
        np.concatenate([data.fixed_effect_tensor for data in prepared])
    )
    padded_tensors = jnp.asarray(padded_tensors)
    padded_masks = jnp.asarray(padded_masks)
    sample_blocks = jnp.asarray(sample_blocks)
    sample_sizes_jax = jnp.asarray(sample_sizes)

    def row_likelihood(observed, mass):
        total = jnp.sum(observed, axis=-1)
        probability = jnp.maximum(mass / jnp.sum(mass, axis=-1, keepdims=True), 1e-300)
        constant = jsp.special.gammaln(total + 1.0) - jnp.sum(
            jsp.special.gammaln(observed + 1.0), axis=-1
        )
        return constant + jnp.sum(observed * jnp.log(probability), axis=-1)

    def padded_objective(parameters):
        free = jnp.einsum("bnip,bp->bni", padded_tensors, parameters)
        logits = jnp.concatenate(
            (
                free,
                jnp.zeros((batch_size, maximum_samples, 1), dtype=free.dtype),
            ),
            axis=2,
        )
        abundance = jnp.exp(logits - jnp.max(logits, axis=2, keepdims=True))
        value = jnp.asarray(0.0)
        for observed, mapping in zip(padded_counts, mappings):
            mass = jnp.einsum("bkt,bnt->bnk", mapping, abundance)
            value += jnp.sum(padded_masks * row_likelihood(observed, mass))
        return -value

    def ragged_objective(parameters):
        free = jnp.einsum("nip,np->ni", ragged_tensors, parameters[sample_blocks])
        logits = jnp.concatenate(
            (free, jnp.zeros((len(free), 1), dtype=free.dtype)), axis=1
        )
        abundance = jnp.exp(logits - jnp.max(logits, axis=1, keepdims=True))
        value = jnp.asarray(0.0)
        for observed, mapping in zip(ragged_counts, mappings):
            mass = jax.lax.ragged_dot(
                abundance,
                jnp.swapaxes(mapping, 1, 2),
                sample_sizes_jax,
            )
            value += jnp.sum(row_likelihood(observed, mass))
        return -value

    results = []
    compiled = {}
    for name, objective in (
        ("padded", padded_objective),
        ("ragged_dot", ragged_objective),
    ):
        function = jax.jit(jax.value_and_grad(objective))
        started = time.perf_counter()
        value, gradient = function(jnp.asarray(initial))
        jax.block_until_ready((value, gradient))
        compile_seconds = time.perf_counter() - started
        timings = []
        for _ in range(int(repetitions)):
            started = time.perf_counter()
            value, gradient = function(jnp.asarray(initial))
            jax.block_until_ready((value, gradient))
            timings.append(time.perf_counter() - started)
        compiled[name] = (float(value), np.asarray(gradient))
        results.append(
            {
                "kernel": name,
                "blocks": batch_size,
                "samples": int(sample_sizes.sum()),
                "maximum_samples": maximum_samples,
                "sample_padding_ratio": (
                    batch_size * maximum_samples / sample_sizes.sum()
                ),
                "compile_seconds": compile_seconds,
                "median_evaluation_seconds": float(np.median(timings)),
                "minimum_evaluation_seconds": float(np.min(timings)),
            }
        )
    np.testing.assert_allclose(
        compiled["ragged_dot"][0], compiled["padded"][0], rtol=1e-10, atol=1e-8
    )
    np.testing.assert_allclose(
        compiled["ragged_dot"][1], compiled["padded"][1], rtol=1e-9, atol=1e-8
    )
    return pd.DataFrame(results)


def main():
    args = parse_args()
    key, prepared = prepare_candidates(args)
    result = benchmark(prepared, args.repetitions, args.seed)
    result.insert(0, "ragged_dot_key", key)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
