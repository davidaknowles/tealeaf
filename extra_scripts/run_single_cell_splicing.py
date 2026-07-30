#!/usr/bin/env python3
"""Fit splice paths with cell-specific baselines and mouse-level tests."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy.special

from extra_scripts import run_differential_splicing as pseudobulk
from tealeaf.sc import differential, glm_cv


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--transcript-to-gene", required=True, type=Path)
    parser.add_argument("--fit-prefix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--barcode-groups", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--block-cache", type=Path)
    parser.add_argument("--target-tests", type=Path)
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-pseudobulk-umis", type=float, default=100_000)
    parser.add_argument("--min-gene-umis", type=float, default=25)
    parser.add_argument("--min-condition-replicates", type=int, default=3)
    parser.add_argument("--max-paths", type=int, default=8)
    parser.add_argument("--max-logratio-variance", type=float, default=0.125)
    parser.add_argument("--min-path-proportion", type=float, default=0.15)
    parser.add_argument("--fit-max-iter", type=int, default=100)
    parser.add_argument("--permutations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def load_cell_state(args, prepared, group_index, cell_totals):
    fit_rows = np.loadtxt(f"{args.fit_prefix}glm_rows.txt", dtype=str)
    fit_transcripts = np.loadtxt(f"{args.fit_prefix}glm_cols.txt", dtype=str)
    if not np.array_equal(fit_transcripts, np.asarray(prepared.features, dtype=str)):
        raise ValueError("hierarchical fit transcripts do not match prepared data")
    prepared_position = {
        str(cell): index for index, cell in enumerate(prepared.barcodes)
    }
    positions = np.array(
        [prepared_position.get(str(cell), -1) for cell in fit_rows],
        dtype=np.int64,
    )
    if np.any(positions < 0):
        raise ValueError("hierarchical fit contains cells absent from preparation")
    with np.load(f"{args.fit_prefix}glm_factors.npz") as saved:
        state = saved["right"].astype(np.float64)
    with np.load(f"{args.fit_prefix}hierarchical_parameters.npz") as saved:
        parameters = {
            name: saved[name].astype(np.float64)
            if saved[name].dtype.kind == "f"
            else saved[name]
            for name in saved.files
        }
    fitted_groups = group_index[positions]
    retained = fitted_groups >= 0
    raw = prepared.cv_raw_counts[positions[retained]].tocsr()
    return (
        state[retained],
        parameters,
        fitted_groups[retained],
        cell_totals[positions[retained]],
        raw,
        fit_transcripts,
    )


def gene_log_normalizers(state, parameters, batch_cells=256):
    intercept = parameters["gene_intercept"]
    loadings = parameters["gene_loadings"]
    normalizers = np.empty(len(state), dtype=np.float64)
    for start in range(0, len(state), int(batch_cells)):
        stop = min(start + int(batch_cells), len(state))
        logits = np.clip(
            intercept[:, None] + loadings @ state[start:stop].T,
            -20,
            20,
        )
        normalizers[start:stop] = scipy.special.logsumexp(logits, axis=0)
    return normalizers


def gene_cell_baselines(
    state,
    parameters,
    transcript_indices,
    gene_cells,
    cell_totals,
    gene_normalizers,
):
    model_genes = np.unique(
        parameters["transcript_gene"][transcript_indices]
    )
    if len(model_genes) != 1:
        raise ValueError("local transcripts do not map to one model gene")
    model_gene = int(model_genes[0])
    logits = (
        parameters["isoform_intercept"][transcript_indices, None]
        + parameters["isoform_loadings"][transcript_indices]
        @ state[gene_cells].T
    ).T
    logits -= logits.max(axis=1, keepdims=True)
    baselines = np.exp(logits)
    baselines /= baselines.sum(axis=1, keepdims=True)
    gene_logits = np.clip(
        parameters["gene_intercept"][model_gene]
        + state[gene_cells] @ parameters["gene_loadings"][model_gene],
        -20,
        20,
    )
    weights = cell_totals[gene_cells] * np.exp(
        gene_logits - gene_normalizers[gene_cells]
    )
    return baselines, weights


def estimate_blocks(
    args,
    blocks,
    groups,
    fitted_groups,
    state,
    parameters,
    cell_totals,
    raw,
    feature_names,
    genes,
    gene_transcripts,
    gene_ecs,
    designs,
):
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    group_metadata = [pseudobulk.parse_group(group) for group in groups]
    gene_normalizers = gene_log_normalizers(state, parameters)
    n_ec = raw.shape[1] // 2
    inputs = defaultdict(list)
    records_path = args.output_dir / "single_cell_estimates.jsonl.gz"
    estimated = 0
    converged = 0
    reliable = 0
    cached_gene = None
    cached = None

    with gzip.open(records_path, "wt") as output:
        for block in blocks:
            gene = gene_to_index.get(block.gene_id)
            if gene is None:
                continue
            if gene != cached_gene:
                cached_gene = gene
                cached = None
                transcript_indices = gene_transcripts[gene]
                ec_rows = gene_ecs[gene]
                if len(transcript_indices) < 2 or not len(ec_rows):
                    continue
                local_designs = tuple(
                    design[ec_rows][:, transcript_indices].tocsr()
                    for design in designs
                )
                local_counts = (
                    raw[:, ec_rows].tocsr(),
                    raw[:, n_ec + ec_rows].tocsr(),
                )
                gene_umis_by_cell = sum(
                    np.asarray(value.sum(axis=1)).ravel()
                    for value in local_counts
                )
                gene_umis_by_group = np.bincount(
                    fitted_groups,
                    weights=gene_umis_by_cell,
                    minlength=len(groups),
                )
                eligible_groups = np.flatnonzero(
                    gene_umis_by_group >= float(args.min_gene_umis)
                )
                if len(eligible_groups) < 2:
                    continue
                gene_cells = np.flatnonzero(
                    np.isin(fitted_groups, eligible_groups)
                )
                baselines, gene_weights = gene_cell_baselines(
                    state,
                    parameters,
                    transcript_indices,
                    gene_cells,
                    cell_totals,
                    gene_normalizers,
                )
                cached = (
                    transcript_indices,
                    local_designs,
                    local_counts,
                    gene_umis_by_group,
                    eligible_groups,
                    gene_cells,
                    baselines,
                    gene_weights,
                )
            if cached is None:
                continue
            (
                transcript_indices,
                local_designs,
                local_counts,
                gene_umis_by_group,
                eligible_groups,
                gene_cells,
                baselines,
                gene_weights,
            ) = cached
            mapping_result = pseudobulk.block_mapping(
                block,
                feature_names[transcript_indices],
            )
            if mapping_result is None:
                continue
            path_index, _ = mapping_result
            n_paths = int(path_index.max()) + 1
            if n_paths > int(args.max_paths):
                continue

            for group in eligible_groups:
                local_cell_positions = np.flatnonzero(
                    fitted_groups[gene_cells] == group
                )
                if not len(local_cell_positions):
                    continue
                cell_rows = gene_cells[local_cell_positions]
                observed = tuple(
                    value[cell_rows].tocsr() for value in local_counts
                )
                fit = differential.fit_shared_path_perturbation(
                    observed,
                    local_designs,
                    baselines[local_cell_positions],
                    path_index,
                    gene_weights[local_cell_positions],
                    max_iter=args.fit_max_iter,
                )
                estimated += 1
                converged += int(fit.converged)
                is_reliable = bool(
                    fit.covariance.identifiable
                    and np.min(fit.path_proportions)
                    >= args.min_path_proportion
                    and np.linalg.eigvalsh(
                        fit.covariance.covariance
                    ).max() <= args.max_logratio_variance
                )
                reliable += int(is_reliable)
                cluster, condition, mouse = group_metadata[group]
                record = {
                    "block_id": block.block_id,
                    "gene_id": block.gene_id,
                    "group": groups[group],
                    "cluster": cluster,
                    "condition": condition,
                    "mouse": mouse,
                    "cells": len(cell_rows),
                    "gene_umis": float(gene_umis_by_group[group]),
                    "path_proportions": fit.path_proportions.tolist(),
                    "path_logratios": fit.path_logratios.tolist(),
                    "covariance": fit.covariance.covariance.tolist(),
                    "identifiable": fit.covariance.identifiable,
                    "reliable": is_reliable,
                    "fit_converged": fit.converged,
                    "fit_iterations": fit.iterations,
                }
                output.write(json.dumps(record) + "\n")
                if fit.converged and is_reliable:
                    inputs[
                        (block.block_id, cluster, "single_cell")
                    ].append(
                        (
                            fit.path_logratios,
                            fit.covariance.covariance,
                            condition,
                        )
                    )
            pseudobulk.emit(
                "single_cell_block_complete",
                block_id=block.block_id,
                estimates=estimated,
            )
    return inputs, {
        "blocks": len(blocks),
        "estimates": estimated,
        "converged": converged,
        "converged_fraction": converged / max(estimated, 1),
        "reliable": reliable,
        "reliable_fraction": reliable / max(estimated, 1),
    }


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    prepared = glm_cv.prepare_paired_primer_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        args.primer_pairs,
        ec_design="weighted",
        regularization_target="theta",
        min_eq=5,
        min_half_umis=args.min_half_umis,
        primer_sampling_model="oligodt_tpm",
    )
    groups, group_index, cell_totals, _ = pseudobulk.aggregate_inputs(
        args,
        prepared,
    )
    (
        state,
        parameters,
        fitted_groups,
        fitted_totals,
        raw,
        feature_names,
    ) = load_cell_state(args, prepared, group_index, cell_totals)
    genes, gene_transcripts, gene_ecs, designs = pseudobulk.gene_structures(
        prepared,
        args.transcript_to_gene,
    )
    blocks = pseudobulk.load_blocks(args.block_cache, args.gtf)
    if args.target_tests is not None:
        target_table = pd.read_csv(args.target_tests, sep="\t")
        target_blocks = set(target_table["block_id"].astype(str))
        blocks = [block for block in blocks if block.block_id in target_blocks]
    if args.max_blocks is not None:
        blocks = blocks[: int(args.max_blocks)]
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must be between zero and shard count")
    blocks = blocks[args.shard_index :: args.shard_count]
    pseudobulk.emit(
        "single_cell_preparation_complete",
        seconds=time.perf_counter() - started,
        cells=len(state),
        blocks=len(blocks),
    )
    inputs, estimate_summary = estimate_blocks(
        args,
        blocks,
        groups,
        fitted_groups,
        state,
        parameters,
        fitted_totals,
        raw,
        feature_names,
        genes,
        gene_transcripts,
        gene_ecs,
        designs,
    )
    test_summary = pseudobulk.differential_tests(
        args,
        inputs,
        np.random.default_rng(args.seed),
    )
    summary = {
        "seconds": time.perf_counter() - started,
        "cells": len(state),
        **estimate_summary,
        **test_summary,
    }
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    pseudobulk.emit("single_cell_analysis_complete", **summary)


if __name__ == "__main__":
    main()
