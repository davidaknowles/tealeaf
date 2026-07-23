"""Utilities for combining independently quantified alevin-fry batches."""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

from tealeaf.sc import sc_utils


def _read_lines(path):
    with open(path) as handle:
        return [line.rstrip("\n") for line in handle]


def _load_quantification(path):
    path = Path(path)
    features = _read_lines(path / "quants_mat_cols.txt")
    barcodes = _read_lines(path / "quants_mat_rows.txt")
    count_cache = path / "geqc_counts.npz"
    counts = (
        sp.load_npz(count_cache).tocsr()
        if count_cache.is_file()
        else scipy.io.mmread(path / "geqc_counts.mtx").tocsr()
    )
    map_cache = path / "gene_eqclass.npz"
    if map_cache.is_file():
        membership = sp.load_npz(map_cache).tocsr()
    else:
        _, _, ecs = sc_utils.read_alevin_ec(path / "gene_eqclass.txt.gz")
        membership = sc_utils.to_coo(
            [ecs[index] for index in range(len(ecs))],
            shape=(len(ecs), len(features)),
        ).tocsr()
    if counts.shape != (len(barcodes), membership.shape[0]):
        raise ValueError(f"inconsistent count dimensions in {path}")
    if membership.shape[1] != len(features):
        raise ValueError(f"inconsistent feature dimensions in {path}")
    return features, barcodes, counts, membership


def merge_alevin_quantifications(inputs, output_dir) -> dict:
    """Merge batches by transcript-set EC identity and prefix cell barcodes.

    Parameters
    ----------
    inputs
        Sequence of ``(batch_name, quantification_directory)`` pairs.
    output_dir
        New directory receiving sparse caches and row/column labels.
    """
    inputs = [(str(name), Path(path)) for name, path in inputs]
    if not inputs:
        raise ValueError("at least one alevin quantification is required")
    if len({name for name, _ in inputs}) != len(inputs):
        raise ValueError("batch names must be unique")

    loaded = [_load_quantification(path) for _, path in inputs]
    feature_order = []
    feature_seen = set()
    for features, _, _, _ in loaded:
        for feature in features:
            if feature not in feature_seen:
                feature_seen.add(feature)
                feature_order.append(feature)
    feature_to_global = {feature: i for i, feature in enumerate(feature_order)}

    ec_to_global = {}
    ec_features = []
    remapped = []
    probability_remaps = []
    all_barcodes = []
    cell_offset = 0
    for (batch, _), (features, barcodes, counts, membership) in zip(inputs, loaded):
        local_to_global_feature = np.asarray(
            [feature_to_global[feature] for feature in features], dtype=np.int64
        )
        local_ec_to_global = np.empty(membership.shape[0], dtype=np.int64)
        local_probability_orders = []
        for ecid in range(membership.shape[0]):
            start, stop = membership.indptr[ecid : ecid + 2]
            local_global_features = local_to_global_feature[
                membership.indices[start:stop]
            ]
            order = np.argsort(local_global_features)
            key = tuple(local_global_features[order])
            global_ecid = ec_to_global.get(key)
            if global_ecid is None:
                global_ecid = len(ec_features)
                ec_to_global[key] = global_ecid
                ec_features.append(key)
            local_ec_to_global[ecid] = global_ecid
            local_probability_orders.append(order)
        coo = counts.tocoo()
        remapped.append(
            sp.csr_matrix(
                (coo.data, (coo.row, local_ec_to_global[coo.col])),
                shape=(counts.shape[0], len(ec_features)),
            )
        )
        all_barcodes.extend(f"{batch}:{barcode}" for barcode in barcodes)
        probability_remaps.append(
            (local_ec_to_global, local_probability_orders, cell_offset)
        )
        cell_offset += len(barcodes)

    n_ec = len(ec_features)
    remapped = [
        matrix if matrix.shape[1] == n_ec else
        sp.csr_matrix(
            (matrix.data, matrix.indices, matrix.indptr),
            shape=(matrix.shape[0], n_ec),
        )
        for matrix in remapped
    ]
    counts = sp.vstack(remapped, format="csr")
    membership = sc_utils.to_coo(
        ec_features, shape=(n_ec, len(feature_order))
    ).tocsr()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    sp.save_npz(output_dir / "geqc_counts.npz", counts)
    sp.save_npz(output_dir / "gene_eqclass.npz", membership)
    with open(output_dir / "quants_mat_rows.txt", "w") as handle:
        handle.write("\n".join(all_barcodes) + "\n")
    with open(output_dir / "quants_mat_cols.txt", "w") as handle:
        handle.write("\n".join(feature_order) + "\n")
    probability_inputs = [
        path / "gene_eqclass_probs.tsv.gz" for _, path in inputs
    ]
    if all(path.is_file() for path in probability_inputs):
        probability_output = output_dir / "gene_eqclass_probs.tsv.gz"
        with gzip.open(probability_output, "wt") as output:
            output.write("cell_idx\teqid\tumi_rank\tprobs\n")
            for probability_path, (
                ec_remap,
                probability_orders,
                offset,
            ) in zip(probability_inputs, probability_remaps):
                with gzip.open(probability_path, "rt") as source:
                    header = next(source, "").rstrip("\n")
                    if header != "cell_idx\teqid\tumi_rank\tprobs":
                        raise ValueError(
                            f"unexpected probability header in {probability_path}"
                        )
                    for line in source:
                        cell, ecid, rank, values = line.rstrip("\n").split("\t", 3)
                        ecid = int(ecid)
                        probabilities = np.fromstring(values, sep=",")
                        order = probability_orders[ecid]
                        if probabilities.size != order.size:
                            raise ValueError(
                                f"probability length mismatch for EC {ecid} "
                                f"in {probability_path}"
                            )
                        reordered = probabilities[order]
                        encoded = ",".join(f"{value:.9g}" for value in reordered)
                        output.write(
                            f"{int(cell) + offset}\t{ec_remap[ecid]}\t"
                            f"{rank}\t{encoded}\n"
                        )
    return {
        "batches": len(inputs),
        "cells": counts.shape[0],
        "equivalence_classes": counts.shape[1],
        "features": membership.shape[1],
        "molecules": float(counts.sum()),
        "probabilities_merged": all(path.is_file() for path in probability_inputs),
    }
