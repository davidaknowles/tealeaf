"""Utilities for combining independent alevin-fry quantifications."""

from __future__ import annotations

import csv
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


def _resolve_quantification_dir(path):
    """Return the directory containing alevin matrix and EC files."""
    path = Path(path)
    if (path / "quants_mat_rows.txt").is_file():
        return path
    nested = path / "alevin"
    if (nested / "quants_mat_rows.txt").is_file():
        return nested
    raise FileNotFoundError(
        f"no alevin quantification found in {path} or {nested}"
    )


def _load_structure(path):
    path = Path(path)
    features = _read_lines(path / "quants_mat_cols.txt")
    map_cache = path / "gene_eqclass.npz"
    if map_cache.is_file():
        membership = sp.load_npz(map_cache).tocsr()
    else:
        _, _, ecs = sc_utils.read_alevin_ec(path / "gene_eqclass.txt.gz")
        membership = sc_utils.to_coo(
            [ecs[index] for index in range(len(ecs))],
            shape=(len(ecs), len(features)),
        ).tocsr()
    if membership.shape[1] != len(features):
        raise ValueError(f"inconsistent feature dimensions in {path}")
    return features, membership


def _load_counts(path):
    path = Path(path)
    barcodes = _read_lines(path / "quants_mat_rows.txt")
    count_cache = path / "geqc_counts.npz"
    counts = (
        sp.load_npz(count_cache).tocsr()
        if count_cache.is_file()
        else scipy.io.mmread(path / "geqc_counts.mtx").tocsr()
    )
    if counts.shape[0] != len(barcodes):
        raise ValueError(f"inconsistent count dimensions in {path}")
    return barcodes, counts


def validate_alevin_quantification(
    path,
    *,
    expected_prefixes=None,
    reference_ids=None,
    primer_pair_file=None,
    min_cell_umis=500,
    min_reference_overlap=1,
    max_total_molecules=None,
):
    """Validate a merged quantification before downstream model fitting."""
    path = _resolve_quantification_dir(path)
    features, membership = _load_structure(path)
    barcodes, counts = _load_counts(path)
    if counts.shape[1] != membership.shape[0]:
        raise ValueError("count columns do not match equivalence classes")
    if len(features) != membership.shape[1]:
        raise ValueError("feature rows do not match the compatibility matrix")
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("cell identifiers are not unique")
    if len(features) != len(set(features)):
        raise ValueError("feature identifiers are not unique")
    if np.any(counts.data < 0) or np.any(membership.data < 0):
        raise ValueError("quantification matrices contain negative values")
    if counts.nnz == 0 or membership.nnz == 0:
        raise ValueError("quantification matrices are empty")

    barcode_set = set(barcodes)
    observed_prefixes = {
        barcode.split(":", 1)[0] for barcode in barcodes if ":" in barcode
    }
    expected_prefixes = set(expected_prefixes or ())
    missing_prefixes = sorted(expected_prefixes - observed_prefixes)
    if missing_prefixes:
        raise ValueError(f"missing expected cell prefixes: {missing_prefixes}")

    totals = np.asarray(counts.sum(axis=1)).ravel()
    total_molecules = float(totals.sum())
    if (
        max_total_molecules is not None
        and total_molecules > float(max_total_molecules)
    ):
        raise ValueError(
            f"quantification has {total_molecules:g} molecules, exceeding "
            f"the limit {float(max_total_molecules):g}"
        )
    eligible = int(np.count_nonzero(totals >= float(min_cell_umis)))
    if eligible == 0:
        raise ValueError(f"no cells have at least {min_cell_umis:g} UMIs")

    reference_ids = set(reference_ids or ())
    reference_overlap = len(barcode_set & reference_ids)
    if reference_ids and reference_overlap < int(min_reference_overlap):
        raise ValueError(
            f"only {reference_overlap} reference IDs overlap; "
            f"minimum is {min_reference_overlap}"
        )

    primer_pairs = complete_primer_pairs = 0
    if primer_pair_file is not None:
        with open(primer_pair_file, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                primer_pairs += 1
                if (
                    row["polydt_barcode"] in barcode_set
                    and row["ranhex_barcode"] in barcode_set
                ):
                    complete_primer_pairs += 1
        if primer_pairs == 0 or complete_primer_pairs == 0:
            raise ValueError("no complete primer pairs are represented")

    probability_file = path / "gene_eqclass_probs.tsv.gz"
    if not probability_file.is_file():
        raise FileNotFoundError(
            f"missing weighted probability sidecar: {probability_file}"
        )
    with gzip.open(probability_file, "rt") as handle:
        header = next(handle, "").rstrip("\n")
        first_probability = next(handle, "").rstrip("\n")
    if header != "cell_idx\teqid\tumi_rank\tprobs":
        raise ValueError(f"unexpected probability sidecar header: {header!r}")
    if not first_probability:
        raise ValueError("weighted probability sidecar has no data rows")

    return {
        "cells": int(counts.shape[0]),
        "equivalence_classes": int(counts.shape[1]),
        "features": int(membership.shape[1]),
        "count_nonzeros": int(counts.nnz),
        "compatibility_nonzeros": int(membership.nnz),
        "molecules": total_molecules,
        "min_cell_umis": float(min_cell_umis),
        "eligible_cells": eligible,
        "observed_prefixes": len(observed_prefixes),
        "expected_prefixes": len(expected_prefixes),
        "reference_ids": len(reference_ids),
        "reference_overlap": reference_overlap,
        "primer_pairs": primer_pairs,
        "complete_primer_pairs": complete_primer_pairs,
    }


def merge_alevin_quantifications(inputs, output_dir) -> dict:
    """Merge runs by transcript-set EC identity and prefix cell barcodes.

    Parameters
    ----------
    inputs
        Sequence of ``(run_name, quantification_directory)`` pairs.
    output_dir
        New directory receiving sparse caches and row/column labels.
    """
    inputs = [
        (str(name), _resolve_quantification_dir(path)) for name, path in inputs
    ]
    if not inputs:
        raise ValueError("at least one alevin quantification is required")
    if len({name for name, _ in inputs}) != len(inputs):
        raise ValueError("batch names must be unique")

    feature_order = []
    feature_to_global = {}
    ec_to_global = {}
    ec_features = []
    ec_remaps = []
    for _, path in inputs:
        features, membership = _load_structure(path)
        for feature in features:
            if feature not in feature_to_global:
                feature_to_global[feature] = len(feature_order)
                feature_order.append(feature)
        local_to_global_feature = np.asarray(
            [feature_to_global[feature] for feature in features], dtype=np.int64
        )
        local_ec_to_global = np.empty(membership.shape[0], dtype=np.int64)
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
        ec_remaps.append(local_ec_to_global)

    remapped = []
    all_barcodes = []
    cell_offsets = []
    cell_offset = 0
    n_ec = len(ec_features)
    for (batch, path), local_ec_to_global in zip(inputs, ec_remaps):
        barcodes, counts = _load_counts(path)
        if counts.shape[1] != local_ec_to_global.size:
            raise ValueError(f"inconsistent count/EC dimensions in {path}")
        coo = counts.tocoo()
        remapped.append(
            sp.csr_matrix(
                (coo.data, (coo.row, local_ec_to_global[coo.col])),
                shape=(counts.shape[0], n_ec),
            )
        )
        all_barcodes.extend(f"{batch}:{barcode}" for barcode in barcodes)
        cell_offsets.append(cell_offset)
        cell_offset += len(barcodes)

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
        with gzip.open(probability_output, "wt", compresslevel=1) as output:
            output.write("cell_idx\teqid\tumi_rank\tprobs\n")
            for probability_path, ec_remap, offset, (_, input_path) in zip(
                probability_inputs, ec_remaps, cell_offsets, inputs
            ):
                features, local_membership = _load_structure(input_path)
                local_to_global_feature = np.asarray(
                    [feature_to_global[feature] for feature in features],
                    dtype=np.int64,
                )
                probability_orders = {}
                probability_sizes = np.diff(local_membership.indptr)
                for ecid in range(local_membership.shape[0]):
                    start, stop = local_membership.indptr[ecid : ecid + 2]
                    global_features = local_to_global_feature[
                        local_membership.indices[start:stop]
                    ]
                    identity = bool(
                        global_features.size < 2
                        or np.all(global_features[:-1] <= global_features[1:])
                    )
                    if not identity:
                        probability_orders[ecid] = np.argsort(global_features)
                with gzip.open(probability_path, "rt") as source:
                    header = next(source, "").rstrip("\n")
                    if header != "cell_idx\teqid\tumi_rank\tprobs":
                        raise ValueError(
                            f"unexpected probability header in {probability_path}"
                        )
                    for line in source:
                        cell, ecid, rank, values = line.rstrip("\n").split("\t", 3)
                        ecid = int(ecid)
                        order = probability_orders.get(ecid)
                        probability_count = values.count(",") + 1
                        if probability_count != probability_sizes[ecid]:
                            raise ValueError(
                                f"probability length mismatch for EC {ecid} "
                                f"in {probability_path}"
                            )
                        if order is None:
                            encoded = values
                        else:
                            probabilities = np.fromstring(values, sep=",")
                            reordered = probabilities[order]
                            encoded = ",".join(
                                f"{value:.17g}" for value in reordered
                            )
                        output.write(
                            f"{int(cell) + offset}\t{ec_remap[ecid]}\t"
                            f"{rank}\t{encoded}\n"
                        )
    return {
        "quantifications": len(inputs),
        "cells": counts.shape[0],
        "equivalence_classes": counts.shape[1],
        "features": membership.shape[1],
        "molecules": float(counts.sum()),
        "probabilities_merged": all(path.is_file() for path in probability_inputs),
    }
