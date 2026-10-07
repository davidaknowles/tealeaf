"""Trace event transcript loss to input EC counts, weighting or indexing.

Read-only input audit. Compatibility support is not a transcript abundance
estimate: one ambiguous molecule can support several transcripts. No external
agreement is used for input preparation or event selection.
"""

import argparse
import gc
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.run_suppa2_tealeaf_hybrid import canonical
from tealeaf.data.alevin import load_alevin_counts, load_alevin_structure
from tealeaf.sc.glm_cv import _transcript_gene_assignment


def binary_transcript_support(membership, ec_totals, minimum_ec=5):
    """Return observed and retained binary support, each length T, not TPM."""
    membership = membership.tocsr(copy=True)
    membership.data[:] = 1
    totals = np.asarray(ec_totals, dtype=float)
    if totals.shape != (membership.shape[0],) or np.any(totals < 0):
        raise ValueError("nonnegative EC totals must align with membership")
    all_support = np.asarray(totals @ membership).ravel()
    retained = np.asarray(np.where(totals >= minimum_ec, totals, 0) @ membership).ravel()
    return all_support, retained


def count_fingerprint(values):
    """Exact CSR shape/ordering/value fingerprint, not merely matching totals."""
    values = values.tocsr()
    digest = hashlib.sha256()
    digest.update(str(values.shape).encode())
    for array in (values.indptr, values.indices, values.data):
        digest.update(str(array.dtype).encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", action="append", required=True, help="label=alevin-directory")
    parser.add_argument("--prepared-control", action="append", default=[], help="label=control-directory containing prepared.pkl and features.txt")
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--transcript-to-gene", type=Path, required=True)
    parser.add_argument("--event-support", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    features = args.features.read_text().splitlines()
    with args.data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(handle)
    assignment, mapped_genes = _transcript_gene_assignment(features, args.transcript_to_gene)
    index_mismatches = []
    for gene, indices in zip(genes, gene_transcripts):
        for index in indices:
            mapped = assignment.indices[assignment.indptr[index]:assignment.indptr[index + 1]]
            if len(mapped) != 1 or mapped_genes[mapped[0]] != gene:
                index_mismatches.append({"feature_index": int(index), "feature_id": features[index], "cache_gene": gene, "mapped_gene": mapped_genes[mapped[0]] if len(mapped) == 1 else None})
    prepared_mass = np.zeros(len(features))
    for design in designs:
        prepared_mass += np.asarray(design.sum(axis=0)).ravel()
    table = pd.DataFrame({"transcript_id": [canonical(value) for value in features], "prepared_EC_column_mass": prepared_mass})
    table = table.groupby("transcript_id", as_index=False).sum()
    production_count_fingerprints = [count_fingerprint(values) for values in counts]
    manifest = {"prepared_features": len(features), "prepared_groups": len(groups), "prepared_nonzero_columns": int(np.sum(prepared_mass > 0)), "index_mismatch_count": len(index_mismatches), "production_count_fingerprints": production_count_fingerprints, "support_interpretation": "binary compatibility duplicated across transcripts, not assigned counts or TPM", "sources": []}
    del counts, designs, assignment
    gc.collect()
    for value in args.source:
        label, directory = value.split("=", 1)
        names, membership = load_alevin_structure(directory)
        barcodes, source_counts = load_alevin_counts(directory)
        totals = np.asarray(source_counts.sum(axis=0)).ravel()
        support, retained = binary_transcript_support(membership, totals)
        manifest["sources"].append({"label": label, "directory": directory, "features": len(names), "cells": len(barcodes), "ECs": membership.shape[0], "counts_sum": float(totals.sum()), "membership_nnz": membership.nnz, "positive_support_transcripts": int(np.sum(support > 0)), "positive_support_after_EC5": int(np.sum(retained > 0))})
        local = pd.DataFrame({"transcript_id": [canonical(value) for value in names], f"{label}_all_support": support, f"{label}_EC5_support": retained}).groupby("transcript_id", as_index=False).sum()
        table = table.merge(local, on="transcript_id", how="outer", validate="one_to_one")
        del membership, source_counts, totals, support, retained, local
        gc.collect()
        print(json.dumps(manifest["sources"][-1]), flush=True)
    manifest["prepared_controls"] = []
    for value in args.prepared_control:
        label, directory = value.split("=", 1)
        directory = Path(directory)
        names = (directory / "features.txt").read_text().splitlines()
        with (directory / "prepared.pkl").open("rb") as handle:
            control_groups, control_counts, control_genes, control_tx, control_ecs, control_designs = pickle.load(handle)
        global_mass = np.sum([np.asarray(design.sum(axis=0)).ravel() for design in control_designs], axis=0)
        gene_mass = np.zeros(len(names))
        for indices, ecs in zip(control_tx, control_ecs):
            for design in control_designs:
                gene_mass[indices] += np.asarray(design[np.asarray(ecs, dtype=int)][:, np.asarray(indices, dtype=int)].sum(axis=0)).ravel()
        local = pd.DataFrame({"transcript_id": [canonical(value) for value in names], f"{label}_global_mass": global_mass, f"{label}_gene_mass": gene_mass}).groupby("transcript_id", as_index=False).sum()
        table = table.merge(local, on="transcript_id", how="outer", validate="one_to_one")
        fingerprints = [count_fingerprint(values) for values in control_counts]
        manifest["prepared_controls"].append({"label": label, "directory": str(directory), "same_feature_order_as_production": names == features, "pseudobulks": len(control_groups), "same_groups_as_production": list(control_groups) == list(groups), "same_exact_count_arrays_as_production": fingerprints == production_count_fingerprints, "count_fingerprints": fingerprints, "primer_UMIs": [float(values.sum()) for values in control_counts], "globally_supported": int(np.sum(global_mass > 0)), "gene_supported": int(np.sum(gene_mass > 0))})
        del control_counts, control_designs, global_mass, gene_mass, local
        gc.collect()
    table = table.fillna(0)
    table.to_csv(args.output_dir / "transcript_support.tsv.gz", sep="\t", index=False)
    lookup = table.set_index("transcript_id")
    events = pd.read_csv(args.event_support, sep="\t").fillna({"missing_transcripts": ""})
    rows = []
    for event in events.itertuples(index=False):
        missing = [value for value in event.missing_transcripts.split(",") if value]
        values = lookup.reindex(missing).fillna(0)
        rows.append({"native_rank": event.native_rank, "contrast_id": event.contrast_id, "feature_id": event.feature_id, "n_missing": len(missing), **{f"n_missing_{column}": int(values[column].gt(0).sum()) for column in values}})
    rows = pd.DataFrame(rows)
    rows.to_csv(args.output_dir / "event_missing_support.tsv.gz", sep="\t", index=False)
    pd.concat([rows.loc[rows.native_rank.le(top)].assign(top=top) for top in (40, 100, 200)]).groupby("top").sum(numeric_only=True).to_csv(args.output_dir / "top_missing_support.tsv", sep="\t")
    (args.output_dir / "index_mismatches.json").write_text(json.dumps(index_mismatches, indent=2) + "\n")
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({key: value for key, value in manifest.items() if key != "sources"}), flush=True)


if __name__ == "__main__":
    main()
