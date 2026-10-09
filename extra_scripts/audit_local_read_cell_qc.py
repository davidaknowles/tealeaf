#!/usr/bin/env python3
"""Reconstruct production retained cells before comparing read and EC models."""

import argparse
import json
from pathlib import Path
import pickle
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.run_differential_splicing import aggregate_inputs, parse_group
from tealeaf.data.alevin import load_alevin_counts
from tealeaf.sc.glm_cv import _read_primer_pairs, paired_primer_row_selection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--read-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new audit output required")
    control_manifest = args.control / "manifest.json"
    settings = json.loads(control_manifest.read_text())
    source = settings["settings"]
    pairs = _read_primer_pairs(source["primer_pairs"])
    barcodes, raw = load_alevin_counts(source["alevin_dir"])
    print(f"Loaded {len(barcodes)} source half-cells and {raw.nnz} EC entries", flush=True)
    complete, _ = paired_primer_row_selection(barcodes, np.asarray(raw.sum(axis=1)).ravel(), pairs, min_half_umis=source["min_half_umis"])
    if len(complete) != settings["paired_cells"]:
        raise ValueError("reconstructed complete-pair count differs from frozen preparation")
    first, second = [np.asarray([row[index] for row in complete], dtype=int) for index in (1, 2)]
    aggregate = np.asarray(raw[first].sum(axis=0) + raw[second].sum(axis=0)).ravel()
    ec_keep = aggregate >= 5.
    totals = np.asarray(raw[:, ec_keep].sum(axis=1)).ravel()
    del raw, aggregate
    # A two-column summary preserves both primer totals and every cell needed
    # by production group QC. It is not a reconstructed EC response matrix.
    summary = sparse.csr_matrix(np.column_stack([totals[first], totals[second]]))
    prepared = SimpleNamespace(barcodes=np.asarray([row[0] for row in complete]), cv_raw_counts=summary)
    group_args = SimpleNamespace(primer_pairs=source["primer_pairs"], barcode_groups=source["barcode_groups"], min_cells=source["min_cells"], min_pseudobulk_umis=source["min_pseudobulk_umis"])
    groups, group_index, _, group_counts = aggregate_inputs(group_args, prepared)
    with (args.control / "prepared.pkl").open("rb") as handle:
        cached = pickle.load(handle)
    primer_totals = [float(value.sum()) for value in group_counts]
    cached_totals = [float(value.sum()) for value in cached[1]]
    if list(groups) != list(cached[0]) or primer_totals != cached_totals or primer_totals != settings["primer_UMIs"]:
        raise ValueError("production pseudobulk identities or exact primer totals changed")
    del cached
    recipe_path = args.read_root / "recipe.json"
    recipe = json.loads(recipe_path.read_text())
    kept, rows = {}, []
    for index, (_, poly_row, hex_row) in enumerate(complete):
        if group_index[index] < 0:
            continue
        cell_type, _, subject = parse_group(groups[group_index[index]])
        for row, primer in ((poly_row, "poly(dT)"), (hex_row, "random hexamer")):
            barcode = barcodes[row]
            group = (subject, cell_type, primer)
            if barcode not in recipe["barcode_groups"] or tuple(recipe["barcode_groups"][barcode]) != group:
                raise ValueError("retained production cell lacks the same frozen read annotation")
            kept[barcode] = group
            rows.append(dict(barcode=barcode, subject=subject, cell_type=cell_type, primer=primer, retained_half_umis=float(totals[row])))
    comparison = []
    for barcode, group in recipe["barcode_groups"].items():
        comparison.append(dict(barcode=barcode, subject=group[0], cell_type=group[1], primer=group[2], production_cell_retained=barcode in kept))
    args.output_dir.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "retained_barcodes.tsv.gz", sep="\t", index=False)
    pd.DataFrame(comparison).to_csv(args.output_dir / "diagnostic_barcode_qc.tsv.gz", sep="\t", index=False)
    report = dict(control_manifest_sha256=file_hash(control_manifest), source_read_recipe_sha256=file_hash(recipe_path), source_barcode_groups_sha256=file_hash(source["barcode_groups"]), source_primer_pairs_sha256=file_hash(source["primer_pairs"]), count_input=dict(path=source["alevin_dir"], rows=len(barcodes), min_eq=5, complete_pairs=len(complete), min_half_umis=source["min_half_umis"]), groups=len(groups), retained_production_barcodes=len(kept), requested_diagnostic_barcodes=len(recipe["barcode_groups"]), diagnostic_barcodes_not_retained=len(recipe["barcode_groups"]) - len(kept), exact_cached_group_and_primer_total_match=True, primer_UMIs=primer_totals, scope="production cell eligibility reconstruction only; no new counts, effect tests, discovery or LR claim", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
