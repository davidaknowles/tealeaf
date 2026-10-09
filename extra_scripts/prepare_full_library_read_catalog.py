#!/usr/bin/env python3
"""Freeze a whole outcome-independent read-marker family and production cell QC."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import build_event_features, file_hash
from tealeaf.sc.differential import read_gtf_exons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library-root", type=Path, required=True)
    parser.add_argument("--cell-qc", type=Path, required=True)
    parser.add_argument("--reference-family", type=Path, required=True, help="Earlier coverage audit, for annotation only, never event selection.")
    parser.add_argument("--reference-family-manifest", type=Path, required=True)
    parser.add_argument("--event-catalog", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new full-family output required")
    library_recipe_path = args.library_root / "recipe.json"
    library_recipe = json.loads(library_recipe_path.read_text())
    qc_path = args.cell_qc / "manifest.json"
    qc = json.loads(qc_path.read_text())
    if not library_recipe.get("library_union") or len(library_recipe["bams"]) != 8 or qc["exact_cached_group_and_primer_total_match"] is not True:
        raise ValueError("complete eight-library source and exact production QC required")
    if qc["source_read_recipe_sha256"] not in library_recipe["diagnostic_inputs"].values():
        raise ValueError("cell QC must refer to the same original alignment metadata recipe")
    settings = json.loads(args.reference_family_manifest.read_text())["candidate_settings"]
    if settings["subject_fold"] is not None or settings["min_gene_umis"] != 25 or settings["min_celltype_mice"] != 4:
        raise ValueError("original full-data catalog coverage recipe required")
    reference = pd.read_csv(args.reference_family, sep="\t")
    reference = reference.loc[reference.source.eq("original_binary")].copy()
    if reference.empty or reference.feature_id.duplicated().any():
        raise ValueError("complete unique original reference audit required")
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    if not set(reference.feature_id) <= set(catalog.index):
        raise ValueError("source catalog must contain the entire earlier audited family")
    # Direct read-marker inference must not inherit a global-isoform-support
    # screen or an EC-dimension limit. Preserve every source event definition.
    selected = catalog.sort_index()
    family = selected[["gene_id", "event_type"]].reset_index()
    family["prior_global_isoform_supported"] = family.feature_id.isin(reference.feature_id)
    family["prior_requested_contrasts"] = family.feature_id.map(reference.set_index("feature_id").n_contrasts).fillna(0).astype(int)
    features, definitions = build_event_features(selected, read_gtf_exons(args.gtf))
    retained_path = args.cell_qc / "retained_barcodes.tsv.gz"
    retained = pd.read_csv(retained_path, sep="\t", dtype=str).set_index("barcode", verify_integrity=True)
    if len(retained) != qc["retained_production_barcodes"] or not set(retained.index) <= set(library_recipe["barcode_groups"]):
        raise ValueError("complete retained production barcode family required")
    groups = {barcode: list(retained.loc[barcode, ["subject", "cell_type", "primer"]]) for barcode in retained.index}
    if any(library_recipe["barcode_groups"][barcode] != group for barcode, group in groups.items()):
        raise ValueError("cell annotations must remain identical")
    packets = []
    seen = set()
    for packet in library_recipe["bams"]:
        local = {barcode: group for barcode, group in packet["barcode_groups"].items() if barcode in groups}
        if seen.intersection(local):
            raise ValueError("physical-library cell scopes must be disjoint")
        seen.update(local)
        packets.append(dict(packet, barcode_groups=local))
    if seen != set(groups):
        raise ValueError("every retained production cell must retain its physical library")
    inputs = {str(path): file_hash(path) for path in (library_recipe_path, qc_path, retained_path, args.reference_family, args.reference_family_manifest, args.event_catalog, args.gtf)}
    recipe = dict(events=features, barcode_groups=groups, bams=packets, library_union=True, full_catalog=True, declared_events=family.feature_id.tolist(), production_cell_qc=qc, source_hashes=inputs, selection="entire original annotation event catalog and production cell QC; no global-isoform support/EC-dimension screen, statistical fits, p-values, significance or LR outcomes", scope="whole local-marker count preparation only; unavailable/markerless definitions remain in declared family; no power, calibration or replication claim", production_changes=False)
    args.output_dir.mkdir(parents=True)
    definitions.to_csv(args.output_dir / "feature_definitions.tsv", sep="\t", index=False)
    family.to_csv(args.output_dir / "declared_family.tsv.gz", sep="\t", index=False)
    (args.output_dir / "recipe.json").write_text(json.dumps(recipe) + "\n")
    summary = dict(declared_events=len(family), marker_geometry_available=len(features), local_marker_status=definitions.status.value_counts().to_dict(), prior_global_isoform_supported_events=int(family.prior_global_isoform_supported.sum()), prior_requested_full_contrasts=int(family.prior_requested_contrasts.sum()), production_barcodes=len(groups), physical_libraries=len(packets), source_hashes=inputs, scope=recipe["scope"], selection=recipe["selection"], production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
