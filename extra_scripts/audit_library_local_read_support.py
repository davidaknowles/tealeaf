#!/usr/bin/env python3
"""Physical-library molecule union for the frozen local-support diagnostic."""

import argparse
import json
from pathlib import Path
import shutil

import pandas as pd

from extra_scripts.audit_event_local_read_support import collate, decode_event, file_hash
from tealeaf.sc.local_read_support import collect_indexed_library_read_support


def prepare(args):
    if args.output_dir.exists():
        raise ValueError("preserve previous outputs, use a new directory")
    source_path = args.source_root / "recipe.json"
    source = json.loads(source_path.read_text())
    provenance_path = args.origins / "manifest.json"
    provenance = json.loads(provenance_path.read_text())
    if not provenance["complete_unique_library_assignment"] or provenance["source_recipe_sha256"] != file_hash(source_path):
        raise ValueError("complete unique physical origins for this frozen recipe required")
    origins_path, runs_path = args.origins / "barcode_library.tsv.gz", args.origins / "library_runs.tsv"
    origins = pd.read_csv(origins_path, sep="\t", dtype=str).set_index("barcode", verify_integrity=True)
    runs = pd.read_csv(runs_path, sep="\t", dtype=str)
    if not origins.status.eq("unique").all() or set(origins.index) != set(source["barcode_groups"]):
        raise ValueError("all original barcodes must have unique physical origins")
    for barcode, group in source["barcode_groups"].items():
        if list(origins.loc[barcode, ["subject", "cell_type", "primer"]]) != list(group):
            raise ValueError("original cell/subject/primer annotations must not change")
    groups = source["barcode_groups"]
    qc_inputs = []
    qc = None
    if getattr(args, "cell_qc", None) is not None:
        qc_manifest = args.cell_qc / "manifest.json"
        qc_barcodes = args.cell_qc / "retained_barcodes.tsv.gz"
        qc = json.loads(qc_manifest.read_text())
        retained = pd.read_csv(qc_barcodes, sep="\t", dtype=str).set_index("barcode", verify_integrity=True)
        if qc["exact_cached_group_and_primer_total_match"] is not True or qc["source_read_recipe_sha256"] != file_hash(source_path) or len(retained) != qc["retained_production_barcodes"] or not set(retained.index) <= set(groups):
            raise ValueError("exact matching production cell QC required")
        restricted = {barcode: list(retained.loc[barcode, ["subject", "cell_type", "primer"]]) for barcode in retained.index}
        if any(groups[barcode] != group for barcode, group in restricted.items()):
            raise ValueError("QC may restrict cells but must not change their annotations")
        groups = restricted
        origins = origins.loc[list(groups)]
        qc_inputs = [qc_manifest, qc_barcodes]
    packets = {Path(row["path"]).parent.name: (index, row) for index, row in enumerate(source["bams"])}
    declared = [value for row in runs.itertuples(index=False) for value in (row.run_a, row.run_b)]
    if runs.library.duplicated().any() or len(set(declared)) != len(declared) or set(declared) != set(packets) or set(origins.library) != set(runs.library):
        raise ValueError("complete disjoint alignment-run partition required")
    libraries = []
    for row in runs.itertuples(index=False):
        inputs = []
        for run in (row.run_a, row.run_b):
            index, packet = packets[run]
            shard = args.source_root / f"shard_{index}"
            receipt_path = shard / "manifest.json"
            receipt = json.loads(receipt_path.read_text())
            index_path = shard / "alignments.bai"
            if receipt["complete"] is not True or receipt["input"] != packet or receipt["recipe_sha256"] != file_hash(source_path) or not index_path.is_file():
                raise ValueError("complete original scan and index required for each run")
            inputs.append(dict(**packet, index=str(index_path.resolve()), index_sha256=file_hash(index_path), scan_receipt=str(receipt_path.resolve()), scan_receipt_sha256=file_hash(receipt_path)))
        local_groups = {barcode: groups[barcode] for barcode in origins.index[origins.library.eq(row.library)]}
        libraries.append(dict(library=row.library, runs=inputs, barcode_groups=local_groups))
    recipe = dict(source)
    recipe.update(bams=libraries, barcode_groups=groups, library_union=True, production_cell_qc=qc, scope="same-read local-marker diagnostic, exact CB/UB/gene union within physical library across its runs; " + ("exact production retained-cell QC; " if qc is not None else "broader barcode-metadata cell scope; ") + "not upstream EC UMI counts, RNA PSI, independent validation or a statistical eligibility rule")
    recipe["diagnostic_inputs"] = dict(source["diagnostic_inputs"])
    recipe["diagnostic_inputs"].update({str(path.resolve()): file_hash(path) for path in (source_path, provenance_path, origins_path, runs_path, *qc_inputs)})
    args.output_dir.mkdir(parents=True)
    for name in ("selected_cases.tsv.gz", "feature_definitions.tsv"):
        shutil.copyfile(args.source_root / name, args.output_dir / name)
    if file_hash(args.output_dir / "selected_cases.tsv.gz") != source["selected_sha256"]:
        raise ValueError("frozen selected-case identities changed")
    (args.output_dir / "recipe.json").write_text(json.dumps(recipe) + "\n")
    print(f"Prepared {len(libraries)} libraries, {len(origins)} unique barcodes, {len(source['events'])} frozen events", flush=True)


def collect(args):
    recipe_path = args.output_dir / "recipe.json"
    recipe = json.loads(recipe_path.read_text())
    if not recipe.get("library_union") or not 0 <= args.shard_index < len(recipe["bams"]):
        raise ValueError("declared library-union shard required")
    packet = recipe["bams"][args.shard_index]
    for source in packet["runs"]:
        bam = Path(source["path"])
        if bam.stat().st_size != source["size"] or bam.stat().st_mtime_ns != source["mtime_ns"] or file_hash(source["index"]) != source["index_sha256"] or file_hash(source["scan_receipt"]) != source["scan_receipt_sha256"]:
            raise ValueError("frozen BAM, index or source receipt changed")
    shard = args.output_dir / f"shard_{args.shard_index}"
    if shard.exists():
        raise ValueError("preserve completed or partial shards, never overwrite")
    events = [decode_event(row) for row in recipe["events"]]
    counts, filters = collect_indexed_library_read_support([row["path"] for row in packet["runs"]], [row["index"] for row in packet["runs"]], packet["barcode_groups"], events)
    shard.mkdir()
    rows = [dict(feature_id=key[0], subject=key[1], cell_type=key[2], primer=key[3], signature=key[4], count=value) for key, value in counts.items()]
    pd.DataFrame(rows, columns=["feature_id", "subject", "cell_type", "primer", "signature", "count"]).to_csv(shard / "support.tsv.gz", sep="\t", index=False)
    receipt = dict(input=packet, recipe_sha256=file_hash(recipe_path), shard=args.shard_index, requested_events=len(events), filters=dict(filters), complete=True)
    (shard / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(dict(library=packet["library"], filters=dict(filters)), indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "collect", "collate"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--origins", type=Path)
    parser.add_argument("--cell-qc", type=Path, help="Exact production retained-cell reconstruction, without changing case selection.")
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--public-dir", type=Path)
    parser.add_argument("--bound-root", type=Path, default=Path("analyses/event_sequence_variance_limits"))
    args = parser.parse_args()
    if args.mode == "prepare":
        if args.source_root is None or args.origins is None:
            parser.error("source counts and physical origins required")
        prepare(args)
    elif args.mode == "collect":
        if args.shard_index is None:
            parser.error("library shard required")
        collect(args)
    else:
        collate(args)


if __name__ == "__main__":
    main()
