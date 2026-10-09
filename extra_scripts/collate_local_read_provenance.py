#!/usr/bin/env python3
"""Archive complete read-key duplication and exact retained-cell QC evidence."""

import argparse
import json
from pathlib import Path
import shutil

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--cell-qc", type=Path, required=True)
    parser.add_argument("--full-catalog", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new provenance archive required")
    support_path = args.support / "manifest.json"
    qc_path = args.cell_qc / "manifest.json"
    catalog_path = args.full_catalog / "manifest.json"
    support, qc, catalog = [json.loads(path.read_text()) for path in (support_path, qc_path, catalog_path)]
    if len(support["shards"]) != 8 or qc["exact_cached_group_and_primer_total_match"] is not True or catalog["production_barcodes"] != qc["retained_production_barcodes"]:
        raise ValueError("whole physical-library family and matching production QC required")
    records = []
    for shard in support["shards"]:
        if shard["complete"] is not True or "library" not in shard["input"]:
            raise ValueError("complete library-union receipt required")
        filters = shard["filters"]
        summed, unioned, duplicated = [int(filters[name]) for name in ("sum_per_run_barcode_umi_event_keys", "barcode_umi_event_keys", "cross_run_duplicate_event_keys")]
        if not 0 <= duplicated == summed - unioned:
            raise ValueError("duplicate-key accounting must match the exact union")
        records.append(dict(library=shard["input"]["library"], run_summed_event_keys=summed, library_union_event_keys=unioned, duplicated_event_keys=duplicated, duplicate_fraction_of_run_sum=duplicated / summed if summed else 0.))
    table = pd.DataFrame(records)
    if table.library.duplicated().any():
        raise ValueError("eight distinct physical libraries required")
    totals = table[["run_summed_event_keys", "library_union_event_keys", "duplicated_event_keys"]].sum().to_dict()
    totals["duplicate_fraction_of_run_sum"] = totals["duplicated_event_keys"] / totals["run_summed_event_keys"]
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "library_key_duplicates.tsv", sep="\t", index=False)
    for name in ("retained_barcodes.tsv.gz", "diagnostic_barcode_qc.tsv.gz"):
        shutil.copyfile(args.cell_qc / name, args.output_dir / name)
    shutil.copyfile(qc_path, args.output_dir / "cell_qc_manifest.json")
    shutil.copyfile(catalog_path, args.output_dir / "full_catalog_manifest.json")
    inputs = (support_path, qc_path, catalog_path, *(args.cell_qc / name for name in ("retained_barcodes.tsv.gz", "diagnostic_barcode_qc.tsv.gz")))
    manifest = dict(input_hashes={str(path): file_hash(path) for path in inputs}, whole_library_key_totals=totals, scope="completed broader-cell diagnostic molecule-key provenance and production cell-QC reconstruction; not statistical or replication results", caveats="Keys are event/barcode/UMI units, not unique RNAs across events. Broader-cell read counts are not the production-cell model inputs. Exact STAR UB union does not recreate joint UMI error correction. No established comparator result is invalidated by this local diagnostic.", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(totals=totals, cell_qc=qc), indent=2), flush=True)


if __name__ == "__main__":
    main()
