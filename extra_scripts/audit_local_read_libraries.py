#!/usr/bin/env python3
"""Resolve frozen alignment barcodes to physical libraries before molecule tests."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import h5py
import pandas as pd

from extra_scripts.export_microglia_primer_pairs import decode_column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--library-runs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new provenance output required")
    recipe_path = args.source_root / "recipe.json"
    recipe = json.loads(recipe_path.read_text())
    run_map = pd.read_csv(args.library_runs, sep="\t", header=None, names=["library", "run_a", "run_b"], dtype=str)
    if run_map.isna().any().any() or run_map.library.duplicated().any():
        raise ValueError("unique nonempty physical library/run mapping required")
    runs = [value for row in run_map.itertuples(index=False) for value in (row.run_a, row.run_b)]
    if len(set(runs)) != len(runs) or set(runs) != {Path(row["path"]).parent.name for row in recipe["bams"]}:
        raise ValueError("each declared alignment run must belong to one library")
    with h5py.File(args.h5ad) as handle:
        obs = handle["obs"]
        poly, hexamer, libraries = [decode_column(obs, field) for field in ("CB_polydT", "CB_ranhex", "sublibrary")]
    if len(poly) != len(hexamer) or len(poly) != len(libraries):
        raise ValueError("aligned source barcode/library observations required")
    origins, occurrences = defaultdict(set), defaultdict(int)
    for first, second, library in zip(poly, hexamer, libraries):
        for barcode in (first, second):
            if barcode and pd.notna(barcode) and library and pd.notna(library):
                origins[str(barcode)].add(str(library))
                occurrences[str(barcode)] += 1
    rows = []
    for barcode, (subject, cell_type, primer) in sorted(recipe["barcode_groups"].items()):
        candidates = origins[barcode]
        library = next(iter(candidates)) if len(candidates) == 1 else ""
        status = "unique" if len(candidates) == 1 and occurrences[barcode] == 1 and library in set(run_map.library) else "missing" if not candidates else "ambiguous_or_reused"
        rows.append(dict(barcode=barcode, subject=subject, cell_type=cell_type, primer=primer, library=library, status=status, source_occurrences=occurrences[barcode], n_libraries=len(candidates)))
    table = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "barcode_library.tsv.gz", sep="\t", index=False)
    run_map.to_csv(args.output_dir / "library_runs.tsv", sep="\t", index=False)
    manifest = dict(source_recipe=str(recipe_path), source_recipe_sha256=hashlib.sha256(recipe_path.read_bytes()).hexdigest(), requested_barcodes=len(table), source_h5ad=dict(path=str(args.h5ad), size=args.h5ad.stat().st_size, mtime_ns=args.h5ad.stat().st_mtime_ns), run_map_sha256=hashlib.sha256(args.library_runs.read_bytes()).hexdigest(), counts=table.status.value_counts().to_dict(), libraries=run_map.library.tolist(), complete_unique_library_assignment=bool(table.status.eq("unique").all()), scope="physical library assignment only; no count unions, modified cell annotations, statistical tests or observed-power claim", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
