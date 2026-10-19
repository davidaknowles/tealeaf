#!/usr/bin/env python3
"""Cell set for native collection from the dataset's cell annotation.

Cells are the annotated cells (combined AnnData, one record per cell with its
poly(dT) and random-hexamer half-cell barcodes, cell type and sublibrary)
whose barcode pair belongs to one record only; a pair used by cells in several
sublibraries cannot be separated by barcode and is left out. Subject and cell
type come from the junction benchmark's barcode labels (the cell set of the
comparator methods), primer from which half the barcode is. Both halves of a
cell are kept, with no count threshold. The output recipe copies --recipe
(sequencing runs per sublibrary) with barcode_groups replaced.
"""

import argparse
import collections
import json
from pathlib import Path

import h5py
import pandas as pd

from extra_scripts.export_microglia_primer_pairs import decode_column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True, help="barcode,cell_type__condition__subject (junction benchmark barcode_groups.csv)")
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with h5py.File(args.h5ad) as handle:
        obs = handle["obs"]
        polydt, ranhex, sublibrary = (decode_column(obs, name) for name in ("CB_polydT", "CB_ranhex", "sublibrary"))
    uses = collections.Counter(zip(polydt, ranhex))
    labels = pd.read_csv(args.labels, header=None, names=["barcode", "group"]).set_index("barcode").group
    recipe = json.loads(args.recipe.read_text())
    libraries = {library["library"]: library for library in recipe["bams"]}
    groups = {name: {} for name in libraries}
    counts = collections.Counter()
    for poly, hexa, library in zip(polydt, ranhex, sublibrary):
        if uses[(poly, hexa)] != 1:
            counts["reused pair"] += 1
            continue
        if poly not in labels.index or hexa not in labels.index:
            counts["unlabelled"] += 1
            continue
        if labels[poly] != labels[hexa]:
            counts["inconsistent labels"] += 1
            continue
        library = str(library)
        if library not in groups:
            counts["unknown sublibrary"] += 1
            continue
        cell_type, _, subject = labels[poly].split("__")
        groups[library][poly] = [subject, cell_type, "poly(dT)"]
        groups[library][hexa] = [subject, cell_type, "random hexamer"]
        counts["cells"] += 1
    for name, library in libraries.items():
        library["barcode_groups"] = groups[name]
    recipe["cell_set"] = {"source": str(args.h5ad), "labels": str(args.labels), "rule": "annotated cells with a barcode pair used by one cell; both halves; no count threshold", **counts}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(recipe))
    print(dict(counts), {name: len(value) for name, value in groups.items()}, flush=True)


if __name__ == "__main__":
    main()
