#!/usr/bin/env python3
"""Build shared pseudobulk junction inputs from STARsolo SJ outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tealeaf.sc.junction_benchmark import (
    aggregate_starsolo_runs,
    annotate_scquint_groups,
    prepare_splitseq_whitelists,
    write_leafcutter_junctions,
)


def parse_run(value: str):
    run, separator, path = value.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("runs must use RUN=PATH")
    return run, Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    whitelist = subparsers.add_parser("whitelists")
    whitelist.add_argument("--barcodes", type=Path, required=True)
    whitelist.add_argument("--output-dir", type=Path, required=True)

    aggregate = subparsers.add_parser("aggregate")
    aggregate.add_argument("--run", action="append", type=parse_run, required=True)
    aggregate.add_argument("--cell-metadata", type=Path, required=True)
    aggregate.add_argument("--output-prefix", type=Path, required=True)
    aggregate.add_argument("--leafcutter-dir", type=Path)
    aggregate.add_argument("--gtf", type=Path)
    args = parser.parse_args()

    if args.command == "whitelists":
        with open(args.barcodes) as handle:
            paths = prepare_splitseq_whitelists(list(handle), args.output_dir)
        print(" ".join(map(str, paths)))
        return

    metadata = pd.read_csv(args.cell_metadata, sep="\t", dtype=str)
    bundle = aggregate_starsolo_runs(dict(args.run), metadata)
    if args.gtf is not None:
        bundle = annotate_scquint_groups(bundle, args.gtf)
    bundle.save(args.output_prefix)
    if args.leafcutter_dir is not None:
        write_leafcutter_junctions(bundle, args.leafcutter_dir)
    print(
        f"samples={len(bundle.samples)} junctions={len(bundle.junctions)} "
        f"umis={int(bundle.counts.sum())}"
    )


if __name__ == "__main__":
    main()
