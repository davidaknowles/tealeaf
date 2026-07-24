#!/usr/bin/env python3
"""Build primer-specific GLM designs from Salmon positional-bias weights."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from tealeaf.data.alevin import load_alevin_counts, load_alevin_structure
from tealeaf.data.salmon import (
    build_positional_ec_design,
    salmon_eqclass_paths,
    summarize_positional_bias_models,
)
from tealeaf.sc.sc_utils import get_feature_weights, get_transcript_lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument(
        "--polydt-quant", action="append", required=True, type=Path
    )
    parser.add_argument(
        "--ranhex-quant", action="append", required=True, type=Path
    )
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--output-report", required=True, type=Path)
    args = parser.parse_args()

    features, membership = load_alevin_structure(args.alevin_dir)
    barcodes, counts = load_alevin_counts(args.alevin_dir)
    ec_counts = {
        "polydt": np.asarray(counts.sum(axis=0)).ravel(),
        "ranhex": np.asarray(counts.sum(axis=0)).ravel(),
    }
    if args.primer_pairs is not None:
        barcode_to_row = {barcode: index for index, barcode in enumerate(barcodes)}
        rows = {"polydt": [], "ranhex": []}
        with open(args.primer_pairs, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for primer, column in (
                    ("polydt", "polydt_barcode"),
                    ("ranhex", "ranhex_barcode"),
                ):
                    index = barcode_to_row.get(row[column])
                    if index is not None:
                        rows[primer].append(index)
        ec_counts = {
            primer: np.asarray(counts[index].sum(axis=0)).ravel()
            for primer, index in rows.items()
        }
    _, fallback = get_feature_weights(
        features, get_transcript_lengths(args.salmon_ref)
    )
    report = {}
    for primer, quant in (
        ("polydt", args.polydt_quant),
        ("ranhex", args.ranhex_quant),
    ):
        eqclasses = salmon_eqclass_paths(quant)
        design, stats, effective_lengths = build_positional_ec_design(
            membership,
            features,
            eqclasses,
            quant,
            ec_counts=ec_counts[primer],
            fallback_inverse_lengths=fallback,
        )
        output = args.alevin_dir / f"gene_eqclass_posbias_{primer}.npz"
        temporary = output.with_name(f".{output.name}.tmp.npz")
        sp.save_npz(temporary, design)
        temporary.replace(output)
        missing_lengths = ~np.isfinite(effective_lengths)
        effective_lengths[missing_lengths] = 1.0 / fallback[missing_lengths]
        length_output = (
            args.alevin_dir / f"salmon_effective_lengths_{primer}.npy"
        )
        temporary_length = length_output.with_name(
            f".{length_output.name}.tmp.npy"
        )
        np.save(temporary_length, effective_lengths)
        temporary_length.replace(length_output)
        stats["output"] = str(output)
        stats["effective_length_output"] = str(length_output)
        stats["positional_models"] = {
            str(path): {
                name: summarize_positional_bias_models(path / "aux_info" / name)
                for name in (
                    "obs5_pos.gz", "obs3_pos.gz",
                    "exp5_pos.gz", "exp3_pos.gz",
                )
            }
            for path in quant
        }
        report[primer] = stats
        print(json.dumps({primer: stats}))
    args.output_report.parent.mkdir(parents=True, exist_ok=True)
    args.output_report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
