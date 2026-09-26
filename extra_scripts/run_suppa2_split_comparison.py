#!/usr/bin/env python3
"""Regenerate native SUPPA2 event tests for a subject-fold manifest.

The native SUPPA2 event catalogue is supplied by ``event-catalog`` so the
rerun does not depend on retaining the external SUPPA2 installation. PSI is
the native transcript-TPM ratio over the catalogue's included and excluded
transcript sets, followed by the same paired Wilcoxon normal approximation
used by ``run_suppa2_comparison.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from run_suppa2_comparison import bh, fast_paired_wilcoxon


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--columns", required=True, type=Path)
    parser.add_argument("--event-catalog", required=True, type=Path)
    parser.add_argument("--contrasts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--minimum-pairs", type=int, default=8)
    parser.add_argument("--fold", type=int, required=True)
    return parser.parse_args()


def load_catalog(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=str).fillna("")


def event_psi(catalog: pd.DataFrame, matrix: sparse.spmatrix, columns: list[str]):
    transcript_index = {name: index for index, name in enumerate(columns)}
    included_rows, included_cols = [], []
    excluded_rows, excluded_cols = [], []
    for event_index, record in enumerate(catalog.itertuples(index=False)):
        for transcript in record.included.split(","):
            if transcript in transcript_index:
                included_rows.append(event_index)
                included_cols.append(transcript_index[transcript])
        for transcript in record.excluded.split(","):
            if transcript in transcript_index:
                excluded_rows.append(event_index)
                excluded_cols.append(transcript_index[transcript])
    shape = (len(catalog), len(columns))
    included = sparse.csr_matrix(
        (np.ones(len(included_rows)), (included_rows, included_cols)), shape=shape
    )
    excluded = sparse.csr_matrix(
        (np.ones(len(excluded_rows)), (excluded_rows, excluded_cols)), shape=shape
    )
    values = matrix.T.tocsr()
    included_values = included @ values
    excluded_values = excluded @ values
    denominator = included_values + excluded_values
    denominator = denominator.toarray()
    numerator = included_values.toarray()
    return np.divide(
        numerator,
        denominator,
        out=np.full(numerator.shape, np.nan, dtype=float),
        where=denominator > 0,
    )


def main():
    args = parse_args()
    matrix = sparse.load_npz(args.matrix).tocsr()
    rows = args.rows.read_text().splitlines()
    columns = args.columns.read_text().splitlines()
    if matrix.shape != (len(rows), len(columns)):
        raise ValueError(f"matrix shape {matrix.shape} does not match rows/columns")
    row_lookup = {
        row.split("__")[-1] + "__" + row.split("__")[0]: index
        for index, row in enumerate(rows)
    }
    catalog = load_catalog(args.event_catalog)
    psi = event_psi(catalog, matrix, columns)
    contrasts = json.loads(args.contrasts.read_text())
    output_rows = []
    for contrast in contrasts:
        level_a, level_b = contrast["level_a"], contrast["level_b"]
        pairs = [
            (row_lookup.get(a), row_lookup.get(b))
            for a, b in zip(contrast["samples_a"], contrast["samples_b"])
        ]
        pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
        if len(pairs) < args.minimum_pairs:
            continue
        first = psi[:, [a for a, _ in pairs]]
        second = psi[:, [b for _, b in pairs]]
        valid = np.isfinite(first) & np.isfinite(second)
        enough = valid.sum(axis=1) >= args.minimum_pairs
        differences = second - first
        counts = valid.sum(axis=1)
        effects = np.divide(
            np.nansum(np.where(valid, differences, np.nan), axis=1),
            counts,
            out=np.full(len(catalog), np.nan),
            where=counts > 0,
        )
        p_values = fast_paired_wilcoxon(differences, valid)
        q_values = bh(p_values)
        for event_index in np.flatnonzero(enough):
            event = catalog.iloc[event_index]
            output_rows.append(
                {
                    "method": "SUPPA2 native PSI",
                    "contrast_id": contrast["contrast_id"],
                    "effect": "cell_type",
                    "stratum": contrast.get("stratum", "all"),
                    "level_a": level_a,
                    "level_b": level_b,
                    "feature_id": event.feature_id,
                    "event_type": event.event_type,
                    "event_id": event.event_id,
                    "p_value": float(p_values[event_index]),
                    "q_value": float(q_values[event_index]),
                    "effect_size": float(effects[event_index]),
                    "n_subjects": int(counts[event_index]),
                    "fold": args.fold,
                    "gene_id": event.gene_id,
                    "gene_name": event.gene_name,
                    "significant": bool(q_values[event_index] < 0.05),
                    "criterion": "SUPPA2 classical paired Wilcoxon, BH q < 0.05",
                }
            )
    result = pd.DataFrame(output_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(result):,} native SUPPA2 split tests from {len(catalog):,} events")


if __name__ == "__main__":
    main()
