#!/usr/bin/env python3
"""Audit SUPPA2 paired event statistics with synchronized sign flips."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.run_suppa2_full_data_comparison import (
    TEST_METHODS,
    bh,
    event_test_differences,
    hybrid_exact_paired_wilcoxon,
)
from extra_scripts.run_suppa2_split_data_comparison import (
    event_psi,
    load_catalog,
)
from tealeaf.sc import differential


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--columns", required=True, type=Path)
    parser.add_argument("--event-catalog", required=True, type=Path)
    parser.add_argument("--contrasts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--test-method", action="append", choices=TEST_METHODS, required=True)
    parser.add_argument("--minimum-pairs", type=int, default=8)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--seed", type=int, default=290927)
    return parser.parse_args()


def synchronized_sign(seed, replicate, subject):
    token = f"{seed}|{replicate}|{subject}".encode()
    return 1.0 if hashlib.sha256(token).digest()[0] & 1 else -1.0


def summarize(p_values, enough, method, contrast_id, replicate):
    values = p_values[enough & np.isfinite(p_values)]
    q_values = bh(p_values)
    return {
        "method": method,
        "contrast_id": contrast_id,
        "replicate": replicate,
        "tests": len(values),
        "reject_0_05": int(np.sum(values < 0.05)),
        "reject_0_01": int(np.sum(values < 0.01)),
        "reject_0_001": int(np.sum(values < 0.001)),
        "bh_0_05": int(np.sum(q_values < 0.05)),
        "minimum_p_value": float(np.min(values)) if len(values) else np.nan,
    }


def main():
    args = parse_args()
    matrix = sparse.load_npz(args.matrix).tocsr()
    rows = args.rows.read_text().splitlines()
    columns = args.columns.read_text().splitlines()
    row_lookup = {
        row.split("__")[-1] + "__" + row.split("__")[0]: index
        for index, row in enumerate(rows)
    }
    catalog = load_catalog(args.event_catalog)
    psi = event_psi(catalog, matrix, columns)
    contrasts = json.loads(args.contrasts.read_text())
    contrasts = [
        contrast
        for index, contrast in enumerate(contrasts)
        if index % args.shard_count == args.shard_index
    ]
    summaries = []
    for contrast in contrasts:
        pairs = [
            (row_lookup.get(a), row_lookup.get(b), str(a).split("__", 1)[0])
            for a, b in zip(contrast["samples_a"], contrast["samples_b"])
        ]
        pairs = [(a, b, subject) for a, b, subject in pairs if a is not None and b is not None]
        if len(pairs) < args.minimum_pairs:
            continue
        first = psi[:, [a for a, _, _ in pairs]]
        second = psi[:, [b for _, b, _ in pairs]]
        valid = np.isfinite(first) & np.isfinite(second)
        enough = valid.sum(axis=1) >= args.minimum_pairs
        subjects = [subject for _, _, subject in pairs]
        for method in args.test_method:
            if method == "wilcoxon":
                continue
            differences = event_test_differences(first, second, method)
            p_values = (
                hybrid_exact_paired_wilcoxon(differences, valid)
                if method == "wilcoxon_exact"
                else differential.vectorized_paired_t_pvalues(
                    differences,
                    valid,
                    moderate=method.endswith("moderated_t"),
                )
            )
            p_values[~enough] = np.nan
            summaries.append(summarize(p_values, enough, method, contrast["contrast_id"], -1))
            for replicate in range(args.null_replicates):
                signs = np.asarray([
                    synchronized_sign(args.seed, replicate, subject)
                    for subject in subjects
                ])
                null_p = (
                    hybrid_exact_paired_wilcoxon(
                        differences * signs[None, :], valid
                    )
                    if method == "wilcoxon_exact"
                    else differential.vectorized_paired_t_pvalues(
                        differences * signs[None, :],
                        valid,
                        moderate=method.endswith("moderated_t"),
                    )
                )
                null_p[~enough] = np.nan
                summaries.append(summarize(null_p, enough, method, contrast["contrast_id"], replicate))
    result = pd.DataFrame(summaries)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(result):,} observed and sign-flip summaries")


if __name__ == "__main__":
    main()
