#!/usr/bin/env python3
"""Score SUPPA/rMATS event directions against Tilgner long-read usage."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, stable_identifier, wilson_interval


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", action="append", required=True, type=Path)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--tilgner-matrix", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--top-per-contrast", type=int, default=500)
    parser.add_argument("--minimum-depth", type=float, default=20.0)
    return parser.parse_args()


def summarize(table, minimum_depth):
    rows = []
    for method, group in table.groupby("method"):
        eligible = group[group.mapping_complete & (group.minimum_pooled_depth >= minimum_depth) & group.pooled_replicated.notna()]
        strict = eligible[eligible.minimum_replicate_depth >= minimum_depth / 2]
        for endpoint, subset, column in (("pooled direction", eligible, "pooled_replicated"), ("both biological replicates", strict, "both_replicates_replicated")):
            values = subset[column].astype(bool)
            successes = int(values.sum())
            low, high = wilson_interval(successes, len(values))
            rows.append({"method": method, "endpoint": endpoint, "scope": "top-ranked calls per contrast", "minimum_depth": minimum_depth, "n_tests": len(values), "n_replicated": successes, "replication_rate": successes / len(values) if len(values) else np.nan, "ci_low": low, "ci_high": high, "conditional_null_rate": 0.5})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    source, features, columns = read_tilgner_matrix(args.tilgner_matrix, args.gtf)
    event_table = pd.read_csv(args.events, sep="\t", compression="infer")
    transcript_rows = defaultdict(list)
    for index, transcript in features.transcript_id.items():
        if pd.notna(transcript):
            transcript_rows[stable_identifier(transcript)].append(index)
    event_rows = []
    event_indices = {}
    for event_index, event in event_table.iterrows():
        included = [row for transcript in str(event.included).split(",") if transcript for row in transcript_rows.get(stable_identifier(transcript), [])]
        excluded = [row for transcript in str(event.excluded).split(",") if transcript for row in transcript_rows.get(stable_identifier(transcript), [])]
        if included and excluded:
            event_indices[str(event.feature_id)] = event_index
            event_rows.extend((event_index, row, 1) for row in included)
            event_rows.extend((event_index, row, 0) for row in excluded)
    inc_map = sparse.csr_matrix((np.ones(sum(kind == 1 for _, _, kind in event_rows)), ([event for event, _, kind in event_rows if kind == 1], [row for _, row, kind in event_rows if kind == 1])), shape=(len(event_table), len(features)))
    exc_map = sparse.csr_matrix((np.ones(sum(kind == 0 for _, _, kind in event_rows)), ([event for event, _, kind in event_rows if kind == 0], [row for _, row, kind in event_rows if kind == 0])), shape=(len(event_table), len(features)))
    event_counts = (inc_map @ source).toarray(), (exc_map @ source).toarray()
    groups = {(cell_type, int(rep)): values["column"].to_numpy(dtype=int) for (cell_type, rep), values in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"])}
    tests = pd.concat((pd.read_csv(path, sep="\t", compression="infer", low_memory=False) for path in args.tests), ignore_index=True)
    tests = tests[tests.effect.eq("cell_type") & tests.p_value.notna()].copy()
    # Keep a stable, method-specific significance ranking without imposing a method-specific FDR cutoff.
    tests["p_value"] = pd.to_numeric(tests.p_value, errors="coerce")
    tests = tests.sort_values(["method", "contrast_id", "p_value", "feature_id"], kind="stable").groupby(["method", "contrast_id"], sort=False).head(args.top_per_contrast)
    rows = []
    for record in tests.itertuples(index=False):
        feature_id = str(record.feature_id)
        if feature_id.startswith("rMATS:"):
            event_key = "SUPPA:" + feature_id.split(":", 1)[1]
            short_sign = -1.0
        else:
            event_key = feature_id
            short_sign = 1.0
        event_index = event_indices.get(event_key)
        levels = (record.level_a, record.level_b)
        if event_index is None or any((level, replicate) not in groups for level in levels for replicate in (1, 2)):
            continue
        # Pool the two source replicates by summing their inclusion and exclusion counts.
        a_inc = event_counts[0][event_index, groups[(levels[0], 1)]].sum() + event_counts[0][event_index, groups[(levels[0], 2)]].sum()
        a_exc = event_counts[1][event_index, groups[(levels[0], 1)]].sum() + event_counts[1][event_index, groups[(levels[0], 2)]].sum()
        b_inc = event_counts[0][event_index, groups[(levels[1], 1)]].sum() + event_counts[0][event_index, groups[(levels[1], 2)]].sum()
        b_exc = event_counts[1][event_index, groups[(levels[1], 1)]].sum() + event_counts[1][event_index, groups[(levels[1], 2)]].sum()
        delta = b_inc / (b_inc + b_exc) - a_inc / (a_inc + a_exc) if (a_inc + a_exc) and (b_inc + b_exc) else np.nan
        replicate_signs = []
        for replicate in (1, 2):
            ai = event_counts[0][event_index, groups[(levels[0], replicate)]].sum()
            ae = event_counts[1][event_index, groups[(levels[0], replicate)]].sum()
            bi = event_counts[0][event_index, groups[(levels[1], replicate)]].sum()
            be = event_counts[1][event_index, groups[(levels[1], replicate)]].sum()
            replicate_signs.append((bi / (bi + be) - ai / (ai + ae)) if (ai + ae) and (bi + be) else np.nan)
        rows.append({"method": record.method, "contrast_id": record.contrast_id, "feature_id": feature_id, "event_type": record.event_type, "gene_id": getattr(record, "gene_id", ""), "level_a": levels[0], "level_b": levels[1], "p_value": record.p_value, "raw_p_value": record.p_value, "statistic": -np.log10(max(record.p_value, 1e-300)), "mapping_complete": True, "minimum_pooled_depth": float(min(a_inc + a_exc, b_inc + b_exc)), "minimum_replicate_depth": float(min(event_counts[0][event_index, groups[(level, replicate)]].sum() + event_counts[1][event_index, groups[(level, replicate)]].sum() for level in levels for replicate in (1, 2))), "pooled_replicated": bool(short_sign * delta > 0) if np.isfinite(delta) else np.nan, "replicate_1_dot_product": short_sign * replicate_signs[0], "replicate_2_dot_product": short_sign * replicate_signs[1], "both_replicates_replicated": bool(short_sign * replicate_signs[0] > 0 and short_sign * replicate_signs[1] > 0) if np.isfinite(replicate_signs).all() else np.nan})
    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")
    summary = summarize(result, args.minimum_depth)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
