#!/usr/bin/env python3
"""Restore full-precision source path counts, validating unchanged published LR effects."""

import argparse
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, load_blocks, block_feature_rows, normalized_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--tilgner-matrix", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    table = pd.read_csv(args.input, sep="\t", low_memory=False)
    with args.candidate_cache.open("rb") as handle:
        candidates = {row[0]: row for row in pickle.load(handle)["candidates"]}
    blocks = load_blocks(args.block_cache)
    matrix, features, columns = read_tilgner_matrix(args.tilgner_matrix, args.gtf)
    values = {(level, int(replicate)): np.asarray(matrix[:, group.column.to_numpy(dtype=int)].sum(axis=1)).ravel() for (level, replicate), group in columns.dropna(subset=["tealeaf_cell_type", "replicate"]).groupby(["tealeaf_cell_type", "replicate"])}
    mapping_cache = {}
    for index, record in table.iterrows():
        if record.test_id not in candidates:
            raise ValueError(f"missing candidate {record.test_id}")
        candidate = candidates[record.test_id]
        signatures = candidate[6]
        if len(signatures) != record.n_paths:
            raise ValueError("source and candidate path dimensions differ")
        key = (record.block_id, json.dumps(signatures))
        if key not in mapping_cache:
            mapped = block_feature_rows(blocks[record.block_id], signatures, features)
            mapping_cache[key] = {path: group.row.to_numpy(dtype=int) for path, group in mapped.groupby("path_number")}
        path_rows = mapping_cache[key]
        counts = {}
        for side, level in (("a", record.level_a), ("b", record.level_b)):
            for replicate in (1, 2):
                vector = np.array([values[(level, replicate)][path_rows.get(path, np.array([], dtype=int))].sum() if (level, replicate) in values else 0. for path in range(1, len(signatures) + 1)])
                counts[(side, replicate)] = vector
                table.at[index, f"counts_{side}_rep{replicate}"] = json.dumps(vector.tolist())
        pooled_a, pooled_b = counts[("a", 1)] + counts[("a", 2)], counts[("b", 1)] + counts[("b", 2)]
        delta = normalized_difference(pooled_a, pooled_b)
        norm = np.linalg.norm(delta)
        depth = min(pooled_a.sum(), pooled_b.sum())
        if not np.isclose(norm, record.long_read_effect_norm, rtol=1e-8, atol=1e-10, equal_nan=True) or not np.isclose(depth, record.minimum_pooled_depth, rtol=1e-8, atol=1e-10):
            raise ValueError(f"source count restoration changes published LR norm/depth for {record.test_id}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    print(f"Restored {len(table)} path-count records at full source precision; published LR norm/depth unchanged", flush=True)


if __name__ == "__main__":
    main()
