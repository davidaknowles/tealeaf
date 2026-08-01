#!/usr/bin/env python3
"""Merge clustered compositional model shards and recompute summaries."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from extra_scripts.run_clustered_compositional_models import summarize


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tables(shards, name, **kwargs):
    tables = []
    for shard in shards:
        try:
            tables.append(pd.read_csv(shard / name, sep="\t", **kwargs))
        except pd.errors.EmptyDataError:
            continue
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    observed = read_tables(args.shards, "clustered_model_observed.tsv")
    null = read_tables(args.shards, "clustered_model_null.tsv.gz")
    power = read_tables(args.shards, "clustered_model_power.tsv.gz")
    summaries = [
        json.loads((shard / "validation_summary.json").read_text())
        for shard in args.shards
    ]
    observed, power, summary = summarize(observed, null, power)
    observed.to_csv(
        args.output_dir / "clustered_model_observed.tsv",
        sep="\t",
        index=False,
    )
    null.to_csv(
        args.output_dir / "clustered_model_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    power.to_csv(
        args.output_dir / "clustered_model_power.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    with gzip.open(args.output_dir / "failures.jsonl.gz", "wt") as output:
        for shard in args.shards:
            with gzip.open(shard / "failures.jsonl.gz", "rt") as source:
                for line in source:
                    output.write(line)
    summary.update({
        "seconds_max_shard": max(item["seconds"] for item in summaries),
        "shards": len(args.shards),
        "permutations": summaries[0]["permutations"],
        "power_replicates": summaries[0]["power_replicates"],
        "effect_sizes": summaries[0]["effect_sizes"],
        "failures": sum(item["failures"] for item in summaries),
    })
    (args.output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
