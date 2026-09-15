#!/usr/bin/env python3
"""Audit long-read agreement across fixed Tealeaf path concentrations.

The long-read table is generated for one discovery set. Other concentrations
are therefore evaluated on their intersection with that reference set; this
is an audit of ranking and direction, not a replacement for rerunning source
mapping for every discovery set.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replication", required=True, type=Path)
    parser.add_argument("--zeta-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top-k", type=int, default=75)
    args = parser.parse_args()
    replication = pd.read_csv(args.replication, sep="\t", low_memory=False)
    eligible = replication[replication.mapping_complete.astype(bool) & (replication.minimum_pooled_depth >= 20) & replication.pooled_replicated.notna()].copy()
    rows = []
    for path in sorted(args.zeta_root.glob("a*/significant_paired_path.tsv")):
        zeta = path.parent.name.removeprefix("a")
        significant = pd.read_csv(path, sep="\t", usecols=["test_id", "p_value", "raw_p_value", "statistic"]).rename(columns={"p_value": "test_p_value", "raw_p_value": "test_raw_p_value", "statistic": "test_statistic"})
        joined = eligible.merge(significant, on="test_id")
        joined = joined.sort_values(["test_p_value", "test_raw_p_value", "test_statistic", "test_id"], ascending=[True, True, False, True], kind="stable")
        for scope, local in (("all eligible intersection", joined), (f"top {args.top_k} eligible intersection", joined.head(args.top_k))):
            rows.append({"zeta": float(zeta), "scope": scope, "n_tests": len(local), "n_replicated": int(local.pooled_replicated.sum()), "replication_rate": float(local.pooled_replicated.mean()) if len(local) else float("nan")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
