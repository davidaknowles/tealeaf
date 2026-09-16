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
    parser.add_argument("--remapped", type=Path, help="Union source-remapping table from assess_tilgner_long_read_replication.py")
    parser.add_argument("--top-k", type=int, default=75)
    args = parser.parse_args()
    replication = pd.read_csv(args.remapped or args.replication, sep="\t", low_memory=False)
    mapping_complete = replication.mapping_complete.astype(str).str.lower().eq("true")
    eligible = replication[mapping_complete & (replication.minimum_pooled_depth >= 20) & replication.pooled_replicated.notna()].copy()
    rows = []
    detailed = []
    for path in sorted(args.zeta_root.glob("a*/significant_paired_path.tsv")):
        zeta = path.parent.name.removeprefix("a")
        significant = pd.read_csv(path, sep="\t", usecols=["test_id", "p_value", "raw_p_value", "statistic"])
        significant = significant.rename(columns={"p_value": "test_p_value", "raw_p_value": "test_raw_p_value", "statistic": "test_statistic"})
        joined = eligible.merge(significant, on="test_id")
        joined = joined.sort_values(["test_p_value", "test_raw_p_value", "test_statistic", "test_id"], ascending=[True, True, False, True], kind="stable")
        joined.insert(0, "zeta", float(zeta))
        detailed.append(joined)
        for scope, local in (("all eligible intersection", joined), (f"top {args.top_k} eligible intersection", joined.head(args.top_k))):
            if args.remapped is not None:
                scope = scope.replace("intersection", "after full remap")
            rows.append({"zeta": float(zeta), "scope": scope, "n_discoveries": len(significant), "n_tests": len(local), "n_replicated": int(local.pooled_replicated.sum()), "replication_rate": float(local.pooled_replicated.mean()) if len(local) else float("nan")})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    if args.remapped is not None:
        detail_path = args.output.with_name("zeta_full_remap.tsv")
        pd.concat(detailed, ignore_index=True).to_csv(detail_path, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
