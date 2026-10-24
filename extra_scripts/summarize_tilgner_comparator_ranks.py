#!/usr/bin/env python3
"""Long-read top-rank agreement of comparator methods from their all-tested source mappings.

Each input table replaces the whole curve of its methods (as
``plot_tilgner_method_replication.py --rank-replacement``): rows with complete
source mapping, minimum pooled source depth 20, an evaluable direction and a
finite p-value are ranked by p-value with raw-p and statistic tie-breaks.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.replication_audit import rank_direction_table, ranked_direction_summary

EVENT_METHODS = {"rMATS": "rMATS paired JCEC", "Tealeaf EC; SUPPA2 event definitions": "Tealeaf/SUPPA2 hybrid"}


def read_source_table(path):
    table = pd.read_csv(path, sep="\t", low_memory=False)
    table["method"] = table["method"].replace(EVENT_METHODS)
    if "raw_p_value" not in table:
        table["raw_p_value"] = table["p_value"]
    if "statistic" not in table:
        table["statistic"] = np.nan
    return table


def evaluable(table):
    """Rows kept by plot_tilgner_method_replication.replace_rank_methods."""
    valid = table.mapping_complete.astype(str).str.lower().eq("true")
    valid &= pd.to_numeric(table.minimum_pooled_depth, errors="coerce").ge(20)
    valid &= table.pooled_replicated.notna() & np.isfinite(pd.to_numeric(table.p_value, errors="coerce"))
    table = table.loc[valid].copy()
    table["pooled_replicated"] = table.pooled_replicated.astype(str).str.lower().eq("true")
    if table.duplicated(["method", "contrast_id", "feature_id"]).any():
        raise ValueError("duplicate method--contrast--feature tests")
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", action="append", required=True, type=Path, help="All-tested junction or event source-mapping table; repeat.")
    parser.add_argument("--rank-table", type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--max-rank", type=int, default=200)
    args = parser.parse_args()
    tables = [evaluable(read_source_table(path)) for path in args.table]
    counts = pd.concat([table.groupby("method").size() for table in tables])
    if counts.index.duplicated().any():
        raise ValueError("a method appears in more than one table")
    ranked = pd.concat([rank_direction_table(table, args.max_rank) for table in tables], ignore_index=True, sort=False)
    summary = pd.DataFrame(ranked_direction_summary(ranked))
    summary["n_evaluable"] = summary.method.map(counts)
    summary["scope"], summary["ranking"] = "all tested", "continuous raw"
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, sep="\t", index=False, na_rep="NA")
    if args.rank_table:
        ranked.to_csv(args.rank_table, sep="\t", index=False, na_rep="NA")
    print(summary[summary.cutoff.eq(100)].to_string(index=False))


if __name__ == "__main__":
    main()
