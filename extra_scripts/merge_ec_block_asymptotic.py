#!/usr/bin/env python3
"""Merge nonbootstrap EC-block GLMM shards and apply BH correction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.merge_ec_block_glmm import supported_partition_representatives
from extra_scripts.run_differential_splicing import benjamini_hochberg


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--candidate-cache", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    tables = []
    failures = []
    for shard in args.shards:
        result = shard / "ec_block_glmm.tsv"
        if result.is_file() and result.stat().st_size:
            try:
                table = pd.read_csv(result, sep="\t")
            except pd.errors.EmptyDataError:
                table = pd.DataFrame()
            if not table.empty:
                tables.append(table)
        failure_path = shard / "failures.json"
        if failure_path.is_file():
            failures.extend(json.loads(failure_path.read_text()))
    if not tables:
        raise ValueError("no nonempty asymptotic EC-block GLMM shards")
    table = pd.concat(tables, ignore_index=True).drop_duplicates(
        ["test_id", "method"], keep="last"
    )
    if args.candidate_cache is not None:
        representatives = supported_partition_representatives(
            args.candidate_cache
        )
        table = table.loc[
            table["test_id"].map(representatives).eq(table["test_id"])
        ].reset_index(drop=True)

    table["p_value"] = table["lrt_p_value"]
    table["fdr"] = np.nan
    eligible = (
        table["null_converged"]
        & table["alternative_converged"]
        & table["p_value"].notna()
    )
    for _, positions in table.loc[eligible].groupby("method").groups.items():
        table.loc[positions, "fdr"] = benjamini_hochberg(
            table.loc[positions, "p_value"].to_numpy()
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "ec_block_glmm.tsv", sep="\t", index=False)
    table.loc[table["fdr"] <= 0.05].to_csv(
        args.output_dir / "significant_ec_block_glmm.tsv", sep="\t", index=False
    )
    (args.output_dir / "failures.json").write_text(
        json.dumps(failures, indent=2) + "\n"
    )
    summaries = []
    for method, records in table.groupby("method"):
        converged = records["null_converged"] & records["alternative_converged"]
        summaries.append({
            "method": method,
            "events": len(records),
            "convergence": float(converged.mean()),
            "nominal_0.05": int(
                np.sum(converged & (records["p_value"] <= 0.05))
            ),
            "fdr_0.05": int(np.sum(records["fdr"] <= 0.05)),
            "bic_bf_above_10": int(
                np.sum(records["bic_log_bayes_factor"] > np.log(10))
            ),
        })
    pd.DataFrame(summaries).to_csv(
        args.output_dir / "summary.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
