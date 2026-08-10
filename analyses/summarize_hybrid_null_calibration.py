#!/usr/bin/env python3
"""Summarize structured-null families from the hybrid EC-GLMM scheduler."""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import pandas as pd


LEVELS = (0.05, 0.01, 0.001)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def summarize_family(path):
    table = pd.read_csv(path, sep="\t")
    converged = table.loc[table.converged].copy()
    summary = {
        "family": Path(path).parent.name,
        "tests": len(table),
        "converged": len(converged),
        "bh_discoveries": int((table.bh_q_value <= 0.05).sum()),
    }
    for level in LEVELS:
        summary[f"p_below_{level:g}"] = float(
            (converged.lrt_p_value <= level).mean()
        )
    return table, summary


def main():
    args = parse_args()
    paths = sorted(glob.glob(args.results))
    if not paths:
        raise ValueError(f"no results matched {args.results}")
    tables = []
    families = []
    for path in paths:
        table, family = summarize_family(path)
        tables.append(table)
        families.append(family)
    family_table = pd.DataFrame(families)
    converged = pd.concat(tables, ignore_index=True).query("converged")
    summary = {
        "families": len(family_table),
        "tests": int(family_table.tests.sum()),
        "converged": int(family_table.converged.sum()),
        "mean_bh_discoveries": float(family_table.bh_discoveries.mean()),
        "maximum_bh_discoveries": int(family_table.bh_discoveries.max()),
        "families_with_bh_discoveries": int(
            (family_table.bh_discoveries > 0).sum()
        ),
    }
    for level in LEVELS:
        summary[f"pooled_p_below_{level:g}"] = float(
            (converged.lrt_p_value <= level).mean()
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    family_table.to_csv(
        args.output_dir / "family_summary.tsv",
        sep="\t",
        index=False,
    )
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
