#!/usr/bin/env python3
"""Select cross-fitted conditional-uncertainty scales by REML evidence."""

import argparse
import json
from pathlib import Path
import zlib

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--minimum-genes", type=int, default=25)
    args = parser.parse_args()
    paths = [path for path in args.inputs if path.is_file() and path.stat().st_size]
    table = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
    table["gene_fold"] = table.gene.map(
        lambda value: zlib.crc32(str(value).encode()) % args.folds
    )
    per_gene = table.groupby(
        ["gene", "gene_fold", "n_paths", "uncertainty_scale"],
        as_index=False,
    ).agg(restricted_objective=("restricted_objective", "mean"))
    records = []
    scores = []
    for held_out in range(args.folds):
        training = per_gene.loc[per_gene.gene_fold.ne(held_out)]
        for n_paths in sorted(per_gene.n_paths.unique()):
            local = training.loc[training.n_paths.eq(n_paths)]
            if local.gene.nunique() < args.minimum_genes:
                local = training
            local_scores = local.groupby("uncertainty_scale").restricted_objective.mean()
            selected = float(local_scores.idxmin())
            records.append({
                "gene_fold": held_out,
                "n_paths": int(n_paths),
                "uncertainty_scale": selected,
                "training_genes": int(local.gene.nunique()),
            })
            scores.extend({
                "gene_fold": held_out,
                "n_paths": int(n_paths),
                "uncertainty_scale": float(scale),
                "mean_restricted_objective": float(value),
                "selected": float(scale) == selected,
            } for scale, value in local_scores.items())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"folds": args.folds, "records": records}, indent=2) + "\n")
    pd.DataFrame(scores).to_csv(args.scores, sep="\t", index=False)


if __name__ == "__main__":
    main()
