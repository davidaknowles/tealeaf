#!/usr/bin/env python3
"""Split and long-read assessments of run_local_path_tests.py fits restricted by block structure.

Blocks are classified from the annotation alone (tealeaf.sc.local_path_reads.
block_structure on all annotated mature paths of the block): "internal" blocks
have both constitutive anchors, "terminal" blocks have alternative first or
last exons. The fits, p-values and reported usages are unchanged; tests outside
the retained structures are removed from the family before the split (gene
pairs without a retained test have p = 1) and long-read assessments.
"""

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from extra_scripts.assess_paired_inference_audit import fit_table, split_assessment, long_read_assessment
from tealeaf.sc.local_path_reads import block_structure

FILTERS = {"internal": lambda s: s["structure"] == "internal", "distinct_junctions": lambda s: s["structure"] == "internal" or not s["same_junctions"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--filter", choices=sorted(FILTERS), required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    structure = {block["block_id"]: block_structure(block["paths"][:len(block["signatures"])]) for block in blocks.values()}
    keep = {block for block, value in structure.items() if FILTERS[args.filter](value)}
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    name = f"{args.model}, {args.filter} blocks"
    folds = [fit_table(args.cache / f"{args.model}_{fold}/merged/paired_path.tsv") for fold in (0, 1)]
    folds = [fold.loc[fold.block_id.isin(keep)].copy() for fold in folds]
    split_assessment(folds, repo, output, name)
    full = pd.read_csv(args.cache / f"{args.model}_full/merged/paired_path.tsv", sep="\t")
    full = full.loc[full.block_id.isin(keep)].copy()
    long_read_assessment(full, args.cache / f"{args.model}_full", args.matrix_dir, args.gtf, args.block_cache, output, name)
    counts = pd.Series([value["structure"] + ("" if value["structure"] == "internal" else (", same junctions" if value["same_junctions"] else ", distinct junctions")) for value in structure.values()]).value_counts().to_dict()
    (output / "manifest.json").write_text(json.dumps({"model": args.model, "filter": args.filter, "blocks_by_structure": counts, "retained_blocks": len(keep), "selection": "annotation-only block structure; fits and p-values unchanged"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
