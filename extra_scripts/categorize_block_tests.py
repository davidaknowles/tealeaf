#!/usr/bin/env python3
"""Event classes of full-data block tests and their long-read agreement.

Each full-data test (internal blocks from --internal-model, terminal blocks
from --terminal-model, one BH family) is classified from the exon chains of
its retained paths (tealeaf.sc.differential.classify_path_chains); a test that
retains an unannotated path is marked novel. Writes the per-test classes and,
per class: tests, BH discoveries at 0.05, long-read-evaluable discoveries and
their direction agreement, and the class composition of the 100 top-ranked
long-read-evaluable tests with their agreement.
"""

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc.differential import classify_path_chains
from tealeaf.sc.junction_benchmark import benjamini_hochberg
from tealeaf.sc.local_path_reads import block_structure, structure_paths, same_junction_groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--internal-model", required=True)
    parser.add_argument("--terminal-model", required=True)
    parser.add_argument("--assessment", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    blocks = {block["block_id"]: block for block in json.load(gzip.open(args.blocks, "rt")).values()}
    structure = {key: block_structure(structure_paths(block))["structure"] for key, block in blocks.items()}
    tests = []
    for model, kind in ((args.internal_model, "internal"), (args.terminal_model, "terminal")):
        table = pd.read_csv(args.cache / f"{model}_full/merged/paired_path.tsv", sep="\t")
        tests.append(table.loc[table.block_id.map(structure).eq(kind)].assign(structure=kind))
    tests = pd.concat(tests, ignore_index=True)
    tests["fdr"] = benjamini_hochberg(tests.p_value.to_numpy())
    events, components, novel = [], [], []
    for row in tests.itertuples():
        block = blocks[row.block_id]
        n_paths = len(block["signatures"])
        kept = json.loads(row.kept_paths) if isinstance(row.kept_paths, str) else []
        groups = same_junction_groups(block["paths"][:n_paths]) if row.structure == "terminal" else [[index] for index in range(n_paths)]
        members = [groups[position] for position in kept]
        if len(members) < 2:
            events.append("untested"), components.append(""), novel.append(False)
            continue
        event, parts = classify_path_chains([block["paths"][member[0]] for member in members], block["strand"])
        events.append(event), components.append(";".join(parts)), novel.append(any(index >= block.get("n_annotated", n_paths) for member in members for index in member))
    tests["event"], tests["components"], tests["novel_path"] = events, components, novel
    mapping = pd.read_csv(args.assessment / "lr_mapping.tsv.gz", sep="\t", usecols=["test_id", "eligible", "pooled_replicated"])
    tests = tests.merge(mapping, on="test_id", how="left")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tests[["test_id", "block_id", "gene_id", "level_a", "level_b", "structure", "event", "components", "novel_path", "n_paths", "p_value", "fdr", "eligible", "pooled_replicated"]].to_csv(args.output_dir / "test_events.tsv.gz", sep="\t", index=False)
    evaluable = tests.loc[tests.eligible.fillna(False).astype(bool)].sort_values(["p_value", "statistic"], ascending=[True, False])
    top = evaluable.head(100)
    rows = []
    for event, local in tests.loc[tests.event.ne("untested")].groupby("event"):
        discoveries = local.loc[local.fdr.le(.05)]
        checked = discoveries.loc[discoveries.eligible.fillna(False).astype(bool)]
        top_local = top.loc[top.event.eq(event)]
        rows.append({"event": event, "tests": len(local), "discoveries": len(discoveries), "discoveries_with_novel_path": int(discoveries.novel_path.sum()), "long_read_evaluable_discoveries": len(checked), "long_read_agreement": checked.pooled_replicated.astype(float).mean() if len(checked) else float("nan"), "top100": len(top_local), "top100_agreeing": int(top_local.pooled_replicated.astype(float).sum())})
    summary = pd.DataFrame(rows).sort_values("discoveries", ascending=False)
    summary.to_csv(args.output_dir / "event_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
