#!/usr/bin/env python3
"""Tidy tables for example block tests: Parse per-subject path usage and long-read proportions.

Examples are the full-data tests with the largest Parse effect (norm of the
path-usage difference) among long-read-evaluable tests with p below --max-p,
at least --min-depth long-read UMIs per cell type and at most --max-paths
paths, one test per gene, with at least one internal and one terminal block.
Selection does not use long-read agreement. For each example the test's
projected classes are rebuilt (prepare_block_test, the fitting code path) and
every subject/cell-type row gets its molecule-share composition
(row_path_compositions); long-read proportions are per source replicate.
Writes examples.tsv, parse_usage.tsv, long_read_usage.tsv and exons.tsv.
"""

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.differential import classify_path_chains
from tealeaf.sc.local_path_reads import block_structure, same_junction_groups
from tealeaf.sc.local_path_test import local_class_masks, prepare_block_test, row_path_compositions


def gene_names(gtf):
    names = {}
    pattern = re.compile(r'gene_id "([^"]+)".*gene_name "([^"]+)"')
    with open(gtf) as handle:
        for line in handle:
            if "\tgene\t" in line:
                match = pattern.search(line)
                if match:
                    names[match.group(1).split(".")[0]] = match.group(2)
    return names


def compressed(exons, max_intron=300):
    """Plot coordinates with introns longer than max_intron shortened to max_intron."""
    points = sorted({x for a, b in exons for x in (a, b)})
    position, mapping = 0., {points[0]: 0.}
    covered = lambda x: any(a <= x < b for a, b in exons)
    for left, right in zip(points, points[1:]):
        width = right - left
        position += width if covered(left) or width <= max_intron else max_intron
        mapping[right] = position
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", type=Path, required=True, help="assess_block_structure_filter.py output directory (lr_mapping.tsv.gz)")
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--internal-model", required=True)
    parser.add_argument("--terminal-model", required=True)
    parser.add_argument("--internal-reads", type=Path, required=True)
    parser.add_argument("--terminal-reads", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--examples", type=int, default=4)
    parser.add_argument("--max-p", type=float, default=1e-8)
    parser.add_argument("--min-depth", type=float, default=50)
    parser.add_argument("--max-paths", type=int, default=4)
    parser.add_argument("--min-share", type=float, default=.01)
    args = parser.parse_args()
    blocks = json.load(gzip.open(args.internal_reads / "blocks.json.gz", "rt"))
    by_id = {block["block_id"]: key for key, block in blocks.items()}
    structure = {block["block_id"]: block_structure(block["paths"][:block.get("n_annotated", len(block["signatures"]))])["structure"] for block in blocks.values()}
    fits = pd.concat([pd.read_csv(args.cache / f"{model}_full/merged/paired_path.tsv", sep="\t").assign(structure=kind) for model, kind in ((args.internal_model, "internal"), (args.terminal_model, "terminal"))], ignore_index=True)
    fits = fits.loc[fits.block_id.map(structure).eq(fits.structure) & fits.converged.astype(bool)]
    mapping = pd.read_csv(args.assessment / "lr_mapping.tsv.gz", sep="\t")
    mapping = mapping.loc[mapping.eligible.astype(bool), ["test_id", "minimum_pooled_depth", "pooled_replicated", "counts_a_rep1", "counts_a_rep2", "counts_b_rep1", "counts_b_rep2"]]
    candidates = fits.merge(mapping, on="test_id")
    candidates = candidates.loc[candidates.p_value.lt(args.max_p) & candidates.minimum_pooled_depth.ge(args.min_depth) & candidates.n_paths.le(args.max_paths)]
    candidates = candidates.sort_values("mean_difference_norm", ascending=False).drop_duplicates("gene_id")
    chosen = list(candidates.head(args.examples).itertuples())
    for kind in ("internal", "terminal"):
        if not any(row.structure == kind for row in chosen):
            extra = candidates.loc[candidates.structure.eq(kind)].head(1)
            chosen = chosen[:-1] + list(extra.itertuples())
    names = gene_names(args.gtf)
    examples, parse_rows, long_rows, exon_rows = [], [], [], []
    for row in chosen:
        key = by_id[row.block_id]
        block = blocks[key]
        reads = args.internal_reads if row.structure == "internal" else args.terminal_reads
        levels = (row.level_a, row.level_b)
        counts = pd.concat([chunk.loc[chunk.path_key.eq(key) & chunk.cell_type.isin(levels)] for chunk in pd.read_csv(reads / "counts.tsv.gz", sep="\t", chunksize=2_000_000)])
        entries = [(r.subject, r.cell_type, r.primer, int(r.mask), float(r.count)) for r in counts.itertuples(index=False)]
        stored = json.load(gzip.open(reads / "opportunities.json.gz", "rt"))[key]
        opportunities = {int(mask): np.asarray(vector, dtype=float) for mask, vector in stored.items()}
        n_paths = len(block["signatures"])
        groups = same_junction_groups(block["paths"][:n_paths]) if row.structure == "terminal" else None
        prepared = prepare_block_test(key, entries, opportunities, n_paths, levels, args.min_share, False, groups)
        kept = prepared["kept"]
        masks = local_class_masks(prepared["opportunities"], len(kept), True)
        compositions = row_path_compositions(prepared["lookup"], key, prepared["opportunities"], prepared["subjects"], prepared["labels"], levels, masks, len(kept))
        gene = row.gene_id.split(".")[0]
        label = f"{names.get(gene, gene)} {row.block_id.split(':')[-1]}"
        examples.append({"test_id": row.test_id, "label": label, "gene_id": gene, "gene_name": names.get(gene, gene), "block_id": row.block_id, "structure": row.structure, "level_a": row.level_a, "level_b": row.level_b, "n_paths": len(kept), "p_value": row.p_value, "effect_norm": row.mean_difference_norm, "long_read_agrees": bool(row.pooled_replicated), "minimum_long_read_depth": row.minimum_pooled_depth, "chromosome": block["chromosome"], "strand": block["strand"], "novel_paths": int(sum(index >= block.get("n_annotated", n_paths) for position in kept for index in (prepared["groups"][position] if prepared["groups"] is not None else [position])))})
        for (subject, level), (composition, molecules) in compositions.items():
            for path, value in enumerate(composition):
                parse_rows.append({"test_id": row.test_id, "subject": subject, "cell_type": levels[level], "path": path + 1, "proportion": value, "molecules": molecules})
        for level, prefix in ((0, "counts_a"), (1, "counts_b")):
            for replicate in (1, 2):
                vector = np.asarray(json.loads(getattr(row, f"{prefix}_rep{replicate}")), dtype=float)
                for path, value in enumerate(vector / vector.sum() if vector.sum() > 0 else vector):
                    long_rows.append({"test_id": row.test_id, "replicate": replicate, "cell_type": levels[level], "path": path + 1, "proportion": value, "umis": vector.sum()})
        chains = [block["paths"][(prepared["groups"][position] if prepared["groups"] is not None else [position])[0]] for position in kept]
        examples[-1]["event"] = classify_path_chains(chains, block["strand"])[0]
        coordinates = compressed([tuple(exon) for chain in chains for exon in chain])
        for path, chain in enumerate(chains):
            for start, end in chain:
                exon_rows.append({"test_id": row.test_id, "path": path + 1, "start": start, "end": end, "x0": coordinates[start], "x1": coordinates[end]})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, rows in (("examples", examples), ("parse_usage", parse_rows), ("long_read_usage", long_rows), ("exons", exon_rows)):
        pd.DataFrame(rows).to_csv(args.output_dir / f"{name}.tsv", sep="\t", index=False)
    print(pd.DataFrame(examples)[["label", "structure", "level_a", "level_b", "n_paths", "p_value", "effect_norm", "long_read_agrees", "novel_paths", "event"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
