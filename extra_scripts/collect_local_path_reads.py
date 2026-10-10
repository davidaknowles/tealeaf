#!/usr/bin/env python3
"""Collect block-local path-compatibility molecule counts from STARsolo BAMs.

Blocks are the distinct (block, path set) pairs in the given pairwise
candidate caches. Libraries, runs, indexes and production-QC barcode groups
come from the existing local-read recipe. Run one task per library and gene
shard, then --collate.
"""

import argparse
import gzip
import json
import pickle
from pathlib import Path
import zlib

import pandas as pd

from tealeaf.sc.junction_paths import block_path_exons
from tealeaf.sc.local_path_reads import collect_library_path_reads


def anchor(value):
    value = json.loads(value) if isinstance(value, str) else value
    return None if value is None else tuple(value)


def path_key(block_id, signatures):
    return f"{block_id}#{zlib.crc32(json.dumps(signatures).encode()):08x}"


def declared_blocks(candidate_caches, block_cache, max_paths=30):
    """path key -> block description for every candidate path set (2..max_paths paths)."""
    annotation = {row["block_id"]: row for row in json.load(gzip.open(block_cache, "rt"))}
    blocks = {}
    for cache in candidate_caches:
        for candidate in pickle.load(open(cache, "rb"))["candidates"]:
            block_id, gene_id, signatures = candidate[1], candidate[2], candidate[6]
            if not 2 <= len(signatures) <= max_paths:
                continue
            key = path_key(block_id, signatures)
            if key not in blocks:
                block = annotation[block_id]
                blocks[key] = {"block_id": block_id, "gene_id": gene_id, "chromosome": block["chromosome"], "signatures": signatures, "paths": [block_path_exons(anchor(block["left_anchor"]), signature, anchor(block["right_anchor"])) for signature in signatures]}
    return blocks


def annotated_blocks(candidate_caches, block_cache, max_paths=10):
    """Every annotated block (2..max_paths paths) in the candidate caches' genes.

    Paths are all annotated paths of the block, independent of EC support. A
    final precursor path, one exon spanning the block window, marks reads with
    no junction inside the window (intronic and exon-body reads).
    """
    genes = set()
    for cache in candidate_caches:
        genes.update(candidate[2] for candidate in pickle.load(open(cache, "rb"))["candidates"])
    blocks = {}
    for block in json.load(gzip.open(block_cache, "rt")):
        signatures = json.loads(block["path_signatures"]) if isinstance(block["path_signatures"], str) else block["path_signatures"]
        if block["gene_id"] not in genes or not 2 <= len(signatures) <= max_paths:
            continue
        paths = [block_path_exons(anchor(block["left_anchor"]), signature, anchor(block["right_anchor"])) for signature in signatures]
        window = (min(a for path in paths for a, _ in path), max(b for path in paths for _, b in path))
        blocks[path_key(block["block_id"], signatures)] = {"block_id": block["block_id"], "gene_id": block["gene_id"], "chromosome": block["chromosome"], "strand": block["strand"], "signatures": signatures, "paths": paths + [[window]], "precursor": True}
    return blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--candidate-cache", type=Path, action="append", required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--library-index", type=int, default=0)
    parser.add_argument("--gene-shard", type=int, default=0)
    parser.add_argument("--gene-shards", type=int, default=1)
    parser.add_argument("--collate", action="store_true")
    parser.add_argument("--annotated", action="store_true", help="all annotated blocks of the candidate genes, with a precursor path")
    args = parser.parse_args()
    recipe = json.loads(args.recipe.read_text())
    if args.collate:
        expected = len(recipe["bams"]) * args.gene_shards
        parts = sorted(args.output_dir.glob("part_*/counts.tsv.gz"))
        if len(parts) != expected:
            raise ValueError(f"expected {expected} parts, found {len(parts)}")
        table = pd.concat([pd.read_csv(path, sep="\t") for path in parts], ignore_index=True)
        table = table.groupby(["path_key", "subject", "cell_type", "primer", "mask"], as_index=False)["count"].sum()
        table.to_csv(args.output_dir / "counts.tsv.gz", sep="\t", index=False)
        blocks = (annotated_blocks if args.annotated else declared_blocks)(args.candidate_cache, args.block_cache)
        with gzip.open(args.output_dir / "blocks.json.gz", "wt") as handle:
            json.dump(blocks, handle)
        filters = {}
        for path in sorted(args.output_dir.glob("part_*/filters.json")):
            for key, value in json.loads(path.read_text()).items():
                filters[key] = filters.get(key, 0) + value
        (args.output_dir / "filters.json").write_text(json.dumps(filters, indent=2) + "\n")
        print(f"{len(blocks)} blocks, {len(table)} count rows, {int(table['count'].sum())} molecules; {filters}", flush=True)
        return
    blocks = (annotated_blocks if args.annotated else declared_blocks)(args.candidate_cache, args.block_cache)
    genes = sorted({block["gene_id"] for block in blocks.values()})
    selected_genes = set(genes[args.gene_shard::args.gene_shards])
    local = {key: block for key, block in blocks.items() if block["gene_id"] in selected_genes}
    library = recipe["bams"][args.library_index]
    groups = {barcode: tuple(group) for barcode, group in library["barcode_groups"].items()}
    counts, filters = collect_library_path_reads([run["path"] for run in library["runs"]], [run["index"] for run in library["runs"]], groups, local)
    out = args.output_dir / f"part_{args.library_index}_{args.gene_shard}"
    out.mkdir(parents=True, exist_ok=True)
    rows = [{"path_key": key[0], "subject": key[1], "cell_type": key[2], "primer": key[3], "mask": key[4], "count": value} for key, value in counts.items()]
    pd.DataFrame(rows, columns=["path_key", "subject", "cell_type", "primer", "mask", "count"]).to_csv(out / "counts.tsv.gz", sep="\t", index=False)
    (out / "filters.json").write_text(json.dumps({**filters, "blocks": len(local), "genes": len(selected_genes)}, indent=2) + "\n")
    print(f"library {library['library']}, shard {args.gene_shard}: {len(local)} blocks, {sum(counts.values())} molecules, {filters}", flush=True)


if __name__ == "__main__":
    main()
