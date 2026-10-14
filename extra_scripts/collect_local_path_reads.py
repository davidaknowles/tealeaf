#!/usr/bin/env python3
"""Collect block-local path-compatibility molecule counts from STARsolo BAMs.

Blocks are the distinct (block, path set) pairs in the given pairwise
candidate caches. Libraries, runs, indexes and production-QC barcode groups
come from the existing local-read recipe. Run one task per library and gene
shard, then --collate, which also stores every block's read opportunities
(read length READ_LENGTH) for the tests. With --junction-only, molecules and
opportunities are restricted to reads spanning a variable junction of the
block's mature paths; pass it, and --max-paths, to both steps. With
--novel-paths (and --annotated), each block also gets unannotated paths from
label-blind STAR junction support of the recipe's libraries
(tealeaf.sc.novel_paths); they follow the annotated paths, whose count is the
block's n_annotated. With --gtf as well, novel junctions outside the
annotated blocks get new blocks.
"""

import argparse
import gzip
import json
from multiprocessing import Pool
import os
from pathlib import Path

import pandas as pd

from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.local_blocks import annotated_blocks, block_opportunities, declared_blocks, junction_block_opportunities, star_junction_support
from tealeaf.sc.local_path_reads import collect_library_path_reads


def select_blocks(args, recipe=None):
    if args.novel_paths:
        if not args.annotated:
            raise ValueError("--novel-paths requires --annotated")
        options = {"min_reads": args.novel_min_reads, "min_libraries": args.novel_min_libraries, "max_novel": args.novel_max, "terminal_length": args.novel_terminal_length}
        return annotated_blocks(args.candidate_cache, args.block_cache, args.max_paths or 10, star_junction_support(recipe), options, read_gtf_exons(args.gtf))
    function = annotated_blocks if args.annotated else declared_blocks
    return function(args.candidate_cache, args.block_cache) if args.max_paths is None else function(args.candidate_cache, args.block_cache, args.max_paths)


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
    parser.add_argument("--include-unassigned", action="store_true", help='also use alignments with STARsolo GX "-" (intronic reads)')
    parser.add_argument("--junction-only", action="store_true", help="only molecules with a read spanning a variable junction of the block")
    parser.add_argument("--max-paths", type=int, help="largest annotated path count per block (default 10 with --annotated, 30 otherwise)")
    parser.add_argument("--novel-paths", action="store_true", help="add unannotated paths from STAR junction support (SJ.out.tab next to each run BAM)")
    parser.add_argument("--gtf", type=Path, help="annotation for new blocks around novel junctions outside annotated blocks (with --novel-paths)")
    parser.add_argument("--novel-min-reads", type=int, default=20)
    parser.add_argument("--novel-min-libraries", type=int, default=4)
    parser.add_argument("--novel-max", type=int, default=8, help="most novel paths per block")
    parser.add_argument("--novel-terminal-length", type=int, default=0, help="length of unannotated first or last exons from deep intronic novel splice sites (tealeaf.sc.novel_paths); 0 skips them")
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
        blocks = select_blocks(args, recipe)
        with gzip.open(args.output_dir / "blocks.json.gz", "wt") as handle:
            json.dump(blocks, handle)
        with Pool(len(os.sched_getaffinity(0))) as pool:
            opportunities = dict(zip(blocks, pool.map(junction_block_opportunities if args.junction_only else block_opportunities, list(blocks.values()), chunksize=8)))
        with gzip.open(args.output_dir / "opportunities.json.gz", "wt") as handle:
            json.dump(opportunities, handle)
        filters = {}
        for path in sorted(args.output_dir.glob("part_*/filters.json")):
            for key, value in json.loads(path.read_text()).items():
                filters[key] = filters.get(key, 0) + value
        (args.output_dir / "filters.json").write_text(json.dumps(filters, indent=2) + "\n")
        print(f"{len(blocks)} blocks, {len(table)} count rows, {int(table['count'].sum())} molecules; {filters}", flush=True)
        return
    blocks = select_blocks(args, recipe)
    genes = sorted({block["gene_id"] for block in blocks.values()})
    selected_genes = set(genes[args.gene_shard::args.gene_shards])
    local = {key: block for key, block in blocks.items() if block["gene_id"] in selected_genes}
    library = recipe["bams"][args.library_index]
    groups = {barcode: tuple(group) for barcode, group in library["barcode_groups"].items()}
    counts, filters = collect_library_path_reads([run["path"] for run in library["runs"]], [run["index"] for run in library["runs"]], groups, local, include_unassigned=args.include_unassigned, junction_only=args.junction_only)
    out = args.output_dir / f"part_{args.library_index}_{args.gene_shard}"
    out.mkdir(parents=True, exist_ok=True)
    rows = [{"path_key": key[0], "subject": key[1], "cell_type": key[2], "primer": key[3], "mask": key[4], "count": value} for key, value in counts.items()]
    pd.DataFrame(rows, columns=["path_key", "subject", "cell_type", "primer", "mask", "count"]).to_csv(out / "counts.tsv.gz", sep="\t", index=False)
    (out / "filters.json").write_text(json.dumps({**filters, "blocks": len(local), "genes": len(selected_genes)}, indent=2) + "\n")
    print(f"library {library['library']}, shard {args.gene_shard}: {len(local)} blocks, {sum(counts.values())} molecules, {filters}", flush=True)


if __name__ == "__main__":
    main()
