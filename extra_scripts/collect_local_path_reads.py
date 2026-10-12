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
from bisect import bisect_left, bisect_right
from collections import defaultdict
import gzip
import json
from multiprocessing import Pool
import os
import pickle
from pathlib import Path
import zlib

import pandas as pd

from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.junction_paths import block_path_exons
from tealeaf.sc.local_path_reads import collect_library_path_reads, path_read_opportunities, variable_junctions
from tealeaf.sc.novel_paths import novel_block_paths, novel_junction_blocks, read_star_junctions

READ_LENGTH = 151


def block_opportunities(block, junction_only=False):
    variable = variable_junctions(block["paths"][:len(block["signatures"])]) if junction_only else None
    return {str(mask): vector.tolist() for mask, vector in path_read_opportunities(block["paths"], READ_LENGTH, variable).items()}


def junction_block_opportunities(block):
    return block_opportunities(block, True)


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


def star_junction_support(recipe):
    """Chromosome -> (sorted junctions, support) from the SJ.out.tab next to each run BAM."""
    libraries = [[str(Path(run["path"]).parent / "SJ.out.tab") for run in library["runs"]] for library in recipe["bams"]]
    by_chromosome = defaultdict(dict)
    for (chromosome, donor, acceptor), value in read_star_junctions(libraries).items():
        by_chromosome[chromosome][(donor, acceptor)] = value
    return {chromosome: (sorted(values), values) for chromosome, values in by_chromosome.items()}


def local_support(junctions, chromosome, strand, start, end):
    """Junction support with donor in [start, end] on chromosome, intron strand matching strand or unknown."""
    if chromosome not in junctions:
        return {}
    keys, values = junctions[chromosome]
    return {key: values[key] for key in keys[bisect_left(keys, (start,)):bisect_right(keys, (end, float("inf")))] if values[key][2] in (strand, ".")}


def annotated_blocks(candidate_caches, block_cache, max_paths=10, junctions=None, novel_options=None, gene_transcripts=None):
    """Every annotated block (2..max_paths paths) in the candidate caches' genes.

    Paths are all annotated paths of the block, independent of EC support. A
    final precursor path, one exon spanning the block window, marks reads with
    no junction inside the window (intronic and exon-body reads). With
    junctions (star_junction_support), novel_block_paths adds unannotated
    paths after the annotated ones (their chains serve as their signatures),
    and with gene_transcripts (read_gtf_exons) novel_junction_blocks adds
    blocks "<gene>:N<k>" for novel junctions outside the annotated blocks.
    """
    genes = set()
    for cache in candidate_caches:
        genes.update(candidate[2] for candidate in pickle.load(open(cache, "rb"))["candidates"])
    blocks, windows = {}, defaultdict(list)
    for block in json.load(gzip.open(block_cache, "rt")):
        signatures = json.loads(block["path_signatures"]) if isinstance(block["path_signatures"], str) else block["path_signatures"]
        if block["gene_id"] not in genes:
            continue
        paths = [block_path_exons(anchor(block["left_anchor"]), signature, anchor(block["right_anchor"])) for signature in signatures]
        window = (min(a for path in paths for a, _ in path), max(b for path in paths for _, b in path))
        windows[block["gene_id"]].append(window)
        if not 2 <= len(signatures) <= max_paths:
            continue
        novel = []
        if junctions is not None:
            novel = novel_block_paths(paths, local_support(junctions, block["chromosome"], block["strand"], *window), **(novel_options or {}))
        signatures = signatures + [chain for chain, _ in novel]
        paths = paths + [chain for chain, _ in novel]
        blocks[path_key(block["block_id"], signatures)] = {"block_id": block["block_id"], "gene_id": block["gene_id"], "chromosome": block["chromosome"], "strand": block["strand"], "signatures": signatures, "n_annotated": len(signatures) - len(novel), "novel_reads": [reads for _, reads in novel], "paths": paths + [[window]], "precursor": True}
    if junctions is not None and gene_transcripts is not None:
        for gene in sorted(genes):
            if gene not in gene_transcripts:
                continue
            record = gene_transcripts[gene]
            chains = list(record["transcripts"].values())
            start, end = min(chain[0][0] for chain in chains), max(chain[-1][1] for chain in chains)
            support = local_support(junctions, record["chromosome"], record["strand"], start, end)
            for index, (annotated, novel) in enumerate(novel_junction_blocks(chains, support, windows[gene], **(novel_options or {}))):
                signatures = annotated + [chain for chain, _ in novel]
                window = (min(path[0][0] for path in signatures), max(path[-1][1] for path in signatures))
                block_id = f"{gene}:N{index + 1}"
                blocks[path_key(block_id, signatures)] = {"block_id": block_id, "gene_id": gene, "chromosome": record["chromosome"], "strand": record["strand"], "signatures": signatures, "n_annotated": len(annotated), "novel_reads": [reads for _, reads in novel], "paths": signatures + [[window]], "precursor": True}
    return blocks


def select_blocks(args, recipe=None):
    if args.novel_paths:
        if not args.annotated:
            raise ValueError("--novel-paths requires --annotated")
        options = {"min_reads": args.novel_min_reads, "min_libraries": args.novel_min_libraries, "max_novel": args.novel_max}
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
