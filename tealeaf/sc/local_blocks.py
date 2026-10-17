"""Block selection and read opportunities for block-local path collections.

Blocks are keyed by path_key (block id and a hash of its path signatures).
annotated_blocks takes every annotated block of the candidate genes, with a
precursor path spanning the block window, and optionally adds unannotated
paths and new blocks from label-blind junction support (junction_support,
tealeaf.sc.novel_paths). declared_blocks keeps the path sets of EC candidate
caches.
"""

from bisect import bisect_left, bisect_right
from collections import defaultdict
import gzip
import json
from pathlib import Path
import pickle
import zlib

from .junction_paths import block_path_exons
from .local_path_reads import path_read_opportunities, variable_junctions
from .novel_paths import novel_block_paths, novel_junction_blocks, read_star_junctions

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


def junction_support(junctions):
    """Chromosome -> (sorted junctions, support) from {(chromosome, donor, acceptor): (reads, libraries, strand)}."""
    by_chromosome = defaultdict(dict)
    for (chromosome, donor, acceptor), value in junctions.items():
        by_chromosome[chromosome][(donor, acceptor)] = value
    return {chromosome: (sorted(values), values) for chromosome, values in by_chromosome.items()}


def star_junction_support(recipe):
    """junction_support from the STAR SJ.out.tab next to each run BAM of the recipe."""
    libraries = [[str(Path(run["path"]).parent / "SJ.out.tab") for run in library["runs"]] for library in recipe["bams"]]
    return junction_support(read_star_junctions(libraries))


def local_support(junctions, chromosome, strand, start, end):
    """Junction support with donor in [start, end] on chromosome, intron strand matching strand or unknown."""
    if chromosome not in junctions:
        return {}
    keys, values = junctions[chromosome]
    return {key: values[key] for key in keys[bisect_left(keys, (start,)):bisect_right(keys, (end, float("inf")))] if values[key][2] in (strand, ".")}


def block_genes(block_cache, max_paths=30):
    """Genes with at least one annotated splice block of 2..max_paths paths."""
    genes = set()
    for block in json.load(gzip.open(block_cache, "rt")):
        signatures = json.loads(block["path_signatures"]) if isinstance(block["path_signatures"], str) else block["path_signatures"]
        if 2 <= len(signatures) <= max_paths:
            genes.add(block["gene_id"])
    return genes


def annotated_blocks(candidate_caches, block_cache, max_paths=10, junctions=None, novel_options=None, gene_transcripts=None, genes=None):
    """Every annotated block (2..max_paths paths) in the target genes.

    Paths are all annotated paths of the block, independent of EC support. A
    final precursor path, one exon spanning the block window, marks reads with
    no junction inside the window (intronic and exon-body reads). With
    junctions (star_junction_support), novel_block_paths adds unannotated
    paths after the annotated ones (their chains serve as their signatures),
    and with gene_transcripts (read_gtf_exons) novel_junction_blocks adds
    blocks "<gene>:N<k>" for novel junctions outside the annotated blocks.
    genes (a set of gene ids, e.g. block_genes) replaces the candidate
    caches' genes.
    """
    if genes is None:
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
