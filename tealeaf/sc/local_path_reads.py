"""Block-local read compatibility classes for local splice paths.

A read (or UMI molecule) is classified only by the part of its alignment that
falls inside the block window, the span of the block's path exons including
its constitutive anchors. Its class is the bitmask of paths it is compatible
with: every clipped aligned segment lies within one exon of the path, and
every splice junction inside the window is a junction of that path. Intronic
or unannotated-junction reads match no path and get mask 0. Reads elsewhere in
the gene never enter, unlike gene-level ECs.

Read opportunities n[mask][s] count the start positions on path s's local
spliced sequence whose (clipped) read falls in each class, using the same
classifier, so a path's expected class frequencies include its effective
length. Dropping the all-paths class gives the anchor-free variant.
"""

from collections import defaultdict
from contextlib import ExitStack

import numpy as np

from .junction_benchmark import normalize_starsolo_barcode
from .sashimi import _alignment_junctions


def path_structures(paths):
    """Sorted exon tuples and junction sets for S genomic exon chains."""
    exons = [tuple(sorted(tuple(int(x) for x in exon) for exon in path)) for path in paths]
    junctions = [{(a[1], b[0]) for a, b in zip(chain, chain[1:]) if a[1] < b[0]} for chain in exons]
    start = min(chain[0][0] for chain in exons)
    end = max(chain[-1][1] for chain in exons)
    return exons, junctions, (start, end)


def compatibility_mask(blocks, junctions, structures):
    """Bitmask of paths compatible with one aligned read inside the window.

    blocks are aligned (start, end) intervals and junctions (donor end,
    acceptor start) gaps, both 0-based half-open as from pysam. Returns 0 when
    the read does not overlap the window or fits no path.
    """
    exons, path_junctions, (start, end) = structures
    clipped = [(max(a, start), min(b, end)) for a, b in blocks if a < end and b > start]
    if not clipped:
        return 0
    inside = [(a, b) for a, b in junctions if a >= start and b <= end]
    mask = 0
    for index, (chain, known) in enumerate(zip(exons, path_junctions)):
        if all(any(x <= a and b <= y for x, y in chain) for a, b in clipped) and all(junction in known for junction in inside):
            mask |= 1 << index
    return mask


def path_read_opportunities(paths, read_length):
    """n[mask] -> length-S vector of start positions per path yielding mask.

    Each path's local transcript is its exon chain, reads of read_length start
    at every position that overlaps it by at least one base, and the part
    outside the chain is clipped as reads leaving the window are in the data.
    Mask 0 positions are not returned.
    """
    structures = path_structures(paths)
    opportunities = defaultdict(lambda: np.zeros(len(paths)))
    for index, chain in enumerate(structures[0]):
        offsets = np.cumsum([0] + [b - a for a, b in chain])
        length = offsets[-1]
        for position in range(1 - read_length, length):
            first, last = max(position, 0), min(position + read_length, length)
            blocks, previous = [], None
            for exon, (a, b) in enumerate(chain):
                low, high = max(first, offsets[exon]), min(last, offsets[exon + 1])
                if low < high:
                    blocks.append((a + low - offsets[exon], a + high - offsets[exon]))
            junctions = [(blocks[i][1], blocks[i + 1][0]) for i in range(len(blocks) - 1) if blocks[i][1] < blocks[i + 1][0]]
            mask = compatibility_mask(blocks, junctions, structures)
            if mask:
                opportunities[mask][index] += 1
    return dict(opportunities)


def collect_library_path_reads(bam_paths, index_paths, barcode_groups, blocks):
    """UMI-molecule path classes for one physical library.

    blocks maps block_id -> dict(gene_id, chromosome, paths). Keys are exact
    (barcode, UB, gene) within the library, unioned across its runs; a
    molecule's mask is the AND of its reads' masks over reads overlapping the
    block window. Only primary NH=1 alignments with CB, UB and a unique
    matching GX gene are used. Returns counts[(block_id, *group, mask)] and
    filter counters; mask-0 (incompatible or conflicting) molecules are
    counted separately.
    """
    import pysam

    by_gene = defaultdict(list)
    for block_id, block in blocks.items():
        by_gene[(block["chromosome"], block["gene_id"])].append((block_id, path_structures(block["paths"])))
    counts, filters = defaultdict(int), defaultdict(int)
    with ExitStack() as stack:
        sources = [stack.enter_context(pysam.AlignmentFile(str(bam), "rb", index_filename=str(index))) for bam, index in zip(bam_paths, index_paths)]
        for (chromosome, gene), local in sorted(by_gene.items()):
            left = min(structures[2][0] for _, structures in local)
            right = max(structures[2][1] for _, structures in local)
            gene_stable = gene.split(".", 1)[0]
            masks = {}
            for source in sources:
                if chromosome not in source.references:
                    continue
                for alignment in source.fetch(chromosome, left, right):
                    if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                        continue
                    if not all(alignment.has_tag(tag) for tag in ("NH", "CB", "UB", "GX")) or alignment.get_tag("NH") != 1:
                        filters["untagged_or_multimapped"] += 1
                        continue
                    genes = str(alignment.get_tag("GX")).split(";")
                    if len(genes) != 1 or genes[0].split(".", 1)[0] != gene_stable:
                        filters["other_gene"] += 1
                        continue
                    barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
                    if barcode not in barcode_groups:
                        continue
                    aligned = tuple(alignment.get_blocks())
                    junctions = tuple(_alignment_junctions(alignment))
                    for block_id, structures in local:
                        start, end = structures[2]
                        if not any(a < end and b > start for a, b in aligned):
                            continue
                        key = (block_id, barcode, str(alignment.get_tag("UB")))
                        mask = compatibility_mask(aligned, junctions, structures)
                        masks[key] = masks.get(key, -1) & mask if key in masks else mask
            for (block_id, barcode, _), mask in masks.items():
                if mask:
                    counts[(block_id, *barcode_groups[barcode], mask)] += 1
                else:
                    filters["incompatible_molecules"] += 1
            filters["molecules"] += len(masks)
    return dict(counts), dict(filters)
