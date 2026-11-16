"""Prototype of local anchors for splice blocks, and how many genes and blocks they change.

Current rule (tealeaf.sc.differential.build_splice_blocks): an anchor is an interval inside an exon of every annotated
transcript of the gene, so one transcript that starts or ends inside the gene removes every anchor it does not reach.
Local rule: an interval is a local anchor if every transcript whose span (first exon start to last exon end) contains
it has an exon covering it, and at least two transcripts cover it. A local block is the window between consecutive
local anchors; its paths are the distinct signatures (exon chains clipped to the window) of the transcripts that span
the whole window. Reports, per annotation, genes whose anchors change, current blocks that the local rule splits, and
new internal blocks (>= 2 spanning signatures), and prints both rules' blocks for a gene.
Usage: prototype.py <gtf> [gene id to show]"""
import sys
from collections import Counter
sys.path.insert(0, ".")
from tealeaf.sc.differential import _clipped_signature, _constitutive_intervals, _merge_intervals, read_gtf_exons


def local_anchor_intervals(chains):
    boundaries = sorted({c for chain in chains for exon in chain for c in exon})
    spans = [(chain[0][0], chain[-1][1]) for chain in chains]
    out = []
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        if start == end:
            continue
        spanning = [chain for chain, (a, b) in zip(chains, spans) if a <= start and end <= b]
        covered = [chain for chain in spanning if any(s <= start and end <= e for s, e in chain)]
        if len(covered) >= 2 and len(covered) == len(spanning):
            out.append((start, end))
    return _merge_intervals(out)


def blocks(chains, anchors, strand, spanning_only):
    result = []
    for left, right in zip(anchors[:-1], anchors[1:]):
        start, end = left[1], right[0]
        if start >= end:
            continue
        members = [c for c in chains if c[0][0] <= left[0] and right[1] <= c[-1][1]] if spanning_only else chains
        signatures = {_clipped_signature(c, start, end, strand) for c in members}
        if len(signatures) >= 2:
            result.append(((start, end), sorted(signatures)))
    return result


genes = read_gtf_exons(sys.argv[1])
show = sys.argv[2] if len(sys.argv) > 2 else None
stats = Counter()
for gene, record in genes.items():
    chains = list(record["transcripts"].values())
    if len(chains) < 2:
        continue
    stats["multi-transcript genes"] += 1
    glob = _constitutive_intervals(chains)
    loc = local_anchor_intervals(chains)
    old = blocks(chains, glob, record["strand"], False)
    new = blocks(chains, loc, record["strand"], True)
    if not glob:
        stats["genes with no global anchor (one whole-gene block)"] += 1
    if loc != glob:
        stats["genes whose anchors change"] += 1
    stats["internal blocks, global rule"] += len(old)
    stats["internal blocks, local rule"] += len(new)
    old_windows = {w: set(s) for w, s in old}
    fresh = [(w, s) for w, s in new if w not in old_windows]
    same = [(w, s) for w, s in new if w in old_windows and set(s) == old_windows[w]]
    stats["local blocks identical to a current internal block"] += len(same)
    stats["local blocks with a new window (variation now tested only inside a terminal or whole-gene block)"] += len(fresh)
    if fresh:
        stats["genes with a new local block"] += 1
    stats["paths in new local blocks"] += sum(len(s) for _, s in fresh)
    clean = [w for w, _ in fresh if all(c[0][0] <= w[0] - 1 and w[1] + 1 <= c[-1][1] or c[-1][1] <= w[0] or c[0][0] >= w[1] for c in chains)]
    stats["new local blocks with no transcript starting or ending inside the window"] += len(clean)
    if gene.split(".")[0] == show:
        print(f"{gene} {record['strand']}: {len(chains)} transcripts")
        print(" global anchors:", glob)
        print(" local anchors:", loc)
        for label, bl in (("global", old), ("local", new)):
            for w, sigs in bl:
                print(f"  {label} block {w}: {len(sigs)} paths, {[len(s) for s in sigs]} exons each")
for k, v in stats.items():
    print(f"{k}: {v:,}")
