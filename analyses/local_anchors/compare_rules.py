"""Blocks of the consistent local-anchor rule against the blocks already fitted: the current (global-anchor) blocks and
the added local-anchor blocks ("L"), matched by content (gene and the set of path exon chains, anchors merged in).
Usage: compare_rules.py <gtf> <current block cache> <added local block cache> <output local-rule cache>"""
import gzip, json, sys
from collections import Counter
sys.path.insert(0, ".")
from tealeaf.sc.differential import build_splice_blocks
from tealeaf.sc.junction_paths import block_path_exons

gtf, current, added, out = sys.argv[1:5]
anchor = lambda a: tuple(a) if a is not None else None
content = lambda b: (b["gene_id"], frozenset(tuple(map(tuple, block_path_exons(anchor(b["left_anchor"]), [tuple(e) for e in s], anchor(b["right_anchor"])))) for s in b["path_signatures"]))
old = {content(b): b["block_id"] for b in json.load(gzip.open(current, "rt"))}
new_l = {content(b): b["block_id"] for b in json.load(gzip.open(added, "rt"))}
blocks = build_splice_blocks(gtf, anchor_rule="local")
rows = [{"gene_id": b.gene_id, "block_id": b.block_id, "chromosome": b.chromosome, "strand": b.strand, "left_anchor": None, "right_anchor": None, "transcripts": b.transcripts, "path_index": b.path_index, "path_signatures": b.path_signatures} for b in blocks]
with gzip.open(out, "wt") as handle:
    json.dump(rows, handle)
c = Counter()
used_old = set()
for r in rows:
    k = content(r)
    if k in old:
        c["identical to a current block"] += 1; used_old.add(k)
    elif k in new_l:
        c["identical to an added local block"] += 1
    else:
        c["new (fit needed)"] += 1
print(f"local-rule blocks: {len(rows):,}")
for k, v in c.items():
    print(f"  {k}: {v:,}")
print(f"current blocks dropped (replaced by split blocks): {len(set(old) - used_old):,} of {len(old):,}")
