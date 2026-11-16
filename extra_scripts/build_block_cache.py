#!/usr/bin/env python3
"""Write a splice-block cache (tealeaf.sc.differential.build_splice_blocks) in the format read by
tealeaf.sc.local_blocks.annotated_blocks. --local-anchors adds blocks between local anchors ("<gene>:L<k>");
--local-only writes only those."""
import argparse
import gzip
import json
from pathlib import Path

from tealeaf.sc.differential import build_splice_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--local-anchors", action="store_true")
    parser.add_argument("--local-only", action="store_true")
    args = parser.parse_args()
    blocks = build_splice_blocks(args.gtf, local_anchors=args.local_anchors or args.local_only, local_only=args.local_only)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as handle:
        json.dump([{"gene_id": b.gene_id, "block_id": b.block_id, "chromosome": b.chromosome, "strand": b.strand, "left_anchor": b.left_anchor, "right_anchor": b.right_anchor, "transcripts": b.transcripts, "path_index": b.path_index, "path_signatures": b.path_signatures} for b in blocks], handle)
    print(f"{len(blocks)} blocks, {sum(':L' in b.block_id for b in blocks)} local-anchor blocks", flush=True)


if __name__ == "__main__":
    main()
