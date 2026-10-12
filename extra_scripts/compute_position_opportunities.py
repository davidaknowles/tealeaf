#!/usr/bin/env python3
"""Primer-specific position-weighted read opportunities for a local path read collection.

For every block of --local-reads (collect_local_path_reads.py --annotated),
each mature path's annotated transcripts (block cache path_index) give its
spliced distances to the transcript ends, and the read-end profiles of
estimate_read_end_profiles.py weight its read start positions
(tealeaf.sc.position_bias.weighted_path_read_opportunities). The precursor
window has uniform weight. Writes opportunities_position.json.gz,
{path key: {primer: {mask: vector}}}. --junction-only must match the
collection.
"""

import argparse
from collections import defaultdict
import gzip
import json
from multiprocessing import Pool
import os
from pathlib import Path
import re

from tealeaf.sc.local_path_reads import variable_junctions
from tealeaf.sc.position_bias import path_end_distances, weighted_path_read_opportunities

READ_LENGTH = 151
STATE = {}


def transcript_exons(gtf, genes):
    exons = defaultdict(list)
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "exon":
                continue
            gene = re.search('gene_id "([^"]+)"', fields[8]).group(1)
            if gene not in genes:
                continue
            exons[re.search('transcript_id "([^"]+)"', fields[8]).group(1)].append((int(fields[3]) - 1, int(fields[4])))
    return {name: sorted(value) for name, value in exons.items()}


def block_opportunities(item):
    key, block = item
    annotation, exons, profiles, bin_width, junction_only = STATE["annotation"], STATE["exons"], STATE["profiles"], STATE["bin_width"], STATE["junction_only"]
    row = annotation[block["block_id"]]
    names = json.loads(row["transcripts"].replace("'", '"')) if isinstance(row["transcripts"], str) else row["transcripts"]
    index = json.loads(row["path_index"]) if isinstance(row["path_index"], str) else row["path_index"]
    n_paths = len(block["signatures"])
    distances = [[] for _ in block["paths"]]
    for name, path in zip(names, index):
        if name in exons and path < n_paths:
            distances[path].extend(path_end_distances(block["paths"][path], [exons[name]], block["strand"]))
    variable = variable_junctions(block["paths"][:n_paths]) if junction_only else None
    result = weighted_path_read_opportunities(block["paths"], READ_LENGTH, block["strand"], distances, profiles, bin_width, variable)
    return key, {primer: {str(mask): vector.tolist() for mask, vector in values.items()} for primer, values in result.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--junction-only", action="store_true")
    args = parser.parse_args()
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    annotation = {row["block_id"]: row for row in json.load(gzip.open(args.block_cache, "rt"))}
    profiles = json.loads(args.profiles.read_text())
    STATE.update(annotation=annotation, exons=transcript_exons(args.gtf, {block["gene_id"] for block in blocks.values()}), profiles={primer: (value["three_prime"], value["five_prime"]) for primer, value in profiles["profiles"].items()}, bin_width=profiles["bin_width"], junction_only=args.junction_only)
    with Pool(len(os.sched_getaffinity(0))) as pool:
        opportunities = dict(pool.map(block_opportunities, list(blocks.items()), chunksize=8))
    with gzip.open(args.local_reads / "opportunities_position.json.gz", "wt") as handle:
        json.dump(opportunities, handle)
    print(f"{len(opportunities)} blocks written; profiles from {profiles['libraries']} libraries", flush=True)


if __name__ == "__main__":
    main()
