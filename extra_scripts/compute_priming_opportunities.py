#!/usr/bin/env python3
"""Primer-specific read opportunities with internal poly(dT) priming at A-runs, for a local path read collection.

For every block of --local-reads, each path (mature paths and the precursor window) gets its sense sequence from
--genome, and poly(dT) read start positions are weighted by the label-free priming kernel around the A-runs of that
sequence (tealeaf.sc.position_bias.primed_path_read_opportunities); random-hexamer positions keep weight 1. The
kernel is the poly(dT) read-end enrichment around A-runs of the internal priming diagnostic (--kernel, kernel.json
with delta and enrichment per primer). --junction-only must match the collection. Writes
opportunities_priming.json.gz, {path key: {primer: {mask: vector}}}, read by run_local_path_tests.py
--position-opportunities --position-opportunities-file opportunities_priming.json.gz.
"""

import argparse
import gzip
import json
from multiprocessing import Pool
import os
from pathlib import Path

from tealeaf.sc.local_path_reads import variable_junctions
from tealeaf.sc.position_bias import primed_path_read_opportunities, priming_kernel

READ_LENGTH = 151
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
STATE = {}


def initialize(genome, kernel, junction_only):
    import pysam
    STATE.update(fasta=pysam.FastaFile(str(genome)), kernel=kernel, junction_only=junction_only)


def sense_sequence(chromosome, chain, strand):
    sequence = "".join(STATE["fasta"].fetch(chromosome, a, b).upper() for a, b in sorted(chain))
    return sequence if strand == "+" else sequence.translate(COMPLEMENT)[::-1]


def block_opportunities(item):
    key, block = item
    sequences = [sense_sequence(block["chromosome"], path, block["strand"]) for path in block["paths"]]
    n_paths = len(block["signatures"])
    variable = variable_junctions(block["paths"][:n_paths]) if STATE["junction_only"] else None
    result = primed_path_read_opportunities(block["paths"], READ_LENGTH, block["strand"], sequences, STATE["kernel"], variable)
    return key, {primer: {str(mask): vector.tolist() for mask, vector in values.items()} for primer, values in result.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--junction-only", action="store_true")
    parser.add_argument("--low", type=int, default=-300, help="smallest kernel offset (A-run start - read end) used")
    parser.add_argument("--high", type=int, default=400, help="largest kernel offset used")
    args = parser.parse_args()
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    table = json.loads(args.kernel.read_text())
    kernel = priming_kernel(table["enrichment"]["poly(dT)"], table["delta"], low=args.low, high=args.high)
    with Pool(len(os.sched_getaffinity(0)), initializer=initialize, initargs=(args.genome, kernel, args.junction_only)) as pool:
        opportunities = dict(pool.map(block_opportunities, list(blocks.items()), chunksize=16))
    with gzip.open(args.local_reads / "opportunities_priming.json.gz", "wt") as handle:
        json.dump(opportunities, handle)
    print(f"{len(opportunities)} blocks written; kernel offsets {args.low}..{args.high}", flush=True)


if __name__ == "__main__":
    main()
