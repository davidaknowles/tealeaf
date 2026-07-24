#!/usr/bin/env python3
"""Stream Parse read 1 into poly(dT)- and random-hexamer-specific FASTQs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.data.parse import demultiplex_parse_transcript_reads, read_barcodes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read1", action="append", required=True, type=Path)
    parser.add_argument("--read2", action="append", required=True, type=Path)
    parser.add_argument("--rt-barcodes", required=True, type=Path)
    parser.add_argument("--polydt-output", required=True, type=Path)
    parser.add_argument("--ranhex-output", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--rt-start", type=int, default=78)
    parser.add_argument("--pigz", type=Path)
    parser.add_argument("--pigz-threads", type=int, default=4)
    parser.add_argument("--max-reads", type=int)
    parser.add_argument("--balanced-prefix-reads", type=int)
    parser.add_argument("--no-hamming1-correction", action="store_true")
    args = parser.parse_args()

    next_progress = 10_000_000

    def report_progress(stats):
        nonlocal next_progress
        if stats["total"] >= next_progress:
            print(json.dumps({"event": "demultiplex_progress", **stats}), flush=True)
            next_progress += 10_000_000

    stats = demultiplex_parse_transcript_reads(
        args.read1,
        args.read2,
        args.polydt_output,
        args.ranhex_output,
        read_barcodes(args.rt_barcodes),
        rt_start=args.rt_start,
        correct_hamming1=not args.no_hamming1_correction,
        pigz=args.pigz,
        pigz_threads=args.pigz_threads,
        max_reads=args.max_reads,
        balanced_prefix_reads=args.balanced_prefix_reads,
        progress_callback=report_progress,
    )
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    args.stats.write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
