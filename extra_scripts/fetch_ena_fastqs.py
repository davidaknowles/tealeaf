#!/usr/bin/env python3
"""Download all FASTQs in a normalized ENA manifest."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tealeaf.data.ena import download_fastq, read_manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--retries", type=int, default=5)
    args = parser.parse_args()
    rows = read_manifest(args.manifest)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                download_fastq,
                row,
                args.destination,
                retries=args.retries,
            ): row
            for row in rows
        }
        for future in as_completed(futures):
            path = future.result()
            print(f"verified\t{path}", flush=True)


if __name__ == "__main__":
    main()
