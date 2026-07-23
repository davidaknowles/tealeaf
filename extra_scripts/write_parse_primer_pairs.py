#!/usr/bin/env python3
"""Write paired Parse poly(dT)/random-hexamer cell identifiers."""

from __future__ import annotations

import argparse
from pathlib import Path

from tealeaf.data.parse import write_primer_pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-rows", required=True, type=Path)
    parser.add_argument("--rt-barcodes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    count = write_primer_pairs(args.cell_rows, args.rt_barcodes, args.output)
    print(f"wrote {count} primer pairs")


if __name__ == "__main__":
    main()
