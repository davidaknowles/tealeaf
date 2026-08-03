#!/usr/bin/env python3
"""Build the reusable constitutive-anchor splice-block cache."""

from __future__ import annotations

import argparse
from pathlib import Path

from extra_scripts.run_differential_splicing import load_blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    blocks = load_blocks(args.output, args.gtf)
    print(f"blocks={len(blocks)} output={args.output}")


if __name__ == "__main__":
    main()
