#!/usr/bin/env python3
"""Validate all per-run Salmon and alevin-fry outputs before merging."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from tealeaf.data.ena import validate_completed_quantifications


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--processed-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    report = validate_completed_quantifications(
        args.manifest,
        args.processed_root,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f".{args.output.name}.{os.getpid()}.tmp")
    with open(temporary, "w") as handle:
        json.dump(report, handle, indent=2)
    os.replace(temporary, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
