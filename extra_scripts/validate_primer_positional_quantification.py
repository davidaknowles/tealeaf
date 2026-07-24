#!/usr/bin/env python3
"""Validate one primer-separated Parse Salmon quantification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.data.salmon import validate_primer_positional_quantification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--source-meta", required=True, type=Path)
    parser.add_argument("--expected-targets", type=int)
    parser.add_argument("--expected-library-type", default="U")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = validate_primer_positional_quantification(
        args.run_root,
        args.source_meta,
        expected_targets=args.expected_targets,
        expected_library_type=args.expected_library_type,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
