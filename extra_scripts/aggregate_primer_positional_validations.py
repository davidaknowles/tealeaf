#!/usr/bin/env python3
"""Aggregate per-run primer positional validation reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.data.salmon import summarize_primer_positional_validations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summary = summarize_primer_positional_validations(args.validation)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
