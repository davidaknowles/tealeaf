#!/usr/bin/env python3
"""Merge independently quantified alevin-fry batches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.data.alevin import merge_alevin_quantifications


def named_path(value):
    name, separator, path = value.partition("=")
    if not separator or not name or not path:
        raise argparse.ArgumentTypeError("inputs must have the form BATCH=PATH")
    return name, Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=named_path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = merge_alevin_quantifications(args.input, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
