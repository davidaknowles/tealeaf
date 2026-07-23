#!/usr/bin/env python3
"""Precompute fixed weighted EC designs for standard and paired-primer GLMs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from tealeaf.sc import glm_cv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--min-half-umis", type=float, default=500)
    args = parser.parse_args()
    standard = glm_cv.prepare_alevin_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        ec_design="weighted",
        regularization_target="theta",
    )
    report = {
        "standard_cells": int(standard.counts.shape[0]),
        "standard_equivalence_classes": int(standard.counts.shape[1]),
        "standard_transcripts": int(standard.compatibility.shape[1]),
    }
    if args.primer_pairs is not None:
        paired = glm_cv.prepare_paired_primer_glm_data(
            args.alevin_dir,
            args.salmon_ref,
            args.primer_pairs,
            ec_design="weighted",
            regularization_target="theta",
            min_half_umis=args.min_half_umis,
        )
        report.update(
            paired_cells=int(paired.counts.shape[0]),
            paired_equivalence_classes=int(paired.counts.shape[1]),
            paired_transcripts=int(paired.compatibility.shape[1]),
        )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
