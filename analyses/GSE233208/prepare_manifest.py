#!/usr/bin/env python3
"""Build the five-batch Parse snRNA-seq manifest for GSE233208."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re

from tealeaf.data.ena import fetch_ena_runs, select_paired_runs, write_manifest


GSM_BATCHES = (
    (7412790, 7412797, "batch1"),
    (7412799, 7412806, "batch2"),
    (7412807, 7412814, "batch3"),
    (7412815, 7412822, "batch4"),
    (7412823, 7412830, "batch5"),
)


def batch_from_run(run):
    match = re.search(r"GSM(\d+)", run.experiment_title)
    if not match:
        raise ValueError(
            f"no GEO sample accession in experiment title: {run.experiment_title}"
        )
    gsm = int(match.group(1))
    for first, last, batch in GSM_BATCHES:
        if first <= gsm <= last:
            return batch
    raise ValueError(f"GSM{gsm} is not a selected Parse snRNA-seq sublibrary")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    runs = fetch_ena_runs("PRJNA975472")
    selected = select_paired_runs(
        runs,
        title_pattern=r"(?i)(?:^|_)sublibrary_?\d+",
        batch_from_run=batch_from_run,
    )
    counts = Counter(run.batch for run in selected)
    if counts != {f"batch{i}": 8 for i in range(1, 6)}:
        raise ValueError(f"expected eight sublibraries per batch, observed {counts}")
    write_manifest(selected, args.output)
    report = {
        "project": "PRJNA975472",
        "series": "GSE233208",
        "runs": len(selected),
        "fastqs": sum(len(run.fastqs) for run in selected),
        "compressed_bytes": sum(
            fastq.bytes for run in selected for fastq in run.fastqs
        ),
        "runs_by_batch": dict(sorted(counts.items())),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
