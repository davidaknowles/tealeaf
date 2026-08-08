#!/usr/bin/env python3
"""Summarize Slurm task time and peak memory for benchmark job groups."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd


FIELDS = (
    "JobID,State,Start,End,ElapsedRaw,TotalCPU,MaxRSS,AllocCPUS"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def cpu_seconds(value):
    """Convert Slurm's day-hour:minute:second duration to seconds."""
    if not value:
        return 0.0
    days = 0
    if "-" in value:
        day, value = value.split("-", 1)
        days = int(day)
    fields = [float(part) for part in value.split(":")]
    if len(fields) == 3:
        hours, minutes, seconds = fields
    elif len(fields) == 2:
        hours = 0.0
        minutes, seconds = fields
    else:
        return fields[0]
    return 86400.0 * days + 3600.0 * hours + 60.0 * minutes + seconds


def memory_kib(value):
    """Convert a Slurm memory value requested in KiB units to KiB."""
    match = re.fullmatch(r"([0-9.]+)([KMGT]?)", value or "")
    if match is None:
        return 0.0
    scale = {
        "": 1.0 / 1024.0,
        "K": 1.0,
        "M": 1024.0,
        "G": 1024.0**2,
        "T": 1024.0**3,
    }[match.group(2)]
    return float(match.group(1)) * scale


def accounting_rows(job_ids):
    command = [
        "sacct",
        "-j",
        ",".join(job_ids),
        "-P",
        "-n",
        "--units=K",
        f"--format={FIELDS}",
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    rows = []
    for line in result.stdout.splitlines():
        fields = line.split("|")
        if len(fields) >= 8:
            rows.append(fields[:8])
    return rows


def summarize(record):
    rows = accounting_rows(str(record.job_ids).split(","))
    roots = [row for row in rows if "." not in row[0]]
    batches = [row for row in rows if row[0].endswith(".batch")]
    elapsed = np.asarray([float(row[4] or 0.0) for row in roots])
    starts = [
        datetime.fromisoformat(row[2])
        for row in roots
        if row[2] not in ("", "Unknown")
    ]
    ends = [
        datetime.fromisoformat(row[3])
        for row in roots
        if row[3] not in ("", "Unknown")
    ]
    states = pd.Series([row[1].split()[0] for row in roots], dtype=str)
    return {
        "dataset": record.dataset,
        "method": record.method,
        "job_ids": record.job_ids,
        "tasks": len(roots),
        "completed_tasks": int(states.str.startswith("COMPLETED").sum()),
        "failed_tasks": int(
            states.str.startswith(("FAILED", "OUT_OF_MEMORY", "TIMEOUT")).sum()
        ),
        "aggregate_task_hours": float(elapsed.sum() / 3600.0),
        "cpu_hours": float(sum(cpu_seconds(row[5]) for row in roots) / 3600.0),
        "median_task_minutes": float(np.median(elapsed) / 60.0),
        "p90_task_minutes": float(np.quantile(elapsed, 0.9) / 60.0),
        "campaign_hours": (
            float((max(ends) - min(starts)).total_seconds() / 3600.0)
            if starts and ends
            else np.nan
        ),
        "peak_memory_gib": float(
            max([memory_kib(row[6]) for row in batches] or [0.0]) / 1024.0**2
        ),
    }


def main():
    args = parse_args()
    jobs = pd.read_csv(args.jobs, sep="\t", dtype=str)
    required = {"dataset", "method", "job_ids"}
    if not required.issubset(jobs):
        raise ValueError(f"job table requires columns {sorted(required)}")
    result = pd.DataFrame(summarize(row) for row in jobs.itertuples(index=False))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
