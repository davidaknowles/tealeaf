"""ENA run discovery and checksum-verified FASTQ acquisition."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from urllib.parse import urlencode
from urllib.request import urlopen


ENA_REPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession",
    "sample_accession",
    "sample_title",
    "experiment_accession",
    "experiment_title",
    "library_layout",
    "fastq_ftp",
    "fastq_bytes",
    "fastq_md5",
)
MANIFEST_FIELDS = (
    "batch",
    "run_accession",
    "sample_accession",
    "sample_title",
    "experiment_accession",
    "experiment_title",
    "read",
    "url",
    "bytes",
    "md5",
    "filename",
)


@dataclass(frozen=True)
class EnaFastq:
    read: int
    url: str
    bytes: int
    md5: str
    filename: str


@dataclass(frozen=True)
class EnaRun:
    run_accession: str
    sample_accession: str
    sample_title: str
    experiment_accession: str
    experiment_title: str
    library_layout: str
    fastqs: tuple[EnaFastq, ...]
    batch: str = ""


def _split_field(value: str) -> list[str]:
    return [] if not value else value.split(";")


def _read_number(path: str) -> int:
    match = re.search(r"(?:^|[_./])([12])\.f(?:ast)?q\.gz$", path)
    if match:
        return int(match.group(1))
    match = re.search(r"_([12])(?:_001)?\.f(?:ast)?q\.gz$", path)
    if match:
        return int(match.group(1))
    raise ValueError(f"cannot infer read number from ENA FASTQ path: {path}")


def parse_ena_report(text: str) -> list[EnaRun]:
    """Parse an ENA filereport TSV into validated paired-run records."""
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    missing = set(ENA_FIELDS) - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"ENA report is missing fields: {sorted(missing)}")
    runs = []
    for row in reader:
        paths = _split_field(row["fastq_ftp"])
        sizes = _split_field(row["fastq_bytes"])
        checksums = _split_field(row["fastq_md5"])
        if not (len(paths) == len(sizes) == len(checksums)):
            raise ValueError(
                f"inconsistent FASTQ metadata for {row['run_accession']}"
            )
        fastqs = []
        for path, size, checksum in zip(paths, sizes, checksums):
            url = path if "://" in path else f"https://{path}"
            fastqs.append(
                EnaFastq(
                    read=_read_number(path),
                    url=url,
                    bytes=int(size),
                    md5=checksum.lower(),
                    filename=Path(path).name,
                )
            )
        fastqs.sort(key=lambda item: item.read)
        runs.append(
            EnaRun(
                run_accession=row["run_accession"],
                sample_accession=row["sample_accession"],
                sample_title=row["sample_title"],
                experiment_accession=row["experiment_accession"],
                experiment_title=row["experiment_title"],
                library_layout=row["library_layout"],
                fastqs=tuple(fastqs),
            )
        )
    return runs


def fetch_ena_runs(accession: str, *, timeout: int = 120) -> list[EnaRun]:
    """Query ENA for all read runs associated with an accession."""
    query = urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
        }
    )
    with urlopen(f"{ENA_REPORT_URL}?{query}", timeout=timeout) as response:
        return parse_ena_report(response.read().decode("utf-8"))


def select_paired_runs(
    runs: list[EnaRun],
    *,
    title_pattern: str | None = None,
    batch_from_run=None,
) -> list[EnaRun]:
    """Filter paired runs and assign a caller-defined biological batch."""
    pattern = re.compile(title_pattern) if title_pattern else None
    selected = []
    for run in runs:
        if pattern and not pattern.search(run.sample_title):
            continue
        if run.library_layout != "PAIRED":
            raise ValueError(f"{run.run_accession} is not paired-end")
        if [item.read for item in run.fastqs] != [1, 2]:
            raise ValueError(
                f"{run.run_accession} does not have exactly reads 1 and 2"
            )
        batch = batch_from_run(run) if batch_from_run else run.batch
        if not batch:
            raise ValueError(f"no batch assigned to {run.run_accession}")
        selected.append(
            EnaRun(
                run_accession=run.run_accession,
                sample_accession=run.sample_accession,
                sample_title=run.sample_title,
                experiment_accession=run.experiment_accession,
                experiment_title=run.experiment_title,
                library_layout=run.library_layout,
                fastqs=run.fastqs,
                batch=str(batch),
            )
        )
    if not selected:
        raise ValueError("no ENA runs matched the requested selection")
    return sorted(selected, key=lambda run: (run.batch, run.run_accession))


def write_manifest(runs: list[EnaRun], path) -> None:
    """Write one normalized manifest row per FASTQ."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(temporary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        for run in runs:
            for fastq in run.fastqs:
                writer.writerow(
                    {
                        "batch": run.batch,
                        "run_accession": run.run_accession,
                        "sample_accession": run.sample_accession,
                        "sample_title": run.sample_title,
                        "experiment_accession": run.experiment_accession,
                        "experiment_title": run.experiment_title,
                        "read": fastq.read,
                        "url": fastq.url,
                        "bytes": fastq.bytes,
                        "md5": fastq.md5,
                        "filename": fastq.filename,
                    }
                )
    os.replace(temporary, path)


def read_manifest(path) -> list[dict]:
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    missing = set(MANIFEST_FIELDS) - set(rows[0] if rows else ())
    if missing:
        raise ValueError(f"manifest is missing fields: {sorted(missing)}")
    for row in rows:
        row["read"] = int(row["read"])
        row["bytes"] = int(row["bytes"])
    return rows


def validate_completed_quantifications(manifest, processed_root) -> dict:
    """Validate restart markers and count reports for every manifest run."""
    rows = read_manifest(manifest)
    runs = {
        (row["batch"], row["run_accession"])
        for row in rows
        if row["read"] == 1
    }
    if not runs:
        raise ValueError("manifest contains no read-1 runs")

    reports = []
    failures = []
    for batch, accession in sorted(runs):
        run_dir = Path(processed_root) / batch / accession
        salmon_dir = run_dir / "salmon_rad"
        fry_dir = run_dir / "alevin_quant"
        required = (
            salmon_dir / ".tealeaf_complete",
            salmon_dir / "map.rad",
            salmon_dir / "aux_info" / "meta_info.json",
            fry_dir / ".tealeaf_complete",
            fry_dir / "validation.json",
        )
        missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
        if missing:
            failures.append(f"{accession}: missing {', '.join(missing)}")
            continue
        with open(salmon_dir / "aux_info" / "meta_info.json") as handle:
            salmon = json.load(handle)
        with open(fry_dir / "validation.json") as handle:
            fry = json.load(handle)
        mapped = salmon.get("num_mapped")
        molecules = fry.get("molecules")
        cells = fry.get("cells")
        equivalence_classes = fry.get("equivalence_classes")
        compatibility_nonzeros = fry.get("compatibility_nonzeros")
        if not isinstance(mapped, (int, float)) or mapped <= 0:
            failures.append(f"{accession}: invalid Salmon num_mapped")
            continue
        if not isinstance(molecules, (int, float)) or not 0 < molecules <= mapped:
            failures.append(
                f"{accession}: molecule count {molecules!r} exceeds mapped count {mapped!r}"
            )
            continue
        if not isinstance(cells, int) or cells <= 0:
            failures.append(f"{accession}: invalid cell count {cells!r}")
            continue
        if not isinstance(equivalence_classes, int) or equivalence_classes <= 0:
            failures.append(
                f"{accession}: invalid equivalence-class count {equivalence_classes!r}"
            )
            continue
        if not isinstance(compatibility_nonzeros, int) or compatibility_nonzeros <= 0:
            failures.append(
                f"{accession}: invalid compatibility count {compatibility_nonzeros!r}"
            )
            continue
        reports.append(
            {
                "batch": batch,
                "run_accession": accession,
                "mapped_fragments": mapped,
                "molecules": molecules,
                "cells": cells,
                "equivalence_classes": equivalence_classes,
            }
        )
    if failures:
        raise ValueError("incomplete quantifications:\n" + "\n".join(failures))
    return {
        "runs": len(reports),
        "mapped_fragments": sum(item["mapped_fragments"] for item in reports),
        "molecules": sum(item["molecules"] for item in reports),
        "cells": sum(item["cells"] for item in reports),
        "quantifications": reports,
    }


def file_md5(path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_download(path, expected_bytes: int, expected_md5: str) -> bool:
    path = Path(path)
    return (
        path.is_file()
        and path.stat().st_size == int(expected_bytes)
        and file_md5(path) == expected_md5.lower()
    )


def download_fastq(row: dict, destination, *, retries: int = 5) -> Path:
    """Resume one curl transfer, verify it, and atomically publish it."""
    destination = Path(destination)
    batch_dir = destination / row["batch"]
    batch_dir.mkdir(parents=True, exist_ok=True)
    output = batch_dir / row["filename"]
    if verify_download(output, row["bytes"], row["md5"]):
        return output
    partial = output.with_suffix(output.suffix + ".part")
    for attempt in range(1, retries + 1):
        command = [
            "curl",
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "--retry",
            "8",
            "--retry-all-errors",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            row["url"],
        ]
        subprocess.run(command, check=True)
        if verify_download(partial, row["bytes"], row["md5"]):
            os.replace(partial, output)
            return output
        if partial.exists():
            quarantine = partial.with_name(
                f"{partial.name}.invalid-{int(time.time())}-attempt{attempt}"
            )
            shutil.move(partial, quarantine)
        if attempt < retries:
            time.sleep(min(60, 2**attempt))
    raise RuntimeError(f"failed to download and verify {row['url']}")
