import hashlib
import json
from pathlib import Path

from tealeaf.data import ena


REPORT = """run_accession\tsample_accession\tsample_title\texperiment_accession\texperiment_title\tlibrary_layout\tfastq_ftp\tfastq_bytes\tfastq_md5
SRR1\tSAM1\tSublibrary_1\tERX1\tGSM1: x\tPAIRED\tftp.sra.ebi.ac.uk/a/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/a/SRR1_2.fastq.gz\t10;20\taaaa;bbbb
"""


def test_parse_select_and_manifest(tmp_path):
    runs = ena.parse_ena_report(REPORT)
    selected = ena.select_paired_runs(
        runs,
        title_pattern=r"Sublibrary",
        batch_from_run=lambda run: "batch1",
    )
    assert selected[0].batch == "batch1"
    assert [fastq.read for fastq in selected[0].fastqs] == [1, 2]
    manifest = tmp_path / "manifest.tsv"
    ena.write_manifest(selected, manifest)
    rows = ena.read_manifest(manifest)
    assert len(rows) == 2
    assert rows[1]["bytes"] == 20
    assert rows[0]["url"].startswith("https://")


def test_download_fastq_from_file_url(tmp_path, monkeypatch):
    payload = b"small fastq payload"
    source = tmp_path / "source.fastq.gz"
    source.write_bytes(payload)
    row = {
        "batch": "b1",
        "filename": "SRR1_1.fastq.gz",
        "url": source.as_uri(),
        "bytes": len(payload),
        "md5": hashlib.md5(payload).hexdigest(),
    }
    output = ena.download_fastq(row, tmp_path / "downloads")
    assert output.read_bytes() == payload
    monkeypatch.setattr(
        ena.subprocess,
        "run",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError()),
    )
    assert ena.download_fastq(row, tmp_path / "downloads") == output


def test_validate_completed_quantifications(tmp_path):
    runs = ena.select_paired_runs(
        ena.parse_ena_report(REPORT),
        batch_from_run=lambda run: "batch1",
    )
    manifest = tmp_path / "manifest.tsv"
    ena.write_manifest(runs, manifest)
    run_dir = tmp_path / "processed" / "batch1" / "SRR1"
    salmon = run_dir / "salmon_rad"
    fry = run_dir / "alevin_quant"
    (salmon / "aux_info").mkdir(parents=True)
    fry.mkdir()
    for path in (
        salmon / ".tealeaf_complete",
        salmon / "map.rad",
        fry / ".tealeaf_complete",
    ):
        path.write_text("complete\n")
    (salmon / "aux_info" / "meta_info.json").write_text(
        json.dumps({"num_mapped": 100})
    )
    (fry / "validation.json").write_text(
        json.dumps(
            {
                "molecules": 60,
                "cells": 4,
                "equivalence_classes": 8,
                "compatibility_nonzeros": 12,
            }
        )
    )

    report = ena.validate_completed_quantifications(
        manifest,
        tmp_path / "processed",
    )
    assert report["runs"] == 1
    assert report["molecules"] == 60

    (fry / "validation.json").write_text(
        json.dumps(
            {
                "molecules": 101,
                "cells": 4,
                "equivalence_classes": 8,
                "compatibility_nonzeros": 12,
            }
        )
    )
    try:
        ena.validate_completed_quantifications(manifest, tmp_path / "processed")
    except ValueError as error:
        assert "exceeds mapped count" in str(error)
    else:
        raise AssertionError("inflated molecule count was accepted")
