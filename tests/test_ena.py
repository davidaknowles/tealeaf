import hashlib
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
