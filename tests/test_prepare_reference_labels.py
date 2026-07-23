import csv
import subprocess
import sys


def test_prepare_reference_labels(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    with open(metadata, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["cell_barcode", "Batch", "annotation", "CaseNum"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "cell_barcode": "ACGT",
                "Batch": "Batch4",
                "annotation": "MG",
                "CaseNum": "donor1",
            }
        )
    labels = tmp_path / "labels.csv"
    groups = tmp_path / "groups.csv"
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/prepare_reference_labels.py",
            "--metadata",
            str(metadata),
            "--batch-map",
            "batch2=Batch4",
            "--labels-output",
            str(labels),
            "--groups-output",
            str(groups),
        ],
        check=True,
    )
    assert labels.read_text() == "batch2:ACGT,MG\n"
    assert groups.read_text() == "batch2:ACGT,donor1\n"


def test_prepare_reference_labels_expands_parse_primer_halves(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "cell_barcode\tBatch\tannotation\tCaseNum\n"
        "AACCGGTTAACCGGTT_1\tBatch1\tEX1\tdonor1\n"
    )
    rt = [f"{value:08d}" for value in range(96)]
    rt_file = tmp_path / "rt.txt"
    rt_file.write_text("\n".join(rt) + "\n")
    labels = tmp_path / "labels.csv"
    groups = tmp_path / "groups.csv"
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/prepare_reference_labels.py",
            "--metadata",
            str(metadata),
            "--batch-map",
            "batch1=Batch1",
            "--parse-rt-barcodes",
            str(rt_file),
            "--labels-output",
            str(labels),
            "--groups-output",
            str(groups),
        ],
        check=True,
    )
    assert labels.read_text().splitlines() == [
        "batch1:AACCGGTTAACCGGTT00000001,EX1",
        "batch1:AACCGGTTAACCGGTT00000049,EX1",
    ]


def test_prepare_reference_labels_uses_run_for_sublibrary(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "cell_barcode\tSublibrary\tBatch\tannotation\tCaseNum\n"
        "AACCGGTTAACCGGTT_1\tSublibrary_1_S1\tBatch5\tEX1\tdonor1\n"
    )
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "batch\trun_accession\tsample_title\tread\n"
        "batch5\tSRR1\tbatch5_Sublibrary_1_S1_L004\t1\n"
    )
    rt_file = tmp_path / "rt.txt"
    rt_file.write_text("\n".join(f"{value:08d}" for value in range(96)) + "\n")
    labels = tmp_path / "labels.csv"
    groups = tmp_path / "groups.csv"
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/prepare_reference_labels.py",
            "--metadata",
            str(metadata),
            "--batch-map",
            "batch5=Batch5",
            "--ena-manifest",
            str(manifest),
            "--parse-rt-barcodes",
            str(rt_file),
            "--labels-output",
            str(labels),
            "--groups-output",
            str(groups),
        ],
        check=True,
    )
    assert labels.read_text().splitlines()[0].startswith("SRR1:")
