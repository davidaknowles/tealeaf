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
