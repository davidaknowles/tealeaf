import csv
import subprocess
import sys

import pytest


def test_aggregate_glm_scores(tmp_path):
    fits = tmp_path / "fits.tsv"
    fits.write_text("first\t/path/first_\nsecond\t/path/second_\n")
    scores = tmp_path / "scores"
    for name, accuracy in (("first", "0.5"), ("second", "0.75")):
        directory = scores / name
        directory.mkdir(parents=True)
        (directory / "label_score_summary.csv").write_text(
            f"name,status,accuracy\n{name},ok,{accuracy}\n"
        )

    subprocess.run(
        [
            sys.executable,
            "extra_scripts/aggregate_glm_scores.py",
            "--fits-file",
            str(fits),
            "--score-output",
            str(scores),
        ],
        check=True,
    )

    with open(scores / "label_score_summary.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["name"] for row in rows] == ["first", "second"]
    assert [row["accuracy"] for row in rows] == ["0.5", "0.75"]


def test_aggregate_glm_scores_rejects_invalid_fit(tmp_path):
    fits = tmp_path / "fits.tsv"
    fits.write_text("failed\t/path/failed_\n")
    scores = tmp_path / "scores" / "failed"
    scores.mkdir(parents=True)
    (scores / "label_score_summary.csv").write_text(
        "name,status,error\nfailed,invalid,nonfinite factors\n"
    )

    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(
            [
                sys.executable,
                "extra_scripts/aggregate_glm_scores.py",
                "--fits-file",
                str(fits),
                "--score-output",
                str(tmp_path / "scores"),
            ],
            check=True,
        )
