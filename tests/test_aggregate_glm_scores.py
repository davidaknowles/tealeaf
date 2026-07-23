import csv
import subprocess
import sys


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
