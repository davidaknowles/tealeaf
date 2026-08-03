import subprocess
import sys

import pandas as pd


def test_combined_metadata_preserves_cell_type_underscores(tmp_path):
    source = tmp_path / "groups.csv"
    source.write_text("AAA,EX_cell_type__case__mouse_1\n")
    output = tmp_path / "cells.tsv"
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/prepare_junction_cell_metadata.py",
            "--combined-groups",
            str(source),
            "--output",
            str(output),
        ],
        check=True,
    )
    row = pd.read_csv(output, sep="\t").iloc[0]
    assert row.to_dict() == {
        "run": "*",
        "barcode": "AAA",
        "cell_type": "EX_cell_type",
        "condition": "case",
        "subject": "mouse_1",
    }
