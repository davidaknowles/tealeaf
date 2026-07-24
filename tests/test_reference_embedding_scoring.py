import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd


def test_score_reference_embedding_maps_and_filters_cells(tmp_path):
    rng = np.random.default_rng(7)
    rows = []
    maps = []
    labels = []
    groups = []
    eligible = []
    for index in range(60):
        label = f"type{index % 3}"
        source = f"published{index}"
        target = f"cell{index}"
        center = np.zeros(4)
        center[index % 3] = 4
        rows.append([source, *(center + rng.normal(scale=0.1, size=4))])
        maps.append([source, target])
        labels.append([target, label])
        groups.append([target, f"donor{index % 6}"])
        if index < 48:
            eligible.append(target)

    pd.DataFrame(
        rows, columns=["cell_id", "dim1", "dim2", "dim3", "dim4"]
    ).to_csv(tmp_path / "embedding.tsv", sep="\t", index=False)
    pd.DataFrame(
        maps, columns=["source_cell_id", "target_cell_id"]
    ).to_csv(tmp_path / "map.csv", index=False)
    pd.DataFrame(labels).to_csv(tmp_path / "labels.csv", header=False, index=False)
    pd.DataFrame(groups).to_csv(tmp_path / "groups.csv", header=False, index=False)
    (tmp_path / "eligible.txt").write_text("\n".join(eligible) + "\n")

    output = tmp_path / "scores"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path.cwd())
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/score_reference_embedding.py",
            "--embedding",
            str(tmp_path / "embedding.tsv"),
            "--cell-map",
            str(tmp_path / "map.csv"),
            "--labels",
            str(tmp_path / "labels.csv"),
            "--groups",
            str(tmp_path / "groups.csv"),
            "--eligible-cells",
            str(tmp_path / "eligible.txt"),
            "--name",
            "reference",
            "--output-dir",
            str(output),
            "--folds",
            "3",
            "--silhouette-sample-size",
            "40",
        ],
        check=True,
        env=environment,
    )
    report = json.loads((output / "reference_label_scores.json").read_text())
    assert report["status"] == "ok"
    assert report["standard_analysis_baseline"] is True
    assert report["n_published_cells"] == 60
    assert report["n_mapped_eligible_cells"] == 48
    assert report["n_scored_cells"] == 48
