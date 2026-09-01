import json

import pandas as pd

from extra_scripts.estimate_path_smoothing import select_cross_fitted_alpha


def test_spike_slab_empirical_bayes_recovers_finite_slab(tmp_path):
    rows = []
    for gene in range(100):
        variable = gene % 5 == 0
        evidence = (
            {1.0: -1.0, 2.0: 0.0, 8.0: -3.0, 128.0: -10.0}
            if variable
            else {1.0: -6.0, 2.0: -5.0, 8.0: -3.0, 128.0: -4.0}
        )
        for alpha, value in evidence.items():
            rows.append({
                "test_id": f"test_{gene}",
                "gene": gene,
                "n_paths": 2 + gene % 2,
                "alpha": alpha,
                "mean_log_evidence": value,
                "mean_point_null_log_evidence": 0.0,
                "n_aggregates": 8,
                "path_pseudocount_scaling": "total",
                "path_prior_center": "uniform",
            })
    input_path = tmp_path / "evidence.tsv"
    output_path = tmp_path / "map.json"
    pd.DataFrame(rows).to_csv(input_path, sep="\t", index=False)
    select_cross_fitted_alpha(
        [input_path], output_path, folds=5, minimum_genes=3
    )
    records = json.loads(output_path.read_text())["records"]
    assert {record["alpha"] for record in records} == {2.0}
    assert all(0 < record["spike_weight"] < 1 for record in records)
    assert len(records) == 5
    assert json.loads(output_path.read_text())["selection_scope"] == "global"
    assert json.loads(output_path.read_text())["path_prior_center"] == "uniform"


def test_single_fold_selects_one_alpha_from_all_genes(tmp_path):
    rows = []
    for gene in range(30):
        for alpha, value in ((1.0, -2.0), (4.0, 0.0), (16.0, -1.0)):
            rows.append({
                "test_id": f"test_{gene}",
                "gene": f"gene_{gene}",
                "n_paths": 2,
                "alpha": alpha,
                "mean_log_evidence": value,
                "mean_point_null_log_evidence": -3.0,
                "n_aggregates": 4,
                "path_pseudocount_scaling": "per_path",
                "path_prior_center": "uniform",
            })
    input_path = tmp_path / "evidence.tsv"
    output_path = tmp_path / "map.json"
    pd.DataFrame(rows).to_csv(input_path, sep="\t", index=False)
    select_cross_fitted_alpha([input_path], output_path, folds=1)
    result = json.loads(output_path.read_text())
    assert result["folds"] == 1
    assert len(result["records"]) == 1
    assert result["records"][0]["alpha"] == 4.0
    assert result["records"][0]["training_genes"] == 30
