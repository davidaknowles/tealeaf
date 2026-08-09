"""Tests for estimate-once path-Wald diagnostic summaries."""

import json

import numpy as np
import pandas as pd

from extra_scripts.summarize_path_wald_diagnostics import (
    calibration_tables,
    null_split_sign_agreement,
    split_direction_tables,
)


def observed_row(test_id, effects, standard_errors, p_value):
    return {
        "test_id": test_id,
        "gene_id": f"gene_{test_id}",
        "statistic": 5.0,
        "p_value": p_value,
        "null_converged": True,
        "alternative_converged": True,
        "tested_levels_json": json.dumps(["a", "b"]),
        "path_effects_json": json.dumps([effects]),
        "path_standard_errors_json": json.dumps([standard_errors]),
    }


def test_calibration_tables_count_null_rejections_and_empirical_rank():
    observed = pd.DataFrame([
        observed_row("x", [1.0, -1.0], [0.2, 0.2], 0.01),
        observed_row("y", [0.5, -0.5], [0.2, 0.2], 0.02),
    ])
    null = pd.DataFrame([
        {
            "test_id": test_id,
            "null_type": "label_permutation",
            "replicate": replicate,
            "statistic": statistic,
            "p_value": p_value,
        }
        for replicate, p_value in enumerate([0.01, 0.4])
        for test_id, statistic in [("x", 2.0), ("y", 7.0)]
    ])
    _, families, summary, empirical = calibration_tables(observed, null)
    assert len(families) == 2
    assert summary.loc[0, "pooled_p_le_0.05"] == 0.5
    ranked = empirical.set_index("test_id").empirical_p_value
    assert ranked["x"] == 1.0 / 3.0
    assert ranked["y"] == 1.0


def test_split_direction_tables_and_null_baseline_track_path_signs():
    fold0 = pd.DataFrame([
        observed_row("x", [1.0, -1.0], [0.2, 0.2], 0.01),
        observed_row("y", [0.5, -0.5], [0.2, 0.2], 0.2),
    ])
    fold1 = pd.DataFrame([
        observed_row("x", [0.8, -0.7], [0.2, 0.2], 0.02),
        observed_row("y", [-0.4, 0.3], [0.2, 0.2], 0.3),
    ])
    coordinates, _, summary = split_direction_tables(fold0, fold1)
    assert coordinates.sign_agreement.tolist() == [True, True, False, False]
    assert summary.loc[summary.subset.eq("all coordinates"), "sign_agreement"].iloc[0] == 0.5
    null0 = pd.DataFrame({
        "test_id": ["x"],
        "null_type": ["label_permutation"],
        "replicate": [0],
        "tested_levels_json": [json.dumps(["a", "b"])],
        "path_effects_json": [json.dumps([[1.0, -1.0]])],
    })
    null1 = null0.copy()
    replicate, null_summary = null_split_sign_agreement(
        null0, null1, {"x"}
    )
    assert replicate.loc[0, "sign_agreement"] == 1.0
    assert null_summary.loc[0, "mean_sign_agreement"] == 1.0
