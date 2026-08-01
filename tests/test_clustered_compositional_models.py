import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.merge_clustered_compositional_models import read_tables
from extra_scripts.run_clustered_compositional_models import summarize


def test_power_calibration_uses_matching_block_null_statistics():
    observed = pd.DataFrame()
    null = pd.DataFrame([
        {
            "block_id": "a",
            "method": "gee_exchangeable",
            "degrees_of_freedom": 1,
            "statistic": statistic,
            "p_value": 0.5,
            "converged": True,
        }
        for statistic in range(20)
    ] + [
        {
            "block_id": "b",
            "method": "gee_exchangeable",
            "degrees_of_freedom": 1,
            "statistic": 100 + statistic,
            "p_value": 0.5,
            "converged": True,
        }
        for statistic in range(20)
    ])
    power = pd.DataFrame([
        {
            "block_id": block_id,
            "method": "gee_exchangeable",
            "effect_size": 0.5,
            "degrees_of_freedom": 1,
            "statistic": 50.0,
            "p_value": 0.01,
            "converged": True,
        }
        for block_id in ("a", "b")
    ])

    _, calibrated, summary = summarize(observed, null, power)

    np.testing.assert_allclose(
        calibrated["calibrated_p_value"], [1.0 / 21.0, 1.0]
    )
    assert summary["power"]["gee_exchangeable"]["0.5"][
        "null_calibrated_power"
    ] == 0.5


def test_shard_merge_skips_empty_tables():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        empty = root / "empty"
        populated = root / "populated"
        empty.mkdir()
        populated.mkdir()
        (empty / "table.tsv").write_text("")
        pd.DataFrame({"value": [3]}).to_csv(
            populated / "table.tsv", sep="\t", index=False
        )

        result = read_tables([empty, populated], "table.tsv")

    assert result["value"].tolist() == [3]


def test_observed_results_receive_event_matched_permutation_p_values():
    observed = pd.DataFrame([{
        "block_id": "a",
        "method": "glmm_laplace",
        "statistic": 10.0,
        "p_value": 0.001,
        "converged": True,
    }])
    null = pd.DataFrame([
        {
            "block_id": "a",
            "method": "glmm_laplace",
            "statistic": statistic,
            "p_value": 0.5,
            "converged": True,
        }
        for statistic in range(20)
    ])

    calibrated, _, summary = summarize(observed, null, pd.DataFrame())

    assert calibrated.loc[0, "permutation_p_value"] == 11.0 / 21.0
    assert summary["observed"]["glmm_laplace"][
        "permutation_p_lt_0.05"
    ] == 0
