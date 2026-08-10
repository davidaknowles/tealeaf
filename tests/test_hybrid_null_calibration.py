"""Tests for hybrid EC-GLMM structured-null summaries."""

import pandas as pd

from analyses.summarize_hybrid_null_calibration import summarize_family


def test_family_summary_excludes_failed_fits_from_nominal_rates(tmp_path):
    directory = tmp_path / "permutation_1_merged"
    directory.mkdir()
    path = directory / "hybrid_results.tsv"
    pd.DataFrame(
        {
            "converged": [True, True, False],
            "lrt_p_value": [0.01, 0.5, 0.001],
            "bh_q_value": [0.03, 0.5, 1.0],
        }
    ).to_csv(path, sep="\t", index=False)
    _, summary = summarize_family(path)
    assert summary["family"] == "permutation_1_merged"
    assert summary["tests"] == 3
    assert summary["converged"] == 2
    assert summary["p_below_0.05"] == 0.5
    assert summary["bh_discoveries"] == 1
