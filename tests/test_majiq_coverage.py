import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from tealeaf.sc.majiq_coverage import coverage_null_contrast, drop_coverage_null_calls, mann_whitney_rows


def test_mann_whitney_rows_matches_scipy_with_ties_and_missing():
    rng = np.random.default_rng(1)
    values = rng.integers(0, 5, size=(20, 14)).astype(float)
    values[3, [0, 9]] = np.nan
    in_a = np.arange(14) < 6
    p = mann_whitney_rows(values, in_a)
    for row in range(20):
        x, y = values[row, in_a], values[row, ~in_a]
        expected = mannwhitneyu(x[~np.isnan(x)], y[~np.isnan(y)], method="asymptotic", use_continuity=True).pvalue
        assert np.isclose(p[row], expected)


def test_coverage_only_null_is_called_when_coverage_differs():
    # two connections of a binary LSV with PSI .9 in every sample; group a has 3 reads, group b 300
    rng = np.random.default_rng(0)
    samples = 40
    in_a = np.arange(samples) < 20
    total = np.where(in_a, 3., 300.)[None, :].repeat(2, 0)
    cover = np.vstack([np.round(.9 * total[0]), np.round(.1 * total[1])])
    passed = np.ones_like(total, bool)
    alpha = np.full(2, .5)
    table = coverage_null_contrast(total, cover, passed, alpha, in_a, .5, rng)
    assert (table.p_null0 < .01).all() and (table.p_obs < 1e-4).all()
    balanced = coverage_null_contrast(np.full_like(total, 300.), np.vstack([np.full(samples, 270.), np.full(samples, 30.)]), passed, alpha, in_a, .5, rng)
    assert (balanced.p_null0 > .001).all()


def test_drop_coverage_null_calls_keeps_unflagged_and_unmatched():
    tests = pd.DataFrame({"contrast_id": ["c"] * 3, "feature_id": ["f1", "f2", "f3"], "p_value": [.001, .002, .003]})
    null = pd.DataFrame({"contrast_id": ["c", "c"], "feature_id": ["f1", "f2"], "p_null0": [.01, .5], "p_null1": [.02, .6], "p_null2": [.2, .4]})
    kept = drop_coverage_null_calls(tests, null)
    assert kept.feature_id.tolist() == ["f2", "f3"]
