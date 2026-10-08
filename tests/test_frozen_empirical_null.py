import numpy as np
import pytest

from tealeaf.sc.empirical_null import leave_parent_out_cdf


def test_repeated_parent_draws_and_ties_exclude_entire_parent():
    result, counts = leave_parent_out_cdf([.1, .2, .1, .1], ["a", "a", "b", "new"], ["s"] * 4, [.1, .1, .2, .1, .3], ["a", "a", "a", "b", "b"], ["s"] * 5)
    np.testing.assert_array_equal(result, [2 / 3, 2 / 3, 3 / 4, 4 / 6])
    np.testing.assert_array_equal(counts, [2, 2, 3, 5])


def test_target_geometry_does_not_reassign_training_parent():
    result, counts = leave_parent_out_cdf([.01, .01, .01], ["a"] * 3, ["balanced", "dominated", "missing"], [.001, .2, .3], ["a", "b", "c"], ["dominated", "balanced", "balanced"])
    np.testing.assert_array_equal(result, [1 / 3, 1., 1.])
    np.testing.assert_array_equal(counts, [2, 0, 0])


def test_matches_physically_excluded_random_training_without_mutation():
    rng = np.random.default_rng(4172)
    training = rng.uniform(size=700)
    parent = np.array([f"p{value}" for value in rng.integers(0, 19, 700)])
    stratum = np.array([f"s{value}" for value in rng.integers(0, 4, 700)])
    target = rng.uniform(size=40)
    target_parent, target_stratum = parent[:40].copy(), stratum[:40].copy()
    before = [value.copy() for value in (target, target_parent, target_stratum, training, parent, stratum)]
    result, counts = leave_parent_out_cdf(target, target_parent, target_stratum, training, parent, stratum)
    for index in range(len(target)):
        values = training[(parent != target_parent[index]) & (stratum == target_stratum[index])]
        assert counts[index] == len(values)
        assert result[index] == (1 + (values <= target[index]).sum()) / (1 + len(values))
    for value, original in zip((target, target_parent, target_stratum, training, parent, stratum), before):
        np.testing.assert_array_equal(value, original)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1., 2.])
def test_nonprobabilities_rejected(bad):
    with pytest.raises(ValueError):
        leave_parent_out_cdf([bad], ["a"], ["s"], [.2], ["b"], ["s"])
