import numpy as np
import pytest

from tealeaf.sc.local_read_null import local_read_pattern_null


def counts():
    return np.array([[[[10, 0], [5, 0]], [[0, 3], [0, 7]]], [[[3, 7], [6, 4]], [[0, 0], [0, 0]]]])


def test_pattern_null_keeps_every_original_stratum_depth_and_boundary():
    original = counts()
    for sd in (0., .8):
        result = local_read_pattern_null(original, np.random.default_rng(17), sd)
        np.testing.assert_array_equal(result.sum(axis=-1), original.sum(axis=-1))
        np.testing.assert_array_equal(result[0], original[0])
        assert not result[1, 1].any()


def test_pattern_null_is_reproducible_with_fresh_not_fixed_class_margins():
    original = np.full((20, 2, 2, 2), 10)
    first = local_read_pattern_null(original, np.random.default_rng(11), .8)
    second = local_read_pattern_null(original, np.random.default_rng(11), .8)
    np.testing.assert_array_equal(first, second)
    assert not np.array_equal(first.sum(axis=2), original.sum(axis=2))


def test_subject_slope_is_shared_across_primers_and_mean_zero_across_draws():
    original = np.full((4000, 2, 2, 2), 100)
    result = local_read_pattern_null(original, np.random.default_rng(3), .8)
    change = result[:, :, 1, 0] / 200 - result[:, :, 0, 0] / 200
    assert np.corrcoef(change.T)[0, 1] > .9
    assert abs(change.mean()) < .01


@pytest.mark.parametrize('bad', [-1, .5, np.nan])
def test_pattern_null_rejects_invalid_counts(bad):
    original = counts().astype(float)
    original[0, 0, 0, 0] = bad
    with pytest.raises(ValueError):
        local_read_pattern_null(original, np.random.default_rng(1))
