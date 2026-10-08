import numpy as np
import pytest

from tealeaf.sc.path_score_mixed import binary_null_minor_usage


def test_binary_shape_recovers_both_orientations_and_extreme_small_root():
    small = np.array([1e-100, 1e-20, 1e-12, .0001, .1, .25, .5])
    shape = 1 / (small * (1 - small))
    np.testing.assert_allclose(binary_null_minor_usage(shape), small, rtol=1e-13)
    moderate = np.array([.5, .75, .9, .9999])
    np.testing.assert_allclose(binary_null_minor_usage(1 / (moderate * (1 - moderate))), 1 - moderate, rtol=1e-13)
    np.testing.assert_allclose(binary_null_minor_usage(4 * (1 - 1e-14)), .5)


@pytest.mark.parametrize("shape", [3., -1., 0., np.nan, np.inf])
def test_invalid_binary_shapes_fail(shape):
    with pytest.raises(ValueError):
        binary_null_minor_usage(shape)
