import numpy as np
import pytest

from tealeaf.sc.path_score_mixed import binary_information_geometry


def test_equal_information_and_shape_are_balanced():
    result = binary_information_geometry(np.ones(5), np.ones(5), np.ones(5))
    assert result["geometry_class"] == "balanced"
    assert result["grid_maximum_precision_share"] == pytest.approx(.2)
    assert result["grid_minimum_effective_subjects"] == pytest.approx(5.)


def test_endpoint_dominance_and_input_nonmutation():
    information, shape, reference = np.array([1e6, 1., 1., 1., 1.]), np.ones(5), np.full(5, 1e6)
    saved = [value.copy() for value in (information, shape, reference)]
    result = binary_information_geometry(information, shape, reference)
    assert result["geometry_class"] == "dominated"
    assert result["grid_maximum_precision_share"] > .9999
    assert result["biological_limit_maximum_share"] == pytest.approx(.2)
    for value, before in zip((information, shape, reference), saved):
        np.testing.assert_array_equal(value, before)
    result = binary_information_geometry(np.ones(5), [1e-7, 1., 1., 1., 1.], np.ones(5))
    assert result["biological_limit_maximum_share"] > .9999


@pytest.mark.parametrize("scale", [1e-90, 1e90])
def test_geometry_is_invariant_to_global_coordinate_units(scale):
    info, shape, reference = np.array([2., 5., 10., 4., 7.]), np.array([4., 9., 6., 8., 3.]), np.full(5, 20.)
    before = binary_information_geometry(info, shape, reference)
    after = binary_information_geometry(info / scale**2, shape * scale**2, reference / scale**2)
    assert after["geometry_class"] == before["geometry_class"]
    for field in ("grid_maximum_precision_share", "grid_minimum_effective_subjects", "measurement_maximum_share", "biological_limit_maximum_share"):
        assert after[field] == pytest.approx(before[field], rel=1e-12)


def test_rank_exclusions_match_relative_information_rule():
    result = binary_information_geometry([1., 2., 3., 4., 1e-14, 0.], np.ones(6), [4., 4., 4., 4., 1., 0.])
    assert result["n_informative_subjects"] == 4
    assert result["measurement_maximum_share"] == pytest.approx(.4)


@pytest.mark.parametrize("info,shape,reference", [([1.] * 3, [1.] * 3, [1.] * 3), ([1.] * 4, [0.] * 4, [1.] * 4), ([2.] * 4, [1.] * 4, [1.] * 4), ([np.nan] * 4, [1.] * 4, [1.] * 4)])
def test_invalid_geometry_is_rejected(info, shape, reference):
    with pytest.raises(ValueError):
        binary_information_geometry(info, shape, reference)
