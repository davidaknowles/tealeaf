import numpy as np
import pytest

from tealeaf.sc.path_score_mixed import mixed_score_test
from tealeaf.sc.score_variance_pooling import scalar_score_panel, fit_shared_biological_variance, scalar_precision_share_bounds


def example_panel():
    rng = np.random.default_rng(23576)
    groups = np.repeat(np.arange(5), [4, 6, 8, 5, 9])
    information = np.exp(rng.normal(size=len(groups)))
    shape = rng.uniform(4., 30., len(groups))
    scores = information * rng.normal(size=len(groups))
    return scores, information, shape, information * 2, groups


@pytest.mark.parametrize("variance", [0., .0001, .1, 30.])
def test_ragged_fixed_variance_matches_existing_kernel(variance):
    arrays = example_panel()
    before = [value.copy() for value in arrays]
    panel = scalar_score_panel(*arrays)
    result = panel.evaluate(variance)
    scores, information, shape, reference, groups = arrays
    for index in range(len(panel.n_subjects)):
        selected = groups == index
        original = mixed_score_test(scores[selected, None], information[selected, None, None], shape[selected, None, None], reference_information=reference[selected, None, None], biological_variance=variance, scalar_fast=True)
        for name in ("p_value", "statistic", "residual_inflation", "restricted_objective"):
            assert result[name][index] == pytest.approx(original[name], rel=3e-6, abs=1e-8)
        assert result["mean_difference"][index] == pytest.approx(original["mean_difference"][0], rel=3e-6, abs=1e-8)
    for value, saved in zip(arrays, before):
        np.testing.assert_array_equal(value, saved)


def test_rank_exclusion_and_signs_retain_original_source_positions():
    arrays = list(example_panel())
    arrays[1][5] = 1e-14
    arrays[3][5] = 1.
    panel = scalar_score_panel(*arrays)
    assert 5 not in panel.source_positions
    assert panel.n_subjects[1] == 5
    result = panel.evaluate(.1)
    reversed_result = panel.evaluate(.1, values=-panel.values)
    np.testing.assert_array_equal(result["p_value"], reversed_result["p_value"])
    np.testing.assert_array_equal(result["mean_difference"], -reversed_result["mean_difference"])


@pytest.mark.parametrize("scale", [1e-30, 1e30])
def test_pooled_variance_and_test_are_invariant_to_coordinate_units(scale):
    scores, information, shape, reference, groups = example_panel()
    original = scalar_score_panel(scores, information, shape, reference, groups)
    changed = scalar_score_panel(scores / scale, information / scale**2, shape * scale**2, reference / scale**2, groups)
    fitted, refitted = fit_shared_biological_variance(original), fit_shared_biological_variance(changed)
    assert refitted["biological_variance"] == pytest.approx(fitted["biological_variance"], rel=1e-5)
    a, b = original.evaluate(fitted["biological_variance"]), changed.evaluate(fitted["biological_variance"])
    np.testing.assert_allclose(a["p_value"], b["p_value"], rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(a["mean_difference"] * scale, b["mean_difference"], rtol=1e-10)


def test_common_variance_recovers_known_simulation_with_free_event_means():
    rng = np.random.default_rng(50772)
    groups = np.repeat(np.arange(500), 8)
    information, shape = np.full(len(groups), 100.), np.full(len(groups), 4.)
    means = rng.normal(0., 2., 500)
    values = means[groups] + rng.normal(0., np.sqrt(.01 + .05 * 4), len(groups))
    panel = scalar_score_panel(values * information, information, shape, information * 2, groups)
    fit = fit_shared_biological_variance(panel)
    assert fit["biological_variance"] == pytest.approx(.05, abs=.005)
    assert not fit["near_search_boundary"]
    assert fit["mean_restricted_objective"] < fit["zero_variance_objective"]


def test_unidentifiable_or_incomplete_panels_are_rejected():
    arrays = list(example_panel())
    arrays[2] = np.zeros(len(arrays[0]))
    with pytest.raises(ValueError, match="unidentified"):
        fit_shared_biological_variance(scalar_score_panel(*arrays))
    with pytest.raises(ValueError, match="consecutive integer"):
        scalar_score_panel(*arrays[:-1], arrays[-1] + 1)
    arrays = list(example_panel())
    arrays[1][0] = 0.
    with pytest.raises(ValueError, match="four informative"):
        scalar_score_panel(*arrays)


def test_precision_bounds_cover_endpoints_and_every_sampled_variance():
    rng = np.random.default_rng(24771)
    information, shape = np.exp(rng.normal(0., 3., 12)), np.exp(rng.normal(0., 3., 12))
    lower, upper = scalar_precision_share_bounds(information, shape)
    for variance in np.r_[0., np.exp(np.linspace(-40., 40., 300))]:
        weights = information / (1 + variance * information * shape)
        shares = weights / weights.sum()
        assert (shares >= lower - 1e-14).all() and (shares <= upper + 1e-14).all()
    limiting = (1 / shape) / np.sum(1 / shape)
    assert (limiting >= lower - 1e-14).all() and (limiting <= upper + 1e-14).all()


def test_one_donor_can_dominate_at_every_possible_biological_variance():
    information, shape = np.array([1e6, 1., 1., 1., 1.]), np.array([1., 1e6, 1e6, 1e6, 1e6])
    lower, upper = scalar_precision_share_bounds(information, shape)
    assert lower[0] > .9999
    np.testing.assert_allclose(lower, upper, rtol=1e-12)
    scaled = scalar_precision_share_bounds(information / 1e100, shape * 1e100)
    np.testing.assert_allclose(scaled[0], lower, rtol=1e-12)
    np.testing.assert_allclose(scaled[1], upper, rtol=1e-12)
