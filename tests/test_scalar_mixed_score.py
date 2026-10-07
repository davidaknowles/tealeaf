import numpy as np
import pytest

from tealeaf.sc.path_score_mixed import mixed_score_test
from tealeaf.sc.path_score_mixed import PathScoreComponents, aggregate_path_scores, signed_path_score_p_value


@pytest.mark.parametrize("reference", [False, True])
@pytest.mark.parametrize("variance", [None, 0., .3])
def test_scalar_specialization_matches_generic(reference, variance):
    rng = np.random.default_rng(824)
    information = np.exp(rng.normal(0, 3, 12))[:, None, None]
    shapes = np.exp(rng.normal(0, 2, 12))[:, None, None]
    scores = rng.normal(.4, 1, 12)[:, None] * information[:, 0]
    target = information * rng.uniform(1, 10, (12, 1, 1)) if reference else None
    arguments = dict(biological_variance=variance, reference_information=target)
    generic = mixed_score_test(scores, information, shapes, **arguments)
    fast = mixed_score_test(scores, information, shapes, scalar_fast=True, **arguments)
    assert generic.keys() == fast.keys()
    for name in generic:
        np.testing.assert_allclose(fast[name], generic[name], rtol=3e-6, atol=1e-8, err_msg=name)


def test_scalar_preserves_relative_rank_and_no_information_failures():
    scores = np.full((5, 1), 1e-18)
    information = np.full((5, 1, 1), 1e-20)
    target = information * 2
    for fast in (False, True):
        with pytest.raises(ValueError, match="four informative"):
            mixed_score_test(scores, information, scalar_fast=fast)
        result = mixed_score_test(scores, information, reference_information=target, scalar_fast=fast)
        assert result["n_subjects"] == 5
        with pytest.raises(ValueError, match="exceeds"):
            mixed_score_test(scores, information, reference_information=information / 2, scalar_fast=fast)


def test_scalar_wrapper_preserves_null_rng_and_multivariate_fallback():
    rng = np.random.default_rng(195)
    for dimension in (1, 3):
        info = np.broadcast_to(np.eye(dimension), (6, dimension, dimension)).copy()
        scores = rng.normal(size=(6, dimension))
        components = PathScoreComponents(scores.copy(), info.copy(), info.copy(), np.arange(6), (0, 1), [], [], reference_information=2 * info)
        for metric in ("absolute", "reference"):
            generic = aggregate_path_scores(components, information_metric=metric)
            fast = aggregate_path_scores(components, information_metric=metric, scalar_fast=True)
            np.testing.assert_allclose(fast["p_value"], generic["p_value"], rtol=3e-6)
            np.testing.assert_array_equal(fast["differences"], generic["differences"])
            first_rng, second_rng = np.random.default_rng(92), np.random.default_rng(92)
            generic_null = signed_path_score_p_value(components, first_rng, information_metric=metric)
            fast_null = signed_path_score_p_value(components, second_rng, information_metric=metric, scalar_fast=True)
            np.testing.assert_allclose(fast_null, generic_null, rtol=3e-6)
            np.testing.assert_array_equal(first_rng.normal(size=10), second_rng.normal(size=10))
        np.testing.assert_array_equal(components.scores, scores)
