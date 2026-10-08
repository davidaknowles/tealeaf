import numpy as np
import pytest

from tealeaf.sc.ec_diagnostics import pooled_ec_mixture_diagnostics


def test_convex_pooling_and_exact_categorical_model_have_no_mismatch():
    mapping = np.array([[1., 0.], [0., 1.], [1., 1.]])
    counts = np.array([[45., 5., 50.], [5., 45., 50.]])
    result = pooled_ec_mixture_diagnostics(counts, mapping)
    assert result["KL"] < 1e-12
    assert result["KL_lower_bound"] == 0
    assert result["log_probability_bound"] == 0


def test_shared_category_mismatch_cannot_be_fixed_by_transcript_mixture():
    mapping = np.array([[1., 0.], [0., 1.], [1., 1.]])
    result = pooled_ec_mixture_diagnostics([1., 1., 98.], mapping)
    exact = .02 * np.log(.02 / .5) + .98 * np.log(.98 / .5)
    np.testing.assert_allclose(result["KL"], exact, atol=1e-12)
    assert result["KL_lower_bound"] <= exact
    assert result["KL_lower_bound"] > exact - 1e-10
    assert result["log_probability_bound"] < -40


def test_dual_lower_bound_is_valid_even_if_optimizer_is_unfinished():
    mapping = np.array([[1., 0.], [0., 1.], [1., 1.]])
    counts = np.array([70., 20., 10.])
    result = pooled_ec_mixture_diagnostics(counts, mapping, max_iter=1)
    exact = .9 * np.log(.9 / .5) + .1 * np.log(.1 / .5)
    assert result["KL_lower_bound"] <= exact + 1e-12
    assert result["KL"] >= exact - 1e-12


def test_fractional_counts_do_not_receive_an_integer_sampling_bound():
    result = pooled_ec_mixture_diagnostics([1.5, 1.5, 97.], [[1., 0.], [0., 1.], [1., 1.]])
    assert not result["integral_counts"] and np.isnan(result["log_probability_bound"])
    assert result["KL_lower_bound"] > 0
    result = pooled_ec_mixture_diagnostics([[.5, .5, 49.], [.5, .5, 49.]], [[1., 0.], [0., 1.], [1., 1.]])
    assert not result["integral_counts"] and np.isnan(result["log_probability_bound"])


def test_negative_observation_counts_cannot_cancel_when_pooled():
    with pytest.raises(ValueError, match="nonnegative"):
        pooled_ec_mixture_diagnostics([[-1., 3.], [2., 0.]], np.eye(2))


def test_single_supported_transcript_uses_its_fixed_categorical_distribution():
    result = pooled_ec_mixture_diagnostics([30., 15.], [[2.], [1.]])
    assert result["fit_converged"] and result["KL"] == 0
    result = pooled_ec_mixture_diagnostics([30., 0.], [[2.], [1.]])
    np.testing.assert_allclose(result["KL_lower_bound"], np.log(1.5), atol=2e-12)


def test_type_bound_does_not_require_identically_distributed_molecules():
    mapping = np.array([[.1, .3], [.9, .7]])
    distribution = np.ones(1)
    for probability in [.1, .3] * 20:
        distribution = np.convolve(distribution, [1 - probability, probability])
    result = pooled_ec_mixture_diagnostics([32., 8.], mapping)
    # Every molecule has a represented composition, but compositions differ.
    # At this threshold the lower-tail count vectors are not as discrepant
    # as 32/40, so the discrepancy event is exactly this upper tail.
    upper_tail = distribution[32:].sum()
    assert result["log_probability_bound"] < 0
    assert upper_tail <= np.exp(result["log_probability_bound"])


def test_absent_primer_and_invalid_support_are_not_rejections():
    result = pooled_ec_mixture_diagnostics([0., 0.], np.eye(2))
    assert not result["positive_counts"] and np.isnan(result["log_probability_bound"])
    with pytest.raises(ValueError, match="represented"):
        pooled_ec_mixture_diagnostics([1., 0.], np.zeros((2, 2)))
