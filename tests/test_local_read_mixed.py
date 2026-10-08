import numpy as np
import pytest
from scipy import special

from tealeaf.sc.local_read_mixed import LocalReadMixed, local_read_mixed_test


def counts():
    return np.array([[[[2, 3], [4, 2]], [[3, 4], [4, 3]]], [[[1, 4], [2, 4]], [[2, 4], [3, 5]]], [[[3, 5], [4, 4]], [[1, 6], [2, 4]]], [[[4, 3], [2, 3]], [[2, 1], [2, 4]]]])


def brute_objective(observed, means, effect, baseline_sd, slope_sd, nodes=101):
    abscissae, weights = special.roots_hermitenorm(nodes)
    first, second = np.meshgrid(abscissae, abscissae, indexing="ij")
    log_weights = (np.log(weights)[:, None] + np.log(weights)[None, :] - np.log(2 * np.pi)).ravel()
    objective = 0.
    for subject in observed:
        kernel = np.zeros(nodes**2)
        for primer in range(observed.shape[1]):
            for level, sign in enumerate((-.5, .5)):
                predictor = means[primer] + baseline_sd * first.ravel() + sign * (effect + slope_sd * second.ravel())
                inc, exc = subject[primer, level]
                kernel -= inc * np.logaddexp(0., -predictor) + exc * np.logaddexp(0., predictor)
        objective -= special.logsumexp(kernel + log_weights)
    return objective


@pytest.mark.parametrize("first,second", [(0., 0.), (.7, 0.), (0., .8), (.6, .8)])
def test_adaptive_binomial_integral_matches_independent_dense_prior_quadrature(first, second):
    observed = counts()
    likelihood = LocalReadMixed(observed)
    actual = likelihood.marginal_objective([-.2, .1], .4, first, second, nodes=31, gradient=False)
    expected = brute_objective(observed, [-.2, .1], .4, first, second)
    np.testing.assert_allclose(actual, expected, atol=1e-8, rtol=1e-10)


def test_binomial_integral_derivatives_include_both_variances_and_primer_means():
    likelihood = LocalReadMixed(counts())
    parameters = np.array([-.2, .1, .4, np.log(.7), np.log(.8)])
    def objective(values):
        return likelihood.marginal_objective(values[:2], values[2], np.exp(values[3]), np.exp(values[4]), nodes=21, gradient=False)
    _, actual = likelihood.marginal_objective(parameters[:2], parameters[2], np.exp(parameters[3]), np.exp(parameters[4]), nodes=21)
    step = 1e-5
    expected = []
    for index in range(5):
        direction = np.eye(5)[index] * step
        expected.append((objective(parameters + direction) - objective(parameters - direction)) / (2 * step))
    np.testing.assert_allclose(actual, expected, atol=1e-6, rtol=1e-6)


def test_swapping_cell_types_reverses_effect_without_changing_likelihood():
    first = LocalReadMixed(counts())
    second = LocalReadMixed(counts()[:, :, ::-1, :])
    a, ga = first.marginal_objective([-.2, .1], .8, .7, .4, nodes=15)
    b, gb = second.marginal_objective([-.2, .1], -.8, .7, .4, nodes=15)
    np.testing.assert_allclose(a, b, atol=1e-10)
    np.testing.assert_allclose(ga, gb * [1, 1, -1, 1, 1], atol=1e-10)


def test_unobserved_subject_and_primer_add_no_likelihood_information():
    observed = counts()
    extra_subject = np.zeros((1, 2, 2, 2), dtype=int)
    extra_primer = np.zeros((4, 1, 2, 2), dtype=int)
    original = LocalReadMixed(observed)
    for other in (LocalReadMixed(np.concatenate([observed, extra_subject])), LocalReadMixed(np.concatenate([observed, extra_primer], axis=1))):
        first = original.marginal_objective([-.2, .1], .4, .7, .8, nodes=15)
        second = other.marginal_objective([-.2, .1], .4, .7, .8, nodes=15)
        np.testing.assert_allclose(first[0], second[0], atol=1e-12)
        np.testing.assert_allclose(first[1], second[1], atol=1e-12)
        assert other.n_paired_subjects == original.n_paired_subjects == 4
    np.testing.assert_array_equal(observed, counts())


def test_counted_zero_class_subjects_remain_in_unconditional_likelihood():
    observed = counts()
    observed[0, :, :, 0] = 0
    assert LocalReadMixed(observed).n_paired_subjects == 4


def test_primer_cell_type_confounding_is_not_rescued_by_latent_priors():
    observed = counts()
    observed[:, 0, 1] = 0
    observed[:, 1, 0] = 0
    with pytest.raises(ValueError, match="confounded"):
        LocalReadMixed(observed)


@pytest.mark.parametrize("problem", [np.ones((2, 2)), np.full((4, 2, 2, 2), .2), np.full((4, 2, 2, 2), -1), np.full((4, 2, 2, 2), np.inf)])
def test_noninteger_or_invalid_counts_fail(problem):
    with pytest.raises(ValueError):
        LocalReadMixed(problem)


def test_complete_nested_binomial_lrt_fits_signal_and_checks_both_likelihoods():
    observed = np.repeat(np.array([[[[8, 22], [20, 10]], [[5, 15], [12, 8]]]]), 8, axis=0)
    result = local_read_mixed_test(LocalReadMixed(observed), nodes=11)
    assert result["converged"]
    assert result["p_value"] < .01
    assert result["log_odds_effect"] > 0
    assert result["quadrature_error"] < 1e-3
    assert result["alternative_objective"] <= result["null_objective"]
