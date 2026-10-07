import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import betaln, expit, roots_jacobi

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_marginal import BinaryECPathLikelihood, beta_quadrature, binary_marginal_test, prepare_binary_ec_likelihood, random_subject_objective


@pytest.mark.parametrize("a,b", [(1, 1), (.1, .9), (.01, 2), (20, 80), (50000, 50000)])
def test_beta_rule_moments(a, b):
    nodes, weights = beta_quadrature(a, b, 16)
    assert weights.sum() == pytest.approx(1)
    assert np.all(nodes > 0) and np.all(nodes < 1)
    assert weights @ nodes == pytest.approx(a / (a + b), abs=2e-13)
    assert weights @ nodes ** 2 == pytest.approx(a * (a + 1) / ((a + b) * (a + b + 1)), abs=2e-13)
    if a + b < 1000:
        reference, mass = roots_jacobi(16, b - 1, a - 1)
        assert np.allclose(nodes, (reference + 1) / 2, atol=1e-13)
        assert np.allclose(weights, mass / mass.sum(), atol=1e-12)


def binomial_likelihood(counts):
    counts = np.asarray(counts, dtype=float)
    return BinaryECPathLikelihood((counts,), (np.array([[0., 1., 0.], [0., 0., 1.]]),), np.arange(len(counts)), np.zeros(len(counts)), counts + .125, np.zeros(len(counts)))


@pytest.mark.parametrize("counts,mean,precision", [([3, 7], .3, 20), ([0, 20], .08, 20), ([10, 190], .1, 20), ([150, 150], .1, 10000)])
@pytest.mark.parametrize("numerical", [False, True])
def test_binomial_integral_agrees_with_exact_beta_binomial(counts, mean, precision, numerical):
    likelihood = binomial_likelihood([counts])
    if numerical:
        likelihood.binomial_terms = np.full((1, 3), np.nan)
    actual = likelihood.log_integrals(np.array([[np.log(mean / (1 - mean))]]), precision, nodes=64)[0, 0]
    a, b = mean * precision, (1 - mean) * precision
    expected = betaln(a + counts[0], b + counts[1]) - betaln(a, b)
    assert actual == pytest.approx(expected, abs=2e-5)


def test_shared_latent_composition_is_not_integrated_separately_for_primers():
    likelihood = binomial_likelihood([[3, 7]])
    likelihood.counts = (*likelihood.counts, np.array([[8., 2.]]))
    likelihood.components = (*likelihood.components, likelihood.components[0])
    actual = likelihood.log_integrals(np.array([[0.]]), 20)[0, 0]
    expected = betaln(21, 19) - betaln(10, 10)
    separate = betaln(13, 17) + betaln(18, 12) - 2 * betaln(10, 10)
    assert actual == pytest.approx(expected)
    assert abs(actual - separate) > .1


def test_exact_likelihood_with_uninformative_ec_and_outside_mass():
    likelihood = BinaryECPathLikelihood((np.array([[3., 7., 8.]]),), (np.array([[0, .7, 0], [0, 0, .7], [.3, 0, 0]]),), np.array([0]), np.array([0]), np.array([[1., 1.]]), np.zeros(1))
    actual = likelihood.log_integrals(np.array([[0.]]), 20)[0, 0]
    expected = betaln(13, 17) - betaln(10, 10) + 10 * np.log(.7) + 8 * np.log(.3)
    assert actual == pytest.approx(expected)


def test_nontrivial_two_primer_ec_integral_against_adaptive_integral():
    components = (np.array([[.03, .7, .05], [.1, .08, .6], [.02, .2, .25]]), np.array([[.02, .3, .1], [.07, .1, .5]]))
    likelihood = BinaryECPathLikelihood((np.array([[2., 4., 2.]]), np.array([[3., 5.]])), components, np.array([0]), np.array([0]), np.array([[4., 9.]]), np.zeros(1))
    mean, precision = .3, 20
    a, b = precision * mean, precision * (1 - mean)
    value = quad(lambda p: np.exp(likelihood.row_log_likelihood(0, np.array([p]))[0] + (a - 1) * np.log(p) + (b - 1) * np.log1p(-p) - betaln(a, b)), 0, 1, epsabs=1e-25, epsrel=1e-10)[0]
    actual = likelihood.log_integrals(np.array([[np.log(mean / (1 - mean))]]), precision, nodes=32)[0, 0]
    assert actual == pytest.approx(np.log(value), abs=1e-6)


def test_flat_likelihood_integrates_exactly_one():
    likelihood = binomial_likelihood([[0, 0]])
    assert np.allclose(likelihood.log_integrals(np.array([[-8., 0., 8.]]), .1), 0, atol=1e-13)


def test_random_subject_score_against_independent_parameter_differences():
    likelihood = binomial_likelihood([[2, 10], [8, 12], [5, 7], [14, 6]])
    likelihood.subjects = np.array([0, 0, 1, 1])
    design = np.column_stack([np.ones(4), [0, 1, 0, 1]])
    parameters = np.array([-.4, .5, np.log(20), np.log(.4)])
    value, gradient = random_subject_objective(parameters, likelihood, design, subject_nodes=8, path_nodes=24)
    assert np.isfinite(value)
    for index in range(len(parameters)):
        step = np.zeros_like(parameters)
        step[index] = 1e-4
        plus = random_subject_objective(parameters + step, likelihood, design, subject_nodes=8, path_nodes=24, gradient=False)
        minus = random_subject_objective(parameters - step, likelihood, design, subject_nodes=8, path_nodes=24, gradient=False)
        assert gradient[index] == pytest.approx((plus - minus) / 2e-4, abs=2e-5)


def test_preparation_retains_separate_primer_counts_and_aggregates_replicates():
    counts = np.array([[3, 7], [2, 8], [6, 4], [4, 6]])
    data = ECGLMMData((counts, counts * 2), (np.eye(2), np.eye(2)), np.ones((4, 1)), np.array([0, 0, 1, 1]))
    result = prepare_binary_ec_likelihood(data, [0, 1], [0, 0, 1, 1], [0, 0, 1, 1], baseline=np.ones(2) / 2)
    assert np.array_equal(result.counts[0], [[5, 15], [10, 10]])
    assert np.array_equal(result.counts[1], [[10, 30], [20, 20]])
    assert np.array_equal(result.labels, [0, 1])
    assert np.allclose(result.components[0], [[0, 1, 0], [0, 0, 1]])
    empty = ECGLMMData((np.zeros_like(counts),), (np.eye(2),), data.design, data.clusters)
    with pytest.raises(ValueError, match="positive-count"):
        prepare_binary_ec_likelihood(empty, [0, 1], [0, 0, 1, 1], [0, 0, 1, 1], baseline=np.ones(2) / 2)


def test_binary_marginal_test_detects_and_reverses_its_own_usage_effect():
    rng = np.random.default_rng(888)
    labels = np.tile([0, 1], 16)
    subjects = np.repeat(np.arange(16), 2)
    mean = expit(-.4 + rng.normal(0, .3, 16)[subjects] + 1.3 * labels)
    probability = rng.beta(30 * mean, 30 * (1 - mean))
    first = rng.binomial(80, probability)
    likelihood = binomial_likelihood(np.column_stack([first, 80 - first]))
    likelihood.subjects, likelihood.labels = subjects, labels
    result = binary_marginal_test(likelihood, subject_nodes=96)
    reverse = binary_marginal_test(likelihood, labels=1 - labels, subject_nodes=96)
    assert result["converged"] and reverse["converged"]
    assert result["p_value"] < 1e-5
    assert result["degrees_of_freedom"] == 1
    assert result["standardized_means"][1, 0] > result["standardized_means"][0, 0]
    assert np.allclose(result["standardized_means"], reverse["standardized_means"][::-1], atol=2e-4)
    assert result["statistic"] == pytest.approx(reverse["statistic"], abs=1e-3)


def test_invalid_ec_likelihood_is_rejected():
    with pytest.raises(ValueError, match="nonzero model probability"):
        BinaryECPathLikelihood((np.ones((1, 3)),), (np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]]),), np.array([0]), np.array([0]), np.ones((1, 2)), np.zeros(1))


def test_mixed_ec_subject_score_against_parameter_differences():
    likelihood = BinaryECPathLikelihood((np.array([[3., 7.], [8., 2.]]),), (np.array([[.1, .7, .05], [.05, .1, .7]]),), np.array([0, 0]), np.array([0, 1]), np.array([[2., 6.], [6., 2.]]), np.zeros(2))
    design = np.array([[1., 0.], [1., 1.]])
    parameters = np.array([-.3, .6, np.log(20), np.log(.3)])
    _, score = random_subject_objective(parameters, likelihood, design, subject_nodes=6, path_nodes=24)
    for index in range(len(parameters)):
        step = np.zeros_like(parameters)
        step[index] = 1e-4
        plus = random_subject_objective(parameters + step, likelihood, design, subject_nodes=6, path_nodes=24, gradient=False)
        minus = random_subject_objective(parameters - step, likelihood, design, subject_nodes=6, path_nodes=24, gradient=False)
        assert score[index] == pytest.approx((plus - minus) / 2e-4, abs=2e-5)


def test_common_objective_error_cannot_cancel_in_the_lr_gate(monkeypatch):
    from types import SimpleNamespace
    from tealeaf.sc import path_marginal
    likelihood = binomial_likelihood(np.tile([5, 5], (8, 1)))
    likelihood.subjects = np.repeat(np.arange(4), 2)
    likelihood.labels = np.tile([0, 1], 4)
    monkeypatch.setattr(path_marginal, "minimize", lambda objective, initial, **kwargs: SimpleNamespace(x=initial, fun=0., success=True))
    monkeypatch.setattr(path_marginal, "random_subject_objective", lambda *args, **kwargs: 2.)
    result = binary_marginal_test(likelihood)
    assert not result["converged"]
    assert result["p_value"] == 1 and result["statistic"] == 0
    assert np.isnan(result["standardized_means"]).all()
    assert result["null_quadrature_error"] == 2
    assert result["alternative_quadrature_error"] == 2


def test_singular_beta_shapes_do_not_use_unresolved_log_moments_as_scores():
    likelihood = BinaryECPathLikelihood((np.array([[3., 7.]]),), (np.array([[.1, .7, .05], [.05, .1, .7]]),), np.array([0]), np.array([0]), np.array([[2., 6.]]), np.zeros(1))
    logits, concentration = np.array([[np.log(.999 / .001)]]), .1
    _, eta, precision = likelihood.log_integrals(logits, concentration, nodes=32, gradient=True)
    step = 1e-3
    expected_eta = (likelihood.log_integrals(logits + step, concentration, nodes=32) - likelihood.log_integrals(logits - step, concentration, nodes=32)) / (2 * step)
    expected_precision = (likelihood.log_integrals(logits, concentration * np.exp(step), nodes=32) - likelihood.log_integrals(logits, concentration * np.exp(-step), nodes=32)) / (2 * step)
    assert np.allclose(eta, expected_eta, atol=1e-10)
    assert np.allclose(precision, expected_precision, atol=1e-10)


def test_exact_integral_score_respects_numerical_mean_clipping():
    likelihood = binomial_likelihood([[3, 7]])
    _, eta, _ = likelihood.log_integrals(np.array([[-50., 50.]]), 20., gradient=True)
    assert np.array_equal(eta, [[0., 0.]])


def test_primer_with_one_invisible_path_has_constant_conditional_likelihood():
    counts = (np.array([[3., 7.]]), np.array([[5., 5.]]))
    components = (np.array([[0, 1, 0], [0, 2, 0]], float), np.array([[0, 1, 0], [0, 0, 1]], float))
    likelihood = BinaryECPathLikelihood(counts, components, np.array([0]), np.array([0]), np.ones((1, 2)), np.zeros(1))
    expected = betaln(15, 15) - betaln(10, 10) + 3 * np.log(1 / 3) + 7 * np.log(2 / 3)
    assert likelihood.log_integrals(np.array([[0.]]), 20)[0, 0] == pytest.approx(expected, abs=1e-12)
    assert np.isfinite(likelihood.row_log_likelihood(0, np.array([1e-14, .5, 1 - 1e-14]))).all()


def test_zero_map_with_zero_primer_counts_does_not_change_integral():
    from tealeaf.sc.path_marginal_quadrature import shared_prior_quadrature

    counts = (np.array([[3., 7.]]), np.zeros((1, 2)))
    components = (np.array([[0, 1, 0], [0, 0, 1]], float), np.zeros((2, 3)))
    likelihood = BinaryECPathLikelihood(counts, components, np.array([0]), np.array([0]), np.ones((1, 2)), np.zeros(1))
    assert shared_prior_quadrature(likelihood).log_integrals(np.array([[0.]]), 20)[0, 0] == pytest.approx(betaln(13, 17) - betaln(10, 10))


def test_proportional_endpoint_maps_and_empty_primer_are_constant_on_numeric_grid():
    from tealeaf.sc.path_marginal_quadrature import matrix_log_likelihood, shared_prior_quadrature

    counts = (np.array([[3., 7.]]), np.zeros((1, 2)))
    components = (np.array([[0, 1, 4], [0, 2, 8]], float), np.zeros((2, 3)))
    likelihood = BinaryECPathLikelihood(counts, components, np.array([0]), np.array([0]), np.ones((1, 2)), np.zeros(1))
    expected = 3 * np.log(1 / 3) + 7 * np.log(2 / 3)
    points = np.array([1e-14, .4, 1 - 1e-14])
    assert np.allclose(matrix_log_likelihood(likelihood, [0], points), expected)
    assert np.allclose(shared_prior_quadrature(likelihood).log_integrals(np.array([[-2., 0., 2.]]), 20), expected)
