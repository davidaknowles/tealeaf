import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import betaln

from tealeaf.sc.path_marginal import BinaryECPathLikelihood
from tealeaf.sc.path_marginal_quadrature import matrix_log_likelihood, prior_quadrature, shared_prior_quadrature


@pytest.mark.parametrize("mean,precision", [(.3, 20), (.999, .1), (.99, 4), (.01, 4)])
def test_prior_rule_matches_independent_algebraic_weight_integral(mean, precision):
    source = BinaryECPathLikelihood((np.array([[3., 7.]]),), (np.array([[.1, .7, .05], [.05, .1, .7]]),), np.array([0]), np.array([0]), np.array([[2., 6.]]), np.zeros(1))
    likelihood = prior_quadrature(source)
    a, b = mean * precision, (1 - mean) * precision
    # QUADPACK's weighted endpoint rule handles Beta exponents > -1,
    # without evaluating log-density singularities at p=0 or p=1.
    integral = quad(lambda p: np.exp(source.row_log_likelihood(0, np.array([p]))[0] - betaln(a, b)), 0, 1, weight="alg", wvar=(a - 1, b - 1), epsabs=1e-25, epsrel=1e-10)[0]
    logits = np.array([[np.log(mean / (1 - mean))]])
    actual = likelihood.log_integrals(logits, precision, nodes=64)[0, 0]
    assert actual == pytest.approx(np.log(integral), abs=1e-7)
    _, eta, _ = likelihood.log_integrals(logits, precision, nodes=64, gradient=True)
    step = 1e-3
    expected = (likelihood.log_integrals(logits + step, precision, nodes=64) - likelihood.log_integrals(logits - step, precision, nodes=64)) / (2 * step)
    assert np.allclose(eta, expected)


def test_changing_proposal_counts_does_not_change_prior_only_quadrature():
    source = BinaryECPathLikelihood((np.array([[3., 7.]]),), (np.array([[.1, .7, .05], [.05, .1, .7]]),), np.array([0]), np.array([0]), np.array([[2., 6.]]), np.zeros(1))
    first = prior_quadrature(source)
    source.proposal_counts = np.array([[200., 600.]])
    second = prior_quadrature(source)
    logits = np.array([[0.]])
    assert np.array_equal(first.log_integrals(logits, 20), second.log_integrals(logits, 20))


def repeated_likelihood():
    counts = (np.array([[3., 7., 4.], [7., 3., 2.], [2., 4., 6.], [4., 2., 1.]]), np.array([[2., 4.], [8., 1.], [4., 3.], [6., 1.]]))
    components = (np.array([[.1, .7, .05], [.05, .1, .7], [.1, .2, .15]]), np.array([[.02, .8, .03], [.04, .05, .6]]))
    return BinaryECPathLikelihood(counts, components, np.repeat([0, 1], 2), np.tile([0, 1], 2), np.ones((4, 2)), np.array([1., 2., 3., 4.]))


def test_shared_probability_matrix_matches_every_independent_row():
    source = repeated_likelihood()
    points, rows = np.array([.01, .2, .7, .99]), np.array([3, 0, 2])
    actual = matrix_log_likelihood(source, rows, points)
    expected = np.stack([source.row_log_likelihood(row, points) for row in rows])
    assert np.allclose(actual, expected, atol=1e-12)


def test_shared_quadrature_is_numerically_equivalent_and_reuses_rules(monkeypatch):
    from tealeaf.sc import path_marginal_quadrature as module
    source = repeated_likelihood()
    unshared, shared = prior_quadrature(source), shared_prior_quadrature(source)
    logits = np.array([[-.4, .6, 8.], [.2, -.3, 4.], [-.4, .6, 8.], [.2, -.3, 4.]])
    calls = []
    original = module.beta_quadrature
    def counted(*args):
        calls.append(args)
        return original(*args)
    monkeypatch.setattr(module, "beta_quadrature", counted)
    expected = unshared.log_integrals(logits, 20, nodes=32)
    assert len(calls) == 12
    calls.clear()
    actual = shared.log_integrals(logits, 20, nodes=32)
    assert len(calls) == 6
    assert np.allclose(actual, expected, atol=1e-12)
    expected = unshared.log_integrals(logits, 20, nodes=32, gradient=True)
    actual = shared.log_integrals(logits, 20, nodes=32, gradient=True)
    assert all(np.allclose(left, right, atol=1e-8) for left, right in zip(actual, expected))


def test_shared_exact_beta_branch_handles_multiple_rows_and_offsets():
    counts = np.array([[3., 7.], [8., 2.], [4., 6.], [5., 5.]])
    source = BinaryECPathLikelihood((counts,), (np.array([[0, 1, 0], [0, 0, 1]]),), np.repeat([0, 1], 2), np.tile([0, 1], 2), np.ones((4, 2)), np.arange(4.))
    logits = np.array([[-.4, .6], [.2, -.3], [-.4, .6], [.2, -.3]])
    assert np.allclose(shared_prior_quadrature(source).log_integrals(logits, 20), prior_quadrature(source).log_integrals(logits, 20), atol=1e-12)
