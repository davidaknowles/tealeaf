import numpy as np
from scipy.integrate import quad
from scipy.special import betaln, expit
from scipy.stats import beta

from tealeaf.sc.path_marginal import BinaryECPathLikelihood
from tealeaf.sc.path_marginal_grid import GridModel, binary_grid_test, row_grids, row_integral


def mixed_likelihood(counts, subjects, labels):
    # ECs: path 0 only, path 1 only, both paths with unequal weight, outside only.
    components = (np.array([[0., .5, 0.], [0., 0., .4], [0., .3, .6], [.2, 0., 0.]]), np.array([[0., .2, 0.], [0., 0., .5], [.1, .4, .3], [.3, .1, .1]]))
    counts = np.asarray(counts, dtype=float)
    n = len(subjects)
    return BinaryECPathLikelihood((counts[:, :4], counts[:, 4:]), components, np.asarray(subjects), np.asarray(labels), np.ones((n, 2)), np.zeros(n))


def test_row_integral_matches_direct_quadrature():
    like = mixed_likelihood([[3, 1, 7, 2, 0, 4, 5, 1], [40, 2, 30, 9, 1, 25, 60, 4]], ["s", "s"], ["a", "b"])
    grids = row_grids(like)
    for row in range(2):
        for mean, kappa in ((0., 5.), (1.5, 200.), (-2., .4)):
            a, c = kappa * expit(mean), kappa * (1 - expit(mean))
            integrand = lambda p: beta.pdf(p, a, c) * np.exp(like.row_log_likelihood(row, np.array([p]))[0] - grids[row].shift)
            direct = np.log(quad(integrand, 0, 1, limit=500, epsabs=0, epsrel=1e-11, points=[1e-6, 1e-3, .5, 1 - 1e-3])[0])
            value = row_integral(grids[row], np.array([mean]), kappa)[0]
            assert abs(value[0] - direct) < 1e-6


def test_exact_rows_use_beta_binomial():
    components = (np.array([[0., 1., 0.], [0., 0., 1.]]),)
    like = BinaryECPathLikelihood((np.array([[6., 2.]]),), components, np.array(["s"]), np.array(["a"]), np.ones((1, 2)), np.zeros(1))
    grid = row_grids(like)[0]
    value = row_integral(grid, np.array([.3]), 7.)[0]
    a, c = 7 * expit(.3), 7 * (1 - expit(.3))
    assert np.isclose(value[0], betaln(a + 6, c + 2) - betaln(a, c))


def simulated(seed, effect=0.):
    rng = np.random.default_rng(seed)
    counts, subjects, labels = [], [], []
    for subject in range(8):
        offset = rng.normal(0, .5)
        for label, shift in (("a", 0.), ("b", effect)):
            psi = expit(.3 + offset + shift + rng.normal(0, .3))
            probabilities = np.array([.5 * psi, .4 * (1 - psi), .3 * psi + .6 * (1 - psi), .2 * .3])
            first = rng.multinomial(60, probabilities / probabilities.sum())
            probabilities = np.array([.2 * psi, .5 * (1 - psi), .1 * .3 + .4 * psi + .3 * (1 - psi), .3 * .3 + .1])
            second = rng.multinomial(40, probabilities / probabilities.sum())
            counts.append(np.r_[first, second])
            subjects.append(f"s{subject}")
            labels.append(label)
    return mixed_likelihood(counts, subjects, labels)


def test_gradient_matches_finite_differences():
    like = simulated(1, .8)
    design = np.column_stack([np.ones(16), like.labels == "b"]).astype(float)
    model = GridModel(row_grids(like), design, like.subjects)
    parameters = np.array([.2, .5, np.log(15.), np.log(.4)])
    _, gradient = model.objective(parameters)
    numeric = []
    for index in range(4):
        step = np.zeros(4)
        step[index] = 1e-5
        numeric.append((model.objective(parameters + step)[0] - model.objective(parameters - step)[0]) / 2e-5)
    assert np.allclose(gradient, numeric, rtol=1e-3, atol=1e-4)


def test_effect_detected_and_null_accurate():
    effect = binary_grid_test(simulated(2, 1.2))
    null = binary_grid_test(simulated(3, 0.))
    assert effect["converged"] and null["converged"]
    assert effect["p_value"] < 1e-3 and effect["coefficients"][1] > 0
    assert effect["quadrature_error"] < 1e-3


def test_row_integral_derivatives_match_finite_differences():
    like = mixed_likelihood([[3, 1, 7, 2, 0, 4, 5, 1], [40, 2, 30, 9, 1, 25, 60, 4]], ["s", "s"], ["a", "b"])
    grid = row_grids(like)[1]
    exact = row_grids(BinaryECPathLikelihood((np.array([[6., 2.]]),), (np.array([[0., 1., 0.], [0., 0., 1.]]),), np.array(["s"]), np.array(["a"]), np.ones((1, 2)), np.zeros(1)))[0]
    for current in (grid, exact):
        for mean, kappa in ((0., 5.), (1.5, 200.), (-2., .4)):
            m, h = np.array([mean]), 1e-4
            value, slope, curvature, d_kappa = row_integral(current, m, kappa)
            up, down = row_integral(current, m + h, kappa), row_integral(current, m - h, kappa)
            assert np.isclose(slope[0], (up[0][0] - down[0][0]) / (2 * h), rtol=1e-4, atol=1e-6)
            assert np.isclose(curvature[0], (up[1][0] - down[1][0]) / (2 * h), rtol=1e-4, atol=1e-6)
            k_up, k_down = row_integral(current, m, kappa * np.exp(h)), row_integral(current, m, kappa * np.exp(-h))
            assert np.isclose(d_kappa[0], (k_up[0][0] - k_down[0][0]) / (2 * h), rtol=1e-4, atol=1e-6)


def test_primer_offset_recovered():
    from tealeaf.sc.path_marginal_grid import estimate_primer_offset
    rng = np.random.default_rng(5)
    rows = []
    for subject in range(10):
        psi = expit(rng.normal(0, 1))
        shifted = expit(np.log(psi / (1 - psi)) + 1.2)
        rows.append(np.r_[rng.multinomial(400, [.5 * psi / (.5 * psi + .4 * (1 - psi)), .4 * (1 - psi) / (.5 * psi + .4 * (1 - psi))]), rng.multinomial(400, [.5 * shifted / (.5 * shifted + .4 * (1 - shifted)), .4 * (1 - shifted) / (.5 * shifted + .4 * (1 - shifted))])])
    components = np.array([[0., .5, 0.], [0., 0., .4]])
    like = BinaryECPathLikelihood((np.array(rows)[:, :2], np.array(rows)[:, 2:]), (components, components), np.array([f"s{i}" for i in range(10)]), np.array(["a"] * 10), np.ones((10, 2)), np.zeros(10))
    assert abs(estimate_primer_offset(like) - 1.2) < .15
