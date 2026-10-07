import numpy as np
import pytest

from tealeaf.sc.contrast_ops import helmert_multiply, helmert_transpose_multiply
from tealeaf.sc.differential import helmert_basis
from tealeaf.sc.path_bias import SharedPathNullProblem


@pytest.mark.parametrize("size", [1, 2, 3, 31, 200])
def test_linear_helmert_matches_dense_convention(size):
    rng = np.random.default_rng(51)
    coordinates, values = rng.normal(size=(3, size - 1)), rng.normal(size=(3, size))
    basis = helmert_basis(size) if size > 1 else np.zeros((1, 0))
    np.testing.assert_allclose(helmert_multiply(coordinates), coordinates @ basis.T, atol=1e-14)
    np.testing.assert_allclose(helmert_transpose_multiply(values), values @ basis, atol=1e-14)


@pytest.mark.parametrize("outside", [False, True])
@pytest.mark.parametrize("transcripts", [4, 25, 200])
def test_vector_null_objective_matches_dense(outside, transcripts):
    rng = np.random.default_rng(671)
    paths = np.arange(transcripts) % 3
    if outside:
        paths[-1] = -1
    maps = tuple(rng.uniform(.01, 2., (rows, transcripts)) for rows in (7, 5))
    counts = tuple(rng.integers(0, 5000, (3, rows)).astype(float) for rows in (7, 5))
    counts[1][0] = 0.
    problem = SharedPathNullProblem(counts, maps, rng.uniform(.1, 1., transcripts), paths)
    parameters = problem.initial + rng.normal(scale=.2, size=problem.dimension)
    original, first_gradient = problem.objective(parameters)
    vector, second_gradient = problem.objective_vector(parameters)
    np.testing.assert_allclose(original, vector, rtol=1e-13, atol=1e-9)
    np.testing.assert_allclose(first_gradient, second_gradient, rtol=1e-9, atol=1e-9)


def test_vector_null_fit_preserves_optimum_and_bounds():
    truth = np.array([[.36, .04, .1, .5], [.028, .252, .07, .65]])
    problem = SharedPathNullProblem((truth * np.array([[1000], [2000]]),), (np.eye(4),), np.ones(4), [0, 0, 1, -1])
    dense = problem.fit(max_iter=2000, multistart=True)
    vector = problem.fit(max_iter=2000, multistart=True, objective_method="vector")
    assert dense.converged and vector.converged
    np.testing.assert_allclose(vector.theta, dense.theta, atol=3e-7)
    np.testing.assert_allclose(vector.objective, dense.objective, rtol=1e-11, atol=1e-7)
    with pytest.raises(ValueError, match="objective method"):
        problem.fit(objective_method="unknown")
