"""Tests for independent batched L-BFGS optimization."""

import numpy as np

from analyses.benchmark_batched_laplace_lbfgs import (
    batched_laplace_functions,
    fit_batched,
)
from tealeaf.sc import ec_glmm
from tealeaf.sc.batched_lbfgs import minimize_batched_lbfgs


def test_batched_lbfgs_uses_independent_histories_and_bounds():
    targets = np.asarray([[2.0, -3.0], [-0.5, 0.75], [0.25, -0.4]])
    curvature = np.asarray([[1.0, 8.0], [20.0, 0.5], [0.2, 30.0]])

    def objective(parameters):
        residual = parameters - targets
        return 0.5 * np.sum(curvature * residual**2, axis=1), curvature * residual

    result = minimize_batched_lbfgs(
        objective,
        np.zeros_like(targets),
        [(-1.0, 1.0), (-2.0, 2.0)],
        max_iter=100,
    )
    np.testing.assert_allclose(
        result.parameters,
        [[1.0, -2.0], [-0.5, 0.75], [0.25, -0.4]],
        atol=2e-5,
    )
    assert result.converged.all()
    assert not result.failed.any()
    assert len(set(result.iterations)) > 1


def test_batched_lbfgs_masks_a_failed_line_search():
    def objective(parameters):
        values = np.sum(parameters**2, axis=1)
        gradients = 2.0 * parameters
        invalid = parameters[:, 0] < 0.0
        values = np.where(invalid, np.nan, values)
        gradients = np.where(invalid[:, None], np.nan, gradients)
        return values, gradients

    result = minimize_batched_lbfgs(
        objective,
        np.asarray([[1.0], [0.0]]),
        [(None, None)],
        max_iter=20,
    )
    assert result.converged.all()
    assert not result.failed.any()
    np.testing.assert_allclose(result.parameters, 0.0, atol=1e-8)


def test_padded_laplace_batch_handles_empty_primers():
    first = ec_glmm.ECGLMMData(
        counts=(np.asarray([[10, 5], [8, 7]]), np.zeros((2, 0))),
        compatibility=(np.eye(2), np.zeros((0, 2))),
        design=np.ones((2, 1)),
        clusters=np.asarray(["a", "b"]),
        fixed_effect_tensor=np.ones((2, 1, 1)),
    )
    second = ec_glmm.ECGLMMData(
        counts=(
            np.asarray([[7, 5, 2], [5, 8, 1], [9, 4, 3]]),
            np.asarray([[6, 8], [7, 7], [8, 8]]),
        ),
        compatibility=(
            np.asarray([[1, 0], [0, 1], [1, 1]]),
            np.eye(2),
        ),
        design=np.ones((3, 1)),
        clusters=np.asarray(["a", "a", "b"]),
        fixed_effect_tensor=np.ones((3, 1, 1)),
    )
    initial, evaluate, evaluate_score = batched_laplace_functions([first, second], 12)
    values, gradients = evaluate(initial)
    assert np.isfinite(values).all()
    assert np.isfinite(gradients).all()
    result, _, gradient_norm, mode_score, converged = fit_batched(
        [first, second], 80, 12
    )
    assert converged.all()
    assert (gradient_norm <= 1e-4).all()
    assert (mode_score <= 1e-4).all()
    scalar = [
        ec_glmm.fit_laplace(
            data,
            max_iter=80,
            mode_steps=12,
            mode_gradient="implicit",
        )
        for data in (first, second)
    ]
    np.testing.assert_allclose(
        result.values,
        [fit["objective"] for fit in scalar],
        atol=2e-3,
    )
