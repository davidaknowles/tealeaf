import numpy as np
import pytest

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_bias import SharedPathNullProblem, expected_ec_counts, paired_null_corrected_path_test, null_corrected_path_responses


@pytest.mark.parametrize("depth_scale", [1., 1e5])
def test_scaled_shared_null_recovers_categorical_mle_at_different_depths(depth_scale):
    counts = np.array([[104., 52., 44.], [106., 50., 44.], [106., 57., 37.]]) * depth_scale
    problem = SharedPathNullProblem((counts,), (np.eye(3),), np.ones(3), [0, 1, 2])
    fitted = problem.fit()
    assert fitted.converged
    assert fitted.termination_message
    assert np.allclose(fitted.path_proportions, counts.sum(axis=0) / counts.sum(), atol=1e-7)
    assert fitted.objective == pytest.approx(problem.objective(problem.basis.T @ np.log(fitted.path_proportions))[0])


@pytest.mark.parametrize("outside", [False, True])
def test_shared_path_null_gradient_matches_independent_differences(outside):
    rng = np.random.default_rng(714)
    paths = np.array([0, 0, 1, 2, 2] + ([-1, -1] if outside else []))
    maps = tuple(rng.uniform(.1, 2, (k, len(paths))) for k in (6, 4))
    counts = tuple(rng.integers(0, 15, (3, k)).astype(float) for k in (6, 4))
    problem = SharedPathNullProblem(counts, maps, np.ones(len(paths)), paths)
    parameters = problem.initial + rng.normal(0, .2, problem.dimension)
    value, gradient = problem.objective(parameters)
    assert np.isfinite(value)
    for index in range(problem.dimension):
        step = np.eye(problem.dimension)[index] * 1e-5
        numeric = (problem.objective(parameters + step)[0] - problem.objective(parameters - step)[0]) / 2e-5
        assert gradient[index] == pytest.approx(numeric, abs=2e-7)
    theta, _, proportions, _, _ = problem.composition(parameters)
    assert np.allclose(theta.sum(axis=1), 1)
    for cell in range(3):
        masses = np.array([theta[cell, paths == path].sum() for path in range(3)])
        assert np.allclose(masses / masses.sum(), proportions)


def test_shared_path_null_allows_different_nuisance_mixtures_and_outside_mass():
    truth = np.array([[.36, .04, .1, .5], [.028, .252, .07, .65]])
    counts = (truth * np.array([[1000], [2000]]),)
    problem = SharedPathNullProblem(counts, (np.eye(4),), np.ones(4), [0, 0, 1, -1])
    fitted = problem.fit()
    assert fitted.converged
    assert np.allclose(fitted.theta, truth, atol=2e-6)
    assert np.allclose(fitted.path_proportions, [.8, .2], atol=2e-6)
    expected = expected_ec_counts(counts, (np.eye(4),), fitted.theta)
    assert np.allclose(expected[0].sum(axis=1), [1000, 2000])


def test_expected_counts_preserve_separate_primer_totals_and_zero_maps():
    counts = (np.array([[5., 5.], [10., 20.]]), np.zeros((2, 3)))
    theta = np.array([[.8, .2], [.3, .7]])
    expected = expected_ec_counts(counts, (np.diag([2., 3.]), np.zeros((3, 2))), theta)
    assert np.allclose(expected[0].sum(axis=1), [10, 30])
    assert np.array_equal(expected[1], np.zeros((2, 3)))


def test_null_correction_removes_deterministic_unequal_depth_smoothing_contrast():
    subjects = np.repeat(np.arange(6), 2)
    labels = np.tile([0, 1], 6)
    depth = np.tile([10., 100.], 6)
    counts = depth[:, None] * np.array([.8, .2])
    data = ECGLMMData((counts,), (np.eye(2),), np.ones((12, 1)), subjects)
    result = paired_null_corrected_path_test(data, [0, 1], labels, subjects, baseline=np.array([.8, .2]))
    assert result["converged"]
    assert np.allclose(result["differences"], 0, atol=1e-7)
    assert result["p_value"] == pytest.approx(1.)
    assert abs(result["observed_proportions"][1, 0] - result["observed_proportions"][0, 0]) > .1


def test_binary_correction_matches_closed_form_and_reverses_direction():
    rng = np.random.default_rng(21)
    subjects = np.repeat(np.arange(8), 2)
    labels = np.tile([0, 1], 8)
    depth = np.tile([20., 100.], 8)
    count = rng.binomial(depth.astype(int), .3 + .3 * labels)
    counts = np.column_stack([count, depth - count]).astype(float)
    data = ECGLMMData((counts,), (np.eye(2),), np.ones((16, 1)), subjects)
    result = paired_null_corrected_path_test(data, [0, 1], labels, subjects, baseline=np.ones(2) / 2)
    expected = []
    for subject in range(8):
        rows = subjects == subject
        common = counts[rows, 0].sum() / depth[rows].sum()
        residual = (counts[rows, 0] - depth[rows] * common) / (depth[rows] + 32 + 1e-4)
        expected.append(np.sqrt(2) * (residual[1] - residual[0]))
    assert np.allclose(result["differences"][:, 0], expected, atol=2e-6)
    reverse = paired_null_corrected_path_test(data, [0, 1], 1 - labels, subjects, baseline=np.ones(2) / 2)
    assert np.allclose(reverse["differences"], -result["differences"], atol=2e-6)
    assert reverse["p_value"] == pytest.approx(result["p_value"], rel=1e-5)


def test_invalid_null_problem_rejects_missing_path_and_negative_counts():
    with pytest.raises(ValueError, match="consecutive"):
        SharedPathNullProblem((np.ones((2, 2)),), (np.eye(2),), np.ones(2), [0, 2])
    with pytest.raises(ValueError, match="nonnegative"):
        SharedPathNullProblem((-np.ones((2, 2)),), (np.eye(2),), np.ones(2), [0, 1])


def test_general_response_correction_retains_three_paths_and_missing_types():
    subjects = np.repeat(np.arange(5), 3)[:-1]
    labels = np.tile(np.arange(3), 5)[:-1]
    rng = np.random.default_rng(811)
    depths = np.take([40, 100, 200], labels)
    counts = np.array([rng.multinomial(depth, [.5, .3, .2]) for depth in depths])
    data = ECGLMMData((counts,), (np.eye(3),), np.ones((len(labels), 1)), subjects)
    result = null_corrected_path_responses(data, [0, 1, 2], labels, subjects, baseline=np.ones(3) / 3, reporting_concentration=1.)
    assert result["values"].shape == (14, 2)
    assert np.array_equal(result["subjects"], subjects)
    assert np.array_equal(result["encoded_labels"], labels)
    assert np.allclose(result["reporting_proportions"].sum(axis=1), 1)
