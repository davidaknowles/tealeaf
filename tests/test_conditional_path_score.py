import numpy as np
import pytest
from scipy import optimize

from tealeaf.sc.conditional_path_score import ConditionalPathNullProblem, efficient_conditional_path_score, binary_fragment_opportunity_kernels
from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_score_mixed import shared_path_score_components, aggregate_path_scores
from tealeaf.sc.path_simulation import simulate_counts


def problem(row_offsets=None):
    mapping = np.array([[1., .1, .2, .3], [.2, 1., .1, .5], [.3, .2, 1., .2], [.4, .2, .5, 1.], [1., 1., 0., 0.], [0., .5, 1., .4], [.7, .3, .2, .9]])
    if row_offsets is not None:
        mapping *= np.asarray(row_offsets)[:, None]
    return ConditionalPathNullProblem((np.array([[20, 30, 10, 4, 8, 9, 11], [15, 21, 22, 9, 14, 7, 8]]),), (mapping,), np.array([.2, .3, .3, .2]), np.array([0, 0, 1, -1]))


def test_conditional_objective_gradient_and_shared_opportunity_cancellation():
    original = problem()
    transformed = problem(np.exp(np.arange(7) - 3.))
    parameters = original.initial + np.linspace(-.2, .3, original.dimension)
    value, gradient = original.objective(parameters)
    shifted_value, shifted_gradient = transformed.objective(parameters)
    np.testing.assert_allclose(value, shifted_value, rtol=1e-12)
    np.testing.assert_allclose(gradient, shifted_gradient, atol=1e-12)
    numerical = optimize.approx_fprime(parameters, lambda p: original.objective(p)[0], 1e-6)
    np.testing.assert_allclose(gradient, numerical, atol=3e-5)
    theta = original.structure.composition(parameters[:original.structure.dimension])[0]
    offsets = parameters[original.structure.dimension:]
    first = efficient_conditional_path_score(original, theta, offsets)
    second = efficient_conditional_path_score(transformed, theta, offsets)
    for left, right in zip(first, second):
        np.testing.assert_allclose(left, right, rtol=1e-8, atol=1e-12)


def test_two_path_cassette_retains_information_and_detects_a_contrast():
    mapping = np.array([[1., 0.], [0., 1.], [1., 1.]])
    subjects, labels = np.repeat(np.arange(6), 2), np.tile([0, 1], 6)
    counts = np.tile([[20., 40., 100.], [40., 20., 100.]], (6, 1))
    data = ECGLMMData((counts,), (mapping,), np.ones((12, 1)), subjects)
    components = shared_path_score_components(data, [0, 1], labels, subjects, baseline=np.array([.5, .5]), count_likelihood="conditional", null_multistart=True)
    result = aggregate_path_scores(components, information_metric="reference")
    assert result["n_subjects"] == 6 and result["mean_difference"][0] > 0
    assert result["p_value"] < .05
    assert (components.information > 0).all()


def test_conditional_model_rejects_no_shared_primer_and_nonbinary_types():
    with pytest.raises(ValueError, match="primer observed"):
        ConditionalPathNullProblem((np.array([[10., 0], [0., 0]]), np.array([[0., 0], [0., 10]])), (np.eye(2), np.eye(2)), [.5, .5], [0, 1])
    with pytest.raises(ValueError, match="exactly two"):
        ConditionalPathNullProblem((np.ones((3, 2)),), (np.eye(2),), [.5, .5], [0, 1])


def test_opportunity_stress_changes_generation_not_analysis_map():
    counts = np.full((4, 3), 50.)
    mapping = np.array([[1., 0.], [0., 1.], [1., 1.]])
    data = ECGLMMData((counts,), (mapping,), np.ones((4, 1)), np.repeat([0, 1], 2))
    first = simulate_counts(data, [.5, .5], data.clusters, np.random.default_rng(1))
    zero = simulate_counts(data, [.5, .5], data.clusters, np.random.default_rng(1), ec_opportunity_scale=0.)
    np.testing.assert_array_equal(first.counts[0], zero.counts[0])
    stressed, details = simulate_counts(data, [.5, .5], data.clusters, np.random.default_rng(1), ec_opportunity_scale=1., return_details=True)
    np.testing.assert_array_equal(stressed.compatibility[0], mapping)
    np.testing.assert_array_equal(stressed.counts[0].sum(axis=1), counts.sum(axis=1))
    assert not np.array_equal(first.counts[0], stressed.counts[0])
    factors = details["ec_opportunity_factors"][0]
    assert factors.shape == (3,) and (factors > 0).all()
    weights = details["observation_weights"]
    np.testing.assert_array_equal(weights[0], weights[1])
    np.testing.assert_array_equal(weights[2], weights[3])


def test_fragment_kernels_restore_opportunity_units_and_reject_weighted_maps():
    membership = np.array([[1., 0.], [1., 1.], [0., 1.], [0., 1.]])
    degree = membership.sum(axis=0)
    lengths = np.array([2., 4.])
    prepared = (membership / degree, membership / degree * lengths)
    fragment = binary_fragment_opportunity_kernels(prepared)
    np.testing.assert_allclose(fragment[0].toarray(), membership / lengths)
    np.testing.assert_allclose(fragment[1].toarray(), membership)
    weighted = prepared[0].copy()
    weighted[0, 0] += .1
    weighted[1, 0] -= .1
    with pytest.raises(ValueError, match="column-constant"):
        binary_fragment_opportunity_kernels((weighted, prepared[1]))


def test_ec_degree_weighting_can_change_the_biological_path_null():
    rna = np.array([[.3, .3, .4], [.5, .1, .4]])
    degree, lengths = np.array([2., 8., 4.]), np.array([1., 2., 1.])
    prepared = rna * (degree / lengths)
    prepared /= prepared.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(rna[:, :2].sum(axis=1), [.6, .6])
    assert not np.isclose(prepared[0, :2].sum(), prepared[1, :2].sum())
    recovered = prepared * (lengths / degree)
    recovered /= recovered.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(recovered, rna)
