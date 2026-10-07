import numpy as np
import pytest
from scipy import stats

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_bias import SharedPathNullProblem
from tealeaf.sc.path_score_mixed import efficient_shared_path_score, mixed_score_test, mixed_path_score_test, score_contrast_proportions


def test_scalar_mixed_score_matches_inverse_variance_and_modified_kh():
    values = np.array([.4, -.2, .6, .1, .8])
    variances = np.array([.05, .5, .1, .2, 1.])
    information = (1 / variances)[:, None, None]
    result = mixed_score_test((values / variances)[:, None], information, biological_variance=.2)
    weights = 1 / (variances + .2)
    mean = weights @ values / weights.sum()
    inflation = max(1., np.sum(weights * (values - mean) ** 2) / 4)
    statistic = mean ** 2 * weights.sum() / inflation
    assert result["mean_difference"][0] == pytest.approx(mean)
    assert result["statistic"] == pytest.approx(statistic)
    assert result["p_value"] == pytest.approx(stats.f.sf(statistic, 1, 4))


def test_partial_information_is_missing_not_zero_variance():
    values = np.array([[.3, 0], [.4, 0], [0, -.2], [0, -.1], [.2, -.3]])
    information = np.array([np.diag([5, 0]), np.diag([5, 0]), np.diag([0, 5]), np.diag([0, 5]), np.diag([5, 5])], dtype=float)
    scores = np.einsum("mij,mj->mi", information, values)
    result = mixed_score_test(scores, information, biological_variance=0.)
    assert np.allclose(result["mean_difference"], [.3, -.2])
    assert result["residual_degrees_of_freedom"] == 4
    assert result["n_subjects"] == 5
    with pytest.raises(ValueError, match="identify"):
        mixed_score_test(scores[:, :], np.tile(np.diag([1, 0]), (5, 1, 1)), biological_variance=0.)


def test_score_matches_categorical_two_sample_fisher_information():
    theta = np.array([[.7, .3], [.7, .3]])
    counts = (np.array([[6., 4.], [75., 25.]]),)
    score, information, shape = efficient_shared_path_score(counts, (np.eye(2),), theta, [0, 1], [0, 1], 2)
    # Helmert path ILR is log(p0/p1)/sqrt(2); its Fisher is 2Np0p1.
    expected_information = 2 * .7 * .3 * 10 * 100 / 110
    expected_score = np.sqrt(2) * ((75 - 70) - 100 / 110 * ((6 - 7) + (75 - 70)))
    assert information[0, 0] == pytest.approx(expected_information)
    assert score[0] == pytest.approx(expected_score)
    assert shape[0, 0] == pytest.approx(1 / .7 + 1 / .3)


def test_different_within_path_composition_and_outside_mass_have_zero_null_score():
    theta = np.array([[.36, .04, .1, .5], [.028, .252, .07, .65]])
    counts = (theta * np.array([[1000], [2000]]),)
    null = SharedPathNullProblem(counts, (np.eye(4),), np.ones(4), [0, 0, 1, -1]).fit()
    assert null.converged
    score, information, _ = efficient_shared_path_score(counts, (np.eye(4),), null.theta, [0, 0, 1, -1], [0, 1], 2)
    assert np.max(np.abs(score)) < 1e-4
    assert information[0, 0] > 0


def test_expected_unequal_primer_totals_cannot_manufacture_effect():
    subjects = np.repeat(np.arange(6), 2)
    labels = np.tile([0, 1], 6)
    counts = np.tile([10., 100.], 6)[:, None] * np.array([.8, .2])
    data = ECGLMMData((counts,), (np.eye(2),), np.ones((12, 1)), subjects)
    result = mixed_path_score_test(data, [0, 1], labels, subjects, baseline=np.array([.8, .2]))
    assert np.max(np.abs(result["mean_difference"])) < 1e-8
    assert result["p_value"] == pytest.approx(1.)


def test_three_types_three_paths_with_missing_reference_in_one_subject():
    rng = np.random.default_rng(304)
    subjects = np.repeat(np.arange(6), 3)[1:]
    labels = np.tile(np.arange(3), 6)[1:]
    counts = np.array([rng.multinomial(200, [.5, .3, .2]) for _ in labels])
    data = ECGLMMData((counts,), (np.eye(3),), np.ones((len(labels), 1)), subjects)
    result = mixed_path_score_test(data, [0, 1, 2], labels, subjects, baseline=np.ones(3) / 3)
    assert result["degrees_of_freedom"] == 4
    assert result["n_subjects"] == 6
    assert result["components"].information.shape == (6, 4, 4)
    assert np.linalg.matrix_rank(result["components"].information[0]) == 2


def test_score_sign_reversal_and_reml_rotation_invariance():
    rng = np.random.default_rng(413)
    scores = rng.normal(size=(7, 2))
    information = np.tile(np.array([[2., .3], [.3, 1.]]), (7, 1, 1))
    shapes = np.tile(np.array([[1., .2], [.2, 2.]]), (7, 1, 1))
    original = mixed_score_test(scores, information, shapes)
    reverse = mixed_score_test(-scores, information, shapes)
    assert reverse["p_value"] == pytest.approx(original["p_value"])
    assert np.allclose(reverse["mean_difference"], -original["mean_difference"])
    rotation = np.array([[.6, -.8], [.8, .6]])
    rotated = mixed_score_test(scores @ rotation.T, np.einsum("ij,mjk,lk->mil", rotation, information, rotation), np.einsum("ij,mjk,lk->mil", rotation, shapes, rotation))
    assert rotated["p_value"] == pytest.approx(original["p_value"], rel=1e-6)


def test_invalid_shapes_and_negative_information_are_rejected():
    with pytest.raises(ValueError, match="finite"):
        mixed_score_test(np.ones((5, 2)), np.ones((5, 2)))
    with pytest.raises(ValueError, match="semidefinite"):
        mixed_score_test(np.ones((5, 1)), -np.ones((5, 1, 1)))


def test_rare_path_does_not_break_analytic_nuisance_rank():
    theta = np.array([[.4, .6 - 1e-15, 1e-15], [.2, .8 - 1e-15, 1e-15]])
    counts = (theta * 1000,)
    score, information, shape = efficient_shared_path_score(counts, (np.eye(3),), theta, [0, 0, 1], [0, 1], 2)
    assert np.isfinite(score).all()
    assert np.isfinite(information).all()
    assert np.isfinite(shape).all()
    assert np.allclose(score, 0)


def test_score_reporting_is_a_simplex_pair_with_exact_ilr_contrast():
    from tealeaf.sc.differential import helmert_basis
    anchor, contrast = np.array([.7, .2, .1]), np.array([.4, -.3])
    pair = score_contrast_proportions(anchor, contrast)
    assert np.allclose(pair.sum(axis=1), 1)
    assert (pair > 0).all()
    assert np.allclose(helmert_basis(3).T @ (np.log(pair[1]) - np.log(pair[0])), contrast)
    assert np.allclose(score_contrast_proportions(anchor, -contrast), pair[::-1])
    assert np.allclose(score_contrast_proportions(anchor, np.zeros(2)), [anchor, anchor])
