import numpy as np
import pytest
from scipy import stats

from tealeaf.sc.ec_glmm import ECGLMMData
from tealeaf.sc.path_bias import SharedPathNullProblem
from tealeaf.sc.path_score_mixed import efficient_shared_path_score, mixed_score_test, mixed_path_score_test, score_contrast_proportions


@pytest.mark.parametrize("size", [2, 3, 4])
def test_proportion_score_is_the_interior_likelihood_coordinate_transform(size):
    from tealeaf.sc.differential import helmert_basis
    from tealeaf.sc.path_score_mixed import path_proportion_covariance
    rng = np.random.default_rng(1230 + size)
    psi = rng.dirichlet(np.ones(size) * 3)
    theta = np.tile(psi, (3, 1))
    maps = tuple(rng.uniform(.1, 2, (size + 3, size)) for _ in range(2))
    counts = tuple(np.array([rng.multinomial(400, mapping @ psi / (mapping @ psi).sum()) for _ in range(3)]) for mapping in maps)
    ilr = efficient_shared_path_score(counts, maps, theta, np.arange(size), [0, 1, 2], 3)
    prop = efficient_shared_path_score(counts, maps, theta, np.arange(size), [0, 1, 2], 3, score_coordinate="proportion")
    basis = helmert_basis(size)
    transform = np.kron(np.eye(2), basis.T @ path_proportion_covariance(psi) @ basis)
    np.testing.assert_allclose(ilr[0], transform.T @ prop[0], atol=1e-9)
    np.testing.assert_allclose(ilr[1], transform.T @ prop[1] @ transform, atol=1e-9)
    np.testing.assert_allclose(prop[2], transform @ ilr[2] @ transform.T, atol=1e-12)


def test_proportion_covariance_is_stable_at_a_nearly_absent_path():
    from tealeaf.sc.path_score_mixed import path_proportion_covariance
    psi = np.array([1., 1e-18])
    covariance = path_proportion_covariance(psi)
    np.testing.assert_allclose(covariance / 1e-18, [[1., -1.], [-1., 1.]], atol=1e-14)
    assert np.linalg.eigvalsh(covariance).min() >= 0
    with pytest.raises(ValueError, match="proportions"):
        path_proportion_covariance([.2, -.1])


def test_absolute_information_floor_can_depend_on_path_coordinate():
    # Fractional expected counts demonstrate the numerical scale issue only,
    # not finite-count validity or biological power at an actual boundary.
    psi = np.array([1., 1e-18])
    theta = np.tile(psi, (2, 1))
    counts = (theta * 100000.,)
    ilr = efficient_shared_path_score(counts, (np.eye(2),), theta, [0, 1], [0, 1], 2)
    prop = efficient_shared_path_score(counts, (np.eye(2),), theta, [0, 1], [0, 1], 2, score_coordinate="proportion")
    with pytest.raises(ValueError, match="four informative"):
        mixed_score_test(np.tile(ilr[0], (5, 1)), np.tile(ilr[1], (5, 1, 1)), np.tile(ilr[2], (5, 1, 1)))
    result = mixed_score_test(np.tile(prop[0], (5, 1)), np.tile(prop[1], (5, 1, 1)), np.tile(prop[2], (5, 1, 1)))
    assert result["n_subjects"] == 5 and result["p_value"] == pytest.approx(1.)


def test_constant_anchor_coordinate_transform_preserves_the_aggregate_test():
    from tealeaf.sc.differential import helmert_basis
    from tealeaf.sc.path_score_mixed import path_proportion_covariance
    rng = np.random.default_rng(1241)
    psi = np.array([.7, .2, .1])
    basis = helmert_basis(3)
    transform = basis.T @ path_proportion_covariance(psi) @ basis
    scores = rng.normal(size=(7, 2))
    information = np.tile(np.array([[2., .3], [.3, 1.]]), (7, 1, 1))
    shape = basis.T @ np.diag(1 / psi) @ basis
    original = mixed_score_test(scores, information, np.tile(shape, (7, 1, 1)))
    inverse = np.linalg.inv(transform)
    prop = mixed_score_test(scores @ inverse.T, np.einsum("ij,mjk,kl->mil", inverse.T, information, inverse), np.tile(transform @ shape @ transform.T, (7, 1, 1)))
    assert prop["p_value"] == pytest.approx(original["p_value"], rel=1e-6)
    np.testing.assert_allclose(prop["mean_difference"], transform @ original["mean_difference"], rtol=1e-6)


def test_target_reference_rank_retains_the_same_model_and_is_unit_invariant():
    rng = np.random.default_rng(1251)
    scores = rng.normal(size=(6, 1))
    information = np.arange(1., 7.)[:, None, None]
    reference = information * 3
    shapes = np.arange(2., 8.)[:, None, None]
    original = mixed_score_test(scores, information, shapes, biological_variance=.2)
    normalized = mixed_score_test(scores, information, shapes, biological_variance=.2, reference_information=reference)
    assert normalized["p_value"] == pytest.approx(original["p_value"])
    np.testing.assert_allclose(normalized["mean_difference"], original["mean_difference"])
    scale = 1e8
    reexpressed = mixed_score_test(scores / scale, information / scale**2, shapes * scale**2, biological_variance=.2, reference_information=reference / scale**2)
    assert reexpressed["p_value"] == pytest.approx(original["p_value"])
    np.testing.assert_allclose(reexpressed["mean_difference"], original["mean_difference"] * scale)


def test_reference_rank_does_not_fabricate_information_for_aliasing():
    for value in (0., 1e-25):
        with pytest.raises(ValueError, match="four informative"):
            mixed_score_test(np.zeros((5, 1)), np.full((5, 1, 1), value), reference_information=np.ones((5, 1, 1)))


def test_target_reference_retains_multivariate_partial_information():
    values = np.array([[.3, 0], [.4, 0], [0, -.2], [0, -.1], [.2, -.3]])
    information = np.array([np.diag([5, 0]), np.diag([5, 0]), np.diag([0, 5]), np.diag([0, 5]), np.diag([5, 5])], dtype=float)
    scores = np.einsum("mij,mj->mi", information, values)
    original = mixed_score_test(scores, information, biological_variance=0.)
    normalized = mixed_score_test(scores, information, biological_variance=0., reference_information=information * 2)
    assert normalized["p_value"] == pytest.approx(original["p_value"])
    np.testing.assert_allclose(normalized["mean_difference"], original["mean_difference"])
    assert normalized["residual_degrees_of_freedom"] == 4


def test_reference_information_cannot_be_smaller_than_profiled_information():
    with pytest.raises(ValueError, match="exceeds"):
        mixed_score_test(np.ones((5, 1)), np.ones((5, 1, 1)), reference_information=np.ones((5, 1, 1)) * .5)


def test_reference_rank_recovers_a_rare_log_ratio_without_changing_coordinates():
    psi = np.array([1., 1e-18])
    theta = np.tile(psi, (2, 1))
    result = efficient_shared_path_score((theta * 100000.,), (np.eye(2),), theta, [0, 1], [0, 1], 2, return_reference=True)
    score, information, shape, reference = (np.tile(value, (5,) + (1,) * value.ndim) for value in result)
    tested = mixed_score_test(score, information, shape, reference_information=reference)
    assert tested["n_subjects"] == 5 and tested["p_value"] == pytest.approx(1.)


def test_binary_score_records_keep_missing_reports_and_raw_information():
    from tealeaf.sc.path_score_mixed import PathScoreComponents, binary_subject_score_records
    reports = [[(0, [.2, .8]), (1, [np.nan, np.nan])], [(0, [.3, .7]), (1, [.6, .4])]]
    information = np.ones((2, 1, 1))
    components = PathScoreComponents(np.array([[.1], [.2]]), information, information * 3, np.array(["A", "B"]), (0, 1), [], reports, reference_information=information * 2)
    rows = binary_subject_score_records(components, "test")
    assert len(rows) == 2 and np.isnan(rows[0]["report_inclusion_b"])
    assert rows[1]["report_inclusion_b"] == .6
    assert rows[0]["reference_information"] == 2 and rows[0]["information"] == 1


def test_paired_reporting_preserves_type_orientation_and_missing_subjects():
    from tealeaf.sc.path_score_mixed import PathScoreComponents, paired_score_reporting
    reports = [[(1, [.6, .4]), (0, [.2, .8])], [(0, [.3, .7]), (1, [.5, .5])]]
    components = PathScoreComponents(np.zeros((2, 1)), np.zeros((2, 1, 1)), np.ones((2, 1, 1)), np.arange(2), (0, 1), [], reports)
    result = paired_score_reporting(components)
    assert np.allclose(result["effect"], [.3, -.3])
    assert result["complete"] and result["n_reported_subjects"] == 2
    components.reporting_proportions[1][1] = (1, [np.nan, np.nan])
    failed = paired_score_reporting(components)
    assert np.isnan(failed["effect"]).all()
    assert not failed["complete"] and failed["n_reported_subjects"] == 1


def test_signed_score_null_refits_without_modifying_components():
    from tealeaf.sc.path_score_mixed import PathScoreComponents, signed_path_score_p_value
    scores = np.arange(1., 7.)[:, None]
    information = np.ones((6, 1, 1))
    components = PathScoreComponents(scores.copy(), information.copy(), information.copy(), np.arange(6), (0, 1), [], [])
    rng = np.random.default_rng(719)
    signs = rng.choice((-1., 1.), size=6)
    expected = mixed_score_test(scores * signs[:, None], information, information)["p_value"]
    assert signed_path_score_p_value(components, np.random.default_rng(719)) == pytest.approx(expected)
    np.testing.assert_array_equal(components.scores, scores)
    np.testing.assert_array_equal(components.information, information)
    components.levels = (0, 1, 2)
    with pytest.raises(ValueError, match="exactly two"):
        signed_path_score_p_value(components, np.random.default_rng(719))


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
    assert result["n_fitted_subjects"] == 6
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
