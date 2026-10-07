import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative
from scipy.special import digamma, softmax

from tealeaf.sc.path_dispersion import dm_coefficient_information, cox_reid_profile, corrected_dm_test
from tealeaf.sc.path_pooling import joint_path_dm_test


@pytest.mark.parametrize("paths,predictors", [(2, 2), (3, 2), (4, 3)])
@pytest.mark.parametrize("concentration", [.25, 20., 102400.])
def test_information_matches_independent_score_finite_differences(paths, predictors, concentration):
    rng = np.random.default_rng(138)
    design = np.column_stack([np.ones(9), rng.normal(size=(9, predictors - 1))])
    counts = rng.gamma(2, 4, size=(9, paths))
    coefficients = rng.normal(0, .4, size=(predictors, paths - 1))

    def score(parameters):
        means = softmax(np.column_stack([design @ parameters.reshape(coefficients.shape), np.zeros(len(counts))]), axis=1)
        alpha = concentration * means
        gradient_alpha = concentration * (digamma(alpha + counts) - digamma(alpha))
        gradient_logits = means * (gradient_alpha - (means * gradient_alpha).sum(axis=1, keepdims=True))
        return (design.T @ gradient_logits[:, :-1]).ravel()

    numeric = -approx_derivative(score, coefficients.ravel(), method="3-point", rel_step=1e-3 if concentration > 1e4 else 1e-4)
    analytic = dm_coefficient_information(counts, design, coefficients, concentration)
    assert np.allclose(analytic, numeric, atol=2e-5, rtol=2e-5)
    assert np.allclose(analytic, analytic.T, atol=1e-12)


def regression_data():
    rng = np.random.default_rng(151)
    labels = np.repeat([0, 1], 16)
    design = np.column_stack([np.ones(len(labels)), labels])
    means = np.asarray([[.7, .2, .1], [.3, .35, .35]])[labels]
    counts = np.asarray([rng.multinomial(200, rng.dirichlet(20 * mean)) for mean in means])
    return counts, design


def test_profile_precision_is_shared_between_hypotheses_not_in_lr_penalty():
    from tealeaf.sc import differential
    counts, design = regression_data()
    grid = [5., 10., 20., 40., 80.]
    profile = cox_reid_profile(counts, design, concentrations=grid)
    result = corrected_dm_test(counts, design[:, :1], design, concentrations=grid)
    assert result["null_converged"] and result["alternative_converged"]
    assert result["null_concentration"] == result["alternative_concentration"] == profile["concentration"]
    assert result["p_value"] < .001
    assert result["f_p_value"] >= result["p_value"]
    alternative = differential._dirichlet_multinomial_fit(counts, design, initial=result["alternative_coefficients"].ravel(), fixed_concentration=profile["concentration"])
    assert result["statistic"] == pytest.approx(2 * (result["null_objective"] - alternative["objective"]), abs=1e-5)


def test_fixed_precision_is_path_permutation_equivariant():
    counts, design = regression_data()
    result = corrected_dm_test(counts, design[:, :1], design, concentration=20.)
    reverse = corrected_dm_test(counts[:, [2, 0, 1]], design[:, :1], design, concentration=20.)
    assert result["null_converged"] and result["alternative_converged"]
    assert reverse["null_converged"] and reverse["alternative_converged"]
    assert reverse["statistic"] == pytest.approx(result["statistic"], rel=1e-6)


def test_corrected_precision_and_test_are_path_reference_equivariant():
    counts, design = regression_data()
    grid = [5., 10., 20., 40., 80.]
    result = corrected_dm_test(counts, design[:, :1], design, concentrations=grid)
    reverse = corrected_dm_test(counts[:, [2, 0, 1]], design[:, :1], design, concentrations=grid)
    assert reverse["alternative_concentration"] == result["alternative_concentration"]
    assert reverse["statistic"] == pytest.approx(result["statistic"], rel=1e-6)


def test_profile_adjustment_rejects_invalid_matrix_shapes():
    for counts, design in ((np.ones(12), np.ones((12, 1))), (np.ones((12, 2)), np.ones(12))):
        with pytest.raises(ValueError, match="matrices"):
            cox_reid_profile(counts, design)


def test_invalid_profile_point_is_not_silently_removed(monkeypatch):
    from tealeaf.sc import path_dispersion
    counts, design = regression_data()
    original = path_dispersion.dm_coefficient_information

    def negative_information(counts, design, coefficients, concentration):
        if concentration == 10.:
            # An even-dimensional negative-definite matrix has positive det.
            return -np.eye(coefficients.size)
        return original(counts, design, coefficients, concentration)

    monkeypatch.setattr(path_dispersion, "dm_coefficient_information", negative_information)
    with pytest.raises(ValueError, match="invalid points"):
        cox_reid_profile(counts, design, concentrations=[5., 10., 20.])


def test_failed_corrected_fit_has_no_significant_tail(monkeypatch):
    from tealeaf.sc import path_dispersion
    original = path_dispersion._fit_at_concentration

    def failed(*args, **kwargs):
        return {**original(*args, **kwargs), "converged": False}

    monkeypatch.setattr(path_dispersion, "_fit_at_concentration", failed)
    counts, design = regression_data()
    result = corrected_dm_test(counts, design[:, :1], design, concentration=20.)
    assert result["p_value"] == result["f_p_value"] == 1
    assert result["statistic"] == 0


@pytest.mark.parametrize("grid", [[1, 1, 2], [0, 1, 2], [1, 2], [1, np.inf, 3]])
def test_invalid_profile_grid_is_rejected(grid):
    counts, design = regression_data()
    with pytest.raises(ValueError, match="ordered positive"):
        cox_reid_profile(counts, design, concentrations=grid)


def test_joint_correction_rejects_null_reuse_and_bad_mode():
    counts, _ = regression_data()
    quantified = {"counts": counts, "subjects": np.tile(np.arange(16), 2), "labels": np.repeat([0, 1], 16), "effective_depths": counts.sum(axis=1)}
    with pytest.raises(ValueError, match="fresh fits"):
        joint_path_dm_test(quantified, dispersion_method="cox_reid", fitted_null={})
    with pytest.raises(ValueError, match="diagnostic mode"):
        joint_path_dm_test(quantified, dispersion_method="fixed")
    with pytest.raises(ValueError, match="unsupported"):
        joint_path_dm_test(quantified, dispersion_method="ml", concentration=20)


def test_corrected_permutation_reselects_precision_and_keeps_failed_null(monkeypatch):
    from extra_scripts import audit_path_reporting_omnibus as audit
    from tealeaf.sc.ec_glmm import ECGLMMData
    labels, subjects = np.tile([0, 1], 12), np.repeat(np.arange(12), 2)
    data = ECGLMMData((np.tile([40., 60.], (24, 1)),), (np.eye(2),), np.ones((24, 1)), subjects)
    calls = []
    def fitted(quantified, **options):
        calls.append(options)
        if "labels" in options:
            raise ValueError("permutation profile failure")
        return {"p_value": .02, "f_p_value": .03, "statistic": 5., "degrees_of_freedom": 1, "converged": True, "n_subjects": 12, "n_observations": 24, "levels": np.asarray([0, 1]), "standardized_means": np.asarray([[.7, .3], [.3, .7]])}
    monkeypatch.setattr(audit, "joint_path_dm_test", fitted)
    observed, null, _ = audit.joint_dm_reports(data, np.ones(2) / 2, np.asarray([0, 1]), labels, subjects, 3, 0, dispersion_method="cox_reid")
    assert len(observed) == 2 and len(null) == 6
    assert sorted(value["p_value"] for value in observed.values()) == [.02, .03]
    assert all(value["p_value"] == 1 and not value["converged"] for value in null)
    assert all(call.get("fitted_null") is None for call in calls)
    assert all(call["dispersion_method"] == "cox_reid" for call in calls)
