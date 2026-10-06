import numpy as np
import pytest
from scipy.stats import f_oneway

from tealeaf.sc.omnibus import regression_omnibus, cluster_max_f
from tealeaf.sc.path_reporting import dirichlet_pooling
from tealeaf.sc.differential import fit_profiled_path_perturbation, fit_event_path_perturbation


def test_scalar_regression_matches_anova():
    rng = np.random.default_rng(22)
    labels = np.repeat(np.arange(3), 15)
    values = rng.normal(size=len(labels)) + labels
    design = np.column_stack([np.ones(len(labels)), labels == 1, labels == 2])
    result = regression_omnibus(values[:, None], design, [1, 2])
    expected = f_oneway(*[values[labels == level] for level in range(3)])
    for name in ("trace F", "Pillai", "maximum-coordinate F"):
        assert result[name]["statistic"] == pytest.approx(expected.statistic)
        assert result[name]["p_value"] == pytest.approx(expected.pvalue)


def test_pillai_is_invariant_to_response_coordinates():
    rng = np.random.default_rng(10)
    labels = np.repeat([0, 1, 2], 20)
    values = rng.normal(size=(60, 2)) + labels[:, None] * [.4, .8]
    design = np.column_stack([np.ones(60), labels == 1, labels == 2])
    first = regression_omnibus(values, design, [1, 2])
    second = regression_omnibus(values @ np.array([[2., 1.], [0., .3]]), design, [1, 2])
    assert first["Pillai"]["p_value"] == pytest.approx(second["Pillai"]["p_value"])


def test_pillai_matches_statsmodels_when_available():
    manova = pytest.importorskip("statsmodels.multivariate.manova")
    rng = np.random.default_rng(15)
    labels = np.repeat([0, 1, 2], 18)
    values = rng.normal(size=(54, 3)) + labels[:, None] * [.1, .3, -.4]
    design = np.column_stack([np.ones(54), labels == 1, labels == 2])
    actual = regression_omnibus(values, design, [1, 2])["Pillai"]
    expected = manova.MANOVA(values, design).mv_test([("cell_type", np.array([[0., 1., 0.], [0., 0., 1.]]))]).results["cell_type"]["stat"].loc["Pillai's trace"]
    assert actual["statistic"] == pytest.approx(expected["F Value"])
    assert actual["p_value"] == pytest.approx(expected["Pr > F"])


def test_regression_rejects_rank_deficient_design():
    with pytest.raises(ValueError, match="full-rank"):
        regression_omnibus(np.ones((10, 1)), np.ones((10, 2)), [1])


def test_cluster_statistic_handles_subject_fixed_effects_and_offset():
    from tealeaf.sc.ec_block_glmm import blocked_multilevel_design
    rng = np.random.default_rng(8)
    subjects = np.repeat(np.arange(15), 3)
    labels = np.tile(np.arange(3), 15)
    values = rng.normal(size=(45, 2)) + labels[:, None] * [.2, .4]
    design, tested, _, _ = blocked_multilevel_design(labels, subjects)
    first = cluster_max_f(values, design, tested, subjects)
    shifted = values + rng.normal(size=(15, 2))[subjects]
    second = cluster_max_f(shifted, design, tested, subjects)
    assert first["statistic"] == pytest.approx(second["statistic"])
    assert 0 < first["p_value"] <= 1


def test_dirichlet_pooling_limits_deep_outlier_and_preserves_orientation():
    # Real subject heterogeneity should prevent a million-read outlier from
    # outweighing nine independent subjects in the population mean.
    theta = np.r_[np.repeat(.2, 8), .3, .9, np.linspace(.5, .8, 10)]
    depths = np.r_[np.repeat(100., 9), 1e6, np.repeat(100., 10)]
    proportions = np.column_stack([1 - theta, theta])
    covariances = np.asarray([(np.diag(p) - np.outer(p, p)) / depth for p, depth in zip(proportions, depths)])
    fitted = dirichlet_pooling(proportions, covariances, np.repeat([0, 1], 10), depths)
    assert fitted["converged"]
    assert fitted["means"][0, 1] < .5
    assert fitted["means"][1, 1] > fitted["means"][0, 1]
    assert fitted["subject_precision_weights"][:10].max() / fitted["subject_precision_weights"][:10].sum() < .25
    assert np.allclose(fitted["effective_depth"], depths)


def test_dirichlet_pooling_is_path_permutation_equivariant():
    rng = np.random.default_rng(4)
    proportions = rng.dirichlet([2., 3., 5.], 20)
    covariances = np.asarray([(np.diag(p) - np.outer(p, p)) / 100 for p in proportions])
    labels = np.repeat([0, 1], 10)
    first = dirichlet_pooling(proportions, covariances, labels, np.repeat(100., 20))
    order = [2, 0, 1]
    second = dirichlet_pooling(proportions[:, order], covariances[:, order][:, :, order], labels, np.repeat(100., 20))
    assert np.allclose(first["means"][:, order], second["means"], atol=2e-5)


def test_profiled_path_fit_matches_binary_event_fit():
    mapping = np.array([[1., 0., 0.], [0., 1., 0.], [0., .8, .2]])
    baseline = np.array([.48, .32, .20])
    theta = np.array([.16, .04, .80])
    counts = mapping @ theta
    counts *= 10000 / counts.sum()
    actual = fit_profiled_path_perturbation((counts,), (mapping,), baseline, [0, 1, -1], path_pseudocount=1.)
    expected = fit_event_path_perturbation((counts,), (mapping,), baseline, [0, 1, -1], path_pseudocount=1.)
    assert actual.converged and actual.path_proportions[0] > .75
    assert np.allclose(actual.path_proportions, expected.path_proportions, atol=1e-5)
    assert np.allclose(actual.theta, expected.theta, atol=1e-5)
    assert np.allclose(actual.covariance.covariance, expected.covariance.covariance, rtol=1e-3)


def test_profiled_path_fit_recovers_multivariate_counts_and_nuisance_mass():
    fit = fit_profiled_path_perturbation((np.array([200., 300., 500., 4000.]),), (np.eye(4),), np.array([.02, .08, .4, .5]), [0, 1, 2, -1], path_pseudocount=3.)
    assert fit.converged
    assert np.allclose(fit.path_proportions, np.array([201., 301., 501.]) / 1003, atol=1e-5)
    assert fit.theta[:3].sum() == pytest.approx(.2, abs=1e-4)
    assert fit.covariance.identifiable


@pytest.mark.parametrize("counts", [np.array([-1., 2.]), np.zeros(2)])
def test_profiled_path_fit_rejects_invalid_count_data(counts):
    with pytest.raises(ValueError, match="EC count"):
        fit_profiled_path_perturbation((counts,), (np.eye(2),), np.array([.5, .5]), [0, 1])


def test_held_permutation_families_do_not_affect_calibration():
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import calibrate_omnibus
    rng = np.random.default_rng(12)
    observed = pd.DataFrame({"test_id": ["a", "b", "c"], "strategy": "trace F", "p_value": [.001, .02, .5], "statistic": [20., 4., .1], "converged": True, "n_subjects": 12, "degrees_of_freedom": 1})
    null = pd.DataFrame([{"test_id": test_id, "strategy": "trace F", "replicate": replicate, "p_value": rng.uniform(), "statistic": 1., "n_subjects": 12, "degrees_of_freedom": 1} for test_id in ("a", "b", "c") for replicate in range(64)])
    first, _, _ = calibrate_omnibus(observed, null)
    null.loc[null.replicate.ge(32), "p_value"] = 0.
    second, held, _ = calibrate_omnibus(observed, null)
    assert np.array_equal(first.p_value, second.p_value)
    assert np.all(held.p_value < .05)


def test_reporting_fallback_keeps_failed_events():
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import add_dirichlet_fallback
    table = pd.DataFrame([{"test_id": "a", "strategy": "subject-mean A1", "effect": "[-0.2, 0.2]", "converged": True}, {"test_id": "a", "strategy": "effective-count Dirichlet pooling", "effect": "[NaN, NaN]", "converged": False}])
    fallback = add_dirichlet_fallback(table).query("strategy == 'Dirichlet pooling with subject-mean fallback'").iloc[0]
    assert fallback.effect == "[-0.2, 0.2]"
    assert fallback.report_fallback and fallback.converged


def test_paired_combination_keeps_nonsignificant_family_members():
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import paired_combination_omnibus
    paired = pd.DataFrame({"block_id": ["b", "b", "b"], "p_value": [.01, .6, .9], "raw_p_value": [.005, .6, .9], "n_subjects": 12, "converged": True})
    template = pd.DataFrame({"block_id": ["b"], "strategy": "trace F", "p_value": [.02], "raw_p_value": [.02], "statistic": [1.]})
    combined = paired_combination_omnibus(paired, template)
    selected = combined.loc[combined.strategy.eq("Bonferroni paired omnibus")].iloc[0]
    assert selected.p_value == pytest.approx(.03)
    assert selected.raw_p_value == pytest.approx(.015)
    assert selected.n_paired_contrasts == 3


def test_paired_combination_uses_explicit_reporting_control():
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import paired_combination_omnibus
    paired = pd.DataFrame({"block_id": ["b"], "p_value": [.01], "raw_p_value": [.01], "n_subjects": [12], "converged": [True]})
    template = pd.DataFrame({"block_id": ["b", "b"], "strategy": ["CR2 wild", "null-variance Wald A32"], "adjusted_effects": ["A32 arrays", "A1 arrays"]})
    combined = paired_combination_omnibus(paired, template)
    assert combined.adjusted_effects.eq("A1 arrays").all()
