import numpy as np
import pytest

from tealeaf.sc.path_reporting import paired_reporting


def inputs(size=10):
    baseline = np.linspace(.05, .7, size)
    proportions = np.stack([np.column_stack([1 - baseline, baseline]), np.column_stack([.9 - baseline, baseline + .1])], axis=1)
    depths = np.column_stack([np.linspace(10, 1000, size), np.linspace(1000, 10, size)])
    covariance = np.array([[(np.diag(p) - np.outer(p, p)) / depth for p, depth in zip(pair, local)] for pair, local in zip(proportions, depths)])
    return proportions, covariance, depths


def test_paired_weights_cancel_subject_baseline():
    proportions, covariance, depths = inputs()
    fitted = paired_reporting(proportions, covariance, depths)
    for name in ("subject arithmetic mean", "paired harmonic depth", "paired random effects"):
        assert np.allclose(fitted[name]["effect"], [-.1, .1])
        assert fitted[name]["weights"].sum() == pytest.approx(1.)
    separately_pooled = np.array([(depths[:, c, None] * proportions[:, c]).sum(axis=0) / depths[:, c].sum() for c in (0, 1)])
    assert separately_pooled[1, 1] - separately_pooled[0, 1] < 0


def test_paired_estimators_reverse_cell_types_and_permute_paths():
    rng = np.random.default_rng(42)
    proportions = rng.dirichlet([2., 3., 5.], size=(15, 2))
    depths = rng.uniform(20, 2000, size=(15, 2))
    covariance = np.array([[(np.diag(p) - np.outer(p, p)) / depth for p, depth in zip(pair, local)] for pair, local in zip(proportions, depths)])
    first = paired_reporting(proportions, covariance, depths)
    reverse = paired_reporting(proportions[:, ::-1], covariance[:, ::-1], depths[:, ::-1])
    order = [2, 0, 1]
    permuted = paired_reporting(proportions[:, :, order], covariance[:, :, order][:, :, :, order], depths)
    for name in first:
        assert np.allclose(first[name]["effect"], -reverse[name]["effect"], atol=1e-7)
        assert np.allclose(first[name]["effect"][order], permuted[name]["effect"], atol=1e-7)


def test_random_effects_caps_deep_outlier_in_paired_differences():
    n = 20
    proportions = np.full((n, 2, 2), .5)
    proportions[:-1, 1] = [.4, .6]
    proportions[-1, 1] = [.9, .1]
    depths = np.full((n, 2), 100.)
    depths[-1] = 1e8
    covariance = np.array([[(np.diag(p) - np.outer(p, p)) / depth for p, depth in zip(pair, local)] for pair, local in zip(proportions, depths)])
    fitted = paired_reporting(proportions, covariance, depths)
    assert fitted["paired harmonic depth"]["effect"][1] < 0
    assert fitted["paired random effects"]["effect"][1] > 0
    assert fitted["paired random effects"]["weights"].max() < .15
    assert fitted["paired random effects"]["between_subject_variance"] > 0


def test_unidentifiable_covariance_uses_explicit_equal_pair_fallback():
    proportions, covariance, depths = inputs()
    covariance[0] = np.nan
    fitted = paired_reporting(proportions, covariance, depths)
    assert fitted["paired random effects"]["fallback"]
    assert np.array_equal(fitted["paired random effects"]["effect"], fitted["subject arithmetic mean"]["effect"])


def test_paired_reporting_validates_shapes():
    proportions, covariance, depths = inputs()
    with pytest.raises(ValueError):
        paired_reporting(proportions, covariance, depths[:, :1])


def test_reporting_control_rejects_cohort_drift_but_records_optimizer_roundoff():
    import pandas as pd
    from extra_scripts.summarize_paired_reporting import check_reporting_control
    new = pd.DataFrame([{"strategy": "subject arithmetic mean A1", "test_id": "a", "effect": "[-0.1, 0.1]", "n_subjects": 8}])
    old = new.assign(strategy="subject-mean A1", effect="[-0.10001, 0.10001]")
    assert check_reporting_control(new, old)["maximum_absolute_component_difference"] == pytest.approx(1e-5)
    with pytest.raises(ValueError, match="subject counts"):
        check_reporting_control(new, old.assign(n_subjects=10))
    with pytest.raises(ValueError, match="subject counts"):
        check_reporting_control(new, old.assign(test_id="other"))
    with pytest.raises(ValueError, match="effect differences"):
        check_reporting_control(new, old.assign(effect="[-0.11, 0.11]"))


def test_failed_omnibus_fits_remain_in_the_fixed_universe():
    import pandas as pd
    from extra_scripts.assess_omnibus_inference_audit import complete_failed_tests
    reference = pd.DataFrame({"block_id": ["a", "b"], "p_value": [.001, .002], "raw_p_value": [.001, .002], "fdr": [.001, .002], "statistic": [10., 9.], "adjusted_effects": ["[[0, 0], [-0.1, 0.1]]"] * 2})
    experimental = reference.iloc[:1].assign(strategy="new")
    output = complete_failed_tests(experimental, reference)
    assert len(output) == 2
    assert output.loc[output.block_id.eq("b"), "p_value"].item() == 1
    assert not output.loc[output.block_id.eq("b"), "fit_available"].item()
    assert output.loc[output.block_id.eq("a"), "fdr"].item() == pytest.approx(.002)


def test_omnibus_shard_loading_accepts_an_explicit_partition_count(tmp_path):
    import pandas as pd
    from extra_scripts.summarize_path_reporting_omnibus import load_shards
    for index in range(2):
        folder = tmp_path / f"shard_{index}"
        folder.mkdir()
        pd.DataFrame({"test_id": [str(index)]}).to_csv(folder / "observed.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="expected 16"):
        load_shards(tmp_path)
    assert len(load_shards(tmp_path, expected=2)) == 2


def test_ec_count_null_preserves_depth_and_removes_cell_type_effect():
    from extra_scripts.audit_ec_count_null import simulate_counts
    from tealeaf.sc.ec_glmm import ECGLMMData
    counts = np.array([[20000., 0.], [0., 30000.], [40000., 0.], [0., 50000.]])
    subjects = np.repeat(["s1", "s2"], 2)
    base = ECGLMMData((counts,), (np.eye(2),), np.ones((4, 1)), subjects)
    simulated = simulate_counts(base, np.array([.8, .2]), subjects, np.random.default_rng(13), subject_scale=0.)
    observed = simulated.counts[0]
    assert np.array_equal(observed.sum(axis=1), counts.sum(axis=1))
    assert np.allclose(observed[:, 1] / observed.sum(axis=1), .2, atol=.01)


def test_omnibus_score_has_no_uniform_prior_depth_contrast():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import omnibus_path_score_test
    depths = np.tile([100., 1000., 10000.], 8)
    labels = np.tile([0, 1, 2], 8)
    subjects = np.repeat(np.arange(8), 3)
    counts = depths[:, None] * [.8, .2]
    data = ECGLMMData((counts,), (np.eye(2),), np.ones((len(counts), 1)), subjects)
    result = omnibus_path_score_test(data, [0, 1], labels, subjects, baseline=np.array([.8, .2]), replicates=16, bootstrap_draws=255)
    assert result["p_value"] == 1
    assert result["statistic"] == 0
    assert result["n_subjects"] == 8
    assert len(result["null"]) == 16
    assert np.allclose(result["cluster_scores"], 0)


def test_omnibus_score_detects_cell_type_effect_and_is_label_invariant():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import omnibus_path_score_test
    labels = np.tile([0, 1, 2], 16)
    subjects = np.repeat(np.arange(16), 3)
    proportions = np.array([[.8, .2], [.5, .5], [.2, .8]])
    depths = np.tile([100., 1000., 10000.], 16)
    counts = depths[:, None] * proportions[labels]
    data = ECGLMMData((counts,), (np.eye(2),), np.ones((len(counts), 1)), subjects)
    kwargs = dict(baseline=np.array([.5, .5]), replicates=16, bootstrap_draws=511, seed=12)
    result = omnibus_path_score_test(data, [0, 1], labels, subjects, **kwargs)
    relabeled = omnibus_path_score_test(data, [0, 1], np.array([2, 0, 1])[labels], subjects, **kwargs)
    assert result["p_value"] < .01
    assert relabeled["p_value"] == result["p_value"]
    assert relabeled["statistic"] == pytest.approx(result["statistic"])
    assert result["one_step_coefficients"].shape == (2, 1)


def test_omnibus_score_is_path_permutation_invariant_with_missing_types():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import omnibus_path_score_test
    labels = np.tile([0, 1, 2], 12)
    subjects = np.repeat(np.arange(12), 3)
    proportions = np.array([[.6, .3, .1], [.3, .3, .4], [.1, .6, .3]])
    counts = 500 * proportions[labels]
    retained = ~((subjects % 3 == 0) & (labels == 2))
    labels, subjects, counts = labels[retained], subjects[retained], counts[retained]
    data = ECGLMMData((counts,), (np.eye(3),), np.ones((len(counts), 1)), subjects)
    kwargs = dict(baseline=np.ones(3) / 3, replicates=8, bootstrap_draws=255, seed=13, null_fit_tolerance=1e-11)
    first = omnibus_path_score_test(data, [0, 1, 2], labels, subjects, **kwargs)
    permuted = omnibus_path_score_test(data, [2, 0, 1], labels, subjects, **kwargs)
    assert first["n_subjects"] == 12
    assert first["p_value"] == permuted["p_value"]
    assert first["statistic"] == pytest.approx(permuted["statistic"])
    assert first["one_step_coefficients"].shape == (2, 2)


def test_free_isoforms_recover_within_path_changes_without_local_effect():
    from tealeaf.sc.differential import fit_free_isoform_paths, fit_path_perturbation
    # Path 0 contains transcripts 0 and 1. The observed shift is within that
    # path, so the conditional path proportions must remain .5/.5.
    counts = (np.array([50., 450., 500.]),)
    baseline = np.array([.45, .05, .5])
    free = fit_free_isoform_paths(counts, (np.eye(3),), baseline, [0, 0, 1], path_pseudocount=1.)
    assert free.converged and free.covariance.identifiable
    assert np.allclose(free.path_proportions, [.5, .5], atol=1e-5)
    assert np.allclose(free.theta, [.05, .45, .5], atol=1e-5)
    # Incomplete/unequal transcript mappings expose the structural concern.
    mapping = np.array([[1., 0., 0.], [0., .1, 0.], [0., 0., 1.]])
    observed = mapping @ np.array([.05, .45, .5])
    observed *= 10000 / observed.sum()
    free = fit_free_isoform_paths((observed,), (mapping,), baseline, [0, 0, 1], path_pseudocount=1.)
    fixed = fit_path_perturbation((observed,), (mapping,), baseline, [0, 0, 1], path_pseudocount=1., path_pseudocount_scaling="total")
    assert free.converged and fixed.converged
    assert np.allclose(free.path_proportions, [.5, .5], atol=1e-3)
    assert abs(fixed.path_proportions[0] - .5) > .2


def test_free_isoform_fit_matches_local_fit_for_singleton_paths():
    from tealeaf.sc.differential import fit_free_isoform_paths, fit_path_perturbation
    counts = (np.array([100., 300., 600.]),)
    baseline = np.array([.2, .5, .3])
    free = fit_free_isoform_paths(counts, (np.eye(3),), baseline, [0, 1, 2], path_pseudocount=16., isoform_pseudocount=0.)
    local = fit_path_perturbation(counts, (np.eye(3),), baseline, [0, 1, 2], path_pseudocount=16., path_pseudocount_scaling="total")
    assert free.converged and free.covariance.identifiable
    assert np.allclose(free.path_proportions, local.path_proportions, atol=1e-5)
    assert np.allclose(free.covariance.covariance, local.covariance.covariance, rtol=1e-3)


def test_path_score_matches_finite_difference_likelihood():
    from tealeaf.sc.path_score import conditional_path_score
    from tealeaf.sc.differential import helmert_basis, _perturbed_theta
    mapping = np.array([[1., .1, .2], [.1, .8, 0.], [0., .2, 1.]])
    theta = np.array([.2, .3, .5])
    counts = np.array([20., 30., 50.])
    score, information = conditional_path_score(theta, [0, 0, 1], (counts,), (mapping,))
    def likelihood(delta):
        perturbed, _ = _perturbed_theta(theta, np.array([0, 0, 1]), helmert_basis(2), np.array([delta]))
        mass = mapping @ perturbed
        return counts @ np.log(mass / mass.sum())
    assert score[0] == pytest.approx((likelihood(1e-5) - likelihood(-1e-5)) / 2e-5, rel=1e-6)
    assert information.shape == (1, 1) and information[0, 0] > 0


def test_paired_score_has_no_coverage_prior_bias_for_exact_binary_null():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import paired_path_score_test
    n = 8
    depths = np.column_stack([np.linspace(100, 200, n), np.linspace(2000, 3000, n)])
    counts = (depths[:, :, None] * np.array([.98, .02])).reshape(-1, 2)
    base = ECGLMMData((counts,), (np.eye(2),), np.ones((2 * n, 1)), np.repeat(np.arange(n), 2))
    result = paired_path_score_test(base, [0, 1], np.tile([0, 1], n), base.clusters, baseline=np.array([.98, .02]), denominator_concentration=32., null_concentration=32.)
    assert result["n_subjects"] == n
    assert np.max(np.abs(result["differences"])) < 1e-12
    assert result["p_value"] == 1.


def test_paired_score_detects_real_effect_and_reverses_orientation():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import paired_path_score_test
    rng = np.random.default_rng(21)
    n = 12
    proportions = np.array([.2, .4])
    counts = np.array([rng.multinomial(1000, [1 - p, p]) for _ in range(n) for p in proportions])
    base = ECGLMMData((counts,), (np.eye(2),), np.ones((2 * n, 1)), np.repeat(np.arange(n), 2))
    labels = np.tile([0, 1], n)
    first = paired_path_score_test(base, [0, 1], labels, base.clusters, baseline=np.array([.7, .3]))
    reverse = paired_path_score_test(base, [0, 1], 1 - labels, base.clusters, baseline=np.array([.7, .3]))
    assert first["p_value"] < 1e-6
    assert np.allclose(first["differences"], -reverse["differences"], atol=1e-6)


def test_subject_centering_eliminates_uniform_prior_coverage_shift():
    from tealeaf.sc.ec_glmm import ECGLMMData
    from tealeaf.sc.path_score import paired_subject_centered_test
    n = 8
    depths = np.column_stack([np.linspace(100, 200, n), np.linspace(2000, 3000, n)])
    counts = (depths[:, :, None] * np.array([.98, .02])).reshape(-1, 2)
    base = ECGLMMData((counts,), (np.eye(2),), np.ones((2 * n, 1)), np.repeat(np.arange(n), 2))
    result = paired_subject_centered_test(base, [0, 1], np.tile([0, 1], n), base.clusters, baseline=np.array([.98, .02]), concentration=32.)
    assert result["n_subjects"] == n
    # The pooled anchor has concentration 1 and finite-depth noise; centering
    # therefore attenuates rather than exactly eliminates the residual bias.
    assert np.max(np.abs(result["differences"])) < .005
def test_linear_simplex_information_matches_logit_pullback_in_interior():
    from tealeaf.sc.differential import conditional_path_information, helmert_basis, simplex_fisher_information

    theta = np.array([.2, .3, .5])
    basis = helmert_basis(3)
    maps = (np.array([[1., .3, .2], [.2, 1., .1], [.1, .2, 2.]]),)
    totals = (200.,)
    logit_information = conditional_path_information(theta, np.arange(3), basis, maps, totals)
    tangent_to_logit = basis.T @ (basis / theta[:, None])
    expected = tangent_to_logit.T @ logit_information @ tangent_to_logit
    np.testing.assert_allclose(simplex_fisher_information(theta, maps, totals), expected, rtol=1e-12, atol=1e-10)


def test_free_path_covariance_remains_identifiable_at_nuisance_boundary():
    from tealeaf.sc.differential import fit_free_isoform_paths
    from tealeaf.sc.path_simulation import simulate_independent_binary_blocks

    data, truth = simulate_independent_binary_blocks(np.random.default_rng(7314159), b_effect=2.2)
    # Row 1 previously converged but failed the numerical covariance rank gate.
    fitted = fit_free_isoform_paths(tuple(count[1] for count in data.counts), data.compatibility, truth["baseline"], truth["path_index"], path_pseudocount=1.)
    assert fitted.converged
    assert fitted.theta.min() < 3e-5
    assert fitted.covariance.identifiable
    assert np.isfinite(fitted.covariance.covariance).all()
    assert 0 < fitted.covariance.covariance[0, 0] < .1
