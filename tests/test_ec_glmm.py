"""Tests for GLMMs that retain equivalence-class counts as observations."""

import numpy as np
import pytest

from tealeaf.sc import ec_block_glmm, ec_glmm, ec_glmm_full
from extra_scripts import run_ec_glmm


def simulated_data(seed=4, effect=0.9, family="multinomial"):
    rng = np.random.default_rng(seed)
    mice = 10
    clusters = np.repeat(np.arange(mice), 2)
    condition = np.repeat(np.arange(mice) >= mice // 2, 2)
    cell_type = np.tile([0, 1], mice)
    design = np.c_[1 - cell_type, cell_type, condition]
    mappings = (
        np.array([[1.0, 0.0], [0.35, 0.65], [0.0, 1.0]]),
        np.array([[1.0, 0.0], [0.7, 0.3], [0.0, 1.0]]),
    )
    random_effect = rng.normal(0.0, 0.35, mice)
    logits = -0.3 * design[:, 0] + 0.25 * design[:, 1]
    logits += effect * condition + random_effect[clusters]
    abundance = np.c_[np.exp(logits), np.ones(len(logits))]
    counts = []
    for mapping in mappings:
        probability = abundance @ mapping.T
        probability /= probability.sum(axis=1, keepdims=True)
        if family == "multinomial":
            observed = np.asarray([rng.multinomial(180, row) for row in probability])
        else:
            latent = np.asarray([rng.dirichlet(30.0 * row) for row in probability])
            observed = np.asarray([rng.multinomial(180, row) for row in latent])
        counts.append(observed)
    return ec_glmm.ECGLMMData(tuple(counts), mappings, design, clusters)


def test_tilted_bound_is_an_upper_bound():
    jax, jnp, _ = ec_glmm._jax()
    means = jnp.asarray([[0.2, -0.4, 0.0]])
    variances = jnp.asarray([[0.3, 0.7, 0.0]])
    bound = float(ec_glmm._tilted_logsumexp_bound(jnp, means, variances)[0])
    rng = np.random.default_rng(2)
    samples = rng.normal(np.asarray(means), np.sqrt(np.asarray(variances)), (200_000, 3))
    expectation = float(np.mean(jax.scipy.special.logsumexp(samples, axis=1)))
    assert bound >= expectation - 3e-3
    assert bound - expectation < 0.08


def test_full_tilted_bound_handles_correlated_logits():
    jax, jnp, _ = ec_glmm._jax()
    means = jnp.asarray([[0.2, -0.4, 0.0]])
    covariance = jnp.asarray(
        [[[0.3, -0.18, 0.0], [-0.18, 0.7, 0.0], [0.0, 0.0, 0.0]]]
    )
    bound = float(
        ec_glmm_full._tilted_logsumexp_bound(jnp, means, covariance)[0]
    )
    rng = np.random.default_rng(2)
    samples = rng.multivariate_normal(
        np.asarray(means[0]), np.asarray(covariance[0]), size=200_000
    )
    expectation = float(np.mean(jax.scipy.special.logsumexp(samples, axis=1)))
    assert bound >= expectation - 3e-3
    assert bound - expectation < 0.08


def test_multinomial_methods_recover_ambiguous_ec_effect():
    data = simulated_data()
    tilted = ec_glmm.fit_tilted_variational(data, max_iter=120)
    assert tilted["converged"]
    assert tilted["coefficients"][2, 0] == pytest.approx(0.9, abs=0.35)
    scores = ec_glmm.evaluate_variational_objectives(
        data, tilted, samples=64
    )
    assert np.isfinite(list(scores.values())).all()
    assert scores["renyi_0.5"] >= scores["elbo_mc"]
    renyi = ec_glmm.fit_variational(
        data,
        alpha=0.5,
        samples=128,
        seed=9,
        initial=tilted["parameters"],
        max_iter=120,
    )
    assert renyi["converged"]
    assert renyi["coefficients"][2, 0] == pytest.approx(0.9, abs=0.4)


def test_laplace_supports_dirichlet_multinomial_ec_counts():
    data = simulated_data(family="dirichlet_multinomial")
    fit = ec_glmm.fit_laplace(
        data,
        family="dirichlet_multinomial",
        max_iter=80,
        mode_steps=20,
    )
    assert fit["converged"]
    assert np.isfinite(fit["objective"])
    assert fit["concentration"] > 0
    assert fit["coefficients"][2, 0] == pytest.approx(0.9, abs=0.55)
    initial = ec_glmm_full.variational_warm_start(
        fit, data.design.shape[1]
    )
    full = ec_glmm_full.fit_variational(
        data,
        family="dirichlet_multinomial",
        objective="monte_carlo",
        samples=32,
        initial=initial,
        max_iter=3,
    )
    assert np.isfinite(full["objective"])
    assert full["concentration"] > 0


def test_tensor_laplace_nested_lrt_detects_effect():
    data = simulated_data(effect=1.2)
    dimension = data.n_isoforms - 1
    null_tensor = ec_block_glmm.unrestricted_tensor(data.design[:, :2], dimension)
    alternative_tensor = ec_block_glmm.unrestricted_tensor(data.design, dimension)

    def with_tensor(tensor):
        return ec_glmm.ECGLMMData(
            data.counts,
            data.compatibility,
            data.design,
            data.clusters,
            fixed_effect_tensor=tensor,
        )

    tests = ec_block_glmm.nested_laplace_tests(
        with_tensor(null_tensor),
        with_tensor(alternative_tensor),
        max_iter=80,
        mode_steps=20,
    )
    assert tests["null_fit"]["converged"]
    assert tests["alternative_fit"]["converged"]
    assert tests["tested_count"] == 1
    assert tests["lrt_p_value"] < 0.05
    assert tests["bic_log_bayes_factor"] > 0


def test_fixed_effect_design_separates_tested_contrast():
    groups = [
        "a__control__m1",
        "b__control__m1",
        "a__case__m2",
        "b__case__m2",
    ]
    null, alternative, clusters, metadata = run_ec_glmm.fixed_effect_design(
        groups, "condition"
    )
    assert null.shape == (4, 2)
    assert alternative.shape == (4, 3)
    np.testing.assert_array_equal(clusters, ["m1", "m1", "m2", "m2"])


def test_laplace_posterior_warm_starts_variational_fit():
    data = simulated_data()
    laplace = ec_glmm.fit_laplace(data, max_iter=50, mode_steps=15)
    initial = ec_glmm.variational_warm_start(
        laplace, data.design.shape[1]
    )
    fit = ec_glmm.fit_tilted_variational(
        data, initial=initial, max_iter=80
    )
    assert fit["converged"]
    assert fit["coefficients"][2, 0] == pytest.approx(0.9, abs=0.35)


def test_full_covariance_and_observation_noise_fits():
    data = simulated_data()
    laplace = ec_glmm.fit_laplace(
        data, observation_noise=True, max_iter=30, mode_steps=12
    )
    assert laplace["observation_noise_sd"] > 0
    assert laplace["latent_covariance"].shape[1] > data.n_isoforms - 1
    initial = ec_glmm_full.variational_warm_start(
        laplace, data.design.shape[1]
    )
    fit = ec_glmm_full.fit_variational(
        data,
        objective="tilted",
        observation_noise=True,
        initial=initial,
        max_iter=5,
    )
    assert np.isfinite(fit["objective"])
    assert fit["random_effect_covariance"].shape == (10, 1, 1)
    scores = ec_glmm_full.evaluate_objectives(data, fit, samples=32)
    assert np.isfinite(list(scores.values())).all()


def test_laplace_supports_paired_random_slopes():
    rng = np.random.default_rng(17)
    subjects = 18
    clusters = np.repeat(np.arange(subjects), 2)
    labels = np.tile([-0.5, 0.5], subjects)
    random_design = np.column_stack((np.ones(2 * subjects), labels))
    design = np.column_stack((np.ones(2 * subjects), labels))
    mapping = np.eye(2)
    intercepts = rng.normal(0.0, 0.25, subjects)
    slopes = rng.normal(0.0, 0.8, subjects)
    logits = 0.5 * labels + intercepts[clusters] + labels * slopes[clusters]
    probabilities = np.column_stack((np.exp(logits), np.ones(len(logits))))
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    counts = np.asarray([rng.multinomial(400, row) for row in probabilities])
    data = ec_glmm.ECGLMMData(
        (counts,),
        (mapping,),
        design,
        clusters,
        random_effect_design=random_design,
    )
    fit = ec_glmm.fit_laplace(
        data, random_slopes=True, max_iter=100, mode_steps=20
    )
    assert fit["converged"]
    assert fit["random_effect_sds"].shape == (2,)
    assert fit["random_slope_sd"] > 0.2
    assert fit["random_term_mean"].shape == (subjects, 2, 1)
