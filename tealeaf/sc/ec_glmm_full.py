"""Full-covariance variational inference for EC-count mixed models."""

from __future__ import annotations

import numpy as np
import scipy.optimize
import scipy.special

from . import ec_glmm


def _latent_layout(data, observation_noise):
    """Return cluster indices and row-to-latent selection matrices."""
    _, cluster_index = np.unique(data.clusters, return_inverse=True)
    n_clusters = int(cluster_index.max()) + 1
    dimension = data.n_isoforms - 1
    slots = np.zeros(len(cluster_index), dtype=int)
    sizes = np.zeros(n_clusters, dtype=int)
    for row, cluster in enumerate(cluster_index):
        slots[row] = sizes[cluster]
        sizes[cluster] += 1
    blocks = 1 + (int(sizes.max()) if observation_noise else 0)
    latent_dimension = blocks * dimension
    selection = np.zeros((len(cluster_index), dimension, latent_dimension))
    selection[:, :, :dimension] = np.eye(dimension)[None, :, :]
    if observation_noise:
        for row, slot in enumerate(slots):
            start = (1 + slot) * dimension
            selection[row, :, start : start + dimension] = np.eye(dimension)
    return cluster_index, selection, blocks


def _cholesky_from_packed(jnp, packed, rows, columns, diagonal):
    transformed = jnp.where(diagonal[None, :], jnp.exp(packed), packed)
    dimension = int((np.sqrt(8 * packed.shape[1] + 1) - 1) / 2)
    result = jnp.zeros(
        (packed.shape[0], dimension, dimension),
        dtype=packed.dtype,
    )
    return result.at[:, rows, columns].set(transformed)


def _pack_cholesky(covariance):
    covariance = np.asarray(covariance, dtype=float)
    dimension = covariance.shape[-1]
    rows, columns = np.tril_indices(dimension)
    covariance = 0.5 * (covariance + np.swapaxes(covariance, 1, 2))
    minimum_eigenvalue = np.linalg.eigvalsh(covariance)[:, 0]
    jitter = np.maximum(1e-8, 1e-8 - minimum_eigenvalue)
    cholesky = np.linalg.cholesky(
        covariance + jitter[:, None, None] * np.eye(dimension)[None, :, :]
    )
    packed = cholesky[:, rows, columns]
    packed[:, rows == columns] = np.log(packed[:, rows == columns])
    return packed.ravel()


def _tilted_logsumexp_bound(jnp, means, covariances, iterations=25):
    """Upper-bound expected log-sum-exp for correlated Gaussian logits."""
    diagonal = jnp.diagonal(covariances, axis1=1, axis2=2)
    local = __import__("jax").nn.softmax(means + 0.5 * diagonal, axis=1)

    def update(_, value):
        covariance_times_local = jnp.einsum("nij,nj->ni", covariances, value)
        return __import__("jax").nn.softmax(
            means + 0.5 * diagonal - covariance_times_local, axis=1
        )

    local = __import__("jax").lax.fori_loop(0, int(iterations), update, local)
    covariance_times_local = jnp.einsum("nij,nj->ni", covariances, local)
    quadratic = 0.5 * jnp.sum(local * covariance_times_local, axis=1)
    adjusted = means + 0.5 * diagonal - covariance_times_local
    return quadratic + __import__("jax.scipy", fromlist=["special"]).special.logsumexp(
        adjusted, axis=1
    )


def fit_variational(
    data,
    *,
    family="multinomial",
    objective="tilted",
    alpha=1.0,
    observation_noise=False,
    samples=64,
    seed=0,
    initial=None,
    max_iter=300,
):
    """Fit a cluster-factorized, full-covariance Gaussian posterior.

    ``objective`` is either ``"tilted"`` for the analytic multinomial ELBO
    lower bound or ``"monte_carlo"`` for an exact-likelihood ELBO/Renyi
    objective.  With observation noise, each cluster block jointly represents
    its shared random intercept and all observation-level residual effects.
    """
    if objective not in ("tilted", "monte_carlo"):
        raise ValueError("objective must be 'tilted' or 'monte_carlo'")
    if objective == "tilted" and family != "multinomial":
        raise ValueError("the tilted objective is only available for multinomial fits")
    if family not in ("multinomial", "dirichlet_multinomial"):
        raise ValueError(f"unknown EC likelihood family: {family}")
    if observation_noise and family != "multinomial":
        raise ValueError("observation noise is only available for multinomial fits")
    if not 0 < float(alpha) <= 1:
        raise ValueError("alpha must lie in (0, 1]")
    (
        jax,
        jnp,
        jsp,
        counts,
        mappings,
        design,
        cluster_index,
        n_clusters,
    ) = ec_glmm._prepared_arrays(data)
    cluster_index_np, selection_np, latent_blocks = _latent_layout(
        data, observation_noise
    )
    np.testing.assert_array_equal(cluster_index_np, np.asarray(cluster_index))
    selection = jnp.asarray(selection_np)
    dimension = data.n_isoforms - 1
    latent_dimension = selection.shape[2]
    fixed_effect_tensor = (
        None
        if data.fixed_effect_tensor is None
        else jnp.asarray(data.fixed_effect_tensor)
    )
    coefficient_count = (
        design.shape[1] * dimension
        if fixed_effect_tensor is None
        else fixed_effect_tensor.shape[2]
    )
    mean_count = n_clusters * latent_dimension
    triangle_rows_np, triangle_columns_np = np.tril_indices(latent_dimension)
    triangle_rows = jnp.asarray(triangle_rows_np)
    triangle_columns = jnp.asarray(triangle_columns_np)
    triangle_diagonal = jnp.asarray(triangle_rows_np == triangle_columns_np)
    cholesky_count = n_clusters * len(triangle_rows_np)
    membership = jax.nn.one_hot(cluster_index, n_clusters, dtype=jnp.float64)

    if initial is None:
        covariance = np.tile(
            (0.3**2) * np.eye(latent_dimension)[None, :, :],
            (n_clusters, 1, 1),
        )
        parameters = np.r_[
            np.zeros(coefficient_count + mean_count),
            _pack_cholesky(covariance),
            np.log(0.5),
        ]
        if observation_noise:
            parameters = np.r_[parameters, np.log(0.3)]
        if family == "dirichlet_multinomial":
            parameters = np.r_[parameters, np.log(100.0)]
    else:
        parameters = np.asarray(initial, dtype=float)

    def unpack(value):
        coefficients = value[:coefficient_count]
        if fixed_effect_tensor is None:
            coefficients = coefficients.reshape(design.shape[1], dimension)
        cursor = coefficient_count
        means = value[cursor : cursor + mean_count].reshape(
            n_clusters, latent_dimension
        )
        cursor += mean_count
        packed = value[cursor : cursor + cholesky_count].reshape(
            n_clusters, len(triangle_rows_np)
        )
        cursor += cholesky_count
        cholesky = _cholesky_from_packed(
            jnp, packed, triangle_rows, triangle_columns, triangle_diagonal
        )
        log_mouse_sd = value[cursor]
        cursor += 1
        log_noise_sd = value[cursor] if observation_noise else value[0] * 0.0
        cursor += int(observation_noise)
        log_kappa = (
            value[cursor]
            if family == "dirichlet_multinomial"
            else value[0] * 0.0
        )
        return (
            coefficients,
            means,
            cholesky,
            log_mouse_sd,
            log_noise_sd,
            log_kappa,
        )

    def fixed_logits(coefficients):
        if fixed_effect_tensor is None:
            return design @ coefficients
        return jnp.einsum("ndp,p->nd", fixed_effect_tensor, coefficients)

    def prior_log_sds(log_mouse_sd, log_noise_sd):
        values = jnp.full((latent_blocks, dimension), log_mouse_sd)
        if observation_noise:
            values = values.at[1:].set(log_noise_sd)
        return values.ravel()

    def moments(value):
        (
            coefficients,
            latent_means,
            cholesky,
            log_mouse_sd,
            log_noise_sd,
            _,
        ) = unpack(value)
        covariance = cholesky @ jnp.swapaxes(cholesky, 1, 2)
        row_means = jnp.einsum(
            "ndi,ni->nd", selection, latent_means[cluster_index]
        )
        row_covariance = jnp.einsum(
            "nai,nij,nbj->nab",
            selection,
            covariance[cluster_index],
            selection,
        )
        free_means = fixed_logits(coefficients) + row_means
        means = jnp.concatenate(
            (free_means, jnp.zeros((len(design), 1), dtype=free_means.dtype)), axis=1
        )
        covariances = jnp.pad(row_covariance, ((0, 0), (0, 1), (0, 1)))
        return (
            coefficients,
            latent_means,
            cholesky,
            covariance,
            means,
            covariances,
            log_mouse_sd,
            log_noise_sd,
        )

    def tilted_bound(value):
        (
            _,
            latent_means,
            cholesky,
            covariance,
            means,
            covariances,
            log_mouse_sd,
            log_noise_sd,
        ) = moments(value)
        result = jnp.asarray(0.0, dtype=means.dtype)
        for observed, mapping in zip(counts, mappings):
            totals = jnp.sum(observed, axis=1)
            constants = jsp.special.gammaln(totals + 1.0) - jnp.sum(
                jsp.special.gammaln(observed + 1.0), axis=1
            )
            log_mapping = jnp.where(mapping > 0, jnp.log(mapping), -jnp.inf)
            numerator = jsp.special.logsumexp(
                means[:, None, :] + log_mapping[None, :, :], axis=2
            )
            lengths = jnp.sum(mapping, axis=0)
            shifted_means = means + jnp.log(lengths)[None, :]
            denominator = _tilted_logsumexp_bound(
                jnp, shifted_means, covariances
            )
            result += jnp.sum(
                constants + jnp.sum(observed * numerator, axis=1) - totals * denominator
            )
        log_prior_sds = prior_log_sds(log_mouse_sd, log_noise_sd)
        prior_variances = jnp.exp(2.0 * log_prior_sds)
        trace = jnp.sum(
            jnp.diagonal(covariance, axis1=1, axis2=2) / prior_variances[None, :]
        )
        quadratic = jnp.sum(latent_means * latent_means / prior_variances[None, :])
        logdet_prior = 2.0 * n_clusters * jnp.sum(log_prior_sds)
        logdet_q = 2.0 * jnp.sum(
            jnp.log(jnp.diagonal(cholesky, axis1=1, axis2=2))
        )
        kl = 0.5 * (
            trace + quadratic - n_clusters * latent_dimension + logdet_prior - logdet_q
        )
        return result - kl

    rng = np.random.default_rng(seed)
    half = max(1, int(samples) // 2)
    standard_noise = rng.standard_normal((half, n_clusters, latent_dimension))
    standard_noise = np.concatenate((standard_noise, -standard_noise), axis=0)[
        : int(samples)
    ]
    standard_noise = jnp.asarray(standard_noise)

    def log_weights(value):
        (
            coefficients,
            means,
            cholesky,
            log_mouse_sd,
            log_noise_sd,
            log_kappa,
        ) = unpack(value)
        effects = means[None, :, :] + jnp.einsum(
            "mij,smj->smi", cholesky, standard_noise
        )
        fixed = fixed_logits(coefficients)
        log_prior_sds = prior_log_sds(log_mouse_sd, log_noise_sd)

        def one_sample(effect, noise):
            random_logits = jnp.einsum(
                "ndi,ni->nd", selection, effect[cluster_index]
            )
            free_logits = fixed + random_logits
            logits = jnp.concatenate(
                (free_logits, jnp.zeros((len(design), 1), dtype=free_logits.dtype)),
                axis=1,
            )
            likelihood_rows = ec_glmm._log_likelihood_rows(
                jnp, jsp, counts, mappings, logits, family, log_kappa
            )
            likelihood = membership.T @ likelihood_rows
            log_prior = -0.5 * jnp.sum(
                (effect / jnp.exp(log_prior_sds)[None, :]) ** 2, axis=1
            )
            log_prior -= jnp.sum(log_prior_sds)
            log_prior -= 0.5 * latent_dimension * jnp.log(2.0 * jnp.pi)
            log_q = -0.5 * jnp.sum(noise * noise, axis=1)
            log_q -= 0.5 * latent_dimension * jnp.log(2.0 * jnp.pi)
            log_q -= jnp.sum(
                jnp.log(jnp.diagonal(cholesky, axis1=1, axis2=2)), axis=1
            )
            return likelihood + log_prior - log_q

        return jax.vmap(one_sample)(effects, standard_noise)

    def monte_carlo_bound(value):
        weights = log_weights(value)
        if float(alpha) == 1.0:
            return jnp.sum(jnp.mean(weights, axis=0))
        scale = 1.0 - float(alpha)
        cluster_bound = jsp.special.logsumexp(scale * weights, axis=0) / scale
        cluster_bound -= jnp.log(len(weights)) / scale
        return jnp.sum(cluster_bound)

    bound = tilted_bound if objective == "tilted" else monte_carlo_bound
    if int(max_iter) > 0:
        value_and_gradient = jax.jit(
            jax.value_and_grad(lambda value: -bound(value))
        )

        def scipy_objective(value):
            result, gradient = value_and_gradient(jnp.asarray(value))
            return float(result), np.asarray(gradient, dtype=float)
    else:
        value_only = jax.jit(lambda value: -bound(value))

        def scipy_objective(value):
            result = value_only(jnp.asarray(value))
            return float(result), np.full_like(value, np.nan, dtype=float)

    non_cholesky_count = coefficient_count + mean_count
    # Packed entries are ordered within clusters, not by triangular coordinate.
    bounds = [(-20.0, 20.0)] * non_cholesky_count + [
        (-8.0, 3.0) if diagonal else (-20.0, 20.0)
        for _ in range(n_clusters)
        for diagonal in (triangle_rows_np == triangle_columns_np)
    ] + [(-8.0, 3.0)] * (1 + int(observation_noise))
    if family == "dirichlet_multinomial":
        bounds.append((-6.0, np.log(1e7)))
    if len(bounds) != len(parameters):
        raise ValueError("initial full-covariance parameter vector has the wrong size")
    if int(max_iter) > 0:
        result = scipy.optimize.minimize(
            scipy_objective,
            parameters,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": int(max_iter), "ftol": 1e-8, "gtol": 1e-4},
        )
        optimized = np.asarray(result.x)
        success = bool(result.success)
        iterations = int(result.nit)
        message = str(result.message)
    else:
        optimized = np.asarray(parameters)
        success = True
        iterations = 0
        message = "objective evaluation only"
    (
        coefficients,
        means,
        cholesky,
        log_mouse_sd,
        log_noise_sd,
        log_kappa,
    ) = unpack(jnp.asarray(optimized))
    covariance = np.asarray(cholesky @ jnp.swapaxes(cholesky, 1, 2))
    final, gradient = scipy_objective(optimized)
    if objective == "monte_carlo":
        final_log_weights = np.asarray(log_weights(jnp.asarray(optimized)))
        importance_scale = 1.0 if float(alpha) == 1.0 else 1.0 - float(alpha)
        normalized = scipy.special.softmax(importance_scale * final_log_weights, axis=0)
        cluster_ess = 1.0 / np.sum(normalized * normalized, axis=0)
    else:
        cluster_ess = None
    return {
        "method": (
            "tilted_elbo_full"
            if objective == "tilted"
            else ("elbo_mc_full" if float(alpha) == 1.0 else f"renyi_{alpha:g}_full")
        ),
        "family": family,
        "alpha": float(alpha),
        "objective": -final,
        "parameters": optimized,
        "coefficients": np.asarray(coefficients),
        "fixed_effect_count": int(coefficient_count),
        "random_effect_sd": float(np.exp(log_mouse_sd)),
        "random_effect_mean": np.asarray(means)[:, :dimension],
        "random_effect_variance": np.diagonal(
            covariance[:, :dimension, :dimension], axis1=1, axis2=2
        ),
        "random_effect_covariance": covariance[:, :dimension, :dimension],
        "latent_mean": np.asarray(means),
        "latent_covariance": covariance,
        "observation_noise": bool(observation_noise),
        "observation_noise_sd": (
            float(np.exp(log_noise_sd)) if observation_noise else 0.0
        ),
        "concentration": (
            float(np.exp(log_kappa))
            if family == "dirichlet_multinomial"
            else np.inf
        ),
        "importance_ess": (
            float(np.median(cluster_ess)) if cluster_ess is not None else np.nan
        ),
        "minimum_importance_ess": (
            float(np.min(cluster_ess)) if cluster_ess is not None else np.nan
        ),
        "importance_samples": int(samples) if cluster_ess is not None else 0,
        "converged": bool(success or np.linalg.norm(gradient, ord=np.inf) <= 1e-3),
        "iterations": iterations,
        "gradient_norm": float(np.linalg.norm(gradient, ord=np.inf)),
        "message": message,
    }


def variational_warm_start(fit, n_design_columns):
    """Convert a Laplace or diagonal-VI fit to full-covariance VI."""
    coefficients = np.asarray(fit["coefficients"], dtype=float)
    expanded = np.zeros((int(n_design_columns), coefficients.shape[1]))
    expanded[: len(coefficients)] = coefficients
    latent_mean = np.asarray(
        fit.get("latent_mean", fit["random_effect_mean"]), dtype=float
    )
    if "latent_covariance" in fit:
        latent_covariance = np.asarray(fit["latent_covariance"], dtype=float)
    else:
        variance = np.asarray(fit["random_effect_variance"], dtype=float)
        latent_covariance = np.asarray([np.diag(row) for row in variance])
    values = np.r_[
        expanded.ravel(),
        latent_mean.ravel(),
        _pack_cholesky(latent_covariance),
        np.log(float(fit["random_effect_sd"])),
    ]
    if fit.get("observation_noise", False):
        values = np.r_[values, np.log(float(fit["observation_noise_sd"]))]
    if fit["family"] == "dirichlet_multinomial":
        values = np.r_[values, np.log(float(fit["concentration"]))]
    return values


def warm_start(fit, n_design_columns):
    """Expand a full-covariance VI fit into a larger fixed-effect design."""
    coefficients = np.asarray(fit["coefficients"], dtype=float)
    expanded = np.zeros((int(n_design_columns), coefficients.shape[1]))
    expanded[: len(coefficients)] = coefficients
    coefficient_count = coefficients.size
    return np.r_[expanded.ravel(), np.asarray(fit["parameters"])[coefficient_count:]]


def fixed_effect_warm_start(fit, n_fixed_effects):
    """Expand a tensor-design fit by appending zero fixed effects."""
    old_count = int(fit["fixed_effect_count"])
    if int(n_fixed_effects) < old_count:
        raise ValueError("warm-start fixed-effect count cannot decrease")
    expanded = np.zeros(int(n_fixed_effects), dtype=float)
    expanded[:old_count] = np.asarray(fit["coefficients"]).ravel()
    return np.r_[expanded, np.asarray(fit["parameters"])[old_count:]]


def evaluate_objectives(data, fit, *, alpha=0.5, samples=2048, seed=1):
    """Evaluate exact-likelihood ELBO and Renyi bounds at a full posterior."""
    initial = np.asarray(fit["parameters"])
    common = dict(
        data=data,
        family=fit["family"],
        objective="monte_carlo",
        observation_noise=fit.get("observation_noise", False),
        samples=samples,
        seed=seed,
        initial=initial,
        max_iter=0,
    )
    elbo = fit_variational(alpha=1.0, **common)["objective"]
    renyi = fit_variational(alpha=alpha, **common)["objective"]
    return {"elbo_mc": float(elbo), f"renyi_{alpha:g}": float(renyi)}
