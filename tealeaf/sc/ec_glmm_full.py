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


def _tilted_logsumexp_bound(
    jnp, means, covariances, iterations=25, differentiate_local=True
):
    """Upper-bound expected log-sum-exp for correlated Gaussian logits."""
    diagonal = jnp.diagonal(covariances, axis1=1, axis2=2)
    local = __import__("jax").nn.softmax(means + 0.5 * diagonal, axis=1)

    def update(_, value):
        covariance_times_local = jnp.einsum("nij,nj->ni", covariances, value)
        return __import__("jax").nn.softmax(
            means + 0.5 * diagonal - covariance_times_local, axis=1
        )

    local = __import__("jax").lax.fori_loop(0, int(iterations), update, local)
    if not differentiate_local:
        local = __import__("jax").lax.stop_gradient(local)
    covariance_times_local = jnp.einsum("nij,nj->ni", covariances, local)
    quadratic = 0.5 * jnp.sum(local * covariance_times_local, axis=1)
    adjusted = means + 0.5 * diagonal - covariance_times_local
    return quadratic + __import__("jax.scipy", fromlist=["special"]).special.logsumexp(
        adjusted, axis=1
    )


def _bouchard_lambda(xi):
    """Jaakkola--Jordan quadratic coefficient used by Bouchard's bound."""
    xi = np.asarray(xi, dtype=float)
    result = np.full_like(xi, 0.125)
    regular = np.abs(xi) >= 1e-7
    result[regular] = np.tanh(xi[regular] / 2.0) / (4.0 * xi[regular])
    return result


def _bouchard_parameters(means, variances, alpha=None, iterations=25):
    """Optimize Bouchard local parameters for Gaussian softmax logits.

    ``means`` and ``variances`` are observations by supported categories.
    The returned ``alpha`` is one scalar per observation and ``xi`` has the
    same shape as ``means``.
    """
    means = np.asarray(means, dtype=float)
    variances = np.asarray(variances, dtype=float)
    if means.ndim != 2 or variances.shape != means.shape:
        raise ValueError("Bouchard moments must be aligned matrices")
    categories = means.shape[1]
    if categories < 2:
        raise ValueError("Bouchard parameters require at least two categories")
    if alpha is None:
        alpha = scipy.special.logsumexp(means, axis=1)
    else:
        alpha = np.asarray(alpha, dtype=float)
    for _ in range(int(iterations)):
        xi = np.sqrt(
            np.maximum(
                variances + (means - alpha[:, None]) ** 2,
                1e-16,
            )
        )
        weights = 2.0 * _bouchard_lambda(xi)
        alpha = (
            categories / 2.0
            - 1.0
            + np.sum(weights * means, axis=1)
        ) / np.sum(weights, axis=1)
    xi = np.sqrt(
        np.maximum(variances + (means - alpha[:, None]) ** 2, 1e-16)
    )
    return alpha, xi


def _bouchard_expected_logsumexp(means, variances, alpha, xi):
    """Evaluate Bouchard's upper bound on expected log-sum-exp."""
    centered_mean = np.asarray(means, dtype=float) - np.asarray(alpha)[:, None]
    xi = np.asarray(xi, dtype=float)
    return np.asarray(alpha, dtype=float) + np.sum(
        _bouchard_lambda(xi)
        * (np.asarray(variances, dtype=float) + centered_mean**2 - xi**2)
        + 0.5 * (centered_mean - xi)
        + np.logaddexp(0.0, xi),
        axis=1,
    )


def _fixed_effect_maps(data):
    """Return observation-by-logit maps from coefficients to free logits."""
    dimension = data.n_isoforms - 1
    if data.fixed_effect_tensor is not None:
        return np.asarray(data.fixed_effect_tensor, dtype=float)
    design = np.asarray(data.design, dtype=float)
    result = np.zeros(
        (len(design), dimension, design.shape[1] * dimension), dtype=float
    )
    for column in range(design.shape[1]):
        start = column * dimension
        result[:, :, start : start + dimension] = (
            design[:, column, None, None] * np.eye(dimension)[None, :, :]
        )
    return result


def _unpack_full_initial(data, observation_noise, initial):
    """Initialize or unpack parameters shared by CAVI and gradient VI."""
    cluster_index, selection, _ = _latent_layout(data, observation_noise)
    n_clusters = int(cluster_index.max()) + 1
    latent_dimension = selection.shape[2]
    fixed_maps = _fixed_effect_maps(data)
    coefficient_count = fixed_maps.shape[2]
    rows, columns = np.tril_indices(latent_dimension)
    if initial is None:
        coefficients = np.zeros(coefficient_count, dtype=float)
        means = np.zeros((n_clusters, latent_dimension), dtype=float)
        covariances = np.tile(
            0.3**2 * np.eye(latent_dimension)[None, :, :],
            (n_clusters, 1, 1),
        )
        mouse_sd = 0.5
        noise_sd = 0.3 if observation_noise else 0.0
    else:
        initial = np.asarray(initial, dtype=float)
        mean_count = n_clusters * latent_dimension
        triangle_count = n_clusters * len(rows)
        expected = coefficient_count + mean_count + triangle_count + 1
        expected += int(observation_noise)
        if len(initial) != expected:
            raise ValueError("initial CAVI parameter vector has the wrong size")
        cursor = 0
        coefficients = initial[cursor : cursor + coefficient_count].copy()
        cursor += coefficient_count
        means = initial[cursor : cursor + mean_count].reshape(
            n_clusters, latent_dimension
        ).copy()
        cursor += mean_count
        packed = initial[cursor : cursor + triangle_count].reshape(
            n_clusters, len(rows)
        ).copy()
        cursor += triangle_count
        packed[:, rows == columns] = np.exp(packed[:, rows == columns])
        cholesky = np.zeros(
            (n_clusters, latent_dimension, latent_dimension), dtype=float
        )
        cholesky[:, rows, columns] = packed
        covariances = cholesky @ np.swapaxes(cholesky, 1, 2)
        mouse_sd = float(np.exp(initial[cursor]))
        cursor += 1
        noise_sd = float(np.exp(initial[cursor])) if observation_noise else 0.0
    return (
        cluster_index,
        selection,
        fixed_maps,
        coefficients,
        means,
        covariances,
        mouse_sd,
        noise_sd,
    )


def _pack_full_parameters(
    coefficients, means, covariances, mouse_sd, noise_sd, observation_noise
):
    values = np.r_[
        np.asarray(coefficients, dtype=float).ravel(),
        np.asarray(means, dtype=float).ravel(),
        _pack_cholesky(np.asarray(covariances, dtype=float)),
        np.log(float(mouse_sd)),
    ]
    if observation_noise:
        values = np.r_[values, np.log(float(noise_sd))]
    return values


def fit_bouchard_cavi(
    data,
    *,
    observation_noise=False,
    initial=None,
    max_iter=300,
    local_steps=25,
    tolerance=1e-5,
):
    """Fit the multinomial EC GLMM by Bouchard coordinate-ascent VI.

    EC-to-isoform responsibilities lower-bound each EC numerator. Bouchard's
    quadratic softmax bound lower-bounds the negative primer normalizer. The
    resulting Gaussian factors, fixed effects, and variance components have
    closed-form coordinate updates. Variance components are point estimates,
    so this is coordinate-ascent variational EM rather than fully Bayesian
    CAVI for every model parameter.
    """
    if int(max_iter) < 1:
        raise ValueError("CAVI max_iter must be positive")
    (
        cluster_index,
        selection,
        fixed_maps,
        coefficients,
        means,
        covariances,
        mouse_sd,
        noise_sd,
    ) = _unpack_full_initial(data, observation_noise, initial)
    counts = tuple(np.asarray(value, dtype=float) for value in data.counts)
    mappings = tuple(np.asarray(value, dtype=float) for value in data.compatibility)
    n_observations = len(cluster_index)
    n_clusters, latent_dimension = means.shape
    dimension = data.n_isoforms - 1
    latent_blocks = latent_dimension // dimension
    alpha = [None] * len(counts)
    lower_sd = float(np.exp(-8.0))
    upper_sd = float(np.exp(3.0))
    objective_history = []
    maximum_change = np.inf

    def moments():
        fixed = np.einsum("ndp,p->nd", fixed_maps, coefficients)
        random = np.einsum(
            "ndi,ni->nd", selection, means[cluster_index]
        )
        free_means = fixed + random
        free_variances = np.einsum(
            "nai,nij,naj->na",
            selection,
            covariances[cluster_index],
            selection,
        )
        full_means = np.column_stack((free_means, np.zeros(n_observations)))
        full_variances = np.column_stack(
            (free_variances, np.zeros(n_observations))
        )
        return fixed, full_means, full_variances

    def local_update(full_means, full_variances):
        weights = np.zeros((n_observations, dimension), dtype=float)
        natural = np.zeros((n_observations, dimension), dtype=float)
        expected_likelihood = np.zeros(n_observations, dtype=float)
        for primer, (observed, mapping) in enumerate(zip(counts, mappings)):
            totals = observed.sum(axis=1)
            expected_likelihood += scipy.special.gammaln(totals + 1.0)
            expected_likelihood -= np.sum(
                scipy.special.gammaln(observed + 1.0), axis=1
            )
            assigned = np.zeros((n_observations, data.n_isoforms), dtype=float)
            allocation_constant = np.zeros(n_observations, dtype=float)
            for ec, compatibility in enumerate(mapping):
                support = compatibility > 0
                if not np.any(support):
                    continue
                logits = (
                    full_means[:, support]
                    + np.log(compatibility[support])[None, :]
                )
                responsibilities = scipy.special.softmax(logits, axis=1)
                assigned[:, support] += observed[:, ec, None] * responsibilities
                allocation_constant += observed[:, ec] * np.sum(
                    responsibilities
                    * (
                        np.log(compatibility[support])[None, :]
                        - np.log(np.maximum(responsibilities, 1e-300))
                    ),
                    axis=1,
                )
            expected_likelihood += allocation_constant
            expected_likelihood += np.sum(assigned * full_means, axis=1)
            natural += assigned[:, :dimension]

            lengths = mapping.sum(axis=0)
            support = lengths > 0
            supported_count = int(np.sum(support))
            shifted_means = (
                full_means[:, support] + np.log(lengths[support])[None, :]
            )
            shifted_variances = full_variances[:, support]
            if supported_count == 0:
                if np.any(totals > 0):
                    raise ValueError("positive primer totals lack isoform support")
                continue
            if supported_count == 1:
                expected_likelihood -= totals * shifted_means[:, 0]
                only = int(np.flatnonzero(support)[0])
                if only < dimension:
                    natural[:, only] -= totals
                continue
            local_alpha, xi = _bouchard_parameters(
                shifted_means,
                shifted_variances,
                alpha=alpha[primer],
                iterations=local_steps,
            )
            alpha[primer] = local_alpha
            lam = _bouchard_lambda(xi)
            expected_likelihood -= totals * _bouchard_expected_logsumexp(
                shifted_means, shifted_variances, local_alpha, xi
            )
            supported_indices = np.flatnonzero(support)
            for local_column, isoform in enumerate(supported_indices):
                if isoform >= dimension:
                    continue
                quadratic = 2.0 * totals * lam[:, local_column]
                weights[:, isoform] += quadratic
                natural[:, isoform] -= totals * (
                    2.0
                    * lam[:, local_column]
                    * (np.log(lengths[isoform]) - local_alpha)
                    + 0.5
                )
        return weights, natural, float(np.sum(expected_likelihood))

    def kl_divergence():
        prior_sds = np.full(latent_dimension, mouse_sd, dtype=float)
        if observation_noise:
            prior_sds[dimension:] = noise_sd
        prior_variances = prior_sds**2
        trace = np.sum(
            np.diagonal(covariances, axis1=1, axis2=2)
            / prior_variances[None, :]
        )
        quadratic = np.sum(means * means / prior_variances[None, :])
        logdet_prior = n_clusters * np.sum(np.log(prior_variances))
        signs, logdet_covariance = np.linalg.slogdet(covariances)
        if not np.all(signs > 0):
            return np.inf
        return 0.5 * (
            trace
            + quadratic
            - n_clusters * latent_dimension
            + logdet_prior
            - np.sum(logdet_covariance)
        )

    converged = False
    iterations = 0
    for iteration in range(int(max_iter)):
        old_parameters = _pack_full_parameters(
            coefficients,
            means,
            covariances,
            mouse_sd,
            noise_sd,
            observation_noise,
        )
        fixed, full_means, full_variances = moments()
        likelihood_weights, natural, _ = local_update(
            full_means, full_variances
        )

        prior_precision = np.full(
            latent_dimension, 1.0 / mouse_sd**2, dtype=float
        )
        if observation_noise:
            prior_precision[dimension:] = 1.0 / noise_sd**2
        for cluster in range(n_clusters):
            rows = np.flatnonzero(cluster_index == cluster)
            precision = np.diag(prior_precision)
            right = np.zeros(latent_dimension, dtype=float)
            for row in rows:
                weighted_selection = (
                    likelihood_weights[row, :, None] * selection[row]
                )
                precision += selection[row].T @ weighted_selection
                right += selection[row].T @ (
                    natural[row]
                    - likelihood_weights[row] * fixed[row]
                )
            covariance = np.linalg.inv(precision)
            covariance = 0.5 * (covariance + covariance.T)
            covariances[cluster] = covariance
            means[cluster] = np.linalg.solve(precision, right)

        fixed_precision = np.zeros(
            (len(coefficients), len(coefficients)), dtype=float
        )
        fixed_right = np.zeros(len(coefficients), dtype=float)
        latent_means = np.einsum(
            "ndi,ni->nd", selection, means[cluster_index]
        )
        for row in range(n_observations):
            weighted_map = likelihood_weights[row, :, None] * fixed_maps[row]
            fixed_precision += fixed_maps[row].T @ weighted_map
            fixed_right += fixed_maps[row].T @ (
                natural[row] - likelihood_weights[row] * latent_means[row]
            )
        coefficients = np.linalg.lstsq(
            fixed_precision, fixed_right, rcond=1e-10
        )[0]

        mouse_second_moment = (
            means[:, :dimension] ** 2
            + np.diagonal(
                covariances[:, :dimension, :dimension], axis1=1, axis2=2
            )
        )
        mouse_sd = float(
            np.clip(
                np.sqrt(np.mean(mouse_second_moment)), lower_sd, upper_sd
            )
        )
        if observation_noise:
            noise_second_moment = (
                means[:, dimension:] ** 2
                + np.diagonal(
                    covariances[:, dimension:, dimension:], axis1=1, axis2=2
                )
            )
            active = np.zeros((n_clusters, latent_blocks - 1), dtype=bool)
            cluster_sizes = np.bincount(cluster_index, minlength=n_clusters)
            for cluster, size in enumerate(cluster_sizes):
                active[cluster, :size] = True
            active = np.repeat(active, dimension, axis=1)
            noise_sd = float(
                np.clip(
                    np.sqrt(np.mean(noise_second_moment[active])),
                    lower_sd,
                    upper_sd,
                )
            )

        _, updated_means, updated_variances = moments()
        _, _, expected_likelihood = local_update(
            updated_means, updated_variances
        )
        objective = expected_likelihood - kl_divergence()
        objective_history.append(float(objective))
        new_parameters = _pack_full_parameters(
            coefficients,
            means,
            covariances,
            mouse_sd,
            noise_sd,
            observation_noise,
        )
        maximum_change = float(np.max(np.abs(new_parameters - old_parameters)))
        iterations = iteration + 1
        if (
            iteration > 0
            and maximum_change <= np.sqrt(float(tolerance))
            and abs(objective_history[-1] - objective_history[-2])
            <= float(tolerance) * max(1.0, abs(objective_history[-2]))
        ):
            converged = True
            break

    parameters = _pack_full_parameters(
        coefficients,
        means,
        covariances,
        mouse_sd,
        noise_sd,
        observation_noise,
    )
    reported_coefficients = coefficients
    if data.fixed_effect_tensor is None:
        reported_coefficients = coefficients.reshape(
            data.design.shape[1], dimension
        )
    return {
        "method": "bouchard_cavi_full",
        "family": "multinomial",
        "alpha": 1.0,
        "objective": float(objective_history[-1]),
        "parameters": parameters,
        "coefficients": reported_coefficients,
        "fixed_effect_count": int(len(coefficients)),
        "random_effect_sd": mouse_sd,
        "random_effect_mean": means[:, :dimension],
        "random_effect_variance": np.diagonal(
            covariances[:, :dimension, :dimension], axis1=1, axis2=2
        ),
        "random_effect_covariance": covariances[:, :dimension, :dimension],
        "latent_mean": means,
        "latent_covariance": covariances,
        "observation_noise": bool(observation_noise),
        "observation_noise_sd": noise_sd if observation_noise else 0.0,
        "concentration": np.inf,
        "importance_ess": np.nan,
        "minimum_importance_ess": np.nan,
        "importance_samples": 0,
        "converged": bool(converged),
        "iterations": iterations,
        "gradient_norm": maximum_change,
        "coordinate_change": maximum_change,
        "objective_history": np.asarray(objective_history),
        "message": (
            "coordinate tolerance reached"
            if converged
            else "coordinate iteration limit reached"
        ),
    }


def fit_cavi_then_tilted(
    data,
    *,
    observation_noise=False,
    initial=None,
    cavi_max_iter=25,
    max_iter=300,
):
    """Initialize the tilted full-covariance fit with Bouchard CAVI."""
    cavi = fit_bouchard_cavi(
        data,
        observation_noise=observation_noise,
        initial=initial,
        max_iter=cavi_max_iter,
    )
    refined = fit_variational(
        data,
        family="multinomial",
        objective="tilted",
        observation_noise=observation_noise,
        initial=cavi["parameters"],
        max_iter=max_iter,
    )
    refined["method"] = "bouchard_cavi_then_tilted_full"
    refined["cavi_objective"] = cavi["objective"]
    refined["cavi_iterations"] = cavi["iterations"]
    refined["cavi_converged"] = cavi["converged"]
    return refined


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
    tilted_local_steps=25,
    differentiate_tilted_local=True,
    optimizer_maxcor=10,
    optimizer_ftol=1e-8,
    standardize_latent=False,
    latent_standardization_power=1.0,
    compute_dtype="float64",
    objective_cache=None,
    objective_cache_key=None,
):
    """Fit a cluster-factorized, full-covariance Gaussian posterior.

    ``objective`` is either ``"tilted"`` for the analytic multinomial ELBO
    lower bound or ``"monte_carlo"`` for an exact-likelihood ELBO/Renyi
    objective.  With observation noise, each cluster block jointly represents
    its shared random intercept and all observation-level residual effects.
    ``compute_dtype`` controls JAX objective evaluation; SciPy retains its
    L-BFGS state in double precision.
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
    dtype_by_name = {
        "float64": jnp.float64,
        "float32": jnp.float32,
        "bfloat16": jnp.bfloat16,
    }
    if compute_dtype not in dtype_by_name:
        raise ValueError("compute_dtype must be 'float64', 'float32', or 'bfloat16'")
    compute_type = dtype_by_name[compute_dtype]
    counts = tuple(jnp.asarray(value, dtype=compute_type) for value in counts)
    mappings = tuple(jnp.asarray(value, dtype=compute_type) for value in mappings)
    design = jnp.asarray(design, dtype=compute_type)
    cluster_index_np, selection_np, latent_blocks = _latent_layout(
        data, observation_noise
    )
    np.testing.assert_array_equal(cluster_index_np, np.asarray(cluster_index))
    selection = jnp.asarray(selection_np, dtype=compute_type)
    dimension = data.n_isoforms - 1
    latent_dimension = selection.shape[2]
    fixed_effect_tensor = (
        None
        if data.fixed_effect_tensor is None
        else jnp.asarray(data.fixed_effect_tensor, dtype=compute_type)
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
    membership = jax.nn.one_hot(cluster_index, n_clusters, dtype=compute_type)

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

    def prior_log_sds(log_mouse_sd, log_noise_sd):
        values = jnp.full((latent_blocks, dimension), log_mouse_sd)
        if observation_noise:
            values = values.at[1:].set(log_noise_sd)
        return values.ravel()

    if standardize_latent:
        parameters = np.asarray(parameters, dtype=float).copy()
        mean_start = coefficient_count
        cholesky_start = mean_start + mean_count
        scale_start = cholesky_start + cholesky_count
        log_mouse_sd = parameters[scale_start]
        log_noise_sd = (
            parameters[scale_start + 1] if observation_noise else log_mouse_sd
        )
        log_prior_sds_np = np.full(latent_dimension, log_mouse_sd, dtype=float)
        if observation_noise:
            log_prior_sds_np[dimension:] = log_noise_sd
        log_prior_sds_np *= float(latent_standardization_power)
        means = parameters[mean_start:cholesky_start].reshape(
            n_clusters, latent_dimension
        )
        means /= np.exp(log_prior_sds_np)[None, :]
        packed = parameters[cholesky_start:scale_start].reshape(
            n_clusters, len(triangle_rows_np)
        )
        diagonal = triangle_rows_np == triangle_columns_np
        packed[:, diagonal] -= log_prior_sds_np[triangle_rows_np[diagonal]]
        packed[:, ~diagonal] /= np.exp(
            log_prior_sds_np[triangle_rows_np[~diagonal]]
        )[None, :]

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
        log_mouse_sd = value[cursor]
        cursor += 1
        log_noise_sd = value[cursor] if observation_noise else value[0] * 0.0
        cursor += int(observation_noise)
        log_kappa = (
            value[cursor]
            if family == "dirichlet_multinomial"
            else value[0] * 0.0
        )
        cholesky = _cholesky_from_packed(
            jnp, packed, triangle_rows, triangle_columns, triangle_diagonal
        )
        if standardize_latent:
            scales = jnp.exp(
                float(latent_standardization_power)
                * prior_log_sds(log_mouse_sd, log_noise_sd)
            )
            means = means * scales[None, :]
            cholesky = cholesky * scales[None, :, None]
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

    def tilted_bound(value, active_counts=counts):
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
        for observed, mapping in zip(active_counts, mappings):
            if observed.shape[1] == 0:
                continue
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
                jnp,
                shifted_means,
                covariances,
                iterations=tilted_local_steps,
                differentiate_local=differentiate_tilted_local,
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
    standard_noise = jnp.asarray(standard_noise, dtype=compute_type)

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
        cache_key = None
        if objective == "tilted" and objective_cache_key is not None:
            cache_key = (
                objective_cache_key,
                "value_and_gradient",
                str(compute_dtype),
                int(tilted_local_steps),
                bool(differentiate_tilted_local),
                bool(standardize_latent),
                float(latent_standardization_power),
            )
        if objective_cache is not None and cache_key in objective_cache:
            value_and_gradient = objective_cache[cache_key]
        elif objective == "tilted" and cache_key is not None:
            value_and_gradient = jax.jit(
                jax.value_and_grad(
                    lambda value, active_counts: -tilted_bound(
                        value, active_counts
                    ),
                    argnums=0,
                )
            )
            objective_cache[cache_key] = value_and_gradient
        else:
            value_and_gradient = jax.jit(
                jax.value_and_grad(lambda value: -bound(value))
            )

        def scipy_objective(value):
            packed = jnp.asarray(value, dtype=compute_type)
            if cache_key is None:
                result, gradient = value_and_gradient(packed)
            else:
                result, gradient = value_and_gradient(packed, counts)
            return float(result), np.asarray(gradient, dtype=float)
    else:
        value_only = jax.jit(lambda value: -bound(value))

        def scipy_objective(value):
            result = value_only(jnp.asarray(value, dtype=compute_type))
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
        starting_value, starting_gradient = scipy_objective(parameters)
        result = scipy.optimize.minimize(
            scipy_objective,
            parameters,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={
                "maxiter": int(max_iter),
                "ftol": float(optimizer_ftol),
                "gtol": 1e-4,
                "maxcor": int(optimizer_maxcor),
            },
        )
        optimized = np.asarray(result.x)
        success = bool(result.success)
        iterations = int(result.nit)
        message = str(result.message)
        tolerance = 1e-10 * max(1.0, abs(starting_value))
        if float(result.fun) > starting_value + tolerance:
            optimized = np.asarray(parameters)
            success = bool(
                np.linalg.norm(starting_gradient, ord=np.inf) <= 1e-3
            )
            message = "rejected optimizer result worse than warm start"
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
    ) = unpack(jnp.asarray(optimized, dtype=compute_type))
    covariance = np.asarray(cholesky @ jnp.swapaxes(cholesky, 1, 2))
    final, gradient = scipy_objective(optimized)
    external_parameters = _pack_full_parameters(
        np.asarray(coefficients).ravel(),
        np.asarray(means),
        covariance,
        float(np.exp(log_mouse_sd)),
        float(np.exp(log_noise_sd)),
        observation_noise,
    )
    if family == "dirichlet_multinomial":
        external_parameters = np.r_[external_parameters, float(log_kappa)]
    if objective == "monte_carlo":
        final_log_weights = np.asarray(
            log_weights(jnp.asarray(optimized, dtype=compute_type))
        )
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
        "parameters": external_parameters,
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
        "tilted_local_steps": int(tilted_local_steps),
        "differentiate_tilted_local": bool(differentiate_tilted_local),
        "optimizer_maxcor": int(optimizer_maxcor),
        "optimizer_ftol": float(optimizer_ftol),
        "standardize_latent": bool(standardize_latent),
        "latent_standardization_power": float(latent_standardization_power),
        "compute_dtype": str(compute_dtype),
    }


def fit_tilted_variational_robust(
    data,
    *,
    observation_noise=False,
    initial=None,
    max_iter=300,
    fallback_iter=50,
    compute_dtype="float64",
    objective_cache=None,
    objective_cache_key=None,
):
    """Fit tilted VI with envelope gradients and prior-scale preconditioning.

    The latent means and Cholesky rows are optimized after scaling by the
    corresponding prior standard deviation to the 0.75 power.  This removes
    most of the variance-component ridge without changing the variational
    family or returned parameter representation.  A relaxed continuation is
    used only when the strict fit reaches its iteration limit.
    """
    common = {
        "data": data,
        "family": "multinomial",
        "objective": "tilted",
        "observation_noise": observation_noise,
        "tilted_local_steps": 8,
        "differentiate_tilted_local": False,
        "standardize_latent": True,
        "latent_standardization_power": 0.75,
        "compute_dtype": compute_dtype,
        "objective_cache": objective_cache,
        "objective_cache_key": objective_cache_key,
    }
    strict = fit_variational(
        initial=initial,
        max_iter=max_iter,
        **common,
    )
    strict_iterations = int(strict["iterations"])
    if strict["converged"] or int(fallback_iter) <= 0:
        strict["strict_iterations"] = strict_iterations
        strict["fallback_iterations"] = 0
        strict["used_relaxed_fallback"] = False
        return strict
    refined = fit_variational(
        initial=strict["parameters"],
        max_iter=fallback_iter,
        optimizer_ftol=1e-7,
        **common,
    )
    fallback_iterations = int(refined["iterations"])
    refined["iterations"] = strict_iterations + fallback_iterations
    refined["strict_iterations"] = strict_iterations
    refined["fallback_iterations"] = fallback_iterations
    refined["used_relaxed_fallback"] = True
    refined["message"] = (
        f"strict: {strict['message']}; relaxed continuation: {refined['message']}"
    )
    return refined


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
