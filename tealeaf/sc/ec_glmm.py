"""Mixed models for paired-primer equivalence-class counts.

The latent linear predictor is defined on isoforms, while observations remain
equivalence-class (EC) counts.  This keeps ambiguity in the likelihood instead
of first converting EC counts to estimated isoform or splice-path counts.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import scipy.optimize
import scipy.special
import scipy.sparse as sp


def _jax():
    try:
        import jax
        import jax.numpy as jnp
        import jax.scipy as jsp
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("EC GLMM fitting requires JAX") from exc
    jax.config.update("jax_enable_x64", True)
    return jax, jnp, jsp


@dataclass(frozen=True)
class ECGLMMData:
    """Aligned EC-count GLMM inputs.

    ``counts[p]`` is observations by ECs and ``compatibility[p]`` is ECs by
    isoforms for primer ``p``.  Compatibility entries may be binary or fixed
    Salmon conditional weights.
    """

    counts: tuple[np.ndarray, ...]
    compatibility: tuple[np.ndarray, ...]
    design: np.ndarray
    clusters: np.ndarray
    fixed_effect_tensor: np.ndarray | None = None

    def __post_init__(self):
        counts = tuple(np.asarray(value, dtype=float) for value in self.counts)
        compatibility = tuple(
            np.asarray(
                value.toarray() if sp.issparse(value) else value,
                dtype=float,
            )
            for value in self.compatibility
        )
        design = np.asarray(self.design, dtype=float)
        clusters = np.asarray(self.clusters)
        fixed_effect_tensor = (
            None
            if self.fixed_effect_tensor is None
            else np.asarray(self.fixed_effect_tensor, dtype=float)
        )
        if not counts or len(counts) != len(compatibility):
            raise ValueError("counts and compatibility need one entry per primer")
        n_observations = counts[0].shape[0]
        n_isoforms = compatibility[0].shape[1]
        for observed, mapping in zip(counts, compatibility):
            if observed.ndim != 2 or mapping.ndim != 2:
                raise ValueError("counts and compatibility must be matrices")
            if observed.shape != (n_observations, mapping.shape[0]):
                raise ValueError("an EC count matrix does not match its mapping")
            if mapping.shape[1] != n_isoforms:
                raise ValueError("primer mappings have different isoform counts")
            if np.any(observed < 0) or np.any(mapping < 0):
                raise ValueError("counts and compatibility must be nonnegative")
            unsupported = (observed.sum(axis=0) > 0) & (mapping.sum(axis=1) <= 0)
            if unsupported.any():
                raise ValueError("positive EC counts lack compatible isoforms")
        combined_support = sum(mapping.sum(axis=0) for mapping in compatibility)
        if np.any(combined_support <= 0):
            raise ValueError("every isoform needs support in at least one primer mapping")
        if design.ndim != 2 or design.shape[0] != n_observations:
            raise ValueError("the fixed-effect design does not align with counts")
        if len(clusters) != n_observations:
            raise ValueError("clusters do not align with counts")
        if np.linalg.matrix_rank(design) != design.shape[1]:
            raise ValueError("the fixed-effect design is rank deficient")
        if fixed_effect_tensor is not None:
            expected = (n_observations, n_isoforms - 1)
            if (
                fixed_effect_tensor.ndim != 3
                or fixed_effect_tensor.shape[:2] != expected
                or fixed_effect_tensor.shape[2] < 1
            ):
                raise ValueError(
                    "fixed_effect_tensor must be observations by free logits "
                    "by coefficients"
                )
        if n_isoforms < 2:
            raise ValueError("an EC GLMM needs at least two isoforms")
        object.__setattr__(self, "counts", counts)
        object.__setattr__(self, "compatibility", compatibility)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "clusters", clusters)
        object.__setattr__(self, "fixed_effect_tensor", fixed_effect_tensor)

    @property
    def n_isoforms(self):
        return self.compatibility[0].shape[1]


def _prepared_arrays(data):
    jax, jnp, jsp = _jax()
    unique_clusters, cluster_index = np.unique(data.clusters, return_inverse=True)
    counts = tuple(jnp.asarray(value) for value in data.counts)
    mappings = tuple(jnp.asarray(value) for value in data.compatibility)
    return (
        jax,
        jnp,
        jsp,
        counts,
        mappings,
        jnp.asarray(data.design),
        jnp.asarray(cluster_index, dtype=jnp.int32),
        len(unique_clusters),
    )


def _log_likelihood_rows(
    jnp, jsp, counts, mappings, logits, family, log_kappa
):
    """Conditional EC log likelihood by observation."""
    abundance = jnp.exp(logits - jnp.max(logits, axis=1, keepdims=True))
    result = jnp.zeros(len(logits), dtype=logits.dtype)
    for observed, mapping in zip(counts, mappings):
        mass = abundance @ mapping.T
        probability = mass / jnp.sum(mass, axis=1, keepdims=True)
        probability = jnp.maximum(probability, 1e-300)
        totals = jnp.sum(observed, axis=1)
        constant = jsp.special.gammaln(totals + 1.0) - jnp.sum(
            jsp.special.gammaln(observed + 1.0), axis=1
        )
        if family == "multinomial":
            result += constant + jnp.sum(
                observed * jnp.log(probability), axis=1
            )
        elif family == "dirichlet_multinomial":
            kappa = jnp.exp(log_kappa)
            alpha = jnp.maximum(kappa * probability, 1e-12)
            result += (
                constant
                + jsp.special.gammaln(kappa)
                - jsp.special.gammaln(kappa + totals)
                + jnp.sum(
                    jsp.special.gammaln(alpha + observed)
                    - jsp.special.gammaln(alpha),
                    axis=1,
                )
            )
        else:
            raise ValueError(f"unknown EC likelihood family: {family}")
    return result


def _log_likelihood(jnp, jsp, counts, mappings, logits, family, log_kappa):
    """Conditional EC log likelihood, including count-only constants."""
    return jnp.sum(
        _log_likelihood_rows(
            jnp, jsp, counts, mappings, logits, family, log_kappa
        )
    )


def _unpack_outer(parameters, n_design, dimension, family):
    coefficient_count = n_design * dimension
    coefficients = parameters[:coefficient_count].reshape(n_design, dimension)
    log_prior_sd = parameters[coefficient_count]
    log_kappa = (
        parameters[coefficient_count + 1]
        if family == "dirichlet_multinomial"
        else parameters[0] * 0.0
    )
    return coefficients, log_prior_sd, log_kappa


def _initial_outer(data, family):
    dimension = data.n_isoforms - 1
    values = np.zeros(data.design.shape[1] * dimension + 1, dtype=float)
    values[-1] = np.log(0.5)
    if family == "dirichlet_multinomial":
        values = np.r_[values, np.log(100.0)]
    return values


def fit_laplace(
    data,
    *,
    family="multinomial",
    observation_noise=False,
    initial=None,
    max_iter=200,
    mode_steps=30,
):
    """Fit an EC-count random-intercept GLMM by Laplace approximation.

    When ``observation_noise`` is true, the latent block for each cluster also
    contains an independent isotropic logistic-normal effect for every
    observation in that cluster.  This supplies multinomial overdispersion on
    the isoform-logit scale while retaining the EC observation model.
    """
    if observation_noise and family != "multinomial":
        raise ValueError("observation noise is currently supported for multinomial fits")
    (
        jax,
        jnp,
        _,
        counts,
        mappings,
        design,
        cluster_index,
        n_clusters,
    ) = _prepared_arrays(data)
    jsp = __import__("jax.scipy", fromlist=["special"])
    dimension = data.n_isoforms - 1
    if observation_noise:
        base = _initial_outer(data, family)
        initial = np.r_[base, np.log(0.3)] if initial is None else np.asarray(initial)
    else:
        initial = _initial_outer(data, family) if initial is None else np.asarray(initial)

    membership = jax.nn.one_hot(
        cluster_index, n_clusters, dtype=jnp.float64
    )
    cluster_index_np = np.asarray(cluster_index)
    slots = np.zeros(len(cluster_index_np), dtype=int)
    cluster_sizes = np.zeros(n_clusters, dtype=int)
    for row, cluster in enumerate(cluster_index_np):
        slots[row] = cluster_sizes[cluster]
        cluster_sizes[cluster] += 1
    max_cluster_size = int(cluster_sizes.max())
    latent_blocks = 1 + (max_cluster_size if observation_noise else 0)
    latent_dimension = latent_blocks * dimension
    selection = np.zeros((len(design), dimension, latent_dimension), dtype=float)
    selection[:, :, :dimension] = np.eye(dimension)[None, :, :]
    if observation_noise:
        for row, slot in enumerate(slots):
            start = (1 + slot) * dimension
            selection[row, :, start : start + dimension] = np.eye(dimension)
    selection = jnp.asarray(selection)

    def unpack_outer(outer):
        coefficient_count = design.shape[1] * dimension
        coefficients = outer[:coefficient_count].reshape(design.shape[1], dimension)
        log_prior_sd = outer[coefficient_count]
        cursor = coefficient_count + 1
        log_kappa = outer[0] * 0.0
        if family == "dirichlet_multinomial":
            log_kappa = outer[cursor]
            cursor += 1
        log_noise_sd = outer[cursor] if observation_noise else outer[0] * 0.0
        return coefficients, log_prior_sd, log_noise_sd, log_kappa

    def prior_log_sds(outer):
        _, log_prior_sd, log_noise_sd, _ = unpack_outer(outer)
        values = jnp.full((latent_blocks, dimension), log_prior_sd)
        if observation_noise:
            values = values.at[1:].set(log_noise_sd)
        return values.ravel()

    def observation_log_likelihood(free_logits, *observed_rows_and_log_kappa):
        observed_rows = observed_rows_and_log_kappa[:-1]
        log_kappa = observed_rows_and_log_kappa[-1]
        logits = jnp.concatenate((free_logits, jnp.zeros(1, dtype=free_logits.dtype)))
        abundance = jnp.exp(logits - jnp.max(logits))
        value = jnp.asarray(0.0, dtype=logits.dtype)
        for observed, mapping in zip(observed_rows, mappings):
            mass = mapping @ abundance
            probability = jnp.maximum(mass / jnp.sum(mass), 1e-300)
            total = jnp.sum(observed)
            constant = jsp.special.gammaln(total + 1.0) - jnp.sum(
                jsp.special.gammaln(observed + 1.0)
            )
            if family == "multinomial":
                value += constant + jnp.sum(observed * jnp.log(probability))
            else:
                kappa = jnp.exp(log_kappa)
                alpha = jnp.maximum(kappa * probability, 1e-12)
                value += (
                    constant
                    + jsp.special.gammaln(kappa)
                    - jsp.special.gammaln(kappa + total)
                    + jnp.sum(
                        jsp.special.gammaln(alpha + observed)
                        - jsp.special.gammaln(alpha)
                    )
                )
        return value

    observation_gradient = jax.vmap(
        jax.grad(observation_log_likelihood, argnums=0),
        in_axes=(0,) + (0,) * len(counts) + (None,),
    )
    observation_hessian = jax.vmap(
        jax.hessian(observation_log_likelihood, argnums=0),
        in_axes=(0,) + (0,) * len(counts) + (None,),
    )

    def mode_derivatives(modes, outer):
        coefficients, _, _, log_kappa = unpack_outer(outer)
        fixed_logits = design @ coefficients
        random_logits = jnp.einsum(
            "ndi,ni->nd", selection, modes[cluster_index]
        )
        free_logits = fixed_logits + random_logits
        score_rows = observation_gradient(
            free_logits, *counts, log_kappa
        )
        log_sds = prior_log_sds(outer)
        variances = jnp.exp(2.0 * log_sds)
        latent_score_rows = jnp.einsum("ndi,nd->ni", selection, score_rows)
        score = membership.T @ latent_score_rows - modes / variances
        curvature_rows = -observation_hessian(
            free_logits, *counts, log_kappa
        )
        latent_curvature_rows = jnp.einsum(
            "nai,nab,nbj->nij", selection, curvature_rows, selection
        )
        hessian = jnp.einsum(
            "im,ijk->mjk", membership, latent_curvature_rows
        ) + jnp.diag(1.0 / variances)[None, :, :]
        return fixed_logits, score, hessian, variances, log_kappa

    def negative_joint(modes, outer):
        fixed_logits, _, _, variances, log_kappa = mode_derivatives(modes, outer)
        random_logits = jnp.einsum(
            "ndi,ni->nd", selection, modes[cluster_index]
        )
        free_logits = fixed_logits + random_logits
        logits = jnp.concatenate(
            (free_logits, jnp.zeros((len(design), 1), dtype=free_logits.dtype)), axis=1
        )
        log_sds = prior_log_sds(outer)
        log_prior = -0.5 * jnp.sum(modes * modes / variances)
        log_prior -= n_clusters * (
            jnp.sum(log_sds) + 0.5 * latent_dimension * jnp.log(2.0 * jnp.pi)
        )
        return -(
            _log_likelihood(jnp, jsp, counts, mappings, logits, family, log_kappa)
            + log_prior
        )

    def posterior_mode(outer):
        current = jnp.zeros((n_clusters, latent_dimension), dtype=jnp.float64)

        def update(_, value):
            _, score, hessian, _, _ = mode_derivatives(value, outer)
            eigen_floor = jnp.maximum(
                1e-6 - jnp.linalg.eigvalsh(hessian)[:, 0], 0.0
            )
            hessian = hessian + eigen_floor[:, None, None] * jnp.eye(latent_dimension)
            step = jnp.linalg.solve(hessian, score[..., None])[..., 0]
            maximum = jnp.max(jnp.abs(step), axis=1, keepdims=True)
            scale = jnp.minimum(1.0, 2.0 / (maximum + 1e-12))
            return value + scale * step

        return jax.lax.fori_loop(0, int(mode_steps), update, current)

    def objective(outer):
        modes = posterior_mode(outer)
        _, _, hessian, _, _ = mode_derivatives(modes, outer)
        sign, log_determinant = jnp.linalg.slogdet(hessian)
        correction = 0.5 * jnp.sum(log_determinant)
        correction -= 0.5 * modes.size * jnp.log(2.0 * jnp.pi)
        return jnp.where(
            jnp.all(sign > 0), negative_joint(modes, outer) + correction, jnp.inf
        )

    value_and_gradient = jax.jit(jax.value_and_grad(objective))

    def scipy_objective(parameters):
        value, gradient = value_and_gradient(jnp.asarray(parameters))
        return float(value), np.asarray(gradient, dtype=float)

    coefficient_count = data.design.shape[1] * dimension
    bounds = [(-20.0, 20.0)] * coefficient_count + [(-8.0, 3.0)]
    if family == "dirichlet_multinomial":
        bounds.append((-6.0, np.log(1e7)))
    if observation_noise:
        bounds.append((-8.0, 3.0))
    result = scipy.optimize.minimize(
        scipy_objective,
        initial,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": int(max_iter), "ftol": 1e-9, "gtol": 1e-5},
    )
    outer = jnp.asarray(result.x)
    modes = posterior_mode(outer)
    _, mode_score, hessian, _, _ = mode_derivatives(modes, outer)
    covariance = np.asarray(jnp.linalg.inv(hessian))
    marginal_variance = np.diagonal(covariance, axis1=1, axis2=2)
    final_value, final_gradient = scipy_objective(result.x)
    mode_score = np.asarray(mode_score)
    coefficients, log_prior_sd, log_noise_sd, log_kappa = unpack_outer(outer)
    return {
        "method": "laplace",
        "family": family,
        "objective": final_value,
        "parameters": np.asarray(result.x),
        "coefficients": np.asarray(coefficients),
        "random_effect_sd": float(np.exp(log_prior_sd)),
        "random_effect_mean": np.asarray(modes)[:, :dimension],
        "random_effect_variance": marginal_variance[:, :dimension],
        "random_effect_covariance": covariance[:, :dimension, :dimension],
        "latent_mean": np.asarray(modes),
        "latent_covariance": covariance,
        "observation_noise_sd": (
            float(np.exp(log_noise_sd)) if observation_noise else 0.0
        ),
        "observation_noise": bool(observation_noise),
        "concentration": (
            float(np.exp(log_kappa)) if family == "dirichlet_multinomial" else np.inf
        ),
        "converged": bool(
            result.success
            or (
                np.isfinite(final_value)
                and np.linalg.norm(final_gradient, ord=np.inf) <= 1e-3
                and np.linalg.norm(mode_score, ord=np.inf) <= 1e-4
            )
        ),
        "iterations": int(result.nit),
        "gradient_norm": float(np.linalg.norm(final_gradient, ord=np.inf)),
        "mode_score_norm": float(np.linalg.norm(mode_score, ord=np.inf)),
        "message": str(result.message),
    }


def _tilted_logsumexp_bound(jnp, means, variances, iterations=25):
    """Knowles--Minka tilted upper bound for expected log-sum-exp."""
    local = jnp.exp(means + 0.5 * variances)
    local /= jnp.sum(local, axis=1, keepdims=True)

    def update(_, value):
        logits = means + 0.5 * (1.0 - 2.0 * value) * variances
        return jnp.exp(logits - jnp.max(logits, axis=1, keepdims=True)) / jnp.sum(
            jnp.exp(logits - jnp.max(logits, axis=1, keepdims=True)),
            axis=1,
            keepdims=True,
        )

    local = __import__("jax").lax.fori_loop(0, int(iterations), update, local)
    adjusted = means + 0.5 * (1.0 - 2.0 * local) * variances
    return 0.5 * jnp.sum(local * local * variances, axis=1) + __import__(
        "jax.scipy", fromlist=["special"]
    ).special.logsumexp(adjusted, axis=1)


def fit_tilted_variational(
    data,
    *,
    initial=None,
    max_iter=300,
):
    """Fit the multinomial EC GLMM with the tilted softmax ELBO bound.

    EC-to-isoform allocations are optimized analytically.  The local tilted
    parameters are updated to their fixed point inside each global update.
    """
    (
        jax,
        jnp,
        jsp,
        counts,
        mappings,
        design,
        cluster_index,
        n_clusters,
    ) = _prepared_arrays(data)
    dimension = data.n_isoforms - 1
    coefficient_count = design.shape[1] * dimension
    variational_count = n_clusters * dimension
    if initial is None:
        parameters = np.r_[
            np.zeros(coefficient_count + variational_count),
            np.full(variational_count, np.log(0.3)),
            np.log(0.5),
        ]
    else:
        parameters = np.asarray(initial, dtype=float)

    def unpack(value):
        coefficients = value[:coefficient_count].reshape(design.shape[1], dimension)
        start = coefficient_count
        means = value[start : start + variational_count].reshape(n_clusters, dimension)
        start += variational_count
        log_sds = value[start : start + variational_count].reshape(n_clusters, dimension)
        return coefficients, means, log_sds, value[-1]

    def bound(value):
        coefficients, random_means, random_log_sds, log_prior_sd = unpack(value)
        free_means = design @ coefficients + random_means[cluster_index]
        means = jnp.concatenate(
            (free_means, jnp.zeros((len(design), 1), dtype=free_means.dtype)), axis=1
        )
        free_variances = jnp.exp(2.0 * random_log_sds)[cluster_index]
        variances = jnp.concatenate(
            (free_variances, jnp.zeros((len(design), 1), dtype=free_variances.dtype)),
            axis=1,
        )
        value_bound = jnp.asarray(0.0, dtype=means.dtype)
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
            denominator = _tilted_logsumexp_bound(
                jnp,
                means + jnp.log(lengths)[None, :],
                variances,
            )
            value_bound += jnp.sum(
                constants + jnp.sum(observed * numerator, axis=1) - totals * denominator
            )
        prior_variance = jnp.exp(2.0 * log_prior_sd)
        q_variance = jnp.exp(2.0 * random_log_sds)
        kl = 0.5 * jnp.sum(
            (q_variance + random_means * random_means) / prior_variance
            - 1.0
            + 2.0 * (log_prior_sd - random_log_sds)
        )
        return value_bound - kl

    value_and_gradient = jax.jit(jax.value_and_grad(lambda value: -bound(value)))

    def objective(value):
        result, gradient = value_and_gradient(jnp.asarray(value))
        return float(result), np.asarray(gradient, dtype=float)

    bounds = (
        [(-20.0, 20.0)] * (coefficient_count + variational_count)
        + [(-8.0, 3.0)] * variational_count
        + [(-8.0, 3.0)]
    )
    result = scipy.optimize.minimize(
        objective,
        parameters,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxiter": int(max_iter), "ftol": 1e-8, "gtol": 1e-4},
    )
    coefficients, means, log_sds, log_prior_sd = unpack(jnp.asarray(result.x))
    final, gradient = objective(result.x)
    return {
        "method": "tilted_elbo",
        "family": "multinomial",
        "objective": -final,
        "parameters": np.asarray(result.x),
        "coefficients": np.asarray(coefficients),
        "random_effect_sd": float(np.exp(log_prior_sd)),
        "random_effect_mean": np.asarray(means),
        "random_effect_variance": np.exp(2.0 * np.asarray(log_sds)),
        "concentration": np.inf,
        "converged": bool(result.success or np.linalg.norm(gradient, ord=np.inf) <= 1e-3),
        "iterations": int(result.nit),
        "gradient_norm": float(np.linalg.norm(gradient, ord=np.inf)),
        "message": str(result.message),
    }


def fit_variational(
    data,
    *,
    family="multinomial",
    alpha=1.0,
    samples=64,
    seed=0,
    initial=None,
    max_iter=300,
):
    """Fit diagonal-Gaussian VI using an ELBO or variational Renyi bound."""
    if not (0 < float(alpha) <= 1):
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
    ) = _prepared_arrays(data)
    dimension = data.n_isoforms - 1
    coefficient_count = design.shape[1] * dimension
    variational_count = n_clusters * dimension
    rng = np.random.default_rng(seed)
    half = max(1, int(samples) // 2)
    noise = rng.standard_normal((half, n_clusters, dimension))
    noise = np.concatenate((noise, -noise), axis=0)[: int(samples)]
    noise = jnp.asarray(noise)
    membership = jax.nn.one_hot(
        cluster_index, n_clusters, dtype=jnp.float64
    )
    if initial is None:
        parameters = np.r_[
            np.zeros(coefficient_count + variational_count),
            np.full(variational_count, np.log(0.3)),
            np.log(0.5),
        ]
        if family == "dirichlet_multinomial":
            parameters = np.r_[parameters, np.log(100.0)]
    else:
        parameters = np.asarray(initial, dtype=float)

    def unpack(value):
        coefficients = value[:coefficient_count].reshape(design.shape[1], dimension)
        start = coefficient_count
        means = value[start : start + variational_count].reshape(n_clusters, dimension)
        start += variational_count
        log_sds = value[start : start + variational_count].reshape(n_clusters, dimension)
        start += variational_count
        log_prior_sd = value[start]
        log_kappa = (
            value[start + 1]
            if family == "dirichlet_multinomial"
            else value[0] * 0.0
        )
        return coefficients, means, log_sds, log_prior_sd, log_kappa

    def log_weights(value):
        coefficients, means, log_sds, log_prior_sd, log_kappa = unpack(value)
        random_effects = means[None, :, :] + jnp.exp(log_sds)[None, :, :] * noise
        fixed = design @ coefficients

        def one_sample(effects):
            free_logits = fixed + effects[cluster_index]
            logits = jnp.concatenate(
                (free_logits, jnp.zeros((len(design), 1), dtype=free_logits.dtype)), axis=1
            )
            likelihood_rows = _log_likelihood_rows(
                jnp, jsp, counts, mappings, logits, family, log_kappa
            )
            likelihood = membership.T @ likelihood_rows
            prior_sd = jnp.exp(log_prior_sd)
            log_prior = -0.5 * jnp.sum((effects / prior_sd) ** 2, axis=1)
            log_prior -= dimension * (
                log_prior_sd + 0.5 * jnp.log(2.0 * jnp.pi)
            )
            standardized = (effects - means) / jnp.exp(log_sds)
            log_q = -0.5 * jnp.sum(standardized * standardized, axis=1)
            log_q -= dimension * 0.5 * jnp.log(2.0 * jnp.pi)
            log_q -= jnp.sum(log_sds, axis=1)
            return likelihood + log_prior - log_q

        return jax.vmap(one_sample)(random_effects)

    def variational_bound(value):
        weights = log_weights(value)
        if float(alpha) == 1.0:
            return jnp.sum(jnp.mean(weights, axis=0))
        scale = 1.0 - float(alpha)
        cluster_bounds = jsp.special.logsumexp(
            scale * weights, axis=0
        ) / scale - jnp.log(len(weights)) / scale
        return jnp.sum(cluster_bounds)

    value_and_gradient = jax.jit(jax.value_and_grad(lambda value: -variational_bound(value)))

    def objective(value):
        result, gradient = value_and_gradient(jnp.asarray(value))
        return float(result), np.asarray(gradient, dtype=float)

    bounds = (
        [(-20.0, 20.0)] * (coefficient_count + variational_count)
        + [(-8.0, 3.0)] * variational_count
        + [(-8.0, 3.0)]
    )
    if family == "dirichlet_multinomial":
        bounds.append((-6.0, np.log(1e7)))
    if int(max_iter) > 0:
        result = scipy.optimize.minimize(
            objective,
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
    coefficients, means, log_sds, log_prior_sd, log_kappa = unpack(
        jnp.asarray(optimized)
    )
    final, gradient = objective(optimized)
    final_log_weights = np.asarray(log_weights(jnp.asarray(optimized)))
    importance_scale = 1.0 if float(alpha) == 1.0 else 1.0 - float(alpha)
    normalized_weights = scipy.special.softmax(
        importance_scale * final_log_weights, axis=0
    )
    cluster_ess = 1.0 / np.sum(np.square(normalized_weights), axis=0)
    importance_ess = float(np.median(cluster_ess))
    return {
        "method": "elbo_mc" if float(alpha) == 1.0 else f"renyi_{alpha:g}",
        "family": family,
        "alpha": float(alpha),
        "objective": -final,
        "parameters": optimized,
        "coefficients": np.asarray(coefficients),
        "random_effect_sd": float(np.exp(log_prior_sd)),
        "random_effect_mean": np.asarray(means),
        "random_effect_variance": np.exp(2.0 * np.asarray(log_sds)),
        "concentration": (
            float(np.exp(log_kappa)) if family == "dirichlet_multinomial" else np.inf
        ),
        "importance_ess": importance_ess,
        "minimum_importance_ess": float(np.min(cluster_ess)),
        "importance_samples": int(len(final_log_weights)),
        "converged": bool(success or np.linalg.norm(gradient, ord=np.inf) <= 1e-3),
        "iterations": iterations,
        "gradient_norm": float(np.linalg.norm(gradient, ord=np.inf)),
        "message": message,
    }


def evaluate_variational_objectives(
    data,
    fit,
    *,
    alpha=0.5,
    samples=2048,
    seed=1,
):
    """Evaluate exact-likelihood Monte Carlo ELBO and Renyi bounds for a fit."""
    mean = np.asarray(fit["random_effect_mean"], dtype=float)
    variance = np.maximum(np.asarray(fit["random_effect_variance"], dtype=float), 1e-12)
    dimension = data.n_isoforms - 1
    n_clusters = len(np.unique(data.clusters))
    coefficient_count = data.design.shape[1] * dimension
    parameters = np.r_[
        np.asarray(fit["coefficients"]).ravel(),
        mean.ravel(),
        0.5 * np.log(variance).ravel(),
        np.log(float(fit["random_effect_sd"])),
    ]
    family = fit["family"]
    if family == "dirichlet_multinomial":
        parameters = np.r_[parameters, np.log(float(fit["concentration"]))]
    expected = coefficient_count + 2 * n_clusters * dimension + 1 + int(
        family == "dirichlet_multinomial"
    )
    if len(parameters) != expected:
        raise ValueError("fit posterior shape does not match the EC GLMM data")
    elbo = fit_variational(
        data,
        family=family,
        alpha=1.0,
        samples=samples,
        seed=seed,
        initial=parameters,
        max_iter=0,
    )["objective"]
    renyi = fit_variational(
        data,
        family=family,
        alpha=alpha,
        samples=samples,
        seed=seed,
        initial=parameters,
        max_iter=0,
    )["objective"]
    return {"elbo_mc": float(elbo), f"renyi_{alpha:g}": float(renyi)}


def warm_start(fit, n_design_columns):
    """Expand a fitted model into an initial vector for a larger design."""
    coefficients = np.asarray(fit["coefficients"], dtype=float)
    if n_design_columns < len(coefficients):
        raise ValueError("warm-start design cannot have fewer columns")
    expanded = np.zeros((int(n_design_columns), coefficients.shape[1]))
    expanded[: len(coefficients)] = coefficients
    if fit["method"] == "laplace":
        values = np.r_[expanded.ravel(), np.log(fit["random_effect_sd"])]
    else:
        values = np.r_[
            expanded.ravel(),
            np.asarray(fit["random_effect_mean"]).ravel(),
            0.5 * np.log(np.maximum(fit["random_effect_variance"], 1e-12)).ravel(),
            np.log(fit["random_effect_sd"]),
        ]
    if fit["family"] == "dirichlet_multinomial":
        values = np.r_[values, np.log(fit["concentration"])]
    if fit.get("observation_noise", False):
        values = np.r_[values, np.log(fit["observation_noise_sd"])]
    return values


def variational_warm_start(fit, n_design_columns):
    """Convert a Laplace or VI posterior into a Gaussian-VI initial vector."""
    coefficients = np.asarray(fit["coefficients"], dtype=float)
    if n_design_columns < len(coefficients):
        raise ValueError("warm-start design cannot have fewer columns")
    expanded = np.zeros((int(n_design_columns), coefficients.shape[1]))
    expanded[: len(coefficients)] = coefficients
    values = np.r_[
        expanded.ravel(),
        np.asarray(fit["random_effect_mean"]).ravel(),
        0.5
        * np.log(
            np.maximum(fit["random_effect_variance"], 1e-12)
        ).ravel(),
        np.log(fit["random_effect_sd"]),
    ]
    if fit["family"] == "dirichlet_multinomial":
        values = np.r_[values, np.log(fit["concentration"])]
    return values
