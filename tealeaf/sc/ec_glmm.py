"""Mixed models for paired-primer equivalence-class counts.

The latent linear predictor is defined on isoforms, while observations remain
equivalence-class (EC) counts.  This keeps ambiguity in the likelihood instead
of first converting EC counts to estimated isoform or splice-path counts.
"""

from __future__ import annotations

from dataclasses import dataclass
import time

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


def laplace_objective_shape_key(data, *, observation_noise=False):
    """Return the static array-shape signature of a Laplace objective.

    Counts, compatibility values, designs, and cluster assignments are dynamic
    JAX arguments.  Two fits can therefore share a compiled objective when
    these arrays have the same shapes and both use either a matrix design or a
    fixed-effect tensor.  Observation-noise fits additionally depend on the
    largest cluster size because that determines their latent dimension.
    """
    _, cluster_sizes = np.unique(data.clusters, return_counts=True)
    maximum_cluster_size = int(cluster_sizes.max())
    tensor_shape = (
        None
        if data.fixed_effect_tensor is None
        else tuple(data.fixed_effect_tensor.shape)
    )
    return (
        tuple(tuple(value.shape) for value in data.counts),
        tuple(tuple(value.shape) for value in data.compatibility),
        tuple(data.design.shape),
        tensor_shape,
        int(len(cluster_sizes)),
        maximum_cluster_size if observation_noise else 0,
    )


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
        if observed.shape[1] == 0:
            continue
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
    coefficient_count = (
        data.design.shape[1] * (data.n_isoforms - 1)
        if data.fixed_effect_tensor is None
        else data.fixed_effect_tensor.shape[2]
    )
    values = np.zeros(coefficient_count + 1, dtype=float)
    values[-1] = np.log(0.5)
    if family == "dirichlet_multinomial":
        values = np.r_[values, np.log(100.0)]
    return values


def _projected_adam(
    objective,
    initial,
    bounds,
    *,
    max_iter,
    learning_rate=0.03,
    beta1=0.9,
    beta2=0.999,
    epsilon=1e-8,
    gradient_tolerance=1e-5,
):
    """Minimize a differentiable objective with projected Adam.

    This small optimizer is intended for controlled comparisons with the
    existing bounded L-BFGS fit. It retains the best finite iterate because a
    fixed first-order schedule need not decrease the objective at every step.
    """
    parameters = np.asarray(initial, dtype=float).copy()
    lower = np.asarray([
        -np.inf if bound[0] is None else bound[0] for bound in bounds
    ])
    upper = np.asarray([
        np.inf if bound[1] is None else bound[1] for bound in bounds
    ])
    first_moment = np.zeros_like(parameters)
    second_moment = np.zeros_like(parameters)
    best_parameters = parameters.copy()
    best_value = np.inf
    best_gradient = np.full_like(parameters, np.nan)
    evaluations = 0
    converged = False
    message = "Adam reached the iteration limit"
    iteration = 0
    for iteration in range(1, int(max_iter) + 1):
        value, gradient = objective(parameters)
        evaluations += 1
        if not np.isfinite(value) or not np.isfinite(gradient).all():
            message = "Adam encountered a non-finite objective or gradient"
            break
        if value < best_value:
            best_value = float(value)
            best_parameters = parameters.copy()
            best_gradient = np.asarray(gradient, dtype=float).copy()
        if np.linalg.norm(gradient, ord=np.inf) <= float(gradient_tolerance):
            converged = True
            message = "Adam gradient tolerance reached"
            break
        first_moment = beta1 * first_moment + (1.0 - beta1) * gradient
        second_moment = beta2 * second_moment + (1.0 - beta2) * gradient**2
        corrected_first = first_moment / (1.0 - beta1**iteration)
        corrected_second = second_moment / (1.0 - beta2**iteration)
        progress = (iteration - 1.0) / max(float(max_iter), 1.0)
        rate = float(learning_rate) * 0.5 * (1.0 + np.cos(np.pi * progress))
        parameters -= rate * corrected_first / (
            np.sqrt(corrected_second) + float(epsilon)
        )
        parameters = np.clip(parameters, lower, upper)
    if not np.isfinite(best_value):
        best_value, best_gradient = objective(initial)
        evaluations += 1
        best_parameters = np.asarray(initial, dtype=float).copy()
    return {
        "x": best_parameters,
        "fun": float(best_value),
        "jac": np.asarray(best_gradient, dtype=float),
        "nit": int(iteration),
        "nfev": int(evaluations),
        "success": bool(converged),
        "message": message,
    }


def fit_laplace(
    data,
    *,
    family="multinomial",
    observation_noise=False,
    initial=None,
    max_iter=200,
    mode_steps=30,
    mode_gradient="unrolled",
    optimizer="lbfgs",
    adam_learning_rate=0.03,
    adam_steps=100,
    objective_cache=None,
    objective_cache_key=None,
    mode_warm_start=False,
    mode_tolerance=0.0,
    return_outer_hessian=False,
):
    """Fit an EC-count GLMM by Laplace approximation.

    The latent model has an isotropic isoform-logit random intercept per
    cluster.

    When ``observation_noise`` is true, the latent block for each cluster also
    contains an independent isotropic logistic-normal effect for every
    observation in that cluster.  This supplies multinomial overdispersion on
    the isoform-logit scale while retaining the EC observation model.

    ``mode_gradient="implicit"`` differentiates the converged mode equation
    instead of backpropagating through the Newton iterations.  A shared
    ``objective_cache_key`` may be shared when array shapes agree.  Counts,
    mappings, designs, fixed-effect tensors, and cluster layouts are passed
    dynamically and may differ between calls sharing that key.  Use
    :func:`laplace_objective_shape_key` to construct a safe global key.
    """
    if observation_noise and family != "multinomial":
        raise ValueError(
            "observation noise is currently supported for multinomial fits"
        )
    if mode_gradient not in ("unrolled", "implicit"):
        raise ValueError("mode_gradient must be 'unrolled' or 'implicit'")
    if optimizer not in ("lbfgs", "adam", "adam_lbfgs"):
        raise ValueError("optimizer must be 'lbfgs', 'adam', or 'adam_lbfgs'")
    if float(adam_learning_rate) <= 0:
        raise ValueError("adam_learning_rate must be positive")
    if int(adam_steps) < 0:
        raise ValueError("adam_steps must be nonnegative")
    if float(mode_tolerance) < 0:
        raise ValueError("mode_tolerance must be nonnegative")
    if float(mode_tolerance) > 0 and mode_gradient != "implicit":
        raise ValueError("adaptive mode termination requires implicit gradients")
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
    base = _initial_outer(data, family)
    if observation_noise:
        initial = (
            np.r_[base, np.log(0.3)]
            if initial is None
            else np.asarray(initial)
        )
    else:
        initial = base if initial is None else np.asarray(initial)

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
        coefficients = outer[:coefficient_count]
        if fixed_effect_tensor is None:
            coefficients = coefficients.reshape(design.shape[1], dimension)
        log_prior_sd = outer[coefficient_count]
        cursor = coefficient_count + 1
        log_kappa = outer[0] * 0.0
        if family == "dirichlet_multinomial":
            log_kappa = outer[cursor]
            cursor += 1
        log_noise_sd = outer[cursor] if observation_noise else outer[0] * 0.0
        return coefficients, log_prior_sd, log_noise_sd, log_kappa

    def fixed_logits(coefficients, active_design, active_fixed_effect_tensor):
        if fixed_effect_tensor is None:
            return active_design @ coefficients
        return jnp.einsum(
            "ndp,p->nd",
            active_fixed_effect_tensor,
            coefficients,
        )

    def prior_log_sds(outer):
        _, log_prior_sd, log_noise_sd, _ = unpack_outer(outer)
        values = jnp.full((latent_blocks, dimension), log_prior_sd)
        if observation_noise:
            values = values.at[1:].set(log_noise_sd)
        return values.ravel()

    def observation_log_likelihood(
        free_logits,
        observed_rows,
        active_mappings,
        log_kappa,
    ):
        logits = jnp.concatenate((free_logits, jnp.zeros(1, dtype=free_logits.dtype)))
        abundance = jnp.exp(logits - jnp.max(logits))
        value = jnp.asarray(0.0, dtype=logits.dtype)
        for observed, mapping in zip(observed_rows, active_mappings):
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
        in_axes=(0, (0,) * len(counts), (None,) * len(mappings), None),
    )
    observation_hessian = jax.vmap(
        jax.hessian(observation_log_likelihood, argnums=0),
        in_axes=(0, (0,) * len(counts), (None,) * len(mappings), None),
    )

    def mode_derivatives(
        modes,
        outer,
        active_design,
        active_fixed_effect_tensor,
        active_counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    ):
        coefficients, _, _, log_kappa = unpack_outer(outer)
        fixed_values = fixed_logits(
            coefficients,
            active_design,
            active_fixed_effect_tensor,
        )
        random_logits = jnp.einsum(
            "ndi,ni->nd", active_selection, modes[active_cluster_index]
        )
        free_logits = fixed_values + random_logits
        score_rows = observation_gradient(
            free_logits,
            active_counts,
            active_mappings,
            log_kappa,
        )
        log_sds = prior_log_sds(outer)
        variances = jnp.exp(2.0 * log_sds)
        latent_score_rows = jnp.einsum(
            "ndi,nd->ni", active_selection, score_rows
        )
        score = active_membership.T @ latent_score_rows - modes / variances
        curvature_rows = -observation_hessian(
            free_logits,
            active_counts,
            active_mappings,
            log_kappa,
        )
        latent_curvature_rows = jnp.einsum(
            "nai,nab,nbj->nij",
            active_selection,
            curvature_rows,
            active_selection,
        )
        hessian = jnp.einsum(
            "im,ijk->mjk", active_membership, latent_curvature_rows
        ) + jnp.diag(1.0 / variances)[None, :, :]
        return fixed_values, score, hessian, variances, log_kappa

    def negative_joint(
        modes,
        outer,
        active_design,
        active_fixed_effect_tensor,
        active_counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    ):
        fixed_values, _, _, variances, log_kappa = mode_derivatives(
            modes,
            outer,
            active_design,
            active_fixed_effect_tensor,
            active_counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        )
        random_logits = jnp.einsum(
            "ndi,ni->nd", active_selection, modes[active_cluster_index]
        )
        free_logits = fixed_values + random_logits
        logits = jnp.concatenate(
            (
                free_logits,
                jnp.zeros((len(active_design), 1), dtype=free_logits.dtype),
            ),
            axis=1,
        )
        log_sds = prior_log_sds(outer)
        log_prior = -0.5 * jnp.sum(modes * modes / variances)
        log_prior -= active_membership.shape[1] * (
            jnp.sum(log_sds)
            + 0.5 * active_selection.shape[2] * jnp.log(2.0 * jnp.pi)
        )
        return -(
            _log_likelihood(
                jnp,
                jsp,
                active_counts,
                active_mappings,
                logits,
                family,
                log_kappa,
            )
            + log_prior
        )

    def solve_posterior_mode(
        outer,
        initial_modes,
        active_design,
        active_fixed_effect_tensor,
        active_counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    ):
        current = initial_modes

        def update(_, value):
            _, score, hessian, _, _ = mode_derivatives(
                value,
                outer,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )
            eigen_floor = jnp.maximum(
                1e-6 - jnp.linalg.eigvalsh(hessian)[:, 0], 0.0
            )
            hessian = hessian + eigen_floor[:, None, None] * jnp.eye(
                active_selection.shape[2]
            )
            step = jnp.linalg.solve(hessian, score[..., None])[..., 0]
            maximum = jnp.max(jnp.abs(step), axis=1, keepdims=True)
            scale = jnp.minimum(1.0, 2.0 / (maximum + 1e-12))
            return value + scale * step

        if float(mode_tolerance) <= 0:
            return jax.lax.fori_loop(0, int(mode_steps), update, current)

        def condition(state):
            iteration, _, maximum_score = state
            return jnp.logical_and(
                iteration < int(mode_steps),
                maximum_score > float(mode_tolerance),
            )

        def adaptive_update(state):
            iteration, value, _ = state
            updated = update(iteration, value)
            _, score, _, _, _ = mode_derivatives(
                updated,
                outer,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )
            return iteration + 1, updated, jnp.max(jnp.abs(score))

        return jax.lax.while_loop(
            condition,
            adaptive_update,
            (0, current, jnp.asarray(jnp.inf, dtype=current.dtype)),
        )[1]

    if mode_gradient == "implicit":
        @jax.custom_vjp
        def posterior_mode(
            outer,
            initial_modes,
            active_design,
            active_fixed_effect_tensor,
            active_counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        ):
            return solve_posterior_mode(
                outer,
                initial_modes,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )

        def posterior_mode_forward(
            outer,
            initial_modes,
            active_design,
            active_fixed_effect_tensor,
            active_counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        ):
            modes = solve_posterior_mode(
                outer,
                initial_modes,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )
            return modes, (
                modes,
                outer,
                initial_modes,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )

        def posterior_mode_backward(residuals, mode_cotangent):
            (
                modes,
                outer,
                initial_modes,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            ) = residuals
            _, _, hessian, _, _ = mode_derivatives(
                modes,
                outer,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )
            solved = jnp.linalg.solve(
                hessian,
                mode_cotangent[..., None],
            )[..., 0]
            _, score_pullback = jax.vjp(
                lambda value: mode_derivatives(
                    modes,
                    value,
                    active_design,
                    active_fixed_effect_tensor,
                    active_counts,
                    active_mappings,
                    active_membership,
                    active_selection,
                    active_cluster_index,
                )[1],
                outer,
            )
            outer_cotangent = score_pullback(solved)[0]
            return (
                outer_cotangent,
                jnp.zeros_like(initial_modes),
                jnp.zeros_like(active_design),
                jnp.zeros_like(active_fixed_effect_tensor),
                tuple(jnp.zeros_like(value) for value in active_counts),
                tuple(jnp.zeros_like(value) for value in active_mappings),
                jnp.zeros_like(active_membership),
                jnp.zeros_like(active_selection),
                None,
            )

        posterior_mode.defvjp(
            posterior_mode_forward,
            posterior_mode_backward,
        )
    else:
        posterior_mode = solve_posterior_mode

    def objective(
        outer,
        initial_modes,
        active_design,
        active_fixed_effect_tensor,
        active_counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    ):
        modes = posterior_mode(
            outer,
            initial_modes,
            active_design,
            active_fixed_effect_tensor,
            active_counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        )
        _, _, hessian, _, _ = mode_derivatives(
            modes,
            outer,
            active_design,
            active_fixed_effect_tensor,
            active_counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        )
        sign, log_determinant = jnp.linalg.slogdet(hessian)
        correction = 0.5 * jnp.sum(log_determinant)
        correction -= 0.5 * modes.size * jnp.log(2.0 * jnp.pi)
        value = jnp.where(
            jnp.all(sign > 0),
            negative_joint(
                modes,
                outer,
                active_design,
                active_fixed_effect_tensor,
                active_counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            ) + correction,
            jnp.inf,
        )
        return value, modes

    cache_key = None
    if objective_cache is not None and objective_cache_key is not None:
        cache_key = (
            objective_cache_key,
            "laplace_value_and_gradient",
            family,
            bool(observation_noise),
            int(mode_steps),
            mode_gradient,
            float(mode_tolerance),
        )
    if cache_key is not None and cache_key in objective_cache:
        value_and_gradient = objective_cache[cache_key]
    else:
        value_and_gradient = jax.jit(
            jax.value_and_grad(objective, argnums=0, has_aux=True)
        )
        if cache_key is not None:
            objective_cache[cache_key] = value_and_gradient
    hessian_function = None
    if return_outer_hessian:
        hessian_cache_key = None
        if objective_cache is not None and objective_cache_key is not None:
            hessian_cache_key = (
                objective_cache_key,
                "laplace_hessian",
                family,
                bool(observation_noise),
                int(mode_steps),
                mode_gradient,
                float(mode_tolerance),
            )
        if hessian_cache_key is not None and hessian_cache_key in objective_cache:
            hessian_function = objective_cache[hessian_cache_key]
        else:
            hessian_function = jax.jit(
                jax.hessian(lambda *arguments: objective(*arguments)[0], argnums=0)
            )
            if hessian_cache_key is not None:
                objective_cache[hessian_cache_key] = hessian_function

    active_fixed_effect_tensor = (
        jnp.empty((0,), dtype=jnp.float64)
        if fixed_effect_tensor is None
        else fixed_effect_tensor
    )
    active_design = design
    active_mappings = mappings
    active_membership = membership
    active_selection = selection
    active_cluster_index = cluster_index
    zero_modes = jnp.zeros(
        (active_membership.shape[1], active_selection.shape[2]),
        dtype=jnp.float64,
    )
    accepted_modes = zero_modes
    evaluated_modes = {}

    def scipy_objective(parameters):
        (value, modes), gradient = value_and_gradient(
            jnp.asarray(parameters),
            accepted_modes if mode_warm_start else zero_modes,
            active_design,
            active_fixed_effect_tensor,
            counts,
            active_mappings,
            active_membership,
            active_selection,
            active_cluster_index,
        )
        if mode_warm_start:
            evaluated_modes[np.asarray(parameters, dtype=float).tobytes()] = modes
        return float(value), np.asarray(gradient, dtype=float)

    evaluation_count = 0
    evaluation_seconds = 0.0

    def timed_scipy_objective(parameters):
        nonlocal evaluation_count, evaluation_seconds
        started = time.perf_counter()
        result = scipy_objective(parameters)
        evaluation_seconds += time.perf_counter() - started
        evaluation_count += 1
        return result

    initial_value, _ = timed_scipy_objective(initial)
    initial_objective_seconds = evaluation_seconds
    optimizer_scale = max(1.0, abs(initial_value))

    def scaled_scipy_objective(parameters):
        value, gradient = timed_scipy_objective(parameters)
        return value / optimizer_scale, gradient / optimizer_scale

    bounds = [(-20.0, 20.0)] * coefficient_count
    bounds.append((-8.0, 3.0))
    if family == "dirichlet_multinomial":
        bounds.append((-6.0, np.log(1e7)))
    if observation_noise:
        bounds.append((-8.0, 3.0))
    optimizer_started = time.perf_counter()
    adam_result = None
    lbfgs_result = None
    if optimizer in ("adam", "adam_lbfgs"):
        steps = int(max_iter) if optimizer == "adam" else int(adam_steps)
        adam_result = _projected_adam(
            scaled_scipy_objective,
            initial,
            bounds,
            max_iter=steps,
            learning_rate=adam_learning_rate,
            gradient_tolerance=1e-5,
        )
    if optimizer in ("lbfgs", "adam_lbfgs") and int(max_iter) > 0:
        lbfgs_initial = (
            np.asarray(initial, dtype=float)
            if adam_result is None else np.asarray(adam_result["x"], dtype=float)
        )
        def accept_iterate(parameters):
            nonlocal accepted_modes
            if not mode_warm_start:
                return
            key = np.asarray(parameters, dtype=float).tobytes()
            if key in evaluated_modes:
                accepted_modes = evaluated_modes[key]
            evaluated_modes.clear()

        lbfgs_result = scipy.optimize.minimize(
            scaled_scipy_objective,
            lbfgs_initial,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            callback=accept_iterate,
            options={
                "maxiter": int(max_iter),
                "ftol": 1e-12,
                "gtol": 1e-5,
                "maxls": 100,
            },
        )
    optimizer_seconds = time.perf_counter() - optimizer_started
    if lbfgs_result is not None:
        parameters = np.asarray(lbfgs_result.x)
        optimizer_success = bool(lbfgs_result.success)
        optimizer_iterations = int(lbfgs_result.nit) + (
            0 if adam_result is None else int(adam_result["nit"])
        )
        optimizer_message = str(lbfgs_result.message)
    elif adam_result is not None:
        parameters = np.asarray(adam_result["x"])
        optimizer_success = bool(adam_result["success"])
        optimizer_iterations = int(adam_result["nit"])
        optimizer_message = str(adam_result["message"])
    else:
        parameters = np.asarray(initial, dtype=float)
        optimizer_success = False
        optimizer_iterations = 0
        optimizer_message = "optimizer skipped because max_iter is zero"
    outer = jnp.asarray(parameters)
    modes = posterior_mode(
        outer,
        accepted_modes if mode_warm_start else zero_modes,
        active_design,
        active_fixed_effect_tensor,
        counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    )
    _, mode_score, hessian, _, _ = mode_derivatives(
        modes,
        outer,
        active_design,
        active_fixed_effect_tensor,
        counts,
        active_mappings,
        active_membership,
        active_selection,
        active_cluster_index,
    )
    covariance = np.asarray(jnp.linalg.inv(hessian))
    marginal_variance = np.diagonal(covariance, axis1=1, axis2=2)
    final_value, final_gradient = timed_scipy_objective(parameters)
    outer_hessian = None
    if hessian_function is not None:
        outer_hessian = np.asarray(
            hessian_function(
                outer,
                accepted_modes if mode_warm_start else zero_modes,
                active_design,
                active_fixed_effect_tensor,
                counts,
                active_mappings,
                active_membership,
                active_selection,
                active_cluster_index,
            )
        )
    mode_score = np.asarray(mode_score)
    gradient_norm = float(np.linalg.norm(final_gradient, ord=np.inf))
    scaled_gradient_norm = gradient_norm / optimizer_scale
    mode_score_norm = float(np.linalg.norm(mode_score, ord=np.inf))
    coefficients, log_prior_sd, log_noise_sd, log_kappa = unpack_outer(outer)
    return {
        "method": "laplace",
        "family": family,
        "objective": final_value,
        "parameters": parameters,
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
            np.isfinite(final_value)
            and scaled_gradient_norm <= 1e-4
            and mode_score_norm <= 1e-4
        ),
        "optimizer": optimizer,
        "optimizer_success": optimizer_success,
        "iterations": optimizer_iterations,
        "adam_iterations": (
            0 if adam_result is None else int(adam_result["nit"])
        ),
        "lbfgs_iterations": (
            0 if lbfgs_result is None else int(lbfgs_result.nit)
        ),
        "objective_evaluations": int(evaluation_count),
        "objective_evaluation_seconds": float(evaluation_seconds),
        "initial_objective_seconds": float(initial_objective_seconds),
        "optimizer_seconds": float(optimizer_seconds),
        "gradient_norm": gradient_norm,
        "scaled_gradient_norm": scaled_gradient_norm,
        "outer_gradient": np.asarray(final_gradient),
        "outer_hessian": outer_hessian,
        "optimizer_scale": float(optimizer_scale),
        "fixed_effect_count": int(coefficient_count),
        "mode_steps": int(mode_steps),
        "mode_gradient": mode_gradient,
        "mode_warm_start": bool(mode_warm_start),
        "mode_tolerance": float(mode_tolerance),
        "mode_score_norm": mode_score_norm,
        "message": optimizer_message,
    }


def efficient_score_statistic(fit, tested_indices):
    """Return a fixed-effect-nuisance-adjusted score statistic.

    ``fit`` is an evaluation of the alternative objective at its constrained
    null point and ``tested_indices`` selects tested fixed effects.  Variance
    components remain fixed at their null estimates.  If the fixed-effect
    Hessian is :math:`H` and the tested gradient is :math:`g`, the statistic
    is :math:`g^T (H^{-1})_{tt} g`, the fixed-effect Schur-complement score
    quadratic.  Holding variance components fixed avoids unstable observed
    cross-curvature while retaining adjustment for all nuisance fixed effects.
    """
    hessian = fit.get("outer_hessian")
    if hessian is None:
        raise ValueError("fit does not contain an outer Hessian")
    gradient = np.asarray(fit["outer_gradient"], dtype=float)
    tested_indices = np.asarray(tested_indices, dtype=int)
    if tested_indices.ndim != 1 or not len(tested_indices):
        raise ValueError("tested_indices must be a nonempty vector")
    fixed_effect_count = int(fit["fixed_effect_count"])
    if np.any(tested_indices < 0) or np.any(tested_indices >= fixed_effect_count):
        raise ValueError("tested indices must select fixed effects")
    fixed_hessian = np.asarray(hessian, dtype=float)[
        :fixed_effect_count, :fixed_effect_count
    ]
    covariance = np.linalg.pinv(fixed_hessian, hermitian=True)
    tested_gradient = gradient[tested_indices]
    tested_covariance = covariance[np.ix_(tested_indices, tested_indices)]
    statistic = float(tested_gradient @ tested_covariance @ tested_gradient)
    return max(0.0, statistic)


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
