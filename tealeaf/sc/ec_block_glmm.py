"""Fixed-effect constructions for splice-block EC-count GLMMs."""

from __future__ import annotations

import numpy as np
import scipy.optimize
import scipy.stats

from . import differential, ec_glmm


def pooled_isoform_weights(data, max_iter=100):
    """Estimate one label-independent isoform mixture from pooled EC counts.

    The returned vector has dimension ``T``, where ``T`` is the number of
    isoforms. Primer totals remain conditioned upon, as in the fitted EC
    likelihood. The estimate is intended only to define fixed within-path
    mixtures before reducing the latent space.
    """
    totals = tuple(
        np.asarray(counts, dtype=float).sum(axis=0)
        for counts in data.counts
    )
    dimension = data.n_isoforms - 1

    def objective(free_logits):
        logits = np.r_[np.asarray(free_logits, dtype=float), 0.0]
        abundance = np.exp(logits - logits.max())
        value = 0.0
        for counts, mapping in zip(totals, data.compatibility):
            mass = np.asarray(mapping, dtype=float) @ abundance
            if counts.sum() > 0:
                value -= float(
                    counts @ np.log(np.maximum(mass, 1e-300))
                )
                value += float(
                    counts.sum()
                    * np.log(np.maximum(mass.sum(), 1e-300))
                )
        return value

    result = scipy.optimize.minimize(
        objective,
        np.zeros(dimension, dtype=float),
        method="L-BFGS-B",
        options={"maxiter": int(max_iter)},
    )
    logits = np.r_[np.asarray(result.x, dtype=float), 0.0]
    weights = np.exp(logits - logits.max())
    return weights / weights.sum()


def collapse_isoforms_to_paths(data, path_index, weights=None):
    """Collapse represented isoforms to paths while retaining nuisance isoforms.

    ``path_index`` has dimension ``T`` and maps represented isoforms to
    nonnegative path labels. Isoforms labeled -1 remain separate normalization
    nuisances. Each collapsed compatibility column is the fixed within-path
    weighted mean of its original columns. Counts are not collapsed, so the
    primer-specific multinomial likelihood retains covariance among ECs.
    """
    path_index = np.asarray(path_index, dtype=int)
    if path_index.shape != (data.n_isoforms,):
        raise ValueError("path_index must have one entry per isoform")
    paths = np.unique(path_index[path_index >= 0])
    if len(paths) < 2:
        raise ValueError(
            "a path collapse needs at least two represented paths"
        )
    if weights is None:
        weights = np.ones(data.n_isoforms, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (data.n_isoforms,) or np.any(weights < 0):
        raise ValueError(
            "weights must be a nonnegative vector of length T"
        )

    groups = [np.flatnonzero(path_index == path) for path in paths]
    groups.extend([
        np.asarray([index]) for index in np.flatnonzero(path_index < 0)
    ])
    collapsed_mappings = []
    for mapping in data.compatibility:
        mapping = np.asarray(mapping, dtype=float)
        columns = []
        for group in groups:
            local_weights = weights[group]
            if local_weights.sum() <= 0:
                local_weights = np.ones(len(group), dtype=float)
            local_weights = local_weights / local_weights.sum()
            columns.append(mapping[:, group] @ local_weights)
        collapsed_mappings.append(np.column_stack(columns))
    collapsed_path_index = np.r_[
        np.arange(len(paths), dtype=int),
        -np.ones(np.sum(path_index < 0), dtype=int),
    ]
    collapsed = ec_glmm.ECGLMMData(
        data.counts,
        tuple(collapsed_mappings),
        data.design,
        data.clusters,
    )
    return collapsed, collapsed_path_index


def estimate_path_logratios(
    data,
    path_index,
    *,
    baseline=None,
    covariance="conditional",
    max_iter=100,
):
    """Estimate each sample's path ILR vector and EC-derived covariance.

    ``data`` contains ``n`` samples and ``T`` isoforms, while ``path_index``
    has dimension ``T``. The returned values have shape ``n`` by ``S - 1``
    and the returned covariances have shape ``n`` by ``S - 1`` by ``S - 1``,
    where ``S`` is the number of represented paths. A single pooled,
    label-independent isoform baseline is used unless ``baseline`` is supplied.
    """
    path_index = np.asarray(path_index, dtype=int)
    if path_index.shape != (data.n_isoforms,):
        raise ValueError("path_index must have one entry per isoform")
    paths = np.unique(path_index[path_index >= 0])
    if not np.array_equal(paths, np.arange(len(paths))) or len(paths) < 2:
        raise ValueError("represented paths must be consecutive labels from zero")
    if covariance not in {"conditional", "profile"}:
        raise ValueError("covariance must be conditional or profile")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    baseline = np.asarray(baseline, dtype=float)
    if baseline.shape != (data.n_isoforms,) or np.any(baseline <= 0):
        raise ValueError("baseline must be a positive vector of length T")

    values = []
    covariances = []
    proportions = []
    iterations = []
    identifiable = []
    converged = []
    for sample in range(len(data.clusters)):
        observed = tuple(
            np.asarray(counts[sample]).ravel() for counts in data.counts
        )
        fit = differential.fit_path_perturbation(
            observed,
            data.compatibility,
            baseline,
            path_index,
            max_iter=max_iter,
        )
        covariance_result = fit.covariance
        if covariance == "profile":
            covariance_result = differential.profiled_path_covariance(
                fit.theta,
                path_index,
                data.compatibility,
                tuple(float(counts.sum()) for counts in observed),
            )
        values.append(fit.path_logratios)
        covariances.append(covariance_result.covariance)
        proportions.append(fit.path_proportions)
        iterations.append(fit.iterations)
        identifiable.append(covariance_result.identifiable)
        converged.append(fit.converged)
    return {
        "values": np.asarray(values),
        "covariances": np.asarray(covariances),
        "proportions": np.asarray(proportions),
        "iterations": np.asarray(iterations, dtype=int),
        "identifiable": np.asarray(identifiable, dtype=bool),
        "converged": np.asarray(converged, dtype=bool),
        "baseline": baseline,
        "covariance_method": covariance,
    }


def path_wald_test(
    data,
    path_index,
    nuisance_design,
    tested_design,
    *,
    baseline=None,
    covariance="conditional",
    max_iter=100,
    cluster_adjustment="cr1",
):
    """Estimate paths once and test fixed effects by a clustered Wald test.

    The nuisance and tested designs have ``n`` rows. Mouse-level correlation is
    represented by a sandwich covariance over the cluster labels in ``data``.
    The tested coefficients receive one finite-cluster multivariate Wald test.
    This is an ablation API: structured real-data nulls show that its analytic
    p-values are anti-conservative for the fitted multi-level cell-type designs.
    """
    nuisance_design = np.asarray(nuisance_design, dtype=float)
    tested_design = np.asarray(tested_design, dtype=float)
    if (
        nuisance_design.ndim != 2
        or tested_design.ndim != 2
        or len(nuisance_design) != len(data.clusters)
        or len(tested_design) != len(data.clusters)
    ):
        raise ValueError("designs and EC samples must align")
    if tested_design.shape[1] < 1:
        raise ValueError("tested_design needs at least one column")
    estimates = estimate_path_logratios(
        data,
        path_index,
        baseline=baseline,
        covariance=covariance,
        max_iter=max_iter,
    )
    result = path_wald_from_estimates(
        estimates,
        nuisance_design,
        tested_design,
        data.clusters,
        cluster_adjustment=cluster_adjustment,
    )
    result["path_estimates"] = estimates
    return result


def path_wald_from_estimates(
    estimates,
    nuisance_design,
    tested_design,
    clusters,
    *,
    cluster_adjustment="cr1",
):
    """Test fixed effects from precomputed path ILRs and covariances."""
    nuisance_design = np.asarray(nuisance_design, dtype=float)
    tested_design = np.asarray(tested_design, dtype=float)
    clusters = np.asarray(clusters)
    values = np.asarray(estimates["values"], dtype=float)
    covariances = np.asarray(estimates["covariances"], dtype=float)
    if (
        nuisance_design.ndim != 2
        or tested_design.ndim != 2
        or len(nuisance_design) != len(values)
        or len(tested_design) != len(values)
        or len(clusters) != len(values)
    ):
        raise ValueError("designs, clusters, and path estimates must align")
    usable = (
        estimates["converged"]
        & estimates["identifiable"]
        & np.isfinite(values).all(axis=1)
        & np.isfinite(covariances).all(axis=(1, 2))
    )
    design = np.column_stack((nuisance_design, tested_design))
    tested_columns = np.arange(nuisance_design.shape[1], design.shape[1])
    if (
        usable.sum() <= design.shape[1]
        or np.linalg.matrix_rank(design[usable]) < design.shape[1]
        or len(np.unique(clusters[usable])) <= design.shape[1]
    ):
        raise np.linalg.LinAlgError(
            "usable path estimates do not identify the regression design"
        )
    result = differential.clustered_multivariate_wald_test(
        values[usable],
        covariances[usable],
        design[usable],
        tested_columns,
        clusters[usable],
        cluster_adjustment=cluster_adjustment,
    )
    estimates["usable"] = usable
    return result


def block_effect_bases(path_index):
    """Return block-path and orthogonal nuisance bases in reference logits.

    ``path_index[t]`` is the represented block path for isoform ``t`` or -1
    when the isoform is retained only as a likelihood-normalization nuisance.
    """
    path_index = np.asarray(path_index, dtype=int)
    paths = np.unique(path_index[path_index >= 0])
    if len(paths) < 2:
        raise ValueError("a splice-block test needs at least two paths")
    remap = {path: position for position, path in enumerate(paths)}
    basis = differential.helmert_basis(len(paths))
    isoform_effects = np.zeros((len(path_index), len(paths) - 1), dtype=float)
    for isoform, path in enumerate(path_index):
        if path >= 0:
            isoform_effects[isoform] = basis[remap[path]]
    centered = isoform_effects - isoform_effects.mean(axis=0, keepdims=True)
    complete, _ = np.linalg.qr(
        np.column_stack((np.ones(len(path_index)), centered)), mode="complete"
    )
    nuisance_full = complete[:, len(paths) :]

    def reference_logits(values):
        return values[:-1] - values[-1][None, :]

    return reference_logits(centered), reference_logits(nuisance_full)


def path_contrast_matrix(path_index):
    """Map block-path contrasts to reference-isoform logits."""
    return block_effect_bases(path_index)[0]


def unrestricted_tensor(design, dimension):
    """Expand a standard design into an unrestricted free-logit tensor."""
    design = np.asarray(design, dtype=float)
    result = np.zeros(
        (len(design), int(dimension), design.shape[1] * int(dimension)),
        dtype=float,
    )
    for column in range(design.shape[1]):
        start = column * int(dimension)
        result[:, :, start : start + int(dimension)] = (
            design[:, column, None, None]
            * np.eye(int(dimension))[None, :, :]
        )
    return result


def full_fixed_effect_tensor(nuisance_design, tested_design, dimension):
    """Build the canonical unrestricted alternative for a tested effect."""
    return np.concatenate(
        (
            unrestricted_tensor(nuisance_design, dimension),
            unrestricted_tensor(tested_design, dimension),
        ),
        axis=2,
    )


def block_fixed_effect_tensors(nuisance_design, tested_design, path_index):
    """Build nested null and block-path alternative fixed-effect tensors."""
    nuisance_design = np.asarray(nuisance_design, dtype=float)
    tested_design = np.asarray(tested_design, dtype=float)
    if len(nuisance_design) != len(tested_design):
        raise ValueError("nuisance and tested designs do not align")
    dimension = len(path_index) - 1
    condition = unrestricted_tensor(nuisance_design, dimension)
    path_contrasts, orthogonal_nuisance = block_effect_bases(path_index)
    null_additions = []
    for column in range(tested_design.shape[1]):
        for contrast in range(orthogonal_nuisance.shape[1]):
            null_additions.append(
                tested_design[:, column, None]
                * orthogonal_nuisance[:, contrast][None, :]
            )
    null = (
        np.concatenate((condition, np.stack(null_additions, axis=2)), axis=2)
        if null_additions
        else condition
    )
    tested_additions = []
    for column in range(tested_design.shape[1]):
        for contrast in range(path_contrasts.shape[1]):
            tested_additions.append(
                tested_design[:, column, None]
                * path_contrasts[:, contrast][None, :]
            )
    alternative = np.concatenate(
        (null, np.stack(tested_additions, axis=2)), axis=2
    )
    return null, alternative, path_contrasts.shape[1] * tested_design.shape[1]


def nested_laplace_tests(
    null_data,
    alternative_data,
    *,
    family="multinomial",
    observation_noise=False,
    max_iter=200,
    mode_steps=30,
):
    """Fit nested tensor-design GLMMs and return LRT/BIC inference.

    The alternative tensor must begin with all null fixed-effect columns.  The
    likelihood-ratio test uses a chi-square reference distribution.  The Bayes
    factor is the BIC approximation based on the number of independent
    clusters, not a prior-integrated marginal likelihood.
    """
    if (
        null_data.fixed_effect_tensor is None
        or alternative_data.fixed_effect_tensor is None
    ):
        raise ValueError("nested tests require fixed-effect tensors")
    null_count = null_data.fixed_effect_tensor.shape[2]
    alternative_count = alternative_data.fixed_effect_tensor.shape[2]
    tested_count = alternative_count - null_count
    if tested_count < 1:
        raise ValueError(
            "the alternative needs at least one additional coefficient"
        )
    if not np.allclose(
        alternative_data.fixed_effect_tensor[:, :, :null_count],
        null_data.fixed_effect_tensor,
    ):
        raise ValueError("alternative fixed effects must begin with the null tensor")

    fit_options = {
        "family": family,
        "observation_noise": observation_noise,
        "max_iter": max_iter,
        "mode_steps": mode_steps,
    }
    null_fit = ec_glmm.fit_laplace(null_data, **fit_options)
    null_parameters = np.asarray(null_fit["parameters"])
    expanded = np.r_[
        null_parameters[:null_count],
        np.zeros(tested_count),
        null_parameters[null_count:],
    ]
    alternative_fit = ec_glmm.fit_laplace(
        alternative_data,
        initial=expanded,
        **fit_options,
    )

    lrt_statistic = max(
        0.0,
        2.0
        * (float(null_fit["objective"]) - float(alternative_fit["objective"])),
    )
    sample_size = len(np.unique(null_data.clusters))
    log_bayes_factor = 0.5 * (
        lrt_statistic - tested_count * np.log(sample_size)
    )
    return {
        "null_fit": null_fit,
        "alternative_fit": alternative_fit,
        "tested_count": tested_count,
        "lrt_statistic": lrt_statistic,
        "lrt_df": tested_count,
        "lrt_p_value": float(scipy.stats.chi2.sf(lrt_statistic, tested_count)),
        "bic_log_bayes_factor": float(log_bayes_factor),
        "bic_bayes_factor": float(
            np.exp(np.clip(log_bayes_factor, -700.0, 700.0))
        ),
        "independent_clusters": sample_size,
    }


def simulate_null_counts(data, fit, rng, *, family, observation_noise=False):
    """Simulate paired-primer EC counts from a fitted tensor-design null."""
    if data.fixed_effect_tensor is None:
        raise ValueError("null simulation requires a fixed-effect tensor")
    coefficients = np.asarray(fit["coefficients"], dtype=float).ravel()
    fixed = np.einsum("ndp,p->nd", data.fixed_effect_tensor, coefficients)
    _, cluster_index = np.unique(data.clusters, return_inverse=True)
    random_effect = rng.normal(
        0.0,
        float(fit["random_effect_sd"]),
        (
            cluster_index.max() + 1,
            data.n_isoforms - 1,
        ),
    )
    free_logits = fixed + random_effect[cluster_index]
    if observation_noise:
        free_logits += rng.normal(
            0.0,
            float(fit["observation_noise_sd"]),
            free_logits.shape,
        )
    logits = np.column_stack((free_logits, np.zeros(len(free_logits))))
    abundance = np.exp(logits - logits.max(axis=1, keepdims=True))
    simulated = []
    for observed, mapping in zip(data.counts, data.compatibility):
        if mapping.shape[0] == 0:
            simulated.append(np.zeros_like(observed, dtype=float))
            continue
        mass = abundance @ np.asarray(mapping, dtype=float).T
        probability = mass / mass.sum(axis=1, keepdims=True)
        totals = np.asarray(observed.sum(axis=1), dtype=int)
        rows = []
        for total, probabilities in zip(totals, probability):
            if family == "dirichlet_multinomial" and total > 0:
                probabilities = rng.dirichlet(
                    np.maximum(
                        float(fit["concentration"]) * probabilities, 1e-12
                    )
                )
            rows.append(rng.multinomial(int(total), probabilities))
        simulated.append(np.asarray(rows, dtype=float))
    return tuple(simulated)
