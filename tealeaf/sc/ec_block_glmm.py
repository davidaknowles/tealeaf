"""Fixed-effect constructions for splice-block EC-count GLMMs."""

from __future__ import annotations

import numpy as np
import scipy.optimize
import scipy.stats

from . import differential, ec_glmm


def pooled_isoform_weights(counts, compatibility, max_iter=200):
    """Estimate label-independent isoform weights from pooled EC counts."""
    counts = tuple(np.asarray(value, dtype=float).sum(axis=0) for value in counts)
    mappings = tuple(np.asarray(value, dtype=float) for value in compatibility)
    n_isoforms = mappings[0].shape[1]
    if n_isoforms < 2:
        return np.ones(n_isoforms, dtype=float)

    def objective(free_logits):
        logits = np.r_[free_logits, 0.0]
        abundance = np.exp(logits - np.max(logits))
        value = 0.0
        gradient = np.zeros(n_isoforms, dtype=float)
        scale = 0.0
        for observed, mapping in zip(counts, mappings):
            total = float(observed.sum())
            if total <= 0:
                continue
            mass = np.maximum(mapping @ abundance, 1e-300)
            total_mass = float(mass.sum())
            value += float(observed @ np.log(mass)) - total * np.log(total_mass)
            gradient += abundance * (
                mapping.T @ (observed / mass)
                - total * mapping.sum(axis=0) / total_mass
            )
            scale += total
        scale = max(scale, 1.0)
        return -value / scale, -gradient[:-1] / scale

    result = scipy.optimize.minimize(
        objective,
        np.zeros(n_isoforms - 1, dtype=float),
        method="L-BFGS-B",
        jac=True,
        bounds=[(-20.0, 20.0)] * (n_isoforms - 1),
        options={"maxiter": int(max_iter), "ftol": 1e-12, "gtol": 1e-7},
    )
    logits = np.r_[result.x, 0.0]
    weights = np.exp(logits - np.max(logits))
    return weights / weights.sum()


def collapse_within_paths(compatibility, weights, path_index):
    """Collapse represented isoforms within paths using fixed pooled weights.

    Each isoform outside the represented block remains its own nuisance class.
    This reduces within-path nuisance dimensions without merging unspliced or
    otherwise unrepresented isoforms.
    """
    path_index = np.asarray(path_index, dtype=int)
    weights = np.asarray(weights, dtype=float)
    if len(path_index) != len(weights):
        raise ValueError("path assignments and pooled weights do not align")
    represented = [
        ("path", int(path))
        for path in sorted(set(path_index[path_index >= 0]))
    ]
    nuisance = [
        ("nuisance", int(index)) for index in np.flatnonzero(path_index < 0)
    ]
    classes = represented + nuisance
    projection = np.zeros((len(weights), len(classes)), dtype=float)
    collapsed_paths = []
    for column, (kind, value) in enumerate(classes):
        members = (
            np.flatnonzero(path_index == value)
            if kind == "path"
            else np.asarray([value], dtype=int)
        )
        local = weights[members]
        local = (
            local / local.sum()
            if local.sum() > 0
            else np.full(len(members), 1 / len(members))
        )
        projection[members, column] = local
        collapsed_paths.append(value if kind == "path" else -1)
    collapsed = tuple(
        np.asarray(mapping, dtype=float) @ projection for mapping in compatibility
    )
    return collapsed, np.asarray(collapsed_paths, dtype=int), projection


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


def simulate_null_counts(
    data, fit, rng, *, family, observation_noise=False, random_slopes=False
):
    """Simulate paired-primer EC counts from a fitted tensor-design null."""
    if data.fixed_effect_tensor is None:
        raise ValueError("null simulation requires a fixed-effect tensor")
    coefficients = np.asarray(fit["coefficients"], dtype=float).ravel()
    fixed = np.einsum("ndp,p->nd", data.fixed_effect_tensor, coefficients)
    _, cluster_index = np.unique(data.clusters, return_inverse=True)
    random_design = (
        np.asarray(data.random_effect_design, dtype=float)
        if random_slopes
        else np.ones((len(data.clusters), 1), dtype=float)
    )
    random_sds = np.asarray(
        fit.get("random_effect_sds", [fit["random_effect_sd"]]), dtype=float
    )
    if len(random_sds) != random_design.shape[1]:
        raise ValueError("fitted random terms do not match random-effect design")
    random_effect = rng.normal(
        0.0,
        random_sds[None, :, None],
        (
            cluster_index.max() + 1,
            random_design.shape[1],
            data.n_isoforms - 1,
        ),
    )
    free_logits = fixed + np.einsum(
        "nr,nrd->nd", random_design, random_effect[cluster_index]
    )
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
