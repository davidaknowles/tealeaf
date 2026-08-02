"""Fixed-effect constructions for splice-block EC-count GLMMs."""

from __future__ import annotations

import numpy as np

from . import differential


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
        (cluster_index.max() + 1, data.n_isoforms - 1),
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
