"""Fixed-effect constructions for splice-block EC-count GLMMs."""

from __future__ import annotations

import numpy as np

from . import differential


def path_contrast_matrix(path_index):
    """Map orthonormal path contrasts to reference-isoform logits.

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
    reference = isoform_effects[-1]
    return isoform_effects[:-1] - reference[None, :]


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
    null = unrestricted_tensor(nuisance_design, dimension)
    path_contrasts = path_contrast_matrix(path_index)
    additions = []
    for column in range(tested_design.shape[1]):
        for contrast in range(path_contrasts.shape[1]):
            additions.append(
                tested_design[:, column, None]
                * path_contrasts[:, contrast][None, :]
            )
    alternative = np.concatenate(
        (null, np.stack(additions, axis=2)), axis=2
    )
    return null, alternative, path_contrasts.shape[1] * tested_design.shape[1]


def permute_within_subject(values, subjects, rng):
    """Permute labels independently within each paired biological subject."""
    values = np.asarray(values)
    subjects = np.asarray(subjects)
    result = values.copy()
    for subject in np.unique(subjects):
        positions = np.flatnonzero(subjects == subject)
        result[positions] = rng.permutation(result[positions])
    return result
