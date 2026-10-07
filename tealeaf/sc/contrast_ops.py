"""Linear-time operations for the package's Helmert contrast convention."""

import numpy as np


def helmert_multiply(coordinates):
    """Apply H_n to last-axis coordinates of length n-1, without dense H."""
    values = np.asarray(coordinates, dtype=float)
    if values.ndim < 1:
        raise ValueError("last-axis contrast coordinates required")
    dimension = values.shape[-1]
    if dimension == 0:
        return np.zeros((*values.shape[:-1], 1))
    indices = np.arange(1, dimension + 1, dtype=float)
    scaled = values / np.sqrt(indices * (indices + 1))
    zeros = np.zeros((*values.shape[:-1], 1))
    tails = np.concatenate((np.cumsum(scaled[..., ::-1], axis=-1)[..., ::-1], zeros), axis=-1)
    return tails - np.concatenate((zeros, indices * scaled), axis=-1)


def helmert_transpose_multiply(values):
    """Apply H_n.T to last-axis values of length n, without dense H."""
    values = np.asarray(values, dtype=float)
    if values.ndim < 1 or values.shape[-1] < 1:
        raise ValueError("nonempty last-axis category values required")
    indices = np.arange(1, values.shape[-1], dtype=float)
    return (np.cumsum(values, axis=-1)[..., :-1] - indices * values[..., 1:]) / np.sqrt(indices * (indices + 1))
