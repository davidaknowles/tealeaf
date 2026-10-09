"""Observed marker-usage reporting, separate from mixed-model testing."""

import numpy as np


def pooled_local_marker_effect(counts):
    """Equal-primer mean of type-B minus type-A pooled inclusion fractions.

    Preserve all declared subjects and zero strata when pooling keys. A primer
    needs a positive total in both types. This is a read-marker fraction, not
    absolute RNA PSI, and does not correct type-dependent capture or outside
    event sources. No pseudocount, effect cutoff, fitted p-value or LR input.
    """
    values = np.asarray(counts)
    if values.ndim != 4 or values.shape[2:] != (2, 2) or not np.isfinite(values).all() or (values < 0).any() or (values > 2 ** 50).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError('original exact nonnegative M/P/two-type/two-class counts required')
    pooled = values.sum(axis=0, dtype=float)
    totals = pooled.sum(axis=-1)
    fractions = np.divide(pooled[..., 0], totals, out=np.full(totals.shape, np.nan), where=totals > 0)
    usable = (totals > 0).all(axis=1)
    differences = fractions[:, 1] - fractions[:, 0]
    return dict(pooled_marker_effect=float(differences[usable].mean()) if usable.any() else np.nan, pooled_marker_n_primers=int(usable.sum()), pooled_marker_primer_effects=[float(value) if np.isfinite(value) else None for value in differences])
