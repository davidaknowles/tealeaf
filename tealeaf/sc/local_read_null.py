"""Read-depth-pattern nulls, not a claim of EC or biological calibration."""

import numpy as np
from scipy.special import expit, logit


PATTERN_NULL_VERSION = 'local_marker_pooled_baseline_null_v1'


def local_read_pattern_null(counts, rng, slope_sd=0.):
    """Resimulate M x P x 2 types x 2 classes with original total depths.

    Fix each subject/primer baseline to inclusion pooled across its two types.
    Draw an independent mean-zero normal contrast shared across primers, then
    fresh binomials. Structural zero/one baselines remain exact. These arbitrary
    fixed baselines need not obey a fitted normal random-intercept model, making
    this a nuisance-misspecification stress check, not that model's own law.
    """
    observed = np.asarray(counts, dtype=float)
    if observed.ndim != 4 or observed.shape[2:] != (2, 2) or not np.isfinite(observed).all() or (observed < 0).any() or (observed > 2 ** 50).any() or not np.equal(observed, np.floor(observed)).all():
        raise ValueError('exact nonnegative M/P/two-type/two-class counts required')
    if not np.isfinite(slope_sd) or slope_sd < 0:
        raise ValueError('finite nonnegative slope standard deviation required')
    totals = observed.sum(axis=-1).astype(np.int64)
    pooled = totals.sum(axis=-1)
    baseline = np.divide(observed[..., 0].sum(axis=-1), pooled, out=np.zeros_like(pooled, dtype=float), where=pooled > 0)
    slope = rng.normal(0., slope_sd, observed.shape[0])
    probability = np.repeat(baseline[..., None], 2, axis=-1)
    interior = (baseline > 0) & (baseline < 1)
    # Do not replace structural boundaries with pseudocounts or epsilon floors.
    if interior.any():
        offset = slope[:, None, None] * np.array([-.5, .5])[None, None, :]
        probability[interior] = expit(logit(baseline[interior])[:, None] + np.broadcast_to(offset, probability.shape)[interior])
    included = rng.binomial(totals, probability)
    return np.stack([included, totals - included], axis=-1)
