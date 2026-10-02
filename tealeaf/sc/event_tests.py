"""Paired signed-rank tests for arrays of event differences."""

from functools import lru_cache

import numpy as np
from scipy.special import ndtr
from scipy.stats import rankdata


def paired_signed_rank(differences, valid, *, exact=False, tie_method="average"):
    """Test M events by N subjects, dropping missing and zero differences.

    Exact tails condition on the observed absolute differences, including ties.
    Each nonzero subject receives an independent fair sign under the null.
    No continuity correction is used for the normal approximation.
    """
    differences = np.asarray(differences, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    if differences.ndim != 2 or valid.shape != differences.shape:
        raise ValueError("differences and valid must have matching event-by-subject shapes")
    if tie_method not in ("average", "ordinal"):
        raise ValueError("tie_method must be 'average' or 'ordinal'")
    valid = valid & np.isfinite(differences)
    nonzero = valid & (differences != 0)
    values = np.where(nonzero, np.abs(differences), np.inf)
    ranks = rankdata(values, axis=1, method=tie_method)
    ranks[~nonzero] = 0
    positive = np.sum(np.where(differences > 0, ranks, 0), axis=1)
    total = ranks.sum(axis=1)
    variance = np.square(ranks).sum(axis=1) / 4
    p = np.full(len(differences), np.nan)
    active = variance > 0
    p[~active & valid.any(axis=1)] = 1
    if not exact:
        p[active] = 2 * ndtr(-np.abs(positive[active] - total[active] / 2) / np.sqrt(variance[active]))
        return p
    if tie_method != "average":
        raise ValueError("exact signed-rank testing requires average ranks")
    for row in np.flatnonzero(active):
        # Average ranks are integers or half-integers; doubling gives a lattice
        # for the exact conditional sign distribution, even with arbitrary ties.
        weights = tuple(sorted(np.rint(2 * ranks[row, nonzero[row]]).astype(int)))
        tail = int(round(2 * min(positive[row], total[row] - positive[row])))
        p[row] = min(1.0, 2 * signed_rank_weight_cdf(weights)[tail])
    return p


@lru_cache(maxsize=256)
def signed_rank_weight_cdf(weights):
    """CDF of a sum of independently included positive integer weights.

    Normalize probabilities at every step rather than counting sign assignments
    in int64, which overflows for sufficiently many subjects.
    """
    probability = np.array([1.0])
    for weight in weights:
        updated = np.zeros(len(probability) + weight)
        updated[:len(probability)] += 0.5 * probability
        updated[weight:] += 0.5 * probability
        probability = updated
    return np.cumsum(probability)


def signed_rank_cdf(n):
    """Exact untied signed-rank CDF for n nonzero subject differences."""
    return signed_rank_weight_cdf(tuple(range(1, int(n) + 1)))
