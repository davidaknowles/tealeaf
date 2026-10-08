"""Checks of categorical EC-mixture assumptions, not differential tests."""

import numpy as np
from scipy import special

from .ec_glmm import ECGLMMData
from .ec_block_glmm import pooled_isoform_weights


def pooled_ec_mixture_diagnostics(counts, compatibility, *, max_iter=300):
    """Best-fit pooled single-primer mixture and a conservative KL certificate.

    Counts are N_observations by K or length K, mapping is K by T. Each
    normalized mapping column is a categorical probability vector. Pooling
    arbitrary observation-specific mixtures remains in their convex hull.
    For convex KL(w), its tangent gives min KL >= KL(w) - dual_gap, even
    when optimization is unfinished. Integral independent categorical draws
    additionally permit a conservative method-of-types probability bound.
    These bounds diagnose the assumed fixed-map count model, never splicing
    itself. Estimated maps and data-dependent category filtering require an
    additional conditioning/selection analysis for a formal sampling claim.
    """
    counts, mapping = np.asarray(counts, dtype=float), np.asarray(compatibility, dtype=float)
    if counts.ndim not in (1, 2) or mapping.ndim != 2 or mapping.shape[0] != counts.shape[-1] or not np.isfinite(counts).all() or not np.isfinite(mapping).all() or (counts < 0).any() or (mapping < 0).any():
        raise ValueError("aligned finite nonnegative counts and K-by-T mapping required")
    integral = bool(np.allclose(counts, np.rint(counts), rtol=0., atol=1e-8))
    if counts.ndim == 2:
        counts = counts.sum(axis=0)
    total = float(counts.sum())
    result = dict(molecules=total, n_ecs=len(counts), n_transcripts=mapping.shape[1], integral_counts=integral)
    if total <= 0:
        return {**result, "positive_counts": False, "fit_converged": False, "KL": np.nan, "KL_lower_bound": np.nan, "KL_dual_gap": np.nan, "deviance": np.nan, "deviance_lower_bound": np.nan, "log_probability_bound": np.nan}
    columns = mapping.sum(axis=0)
    active = columns > 0
    if not active.any() or ((mapping[:, active].sum(axis=1) <= 0) & (counts > 0)).any():
        raise ValueError("observed ECs lack any represented transcript")
    mapping, columns = mapping[:, active], columns[active]
    categorical = mapping / columns
    if mapping.shape[1] == 1:
        theta, converged = np.ones(1), True
    else:
        data = ECGLMMData((counts[None, :],), (mapping,), np.ones((1, 1)), np.zeros(1))
        theta, converged = pooled_isoform_weights(data, max_iter=max_iter, return_status=True)
    weights = columns * theta
    weights /= weights.sum()
    probability = categorical @ weights
    positive = counts > 0
    if (probability[positive] <= 0).any():
        # A feasible interior point still gives a valid lower-bound certificate.
        weights = np.full(len(columns), 1 / len(columns))
        probability = categorical @ weights
        converged = False
    empirical = counts / total
    KL = max(float(np.sum(empirical[positive] * np.log(empirical[positive] / probability[positive]))), 0.)
    gradient = -categorical[positive].T @ (empirical[positive] / probability[positive])
    gap = max(float(gradient @ weights - gradient.min()), 0.)
    lower = max(KL - gap - 1e-12, 0.) if np.isfinite(gap) else 0.
    log_bound = np.nan
    if result["integral_counts"]:
        molecules = int(round(total))
        # Number of possible K-category count vectors, rather than a Wilks df.
        log_types = special.gammaln(molecules + len(counts)) - special.gammaln(molecules + 1) - special.gammaln(len(counts))
        log_bound = min(0., float(log_types - total * lower))
    return {**result, "positive_counts": True, "fit_converged": converged, "KL": KL, "KL_lower_bound": lower, "KL_dual_gap": gap, "deviance": 2 * total * KL, "deviance_lower_bound": 2 * total * lower, "log_probability_bound": log_bound}
