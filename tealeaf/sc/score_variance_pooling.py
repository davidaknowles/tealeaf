"""Experimental global biological variance for ragged binary EC score panels."""

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats, special


@dataclass
class ScalarScorePanel:
    """K tests and N retained subject scores, with no padded subject arrays."""

    values: np.ndarray
    log_measurement: np.ndarray
    biological_shape: np.ndarray
    groups: np.ndarray
    source_positions: np.ndarray
    n_subjects: np.ndarray
    reference_objective_offset: np.ndarray

    def evaluate(self, variance, *, values=None):
        """Same common-mean Gaussian model and modified KH F, fixed variance.

        Optional values is an N-vector of pseudo-effects (e.g. signed values),
        not raw efficient scores. Precision shares use per-test log scaling;
        residual sums use direct centered residuals to avoid cancellation.
        """
        if not np.isfinite(variance) or variance < 0:
            raise ValueError("finite nonnegative biological variance required")
        values = self.values if values is None else np.asarray(values, dtype=float)
        if values.shape != self.values.shape or not np.isfinite(values).all():
            raise ValueError("finite aligned retained pseudo-effects required")
        log_variance = self.log_measurement.copy()
        if variance > 0:
            positive = self.biological_shape > 0
            log_variance[positive] = np.logaddexp(log_variance[positive], np.log(variance) + np.log(self.biological_shape[positive]))
        size = len(self.n_subjects)
        log_scale = np.full(size, -np.inf)
        np.maximum.at(log_scale, self.groups, -log_variance)
        relative_weights = np.exp(-log_variance - log_scale[self.groups])
        relative_precision = np.bincount(self.groups, weights=relative_weights, minlength=size)
        shares = relative_weights / relative_precision[self.groups]
        means = np.bincount(self.groups, weights=shares * values, minlength=size)
        log_precision = log_scale + np.log(relative_precision)
        centered = (values - means[self.groups]) * np.exp(-.5 * log_variance)
        residual = np.bincount(self.groups, weights=centered**2, minlength=size)
        inflation = np.maximum(1., residual / (self.n_subjects - 1))
        standardized = means * np.exp(.5 * log_precision)
        statistic = standardized**2 / inflation
        objective = np.bincount(self.groups, weights=log_variance, minlength=size) + log_precision + residual + self.reference_objective_offset
        maximum_share = np.zeros(size)
        np.maximum.at(maximum_share, self.groups, shares)
        effective = 1 / np.bincount(self.groups, weights=shares**2, minlength=size)
        if any(not np.isfinite(value).all() for value in (means, statistic, objective, maximum_share, effective)):
            raise ValueError("score panel evaluation exceeds finite range")
        return dict(p_value=stats.f.sf(statistic, 1, self.n_subjects - 1), statistic=statistic, mean_difference=means, residual_inflation=inflation, restricted_objective=objective, maximum_precision_share=maximum_share, effective_weighted_subjects=effective, n_subjects=self.n_subjects.copy())


def scalar_score_panel(scores, information, biological_shape, reference_information, groups):
    """Prepare aligned N-vectors using the existing relative information rank.

    Integer groups cover 0,...,K-1; each retained group has at least four
    subjects. source_positions preserves the original rank mask for signs.
    Source arrays are copied/selected, never edited. No counts are refitted.
    """
    scores, information, shape, reference = (np.asarray(value, dtype=float) for value in (scores, information, biological_shape, reference_information))
    groups = np.asarray(groups)
    if scores.ndim != 1 or not len(scores) or any(value.shape != scores.shape or not np.isfinite(value).all() for value in (scores, information, shape, reference)) or groups.shape != scores.shape or not np.issubdtype(groups.dtype, np.integer) or not np.array_equal(np.unique(groups), np.arange(len(np.unique(groups)))):
        raise ValueError("finite aligned scalar score arrays and consecutive integer groups required")
    if (information < 0).any() or (shape < 0).any() or (reference < 0).any():
        raise ValueError("nonnegative information, biological shape and reference required")
    supported = reference > np.maximum(reference, np.finfo(float).tiny) * 1e-10
    relative = information[supported] / reference[supported]
    if (relative > 1 + 1e-8).any():
        raise ValueError("profiled information exceeds its reference")
    keep = supported.copy()
    keep[supported] = relative > np.maximum(relative, 1.) * 1e-10
    n_subjects = np.bincount(groups[keep], minlength=len(np.unique(groups)))
    if (n_subjects < 4).any():
        raise ValueError("four informative subjects required in every panel test")
    values = scores[keep] / information[keep]
    if not np.isfinite(values).all():
        raise ValueError("retained pseudo-effects exceed finite range")
    retained_groups = groups[keep]
    log_reference = np.log(reference[keep])
    maximum_reference = np.full(len(n_subjects), -np.inf)
    np.maximum.at(maximum_reference, retained_groups, log_reference)
    log_reference_sum = maximum_reference + np.log(np.bincount(retained_groups, weights=np.exp(log_reference - maximum_reference[retained_groups]), minlength=len(n_subjects)))
    # Match the original reference-whitened REML objective, including its
    # variance-independent sum(log R)-log(sum R) coordinate constant.
    offset = np.bincount(retained_groups, weights=log_reference, minlength=len(n_subjects)) - log_reference_sum
    return ScalarScorePanel(values, -np.log(information[keep]), shape[keep].copy(), retained_groups.copy(), np.flatnonzero(keep), n_subjects, offset)


def fit_shared_biological_variance(panel, *, values=None):
    """One global REML multiplier, with one free common mean per panel test.

    This is variance pooling, not a target smoothing prior or independent
    replicated tests. Dataset callers must freeze a gene-balanced training
    panel, report unavailable selected cases, and refit within null families.
    Both the finite-variance optimum and exact zero boundary are evaluated.
    """
    positive = panel.biological_shape > 0
    if not positive.any():
        raise ValueError("biological variance is unidentified when every shape is zero")
    log_scale = float(np.median(panel.log_measurement[positive] - np.log(panel.biological_shape[positive])))
    lower = max(log_scale - 24., np.log(np.finfo(float).tiny))
    upper = min(log_scale + 16., np.log(np.finfo(float).max) - 1.)
    if lower >= upper:
        raise ValueError("no finite biological variance search interval")
    def objective(variance):
        return float(panel.evaluate(variance, values=values)["restricted_objective"].mean())
    result = optimize.minimize_scalar(lambda value: objective(np.exp(value)), bounds=(lower, upper), method="bounded", options={"xatol": 1e-6})
    if not result.success or not np.isfinite(result.fun):
        raise ValueError("pooled biological-variance REML failed")
    variance = min((0., float(np.exp(result.x))), key=objective)
    return dict(biological_variance=variance, mean_restricted_objective=objective(variance), zero_variance_objective=objective(0.), training_tests=len(panel.n_subjects), informative_training_subjects=len(panel.values), log_search_lower=lower, log_search_upper=upper, near_search_boundary=bool(variance > 0 and min(result.x - lower, upper - result.x) < 1e-3), optimizer_evaluations=int(result.nfev), interpretation="one pooled variance, unrestricted common mean per frozen training test; estimated covariance and Gaussian EC-score approximation")


def scalar_precision_share_bounds(information, biological_shape):
    """Rigorous subject share bounds over every variance in [0,infinity].

    For M positive information/shape values, w_i(t)=I_i/(1+t*I_i*B_i).
    The ratio w_j/w_i lies between I_j/I_i and B_i/B_j, so subject i's
    share is bounded by reciprocal sums of their elementwise maxima/minima.
    Returns two M-vectors, lower and upper bounds; these are not claims that
    either bound is attainable, nor a fitted or replacement statistical test.
    """
    information, shape = (np.asarray(value, dtype=float) for value in (information, biological_shape))
    if information.ndim != 1 or not len(information) or shape.shape != information.shape or any(not np.isfinite(value).all() or (value <= 0).any() for value in (information, shape)):
        raise ValueError("aligned finite positive information and shape required")
    log_info, log_shape = np.log(information), np.log(shape)
    zero_ratios = log_info[None, :] - log_info[:, None]
    infinite_ratios = log_shape[:, None] - log_shape[None, :]
    lower = np.exp(-special.logsumexp(np.maximum(zero_ratios, infinite_ratios), axis=1))
    upper = np.exp(-special.logsumexp(np.minimum(zero_ratios, infinite_ratios), axis=1))
    return lower, upper
