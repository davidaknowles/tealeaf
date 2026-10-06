"""Subject-aware descriptive pooling of uncertain local-path proportions."""

import numpy as np
from scipy.optimize import minimize
from scipy.special import digamma, gammaln, softmax

from .differential import helmert_basis


def proportion_covariance(fit):
    """Propagate an S-1 dimensional ILR covariance to S proportions."""
    proportions = np.asarray(fit.path_proportions, dtype=float)
    jacobian = (np.diag(proportions) - np.outer(proportions, proportions)) @ helmert_basis(len(proportions))
    return jacobian @ fit.covariance.covariance @ jacobian.T


def dirichlet_pooling(proportions, covariances, labels, depths, *, max_iter=300):
    """Approximate Dirichlet random-effects reporting, not a count-based test.

    Proportions have shape N by S, conditional covariance N by S by S,
    and labels/depths length N. EC-derived estimates are not observed path
    counts. Moment matching replaces their uncertainty by an effective
    multinomial depth, capped at the actual EC depth. Fractional effective
    counts are fitted with a generalized Dirichlet-multinomial objective.
    One concentration is shared across the C label-specific population means.
    It limits high-depth subject leverage when biological heterogeneity is
    present. This is an approximation: anisotropic EC uncertainty, dependence
    between a subject's labels, and baseline uncertainty are not modeled.
    No biological p-value or posterior confidence interval is returned.
    """
    proportions = np.asarray(proportions, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    depths = np.asarray(depths, dtype=float)
    labels = np.asarray(labels)
    if proportions.ndim != 2 or proportions.shape[1] < 2:
        raise ValueError("proportions must have N rows and at least two paths")
    n, size = proportions.shape
    if covariances.shape != (n, size, size) or labels.shape != (n,) or depths.shape != (n,):
        raise ValueError("pooling inputs must align")
    if not np.isfinite(proportions).all() or np.any(proportions <= 0) or not np.allclose(proportions.sum(axis=1), 1):
        raise ValueError("finite positive simplex proportions required")
    if not np.isfinite(covariances).all() or np.any(depths <= 0) or not np.isfinite(depths).all():
        raise ValueError("finite covariance and positive depths required")
    if not np.allclose(covariances, covariances.transpose(0, 2, 1)) or np.min(np.linalg.eigvalsh(covariances)) < -1e-10:
        raise ValueError("symmetric positive-semidefinite covariance required")
    traces = np.trace(covariances, axis1=1, axis2=2)
    if np.any(traces <= 0):
        raise ValueError("positive conditional uncertainty required")
    effective_depth = np.minimum(depths, (1 - np.square(proportions).sum(axis=1)) / traces)
    levels, encoded = np.unique(labels, return_inverse=True)
    if any(np.sum(encoded == index) < 2 for index in range(len(levels))):
        raise ValueError("at least two subjects per label required")
    counts = effective_depth[:, None] * proportions
    initial = np.array([proportions[encoded == index].mean(axis=0) for index in range(len(levels))])
    logits = np.log(initial[:, :-1]) - np.log(initial[:, -1:])

    def objective(parameters):
        means = softmax(np.column_stack([parameters[:-1].reshape(len(levels), size - 1), np.zeros(len(levels))]), axis=1)
        concentration = np.exp(parameters[-1])
        alpha = concentration * means[encoded]
        score = digamma(alpha + counts) - digamma(alpha)
        likelihood = gammaln(concentration) - gammaln(concentration + effective_depth) + (gammaln(alpha + counts) - gammaln(alpha)).sum(axis=1)
        # Total uniform population-mean pseudocount 1, separate from subject fits.
        likelihood_sum = likelihood.sum() + np.log(means).sum() / size
        mean_gradient = np.zeros_like(means)
        np.add.at(mean_gradient, encoded, concentration * score)
        mean_gradient += 1 / (size * means)
        logit_gradient = means * (mean_gradient - (mean_gradient * means).sum(axis=1, keepdims=True))
        concentration_gradient = concentration * (digamma(concentration) - digamma(concentration + effective_depth) + (means[encoded] * score).sum(axis=1)).sum()
        return -float(likelihood_sum), -np.r_[logit_gradient[:, :-1].ravel(), concentration_gradient]

    fitted = minimize(objective, np.r_[logits.ravel(), np.log(20.)], jac=True, method="L-BFGS-B", bounds=[(-25, 25)] * logits.size + [(-5, 14)], options={"maxiter": max_iter, "ftol": 1e-10})
    means = softmax(np.column_stack([fitted.x[:-1].reshape(len(levels), size - 1), np.zeros(len(levels))]), axis=1)
    concentration = float(np.exp(fitted.x[-1]))
    weights = effective_depth * (concentration + 1) / (effective_depth + concentration)
    return {"means": means, "levels": levels, "concentration": concentration, "effective_depth": effective_depth, "subject_precision_weights": weights, "converged": bool(fitted.success), "objective": float(fitted.fun)}
