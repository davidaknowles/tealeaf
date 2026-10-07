"""Subject-aware descriptive pooling of uncertain local-path proportions."""

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import digamma, gammaln, softmax

from .differential import helmert_basis


def paired_reporting(proportions, covariances, depths):
    """Descriptive estimators using identical weights for both members of a pair.

    Inputs have shapes M x 2 x S, M x 2 x S x S and M x 2. The
    random-effects estimator matches each difference covariance to an
    isotropic variance in the S-1 dimensional simplex tangent space and
    estimates between-subject variance by restricted likelihood. It is a
    scalar-weight reporting approximation, not a statistical test. It omits
    shared-baseline measurement dependence and anisotropic uncertainty.
    No long-read outcomes or other subject folds enter any estimator.
    """
    proportions = np.asarray(proportions, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    depths = np.asarray(depths, dtype=float)
    if proportions.ndim != 3 or proportions.shape[1] != 2 or proportions.shape[2] < 2:
        raise ValueError("proportions must have shape M by 2 by S")
    m, _, s = proportions.shape
    if m < 2 or covariances.shape != (m, 2, s, s) or depths.shape != (m, 2):
        raise ValueError("at least two aligned subject pairs required")
    if not np.isfinite(proportions).all() or np.any(proportions <= 0) or not np.allclose(proportions.sum(axis=2), 1):
        raise ValueError("finite positive simplex proportions required")
    if not np.isfinite(depths).all() or np.any(depths <= 0):
        raise ValueError("positive finite paired depths required")
    differences = proportions[:, 1] - proportions[:, 0]
    harmonic = 1 / np.sum(1 / depths, axis=1)
    harmonic /= harmonic.sum()
    equal = np.full(m, 1 / m)
    log_proportions = np.log(proportions)
    geometric = softmax(log_proportions.mean(axis=0), axis=1)
    # Median CLR differences are path-permutation equivariant, unlike taking
    # coordinate medians in an arbitrary ILR basis.
    clr = log_proportions - log_proportions.mean(axis=2, keepdims=True)
    median = np.median(clr[:, 1] - clr[:, 0], axis=0)
    anchor = clr.mean(axis=(0, 1))
    median_means = softmax(np.stack([anchor - median / 2, anchor + median / 2]), axis=1)
    reports = {
        "subject arithmetic mean": {"effect": differences.mean(axis=0), "weights": equal},
        "subject geometric mean": {"effect": geometric[1] - geometric[0], "weights": equal},
        "paired median CLR": {"effect": median_means[1] - median_means[0], "weights": equal},
        "paired harmonic depth": {"effect": harmonic @ differences, "weights": harmonic},
    }
    try:
        if not np.isfinite(covariances).all() or not np.allclose(covariances, covariances.swapaxes(-1, -2)) or np.min(np.linalg.eigvalsh(covariances)) < -1e-10:
            raise ValueError("finite symmetric positive-semidefinite covariance required")
        variances = np.trace(covariances.sum(axis=1), axis1=1, axis2=2) / (s - 1)
        if np.any(variances <= 0):
            raise ValueError("positive paired measurement uncertainty required")
        tangent = differences @ helmert_basis(s)

        def reml(tau):
            total = variances + tau
            precision = 1 / total
            mean = precision @ tangent / precision.sum()
            return (s - 1) * (np.log(total).sum() + np.log(precision.sum())) + (np.square(tangent - mean) * precision[:, None]).sum()

        upper = max(float(np.square(tangent - tangent.mean(axis=0)).sum() / ((m - 1) * (s - 1))), float(variances.max()), 1e-6) * 10
        optimum = minimize_scalar(reml, bounds=(0., upper), method="bounded", options={"xatol": 1e-10})
        if not optimum.success:
            raise ValueError("paired reporting REML failed")
        tau = float(optimum.x) if reml(optimum.x) < reml(0.) else 0.
        weights = 1 / (variances + tau)
        weights /= weights.sum()
        reports["paired random effects"] = {"effect": weights @ differences, "weights": weights, "between_subject_variance": tau, "fallback": False}
    except (ValueError, np.linalg.LinAlgError) as exception:
        reports["paired random effects"] = {"effect": differences.mean(axis=0), "weights": equal, "between_subject_variance": np.nan, "fallback": True, "error": str(exception)}
    return reports


def proportion_covariance(fit):
    """Propagate an S-1 dimensional ILR covariance to S proportions."""
    proportions = np.asarray(fit.path_proportions, dtype=float)
    if not fit.covariance.identifiable or not np.isfinite(fit.covariance.covariance).all():
        return np.full((len(proportions), len(proportions)), np.nan)
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
