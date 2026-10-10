"""Fast evaluation of the binary-path EC random-subject model.

Same model as ``path_marginal``: for subject u and cell type c the path-1
proportion psi_uc ~ Beta(kappa mu_uc, kappa (1 - mu_uc)), logit mu_uc =
x_uc' beta + b_u, b_u ~ N(0, sigma^2), and both primer EC count vectors are
multinomial with probabilities affine in psi_uc. Only the numerics differ.

Each row's EC log-likelihood l_n(t) depends on the scalar t = logit psi, so it
is computed once on a fixed t-grid. The Beta integral is then a trapezoid sum
in t, where the integrand psi^a (1 - psi)^c L_n is smooth and trapezoid
quadrature converges spectrally; |t| > LIMIT tails are added analytically with
L_n held at its boundary value. The subject integral uses adaptive
Gauss-Hermite quadrature centred on each subject's posterior mode of b_u,
warm-started from the previous evaluation. Gradients are the quadrature of the
analytic score, from the same posterior weights.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, digamma, expit, log_expit, polygamma, roots_hermite
from scipy.stats import chi2

from .ec_block_glmm import collapse_isoforms_to_paths
from .path_marginal import BinaryECPathLikelihood

LIMIT = 18.
KAPPA_BOUNDS = (.1, 2000.)
SIGMA_BOUNDS = (.01, 3.)


def logsumexp(values, axis=None, keepdims=False):
    """Plain numpy log-sum-exp; scipy's version has large per-call overhead."""
    peak = np.max(values, axis=axis, keepdims=True)
    peak = np.where(np.isfinite(peak), peak, 0.)
    total = np.log(np.sum(np.exp(values - peak), axis=axis, keepdims=True)) + peak
    return total if keepdims else np.squeeze(total, axis=axis)


@dataclass
class RowGrid:
    """Per-row trapezoid grid in t = logit psi, or exact Beta-binomial terms."""
    t: np.ndarray = None
    log_weight: np.ndarray = None
    ell: np.ndarray = None
    exact: tuple = None
    shift: float = 0.

    @property
    def log_p(self):
        return log_expit(self.t)

    @property
    def log_q(self):
        return log_expit(-self.t)


def row_grids(likelihood, *, refine=1., coarse_step=.05):
    """Grid each row's EC log-likelihood, spacing from its own curvature.

    Spacing is min(.03, sd/1.5) / refine, where sd is the Laplace width of the
    row likelihood in t; .03 resolves the narrowest Beta prior allowed by
    KAPPA_BOUNDS. Exact pure-path rows keep the Beta-binomial identity.
    Each grid is shifted by the row's coarse-grid maximum, which does not
    depend on refine or on any model parameter, so it cancels from likelihood
    ratios and from refined-grid accuracy checks.
    """
    grids = []
    coarse = np.arange(-LIMIT, LIMIT + coarse_step / 2, coarse_step)
    for row in range(len(likelihood.subjects)):
        terms = likelihood.binomial_terms[row]
        if np.isfinite(terms).all():
            grids.append(RowGrid(exact=tuple(terms)))
            continue
        ell = likelihood.row_log_likelihood(row, expit(coarse))
        peak = np.argmax(ell)
        width = coarse[ell >= ell[peak] - 2.]
        sd = max((width.max() - width.min()) / 4., coarse_step / 4.)
        step = min(.03, sd / 1.5) / refine
        t = np.linspace(-LIMIT, LIMIT, int(np.ceil(2 * LIMIT / step)) + 1)
        shift = float(ell[peak])
        ell = likelihood.row_log_likelihood(row, expit(t)) - shift
        weight = np.full(len(t), t[1] - t[0])
        weight[[0, -1]] /= 2
        grids.append(RowGrid(t=t, log_weight=np.log(weight), ell=ell, shift=shift))
    return grids


def row_integral(grid, means, kappa):
    """log E_Beta[L] and derivatives at mean logits m (any shape).

    Returns the value, its first and second derivatives in m, and kappa times
    its derivative in kappa. With a = kappa s and c = kappa (1 - s), the
    a/c derivatives are posterior means and covariances of the grid scores
    u = log psi, v = log(1 - psi) (analytic tail scores beyond LIMIT) minus
    the Beta normalizer's digamma and trigamma terms.
    """
    s, r = expit(means), expit(-means)
    a, c = np.maximum(kappa * s, 1e-300), np.maximum(kappa * r, 1e-300)
    trigamma_total = polygamma(1, a + c)
    if grid.exact is not None:
        first, second, constant = grid.exact
        value = betaln(a + first, c + second) - betaln(a, c) + constant
        total = digamma(a + c + first + second)
        da = digamma(a + first) - total - digamma(a) + digamma(a + c)
        dc = digamma(c + second) - total - digamma(c) + digamma(a + c)
        both = -polygamma(1, a + c + first + second) + trigamma_total
        daa = polygamma(1, a + first) - polygamma(1, a) + both
        dcc = polygamma(1, c + second) - polygamma(1, c) + both
        dac = both
    else:
        log_p, log_q = grid.log_p, grid.log_q
        flat_a, flat_c = a.reshape(-1, 1), c.reshape(-1, 1)
        interior = grid.ell + grid.log_weight + flat_a * log_p + flat_c * log_q
        left = grid.ell[0] - LIMIT * flat_a - np.log(flat_a)
        right = grid.ell[-1] - LIMIT * flat_c - np.log(flat_c)
        joint = np.concatenate([interior, left, right], axis=1)
        total = logsumexp(joint, axis=1, keepdims=True)
        posterior = np.exp(joint - total)
        interior_posterior, left_posterior, right_posterior = posterior[:, :-2], posterior[:, -2:-1], posterior[:, -1:]
        mean_u = interior_posterior @ log_p + (left_posterior * (-LIMIT - 1 / flat_a)).ravel()
        mean_v = interior_posterior @ log_q + (right_posterior * (-LIMIT - 1 / flat_c)).ravel()
        square_u = interior_posterior @ np.square(log_p) + (left_posterior * np.square(-LIMIT - 1 / flat_a)).ravel()
        square_v = interior_posterior @ np.square(log_q) + (right_posterior * np.square(-LIMIT - 1 / flat_c)).ravel()
        cross = interior_posterior @ (log_p * log_q)
        normalizer = digamma(a + c).ravel()
        value = total.reshape(a.shape) - betaln(a, c)
        da = (mean_u - digamma(a).ravel() + normalizer).reshape(a.shape)
        dc = (mean_v - digamma(c).ravel() + normalizer).reshape(a.shape)
        daa = (square_u - np.square(mean_u) + (left_posterior / np.square(flat_a)).ravel()).reshape(a.shape) - polygamma(1, a) + trigamma_total
        dcc = (square_v - np.square(mean_v) + (right_posterior / np.square(flat_c)).ravel()).reshape(a.shape) - polygamma(1, c) + trigamma_total
        dac = (cross - mean_u * mean_v).reshape(a.shape) + trigamma_total
    scale = kappa * s * r
    slope = scale * (da - dc)
    curvature = scale * (r - s) * (da - dc) + np.square(scale) * (daa - 2 * dac + dcc)
    return value, slope, curvature, a * da + c * dc


@dataclass
class GridModel:
    """Binary EC random-subject model on precomputed row grids.

    design is N x C (rows in likelihood order), subjects length N. Mode
    caches hold each subject's last posterior mode of b_u for warm starts.
    """
    grids: list
    design: np.ndarray
    subjects: np.ndarray
    nodes: int = 9
    modes: dict = field(default_factory=dict)

    def __post_init__(self):
        self.design = np.asarray(self.design, dtype=float)
        self.groups = [np.flatnonzero(self.subjects == subject) for subject in np.unique(self.subjects)]
        self.hermite = roots_hermite(self.nodes)

    def _subject_terms(self, rows, eta, kappa, sigma, b):
        """f_u(b), its first two b-derivatives and row scores, b any 1-D array."""
        value = -np.square(b) / (2 * sigma ** 2) - np.log(sigma) - .5 * np.log(2 * np.pi)
        slope = -b / sigma ** 2
        curvature = np.full_like(b, -1 / sigma ** 2)
        scores = []
        for row in rows:
            integral, d_mean, d2_mean, d_kappa = row_integral(self.grids[row], eta[row] + b, kappa)
            value, slope, curvature = value + integral, slope + d_mean, curvature + d2_mean
            scores.append((d_mean, d_kappa))
        return value, slope, curvature, scores

    def _mode(self, index, rows, eta, kappa, sigma):
        """Damped Newton for the posterior mode of b_u, warm-started."""
        b = self.modes.get(index, 0.)
        value, slope, curvature, _ = self._subject_terms(rows, eta, kappa, sigma, np.array([b]))
        for _ in range(100):
            step = -slope[0] / curvature[0] if curvature[0] < 0 else slope[0] * sigma ** 2
            for _ in range(40):
                trial = np.array([b + step])
                new = self._subject_terms(rows, eta, kappa, sigma, trial)
                if new[0][0] >= value[0] - 1e-12:
                    break
                step /= 2
            b = trial[0]
            value, slope, curvature, _ = new
            if abs(step) < 1e-9 or abs(slope[0]) < 1e-9:
                break
        self.modes[index] = b
        return b, 1 / np.sqrt(-curvature[0]) if curvature[0] < 0 else sigma

    def objective(self, parameters):
        """Negative marginal log-likelihood and its gradient."""
        coefficients, kappa, sigma = parameters[:-2], np.exp(parameters[-2]), np.exp(parameters[-1])
        eta = self.design @ coefficients
        roots, weights = self.hermite
        total, gradient = 0., np.zeros_like(parameters)
        for index, rows in enumerate(self.groups):
            mode, scale = self._mode(index, rows, eta, kappa, sigma)
            points = mode + np.sqrt(2) * scale * roots
            value, _, _, scores = self._subject_terms(rows, eta, kappa, sigma, points)
            joint = np.log(weights) + np.square(roots) + value
            evidence = logsumexp(joint)
            total += evidence + np.log(np.sqrt(2) * scale)
            posterior = np.exp(joint - evidence)
            for row, (d_mean, d_kappa) in zip(rows, scores):
                gradient[:-2] += self.design[row] * (posterior @ d_mean)
                gradient[-2] += posterior @ d_kappa
            gradient[-1] += posterior @ (np.square(points) / sigma ** 2 - 1)
        return -total, -gradient


def prepare_grid_likelihood(data, path_index, labels, subjects, baseline):
    """Subject-by-type EC aggregates with fixed outside and within-path shares.

    Same construction as ``prepare_binary_ec_likelihood`` without its
    per-row proposal fits, which the grid integrator does not need.
    """
    path_index = np.asarray(path_index, dtype=int)
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if not np.array_equal(np.unique(path_index[path_index >= 0]), [0, 1]):
        raise ValueError("binary paths 0 and 1 required")
    baseline = np.maximum(np.asarray(baseline, dtype=float), 1e-12)
    baseline /= baseline.sum()
    collapsed, _ = collapse_isoforms_to_paths(data, path_index, baseline)
    mass = baseline[path_index >= 0].sum()
    components = tuple(np.column_stack([mapping[:, 2:] @ baseline[path_index < 0], mass * mapping[:, 0], mass * mapping[:, 1]]) for mapping in collapsed.compatibility)
    keys, counts = [], []
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            selected = (subjects == subject) & (labels == label)
            observed = tuple(values[selected].sum(axis=0) for values in data.counts)
            if sum(values.sum() for values in observed) > 0:
                keys.append((subject, label))
                counts.append(observed)
    if not counts:
        raise ValueError("at least one positive-count aggregate required")
    aggregated = tuple(np.stack([values[primer] for values in counts]) for primer in range(len(data.counts)))
    n = len(keys)
    return BinaryECPathLikelihood(aggregated, components, np.asarray([key[0] for key in keys]), np.asarray([key[1] for key in keys]), np.ones((n, 2)), np.zeros(n))


def _initial(grids, design):
    """Least-squares logits from each row's grid maximum, clipped to +/-4."""
    peaks = []
    for grid in grids:
        if grid.exact is not None:
            first, second, _ = grid.exact
            peaks.append(np.log((first + .5) / (second + .5)))
        else:
            peaks.append(grid.t[np.argmax(grid.ell)])
    peaks = np.clip(peaks, -4., 4.)
    return np.linalg.lstsq(design, peaks, rcond=None)[0]


def _fit(model, initial, max_iter):
    bounds = [(-12., 12.)] * (len(initial) - 2) + [tuple(np.log(KAPPA_BOUNDS)), tuple(np.log(SIGMA_BOUNDS))]
    return minimize(model.objective, initial, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": max_iter, "ftol": 1e-10, "gtol": 1e-6})


def binary_grid_test(likelihood, *, nodes=9, max_iter=200, tolerance=1e-3, refine=1., retries=1):
    """Random-subject binary EC likelihood-ratio test on two or more types.

    Null and alternative refit kappa and sigma. Accuracy is checked by
    re-evaluating both optima with halved grid spacing and 2 * nodes + 1
    subject nodes. A fit failing only that check is refitted at the finer
    setting (up to retries times); remaining failures return p = 1 with
    converged False.
    """
    labels, subjects = likelihood.labels, likelihood.subjects
    levels = np.unique(labels)
    if len(levels) < 2:
        raise ValueError("two type levels required")
    repeated = sum(len(np.unique(labels[subjects == subject])) >= 2 for subject in np.unique(subjects))
    if repeated < 4:
        raise ValueError("at least four subjects with repeated types required")
    design = np.column_stack([np.ones(len(labels)), *[labels == level for level in levels[1:]]]).astype(float)
    grids = row_grids(likelihood, refine=refine)
    start = _initial(grids, design)
    alternative_model = GridModel(grids, design, subjects, nodes)
    null_model = GridModel(grids, design[:, :1], subjects, nodes)
    null = _fit(null_model, np.r_[_initial(grids, design[:, :1]), np.log(20.), np.log(.5)], max_iter)
    alternative = _fit(alternative_model, np.r_[null.x[:1], start[1:], null.x[-2:]], max_iter)
    if alternative.fun > null.fun:
        alternative = _fit(alternative_model, np.r_[null.x[:1], np.zeros(len(levels) - 1), null.x[-2:]], max_iter)
    fine_grids = row_grids(likelihood, refine=2 * refine)
    fine_null = GridModel(fine_grids, design[:, :1], subjects, 2 * nodes + 1).objective(null.x)[0]
    fine_alternative = GridModel(fine_grids, design, subjects, 2 * nodes + 1).objective(alternative.x)[0]
    statistic, fine_statistic = 2 * (null.fun - alternative.fun), 2 * (fine_null - fine_alternative)
    error = max(abs(fine_null - null.fun), abs(fine_alternative - alternative.fun), abs(statistic - fine_statistic))
    optimized = bool(null.success and alternative.success and np.isfinite([fine_null, fine_alternative]).all() and min(statistic, fine_statistic) >= -1e-6)
    if optimized and error > tolerance and retries > 0:
        return binary_grid_test(likelihood, nodes=2 * nodes + 1, max_iter=max_iter, tolerance=tolerance, refine=2 * refine, retries=retries - 1)
    converged = optimized and error <= tolerance
    statistic = max(float(fine_statistic), 0.) if converged else 0.
    roots, weights = roots_hermite(64)
    sigma = np.exp(alternative.x[-1])
    means = np.asarray([weights @ expit(alternative.x[0] + (alternative.x[index] if index else 0.) + np.sqrt(2) * sigma * roots) / np.sqrt(np.pi) for index in range(len(levels))])
    return {"p_value": float(chi2.sf(statistic, len(levels) - 1)) if converged else 1., "statistic": statistic, "degrees_of_freedom": len(levels) - 1, "converged": converged, "levels": levels, "coefficients": alternative.x[:-2], "standardized_means": np.column_stack([means, 1 - means]) if converged else np.full((len(levels), 2), np.nan), "null_concentration": float(np.exp(null.x[-2])), "alternative_concentration": float(np.exp(alternative.x[-2])), "null_subject_sd": float(np.exp(null.x[-1])), "alternative_subject_sd": float(sigma), "quadrature_error": float(error), "null_fit": null, "alternative_fit": alternative, "n_subjects": len(np.unique(subjects)), "n_observations": len(labels), "nodes": nodes, "refine": refine}
