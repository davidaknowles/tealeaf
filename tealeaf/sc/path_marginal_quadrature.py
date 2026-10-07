"""Alternative numerical integration for the experimental binary EC model.

These evaluators change quadrature, not the latent model or statistical test.
The prior-only rule isolates potential importance-mixture normalization error
near singular Beta endpoints. It is not a production inference backend.
"""

import numpy as np
from scipy.special import betaln, expit, logsumexp

from .path_marginal import BinaryECPathLikelihood, beta_quadrature


class PriorBinaryECPathLikelihood(BinaryECPathLikelihood):
    """Direct normalized Gauss-Jacobi integration under the Beta prior.

    No numerical importance-mixture normalization is used. Strong likelihood
    peaks may need more nodes, so the existing doubled-order validation gate
    remains necessary. Numerical scores differentiate this actual quadrature;
    log-density moments are not assumed accurately integrated at endpoints.
    """

    def log_integrals(self, logits, concentration, *, nodes=16, gradient=False):
        logits = np.asarray(logits, dtype=float)
        if logits.ndim != 2 or len(logits) != len(self.subjects) or concentration <= 0 or not np.isfinite(concentration):
            raise ValueError("aligned mean logits and positive precision required")
        if gradient:
            step = 1e-3
            value = self.log_integrals(logits, concentration, nodes=nodes)
            eta = (self.log_integrals(logits + step, concentration, nodes=nodes) - self.log_integrals(logits - step, concentration, nodes=nodes)) / (2 * step)
            precision = (self.log_integrals(logits, concentration * np.exp(step), nodes=nodes) - self.log_integrals(logits, concentration * np.exp(-step), nodes=nodes)) / (2 * step)
            return value, eta, precision
        means = np.clip(expit(logits), 1e-12, 1 - 1e-12)
        output = np.zeros_like(means)
        for row in range(len(means)):
            if np.isfinite(self.binomial_terms[row]).all():
                first, second, constant = self.binomial_terms[row]
                a, b = concentration * means[row], concentration * (1 - means[row])
                output[row] = betaln(a + first, b + second) - betaln(a, b) + constant
                continue
            for column, mean in enumerate(means[row]):
                points, weights = beta_quadrature(concentration * mean, concentration * (1 - mean), nodes)
                output[row, column] = logsumexp(np.log(np.maximum(weights, 1e-300)) + self.row_log_likelihood(row, points))
        return output


def prior_quadrature(likelihood):
    """Reuse prepared actual EC likelihood, changing only its integrator."""
    return PriorBinaryECPathLikelihood(likelihood.counts, likelihood.components, likelihood.subjects, likelihood.labels, likelihood.proposal_counts, likelihood.offsets)


def matrix_log_likelihood(likelihood, rows, points):
    """Reuse EC probabilities across rows evaluated at the same path grid.

    Rows has length N_selected, points has length L. The returned array is
    N_selected by L. Primer normalizers remain path dependent. This is a
    matrix evaluation of the same per-row likelihood, with no count pooling.
    """
    rows, points = np.asarray(rows, dtype=int), np.asarray(points, dtype=float)
    result = np.zeros((len(rows), len(points)))
    for counts, components in zip(likelihood.counts, likelihood.components):
        if counts[rows].sum() == 0:
            continue
        mass = components[:, 0, None] + points[None, :] * components[:, 1, None] + (1 - points[None, :]) * components[:, 2, None]
        observed = counts[rows]
        result += observed @ np.log(np.maximum(mass, 1e-300)) - observed.sum(axis=1, keepdims=True) * np.log(mass.sum(axis=0))[None, :]
    return result - likelihood.offsets[rows, None]


class SharedPriorBinaryECPathLikelihood(PriorBinaryECPathLikelihood):
    """Share quadrature and EC probability grids for identical type designs.

    Subject likelihoods are still independent before subject-level products.
    In the random-subject model, prior offset grids are identical across all
    subjects; counts, not prior grids, produce their different posteriors.
    Identical full logit rows are grouped exactly, without rounding or padding.
    """

    def log_integrals(self, logits, concentration, *, nodes=16, gradient=False):
        if gradient:
            return super().log_integrals(logits, concentration, nodes=nodes, gradient=True)
        logits = np.asarray(logits, dtype=float)
        if logits.ndim != 2 or len(logits) != len(self.subjects) or concentration <= 0 or not np.isfinite(concentration):
            raise ValueError("aligned mean logits and positive precision required")
        means = np.clip(expit(logits), 1e-12, 1 - 1e-12)
        output = np.zeros_like(means)
        exact = np.isfinite(self.binomial_terms).all(axis=1)
        if exact.any():
            first, second, constant = self.binomial_terms[exact].T
            a, b = concentration * means[exact], concentration * (1 - means[exact])
            output[exact] = betaln(a + first[:, None], b + second[:, None]) - betaln(a, b) + constant[:, None]
        unique, grouping = np.unique(logits, axis=0, return_inverse=True)
        for group in np.unique(grouping[~exact]):
            rows = np.flatnonzero((grouping == group) & ~exact)
            rules = [beta_quadrature(concentration * mean, concentration * (1 - mean), nodes) for mean in np.clip(expit(unique[group]), 1e-12, 1 - 1e-12)]
            points = np.stack([rule[0] for rule in rules])
            weights = np.stack([rule[1] for rule in rules])
            values = matrix_log_likelihood(self, rows, points.ravel()).reshape(len(rows), len(unique[group]), nodes)
            output[rows] = logsumexp(values + np.log(np.maximum(weights, 1e-300))[None, :, :], axis=2)
        return output


def shared_prior_quadrature(likelihood):
    """Change only evaluation batching, preserving prior-only quadrature."""
    return SharedPriorBinaryECPathLikelihood(likelihood.counts, likelihood.components, likelihood.subjects, likelihood.labels, likelihood.proposal_counts, likelihood.offsets)
