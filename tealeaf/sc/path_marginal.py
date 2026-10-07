"""Experimental binary-path EC likelihood with integrated subject effects.

Both primers share a latent Beta path composition. A Gaussian subject logit
offset is shared across types. Counts are never replaced by estimated path
counts. Quadrature accuracy, biological calibration, and real endpoints must
be assessed before use; this is not a production backend.
"""

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize
from scipy.special import betaln, digamma, expit, logsumexp, roots_hermitenorm
from scipy.stats import chi2

from . import differential
from .ec_block_glmm import collapse_isoforms_to_paths, pooled_isoform_weights

MODEL_VERSION = "binary_ec_random_subject_beta_v4_supported_primer_endpoints"


def beta_quadrature(a, b, nodes=16):
    """Normalized Gauss-Jacobi rule on (0,1), avoiding Gamma-mass overflow.

    Positive Beta shapes a,b are scalars; returns two length-Q vectors.
    Golub-Welsch eigenvector squares are probability weights, including at
    high concentration where unnormalized Jacobi weights can overflow.
    """
    a, b, nodes = float(a), float(b), int(nodes)
    if not np.isfinite(a + b) or min(a, b) <= 0 or nodes < 2:
        raise ValueError("positive finite Beta shapes and at least two nodes required")
    total = a + b
    index = np.arange(1, nodes, dtype=float)
    diagonal = np.r_[(a - b) / total, (a - b) * (total - 2) / ((2 * index + total - 2) * (2 * index + total))]
    off_diagonal = np.zeros(nodes - 1)
    off_diagonal[0] = np.sqrt(4 * a * b / (total ** 2 * (total + 1)))
    index = index[1:]
    off_diagonal[1:] = np.sqrt(4 * index * (index + a - 1) * (index + b - 1) * (index + total - 2) / ((2 * index + total - 2) ** 2 * (2 * index + total - 1) * (2 * index + total - 3)))
    values, vectors = eigh_tridiagonal(diagonal, off_diagonal, check_finite=False)
    weights = np.square(vectors[0])
    weights /= weights.sum()
    return np.clip((values + 1) / 2, 1e-14, 1 - 1e-14), weights


@dataclass
class BinaryECPathLikelihood:
    """N aggregates, each with primer-specific K-vectors and affine EC masses.

    Components are K by 3, columns are outside mass, path-0 mass, path-1
    mass. Proposal counts N by 2 are ONLY numerical integration proposals,
    not a response, a prior, or extra biological observations. Row offsets
    subtract a fixed EC log-likelihood anchor for numerical stability.
    """
    counts: tuple
    components: tuple
    subjects: np.ndarray
    labels: np.ndarray
    proposal_counts: np.ndarray
    offsets: np.ndarray

    def __post_init__(self):
        self.counts = tuple(np.asarray(values, dtype=float) for values in self.counts)
        self.components = tuple(np.asarray(values, dtype=float) for values in self.components)
        self.subjects, self.labels = np.asarray(self.subjects), np.asarray(self.labels)
        self.proposal_counts, self.offsets = np.asarray(self.proposal_counts, dtype=float), np.asarray(self.offsets, dtype=float)
        n = len(self.subjects)
        if self.subjects.shape != (n,) or self.labels.shape != (n,) or self.proposal_counts.shape != (n, 2) or self.offsets.shape != (n,) or not np.isfinite(self.offsets).all() or not np.isfinite(self.proposal_counts).all() or (self.proposal_counts < 0).any():
            raise ValueError("aligned subject/type aggregates and nonnegative integration proposals required")
        if not self.counts or len(self.counts) != len(self.components):
            raise ValueError("aligned primer counts and EC components required")
        for counts, components in zip(self.counts, self.components):
            if counts.ndim != 2 or counts.shape[0] != n or components.shape != (counts.shape[1], 3) or not np.isfinite(counts).all() or not np.isfinite(components).all() or (counts < 0).any() or (components < 0).any():
                raise ValueError("finite nonnegative N by K counts and K by 3 EC components required")
            if max((components[:, 0] + components[:, 1]).sum(), (components[:, 0] + components[:, 2]).sum()) <= 0 and counts.sum() > 0:
                raise ValueError("observed primer requires positive mass at an interior path composition")
            impossible = components.sum(axis=1) == 0
            if (counts[:, impossible] > 0).any():
                raise ValueError("positive EC counts require a nonzero model probability")

    @cached_property
    def binomial_terms(self):
        """Exact Beta integral when observed ECs are pure paths or constants.

        Returns N by 3, path-0 exponent, path-1 exponent, log constant.
        Mixed ECs or path-dependent primer normalizers require quadrature;
        their entries are NaN. This is a likelihood identity, not rounded EC
        allocation or a change to the statistical model.
        """
        terms = np.zeros((len(self.subjects), 3))
        for counts, components in zip(self.counts, self.components):
            if counts.sum() == 0:
                continue
            first = components[:, 0] + components[:, 1]
            second = components[:, 0] + components[:, 2]
            # A primer can see just one path, or both in the same EC ratios.
            # Conditional on its gene total it then carries no path-usage
            # information, and its normalized likelihood is constant. A zero
            # endpoint normalizer is not an invalid interior likelihood.
            first_sum, second_sum = first.sum(), second.sum()
            if min(first_sum, second_sum) == 0 or np.allclose(first / first_sum, second / second_sum, rtol=1e-12, atol=0):
                probabilities = first / first_sum if first_sum > 0 else second / second_sum
                terms[:, 2] += counts @ np.log(np.maximum(probabilities, 1e-300))
                continue
            if not np.isclose(first.sum(), second.sum(), rtol=1e-12, atol=1e-15):
                terms[counts.sum(axis=1) > 0] = np.nan
                continue
            pure_first = (second == 0) & (first > 0)
            pure_second = (first == 0) & (second > 0)
            constant = np.isclose(first, second, rtol=1e-12, atol=0)
            unsupported = ~(pure_first | pure_second | constant)
            terms[(counts[:, unsupported] > 0).any(axis=1)] = np.nan
            terms[:, 0] += counts[:, pure_first].sum(axis=1)
            terms[:, 1] += counts[:, pure_second].sum(axis=1)
            scale = np.maximum(first, second) / max(first.sum(), 1e-300)
            terms[:, 2] += counts @ np.log(np.maximum(scale, 1e-300))
        terms[:, 2] -= self.offsets
        return terms

    def row_log_likelihood(self, row, proportions):
        proportions = np.asarray(proportions, dtype=float)
        result = np.zeros(proportions.shape)
        for counts, components in zip(self.counts, self.components):
            observed = counts[row]
            if observed.sum() <= 0:
                continue
            mass = components[:, 0] + proportions[:, None] * components[:, 1] + (1 - proportions[:, None]) * components[:, 2]
            result += np.log(np.maximum(mass, 1e-300)) @ observed - observed.sum() * np.log(mass.sum(axis=1))
        return result - self.offsets[row]

    def log_integrals(self, logits, concentration, *, nodes=16, gradient=False):
        """Log E_Beta[EC likelihood] at N by Q_subject mean logits.

        A mixture of prior and evidence-focused Gauss-Jacobi rules retains
        both endpoints and likelihood peaks. Importance weights are normalized
        numerically, so a flat likelihood integrates to exactly one. Doubling
        both quadrature orders independently assesses approximation error.
        """
        logits = np.asarray(logits, dtype=float)
        if logits.ndim != 2 or len(logits) != len(self.subjects) or concentration <= 0 or not np.isfinite(concentration):
            raise ValueError("aligned mean logits and positive precision required")
        original_means = expit(logits)
        means = np.clip(original_means, 1e-12, 1 - 1e-12)
        output = np.zeros_like(means)
        if gradient:
            eta_score, precision_score = np.zeros_like(means), np.zeros_like(means)
        for row in range(len(means)):
            if np.isfinite(self.binomial_terms[row]).all():
                count_a, count_b, constant = self.binomial_terms[row]
                a, b = concentration * means[row], concentration * (1 - means[row])
                output[row] = betaln(a + count_a, b + count_b) - betaln(a, b) + constant
                if gradient:
                    score_a = digamma(a + count_a) - digamma(a)
                    score_b = digamma(b + count_b) - digamma(b)
                    eta_score[row] = concentration * means[row] * (1 - means[row]) * (score_a - score_b)
                    eta_score[row] *= (original_means[row] > 1e-12) & (original_means[row] < 1 - 1e-12)
                    precision_score[row] = concentration * (means[row] * score_a + (1 - means[row]) * score_b + digamma(concentration) - digamma(concentration + count_a + count_b))
                continue
            proposal = self.proposal_counts[row]
            for column, mean in enumerate(means[row]):
                a, b = concentration * mean, concentration * (1 - mean)
                first, first_weights = beta_quadrature(a, b, nodes)
                second, second_weights = beta_quadrature(a + proposal[0], b + proposal[1], nodes)
                points = np.r_[first, second]
                log_p, log_other = np.log(points), np.log1p(-points)
                prior_density = (a - 1) * log_p + (b - 1) * log_other - betaln(a, b)
                proposal_density = (a + proposal[0] - 1) * log_p + (b + proposal[1] - 1) * log_other - betaln(a + proposal[0], b + proposal[1])
                mixture = np.logaddexp(prior_density, proposal_density) - np.log(2.)
                weights = np.log(np.maximum(np.r_[first_weights, second_weights], 1e-300)) - np.log(2.) + prior_density - mixture
                weights -= logsumexp(weights)
                output[row, column] = logsumexp(weights + self.row_log_likelihood(row, points))
                if gradient:
                    posterior = np.exp(weights + self.row_log_likelihood(row, points) - output[row, column])
                    eta_score[row, column] = posterior @ (concentration * mean * (1 - mean) * (log_p - log_other - digamma(a) + digamma(b)))
                    precision_score[row, column] = posterior @ (concentration * (mean * log_p + (1 - mean) * log_other - mean * digamma(a) - (1 - mean) * digamma(b) + digamma(concentration)))
            if gradient:
                # A smooth EC likelihood can integrate accurately even when
                # its log-Beta score has unresolved endpoint singularities.
                # Differentiate the quadrature itself in this regime rather
                # than treating the inaccurate log-moment rule as a score.
                singular = concentration * np.minimum(means[row], 1 - means[row]) < 2
                if singular.any():
                    local = BinaryECPathLikelihood(tuple(values[row:row + 1] for values in self.counts), self.components, self.subjects[row:row + 1], self.labels[row:row + 1], self.proposal_counts[row:row + 1], self.offsets[row:row + 1])
                    selected = logits[row:row + 1, singular]
                    step = 1e-3
                    eta_score[row, singular] = (local.log_integrals(selected + step, concentration, nodes=nodes)[0] - local.log_integrals(selected - step, concentration, nodes=nodes)[0]) / (2 * step)
                    precision_score[row, singular] = (local.log_integrals(selected, concentration * np.exp(step), nodes=nodes)[0] - local.log_integrals(selected, concentration * np.exp(-step), nodes=nodes)[0]) / (2 * step)
        return (output, eta_score, precision_score) if gradient else output


def prepare_binary_ec_likelihood(data, path_index, labels, subjects, *, baseline=None):
    """Collapse fixed within-path shares, aggregate subject/type EC counts.

    Baseline is length T, path_index length T, and labels/subjects length I.
    Only consecutive binary paths 0 and 1 are supported in this prototype.
    Counts remain in the original primer-specific EC space. Outside mass and
    within-path ratios condition on a label-blind baseline, not on LR outcomes.
    """
    path_index = np.asarray(path_index, dtype=int)
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if path_index.shape != (data.n_isoforms,) or not np.array_equal(np.unique(path_index[path_index >= 0]), [0, 1]) or subjects.shape != labels.shape or subjects.shape != (len(data.counts[0]),):
        raise ValueError("aligned binary paths, labels and subjects required")
    if baseline is None:
        baseline, converged = pooled_isoform_weights(data, max_iter=250, return_status=True)
        if not converged:
            raise ValueError("pooled transcript baseline did not converge")
    else:
        baseline = np.asarray(baseline, dtype=float)
    if baseline.shape != path_index.shape or not np.isfinite(baseline).all() or np.any(baseline < 0) or baseline.sum() <= 0:
        raise ValueError("valid transcript baseline required")
    baseline = np.maximum(baseline, 1e-12)
    baseline /= baseline.sum()
    collapsed, _ = collapse_isoforms_to_paths(data, path_index, baseline)
    mass = baseline[path_index >= 0].sum()
    components = tuple(np.column_stack([mapping[:, 2:] @ baseline[path_index < 0], mass * mapping[:, 0], mass * mapping[:, 1]]) for mapping in collapsed.compatibility)
    selected_labels, selected_subjects, counts, proposals = [], [], [], []
    for subject in np.unique(subjects):
        for label in np.unique(labels[subjects == subject]):
            selected = (subjects == subject) & (labels == label)
            observed = tuple(values[selected].sum(axis=0) for values in data.counts)
            depth = sum(values.sum() for values in observed)
            if depth <= 0:
                continue
            fitted = differential.fit_path_perturbation(observed, data.compatibility, baseline, path_index, path_pseudocount=.25, path_pseudocount_scaling="total")
            if not fitted.converged or not fitted.covariance.identifiable:
                raise ValueError("evidence-proposal fit failed or path contrast is unidentifiable")
            # This depth affects quadrature placement only, not the likelihood.
            effective = differential.effective_multinomial_size(fitted.path_proportions, fitted.covariance.covariance, maximum=depth)
            proposals.append(effective * fitted.path_proportions)
            counts.append(observed)
            selected_labels.append(label)
            selected_subjects.append(subject)
    if not counts:
        raise ValueError("at least one positive-count aggregate required")
    aggregated = tuple(np.stack([values[primer] for values in counts]) for primer in range(len(data.counts)))
    result = BinaryECPathLikelihood(aggregated, components, np.asarray(selected_subjects), np.asarray(selected_labels), np.asarray(proposals), np.zeros(len(counts)))
    for row, proposal in enumerate(proposals):
        result.offsets[row] = result.row_log_likelihood(row, np.asarray([proposal[0] / np.sum(proposal)]))[0]
    return result


def random_subject_objective(parameters, likelihood, design, *, subject_nodes=12, path_nodes=16, gradient=True):
    """Marginal EC negative log-likelihood and integrated analytic score.

    Parameters are C mean-logit coefficients, log kappa, and log subject SD.
    X is N by C, subject offset is shared across all types/primers, and each
    subject/type Beta composition is shared by its two primer likelihoods.
    Beta posterior expectations of its density score differentiate the path
    integral. Pure-path cases use exact digamma differences. This avoids five
    independent quadrature evaluations per regular optimizer step. Beta shapes
    below two use numerical quadrature differences, since score log-moments
    can have unresolved endpoint singularities despite accurate objectives.
    Score/objective agreement still depends on quadrature accuracy for mixed ECs.
    """
    parameters, design = np.asarray(parameters, dtype=float), np.asarray(design, dtype=float)
    if design.ndim != 2 or parameters.shape != (design.shape[1] + 2,) or design.shape[0] != len(likelihood.subjects) or not np.isfinite(design).all() or not np.isfinite(parameters).all():
        raise ValueError("aligned coefficients, variance parameters and design required")
    coefficients, concentration, subject_sd = parameters[:-2], np.exp(parameters[-2]), np.exp(parameters[-1])
    nodes, weights = roots_hermitenorm(subject_nodes)
    weights /= weights.sum()
    offsets = subject_sd * nodes
    logits = (design @ coefficients)[:, None] + offsets
    if gradient:
        values, eta_score, precision_score = likelihood.log_integrals(logits, concentration, nodes=path_nodes, gradient=True)
        score = np.zeros_like(parameters)
    else:
        values = likelihood.log_integrals(logits, concentration, nodes=path_nodes)
    objective = 0.
    for subject in np.unique(likelihood.subjects):
        selected = likelihood.subjects == subject
        joint = values[selected].sum(axis=0) + np.log(weights)
        evidence = logsumexp(joint)
        objective -= evidence
        if gradient:
            posterior = np.exp(joint - evidence)
            score[:-2] += design[selected].T @ (eta_score[selected] @ posterior)
            score[-2] += posterior @ precision_score[selected].sum(axis=0)
            score[-1] += posterior @ (eta_score[selected].sum(axis=0) * offsets)
    return (float(objective), -score) if gradient else float(objective)


def binary_marginal_test(likelihood, *, labels=None, max_iter=100, subject_nodes=12, path_nodes=16, quadrature_tolerance=1e-3):
    """Random-subject binary-path EC marginal LRT, not a fixed-intercept DM.

    No M nuisance means are estimated. Kappa and Gaussian subject SD are
    refitted under null and alternative. Each fit must pass doubled-quadrature
    LR validation; analytic chi-square tails still need count-null assessment.
    Unsupported multi-path blocks must stay p=1 in an eligible full universe.
    """
    labels = np.asarray(likelihood.labels if labels is None else labels)
    levels = np.unique(labels)
    if labels.shape != likelihood.subjects.shape or len(levels) < 2 or len(np.unique(likelihood.subjects)) < 4:
        raise ValueError("aligned labels, two levels and at least four subjects required")
    repeated = [len(np.unique(labels[likelihood.subjects == subject])) >= 2 for subject in np.unique(likelihood.subjects)]
    if sum(repeated) < 4:
        raise ValueError("at least four subjects with repeated types required")
    design = np.column_stack([np.ones(len(labels)), *[labels == level for level in levels[1:]]]).astype(float)
    if len(labels) <= design.shape[1] + 2 or np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("full-rank type design and residual observations required")

    # Explicit wrappers preserve the requested quadrature order in optimization.
    def fit(local_design, initial):
        return minimize(lambda parameters: random_subject_objective(parameters, likelihood, local_design, subject_nodes=subject_nodes, path_nodes=path_nodes), initial, jac=True, method="L-BFGS-B", bounds=[(-8., 8.)] * local_design.shape[1] + [(np.log(.1), np.log(1e4)), (np.log(.02), np.log(2.))], options={"maxiter": max_iter, "ftol": 1e-9, "gtol": 1e-4})

    fractions = likelihood.proposal_counts[:, 0] / likelihood.proposal_counts.sum(axis=1)
    mean = np.clip(fractions.mean(), .01, .99)
    null = fit(design[:, :1], np.r_[np.log(mean / (1 - mean)), np.log(20.), np.log(.5)])
    initial = np.r_[null.x[:1], np.zeros(len(levels) - 1), null.x[-2:]]
    alternative = fit(design, initial)
    null_fine = random_subject_objective(null.x, likelihood, design[:, :1], subject_nodes=2 * subject_nodes, path_nodes=2 * path_nodes, gradient=False)
    alt_fine = random_subject_objective(alternative.x, likelihood, design, subject_nodes=2 * subject_nodes, path_nodes=2 * path_nodes, gradient=False)
    difference = 2 * (null.fun - alternative.fun)
    fine_difference = 2 * (null_fine - alt_fine)
    null_error, alternative_error = abs(null.fun - null_fine), abs(alternative.fun - alt_fine)
    error = max(abs(difference - fine_difference), null_error, alternative_error)
    converged = bool(null.success and alternative.success and np.isfinite([null_fine, alt_fine]).all() and min(difference, fine_difference) >= -1e-6 and error <= quadrature_tolerance)
    statistic = max(float(fine_difference), 0.) if converged else 0.
    mean_nodes, mean_weights = roots_hermitenorm(64)
    mean_weights /= mean_weights.sum()
    means = []
    for level in range(len(levels)):
        coefficient = alternative.x[0] + (alternative.x[level] if level else 0.)
        value = mean_weights @ expit(coefficient + np.exp(alternative.x[-1]) * mean_nodes)
        means.append([value, 1 - value])
    validated_means = np.asarray(means) if converged else np.full((len(levels), 2), np.nan)
    return {"model": "binary_ec_random_subject_beta", "model_version": MODEL_VERSION, "p_value": float(chi2.sf(statistic, len(levels) - 1)) if converged else 1., "statistic": statistic, "degrees_of_freedom": len(levels) - 1, "converged": converged, "n_subjects": len(np.unique(likelihood.subjects)), "n_observations": len(labels), "levels": levels, "standardized_means": validated_means, "null_concentration": float(np.exp(null.x[-2])), "alternative_concentration": float(np.exp(alternative.x[-2])), "null_subject_sd": float(np.exp(null.x[-1])), "alternative_subject_sd": float(np.exp(alternative.x[-1])), "quadrature_error": error, "null_quadrature_error": null_error, "alternative_quadrature_error": alternative_error, "null_fit": null, "alternative_fit": alternative}
