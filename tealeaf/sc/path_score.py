"""Experimental paired EC-score inference without cell-specific prior shifts."""

import numpy as np
from scipy import linalg

from . import differential
from .ec_block_glmm import pooled_isoform_weights


def omnibus_path_score_test(data, path_index, labels, subjects, *, baseline=None, denominator_concentration=32., null_concentration=1., replicates=64, bootstrap_draws=4095, seed=0, null_fit_tolerance=1e-7):
    """Experimental subject-cluster wild EC-score omnibus.

    A common S-path mixture is fitted per subject, aggregating all observed
    cell types. With C types, scores have (C-1)*(S-1) coordinates. Subject
    mixture nuisance scores are projected out before forming any contrast.
    Expected likelihood information plus a zero-centered denominator penalty
    preconditions the aggregate score; it is not a posterior cell-type fit.
    One Rademacher sign per subject resamples the entire score vector, so
    cell-type contrasts remain dependent. Wild p-values require approximately
    symmetric independent subject scores under the null. They are not exact
    count-likelihood p-values; actual count-null checks remain necessary.
    Within-path shares and outside mass condition on the baseline mixture.
    Neither this procedure nor its one-step coefficients are production defaults.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    path_index = np.asarray(path_index, dtype=int)
    if labels.shape != subjects.shape or labels.shape != (data.counts[0].shape[0],):
        raise ValueError("labels and subjects must align with observations")
    levels, encoded = np.unique(labels, return_inverse=True)
    if len(levels) < 2 or denominator_concentration < 0 or null_concentration < 0:
        raise ValueError("at least two cell types and nonnegative concentrations required")
    if replicates < 0 or bootstrap_draws < max(replicates, 1):
        raise ValueError("bootstrap draws must cover the requested null families")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    size = int(path_index[path_index >= 0].max()) + 1
    path_basis, type_basis = differential.helmert_basis(size), differential.helmert_basis(len(levels))
    dimension = (len(levels) - 1) * (size - 1)
    scores, information, penalties, retained, omitted = [], [], [], [], []
    for subject in np.unique(subjects):
        positions = np.flatnonzero(subjects == subject)
        observed_levels = np.unique(encoded[positions])
        if len(observed_levels) < 2:
            continue
        cell_counts = [tuple(np.asarray(matrix[positions[encoded[positions] == level]], dtype=float).sum(axis=0) for matrix in data.counts) for level in observed_levels]
        positive = [index for index, counts in enumerate(cell_counts) if sum(array.sum() for array in counts) > 0]
        if len(positive) < 2:
            continue
        cell_counts = [cell_counts[index] for index in positive]
        observed_levels = observed_levels[positive]
        combined = tuple(sum(pair[primer] for pair in cell_counts) for primer in range(len(data.counts)))
        fitted = differential.fit_path_perturbation(combined, data.compatibility, baseline, path_index, path_pseudocount=null_concentration, path_pseudocount_scaling="total", tolerance=null_fit_tolerance)
        if not fitted.converged:
            omitted.append(str(subject))
            continue
        local_score, local_info = zip(*[conditional_path_score(fitted.theta, path_index, counts, data.compatibility) for counts in cell_counts])
        total_info = sum(local_info)
        inverse = linalg.pinvh(total_info, rtol=1e-10)
        cross = sum(np.kron(type_basis[level, :, None], matrix) for level, matrix in zip(observed_levels, local_info))
        contrast = sum(np.kron(type_basis[level], value) for level, value in zip(observed_levels, local_score))
        efficient_score = contrast - cross @ inverse @ sum(local_score)
        roundoff = 32 * np.finfo(float).eps * max(sum(array.sum() for array in combined), 1.)
        efficient_score[np.abs(efficient_score) < roundoff] = 0.
        efficient_info = sum(np.kron(np.outer(type_basis[level], type_basis[level]), matrix) for level, matrix in zip(observed_levels, local_info)) - cross @ inverse @ cross.T
        proportions = fitted.path_proportions
        path_penalty = path_basis.T @ (np.diag(proportions) - np.outer(proportions, proportions)) @ path_basis
        scores.append(efficient_score)
        information.append((efficient_info + efficient_info.T) / 2)
        penalties.append(np.kron(np.eye(len(levels) - 1), path_penalty))
        retained.append(subject)
    if len(scores) < 4:
        raise ValueError("at least four fitted subject clusters required")
    scores = np.asarray(scores).reshape(-1, dimension)
    total_info = sum(information)
    rank = np.linalg.matrix_rank(total_info, tol=max(float(linalg.norm(total_info, 2)), 1.) * 1e-10)
    if rank == 0:
        raise ValueError("no identifiable omnibus score contrasts")
    weight = linalg.pinvh(total_info + denominator_concentration * sum(penalties), rtol=1e-10)
    gram = scores @ weight @ scores.T
    gram = (gram + gram.T) / 2
    rng = np.random.default_rng(seed)
    signs = np.vstack([np.ones(len(scores)), rng.choice([-1., 1.], size=(bootstrap_draws, len(scores)))])
    statistics = np.maximum(np.einsum("bi,ij,bj->b", signs, gram, signs, optimize=True), 0.)
    # Include the observed sign vector; ties count conservatively in the tail.
    tails = (len(statistics) - np.searchsorted(np.sort(statistics), statistics, side="left")) / len(statistics)
    result = {"statistic": float(statistics[0]), "p_value": float(tails[0]), "degrees_of_freedom": rank, "n_subjects": len(scores), "converged": True}
    null = [{"replicate": index, "statistic": float(statistics[index + 1]), "p_value": float(tails[index + 1]), "degrees_of_freedom": rank} for index in range(replicates)]
    coefficients = (linalg.pinvh(total_info, rtol=1e-10) @ scores.sum(axis=0)).reshape(len(levels) - 1, size - 1)
    return {**result, "null": null, "cluster_scores": scores, "subject_ids": np.asarray(retained), "omitted_subjects": omitted, "levels": tuple(levels), "one_step_coefficients": coefficients, "bootstrap_draws": bootstrap_draws}


def paired_subject_centered_test(data, path_index, labels, subjects, *, baseline=None, concentration=32.):
    """Center both cell-type priors on one subject's common EC mixture.

    The common mixture is fitted from both cell types with all transcript
    weights free, using total path concentration 1. Both separate local fits
    then use its same within-path shares and prior mean. This removes the
    direct uniform-prior/coverage contrast under a fixed-subject null, but
    reuses data for the prior mean and still requires null calibration.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if labels.shape != subjects.shape or len(labels) != data.counts[0].shape[0]:
        raise ValueError("labels and subjects must align with observations")
    levels = np.unique(labels)
    if len(levels) != 2 or concentration <= 0:
        raise ValueError("two cell types and positive concentration required")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    size = int(np.max(np.asarray(path_index)[np.asarray(path_index) >= 0])) + 1
    responses, fits, retained = [], [], []
    for subject in np.unique(subjects):
        positions = [np.flatnonzero((subjects == subject) & (labels == level)) for level in levels]
        if any(len(index) == 0 for index in positions):
            continue
        counts = [tuple(np.asarray(matrix[index], dtype=float).sum(axis=0) for matrix in data.counts) for index in positions]
        if any(sum(values.sum() for values in pair) <= 0 for pair in counts):
            continue
        combined = tuple(first + second for first, second in zip(*counts))
        common = differential.fit_free_isoform_paths(combined, data.compatibility, baseline, path_index, path_pseudocount=1., path_pseudocount_scaling="total")
        if not common.converged:
            continue
        local = [differential.fit_path_perturbation(pair, data.compatibility, common.theta, path_index, path_pseudocount=concentration, path_pseudocount_scaling="total", path_prior_center="baseline") for pair in counts]
        if not all(fit.converged for fit in local):
            continue
        responses.append(local[1].path_logratios - local[0].path_logratios)
        fits.append(local)
        retained.append(subject)
    responses = np.asarray(responses).reshape(-1, size - 1)
    result = differential.paired_mean_test(responses)
    return {**result, "differences": responses, "path_fits": fits, "subject_ids": np.asarray(retained), "levels": tuple(levels), "concentration": concentration}


def conditional_path_score(theta, path_index, counts, designs):
    """Unpenalized EC score (S-1,) and expected information (S-1,S-1)."""
    theta = np.asarray(theta, dtype=float)
    path_index = np.asarray(path_index, dtype=int)
    selected = path_index >= 0
    size = int(path_index[selected].max()) + 1
    basis = differential.helmert_basis(size)
    features = basis[path_index[selected]]
    log_jacobian = features - (theta[selected] / theta[selected].sum()) @ features
    score = np.zeros(size - 1)
    totals = []
    for observed, mapping in zip(counts, designs):
        observed = np.asarray(observed, dtype=float)
        mass = np.asarray(mapping @ theta).ravel()
        total, normalizer = float(observed.sum()), float(mass.sum())
        totals.append(total)
        if total <= 0:
            continue
        if normalizer <= 0 or np.any((observed > 0) & (mass <= 0)):
            raise ValueError("positive counts require positive EC compatibility mass")
        transcript_score = theta * np.asarray(mapping.T @ (observed / np.maximum(mass, 1e-300))).ravel()
        transcript_score -= total * theta * np.asarray(mapping.sum(axis=0)).ravel() / normalizer
        score += log_jacobian.T @ transcript_score[selected]
    information = differential.conditional_path_information(theta, path_index, basis, designs, totals)
    return score, information


def paired_path_score_test(data, path_index, labels, subjects, *, baseline=None, denominator_concentration=32., null_concentration=1., free_null=False):
    """Paired mean test of nuisance-projected, regularized one-step EC effects.

    A common null path mixture is fitted once per subject, not once per cell
    type. Prior gradients never enter the contrast score. Common-null nuisance
    scores are projected out using expected likelihood information. The
    denominator prior regularizes the efficient-score response toward zero,
    rather than shrinking two differently covered cell types toward a common
    uniform usage. Conditional uncertainty and baseline estimation are not
    fully modeled, so count-null and biological-null audits remain necessary.
    This experimental procedure is not a production default.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    path_index = np.asarray(path_index, dtype=int)
    if labels.shape != subjects.shape or len(labels) != data.counts[0].shape[0]:
        raise ValueError("labels and subjects must align with observations")
    if denominator_concentration < 0 or null_concentration < 0:
        raise ValueError("nonnegative score regularization required")
    levels = np.unique(labels)
    if len(levels) != 2:
        raise ValueError("two cell types required")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    size = int(path_index[path_index >= 0].max()) + 1
    basis = differential.helmert_basis(size)
    responses, retained, fits, omitted = [], [], [], []
    for subject in np.unique(subjects):
        positions = [np.flatnonzero((subjects == subject) & (labels == level)) for level in levels]
        if any(len(index) == 0 for index in positions):
            continue
        counts = [tuple(np.asarray(matrix[index], dtype=float).sum(axis=0) for matrix in data.counts) for index in positions]
        if any(sum(values.sum() for values in pair) <= 0 for pair in counts):
            continue
        combined = tuple(first + second for first, second in zip(*counts))
        fitter = differential.fit_free_isoform_paths if free_null else differential.fit_path_perturbation
        fitted = fitter(combined, data.compatibility, baseline, path_index, path_pseudocount=null_concentration, path_pseudocount_scaling="total")
        if not fitted.converged:
            omitted.append(str(subject))
            continue
        scores, information = zip(*[conditional_path_score(fitted.theta, path_index, pair, data.compatibility) for pair in counts])
        common_inverse = linalg.pinvh(information[0] + information[1], rtol=1e-10)
        efficient_score = scores[1] - information[1] @ common_inverse @ (scores[0] + scores[1])
        roundoff = 32 * np.finfo(float).eps * max(sum(values.sum() for pair in counts for values in pair), 1.)
        efficient_score[np.abs(efficient_score) < roundoff] = 0.
        efficient_information = information[1] - information[1] @ common_inverse @ information[1]
        efficient_information = (efficient_information + efficient_information.T) / 2
        eigenvalues = linalg.eigvalsh(efficient_information)
        if eigenvalues[0] <= max(float(eigenvalues[-1]), 1.) * 1e-10:
            omitted.append(str(subject))
            continue
        proportions = fitted.path_proportions
        regularizer = denominator_concentration * basis.T @ (np.diag(proportions) - np.outer(proportions, proportions)) @ basis
        response = linalg.solve(efficient_information + regularizer, efficient_score, assume_a="pos")
        responses.append(response)
        retained.append(subject)
        fits.append(fitted)
    responses = np.asarray(responses).reshape(-1, size - 1)
    result = differential.paired_mean_test(responses)
    return {**result, "differences": responses, "subject_ids": np.asarray(retained), "null_fits": fits, "omitted_subjects": omitted, "levels": tuple(levels), "denominator_concentration": denominator_concentration, "null_concentration": null_concentration, "free_null": free_null}
