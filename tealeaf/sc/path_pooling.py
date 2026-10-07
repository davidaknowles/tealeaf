"""Experimental joint regression over covariance-matched local-path counts."""

import numpy as np
from dataclasses import replace
from scipy.special import softmax

from . import differential
from .ec_block_glmm import pooled_isoform_weights, blocked_multilevel_design
from .ec_glmm import ECGLMMData
from .path_reporting import proportion_covariance


def subject_isoform_baselines(data, subjects, *, max_iter=250):
    """Label-blind T-dimensional nuisance mixture estimated for each subject.

    All types of a subject share within-path transcript ratios and outside
    mass. No cell-type labels or external outcomes enter this estimate.
    Failed positive-depth baseline fits fail the sensitivity, not silently
    reverting to another baseline model. Baseline uncertainty is conditional.
    """
    subjects = np.asarray(subjects)
    if subjects.shape != (len(data.counts[0]),):
        raise ValueError("subjects must align with EC observations")
    result = {}
    for subject in np.unique(subjects):
        selected = subjects == subject
        local = ECGLMMData(tuple(counts[selected] for counts in data.counts), data.compatibility, data.design[selected], data.clusters[selected])
        weights, converged = pooled_isoform_weights(local, max_iter=max_iter, return_status=True)
        if not converged:
            raise ValueError(f"subject baseline fit failed for {subject}")
        result[subject] = weights
    return result


def quantify_effective_paths(data, path_index, labels, subjects, *, baseline=None, concentration=.25, covariance_source="likelihood", max_iter=100, subject_baselines=None):
    """Aggregate primer counts by subject/type and retain path uncertainty.

    Returned proportions/counts are N by S, labels/subjects/depths length N.
    Effective depth matches the conditional ILR covariance to a multinomial,
    capped at observed EC depth. Likelihood covariance excludes quantification
    prior information; posterior covariance is an explicitly named sensitivity.
    All positive-count subject/type aggregates must converge and be identifiable.
    A failure rejects the block instead of silently changing its design.
    Effective counts are fractional measurement approximations, not read counts.
    Within-path transcript shares and outside-block mass condition on baseline.
    Optional subject_baselines maps subject IDs to their fixed T-vector of
    nuisance weights, shared across all labels. The default remains pooled.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    path_index = np.asarray(path_index, dtype=int)
    if labels.shape != subjects.shape or labels.shape != (data.counts[0].shape[0],):
        raise ValueError("labels and subjects must align with EC observations")
    if covariance_source not in ("likelihood", "posterior") or not np.isfinite(concentration) or concentration <= 0:
        raise ValueError("positive concentration and a supported covariance source required")
    if path_index.shape != (data.n_isoforms,) or not np.any(path_index >= 0):
        raise ValueError("path indices must align with transcripts")
    size = int(path_index[path_index >= 0].max()) + 1
    if size < 2 or not np.array_equal(np.unique(path_index[path_index >= 0]), np.arange(size)):
        raise ValueError("at least two consecutive path categories required")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    basis = differential.helmert_basis(size)
    proportions, depths, effective, selected_labels, selected_subjects, covariances = [], [], [], [], [], []
    for subject in np.unique(subjects):
        local_baseline = baseline
        if subject_baselines is not None:
            if subject not in subject_baselines:
                raise ValueError(f"missing subject baseline for {subject}")
            local_baseline = np.asarray(subject_baselines[subject], dtype=float)
            if local_baseline.shape != (data.n_isoforms,) or not np.isfinite(local_baseline).all() or np.any(local_baseline < 0) or local_baseline.sum() <= 0:
                raise ValueError("valid aligned subject baseline required")
            local_baseline = local_baseline / local_baseline.sum()
        for level in np.unique(labels[subjects == subject]):
            selected = (subjects == subject) & (labels == level)
            counts = tuple(np.asarray(matrix[selected], dtype=float).sum(axis=0) for matrix in data.counts)
            totals = [float(array.sum()) for array in counts]
            depth = sum(totals)
            if depth <= 0:
                continue
            fitted = differential.fit_path_perturbation(counts, data.compatibility, local_baseline, path_index, path_pseudocount=concentration, path_pseudocount_scaling="total", max_iter=max_iter)
            if not fitted.converged:
                raise ValueError(f"path quantification failed for subject {subject}, level {level}")
            covariance = fitted.covariance
            if covariance_source == "likelihood":
                information = differential.conditional_path_information(fitted.theta, path_index, basis, data.compatibility, totals)
                covariance = differential.identifiable_covariance(information, np.eye(size - 1), rtol=1e-10)
            if not covariance.identifiable:
                raise ValueError(f"unidentifiable {covariance_source} covariance for subject {subject}, level {level}")
            number = differential.effective_multinomial_size(fitted.path_proportions, covariance.covariance, maximum=depth)
            proportions.append(fitted.path_proportions)
            covariances.append(proportion_covariance(replace(fitted, covariance=covariance)))
            depths.append(depth)
            effective.append(number)
            selected_labels.append(level)
            selected_subjects.append(subject)
    if len(proportions) < 4:
        raise ValueError("fewer than four positive-count subject/type aggregates")
    proportions = np.asarray(proportions)
    effective = np.asarray(effective)
    return {"proportions": proportions, "counts": effective[:, None] * proportions, "depths": np.asarray(depths), "effective_depths": effective, "labels": np.asarray(selected_labels), "subjects": np.asarray(selected_subjects), "proportion_covariances": np.asarray(covariances), "scalar_proportion_covariances": (proportions[:, :, None] * np.eye(size) - proportions[:, :, None] * proportions[:, None, :]) / effective[:, None, None], "concentration": concentration, "covariance_source": covariance_source}


def joint_path_dm_test(quantified, *, labels=None, fitted_null=None, max_iter=250, dispersion_method="ml", concentration=None):
    """Subject-blocked DM likelihood-ratio test, re-estimating dispersion.

    Conditional means are softmax(X B), with B of shape P by (S-1).
    Subject intercepts occur in null and alternative; only C-1 cell-type
    coefficients are tested. One biological precision is estimated per block,
    separately under null and alternative. It is not the quantification prior.
    Gamma-function likelihood on fractional effective counts is approximate.
    The analytic chi-square tail requires calibration with count-level nulls;
    label exchangeability alone does not establish biological calibration.
    Experimental Cox-Reid selection refits precision for each tested design;
    fixed precision is a diagnostic. Neither permits reuse of fitted_null.
    """
    counts = np.asarray(quantified["counts"], dtype=float)
    subjects = np.asarray(quantified["subjects"])
    labels = np.asarray(quantified["labels"] if labels is None else labels)
    if len(np.unique(subjects)) < 4:
        raise ValueError("at least four subject clusters required")
    design, tested, levels, _ = blocked_multilevel_design(labels, subjects)
    if len(counts) <= design.shape[1]:
        raise ValueError("subject-blocked regression needs residual observations")
    null_design = design[:, :tested[0]]
    if dispersion_method == "ml" and concentration is None:
        result = differential.dirichlet_multinomial_test(counts, null_design, design, max_iter=max_iter, fix_null_concentration=False, fitted_null=fitted_null)
    elif dispersion_method in ("cox_reid", "fixed"):
        from .path_dispersion import corrected_dm_test
        if fitted_null is not None or (dispersion_method == "fixed") != (concentration is not None):
            raise ValueError("corrected dispersion needs fresh fits and fixed precision only in diagnostic mode")
        result = corrected_dm_test(counts, null_design, design, concentration=concentration, max_iter=max_iter)
    else:
        raise ValueError("unsupported dispersion method or precision")
    converged = bool(result["null_converged"] and result["alternative_converged"])
    if not converged:
        result["p_value"] = 1.
    coefficients = result["alternative_coefficients"]
    # Average each modeled cell-type mean over the SAME represented subjects.
    # This is a descriptive standardization, not a posterior usage interval.
    subject_design = null_design[np.unique(subjects, return_index=True)[1]]
    subject_logits = subject_design @ coefficients[:tested[0]]
    means = []
    for index in range(len(levels)):
        logits = subject_logits if index == 0 else subject_logits + coefficients[tested[index - 1]]
        means.append(softmax(np.column_stack([logits, np.zeros(len(logits))]), axis=1).mean(axis=0))
    return {**result, "converged": converged, "n_subjects": len(subject_design), "n_observations": len(counts), "levels": levels, "standardized_means": np.asarray(means), "effective_depth_median": float(np.median(quantified["effective_depths"]))}
