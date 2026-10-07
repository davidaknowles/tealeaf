"""Experimental uncertainty-weighted EC scores with biological heterogeneity.

No target smoothing prior enters either the null fit or its efficient score.
The Gaussian approximation, estimated nuisance parameters and small-sample F
reference require independent EC count-null checks, not just sign-flip checks.
"""

from dataclasses import dataclass

import numpy as np
from scipy import linalg, optimize, sparse, stats, special

from . import differential
from .ec_block_glmm import pooled_isoform_weights
from .path_bias import SharedPathNullProblem

MODEL_VERSION = "v3_scaled_shared_null"


def score_contrast_proportions(anchor, contrast):
    """Approximate A/B usages from a shared S-path anchor and fitted ILR delta.

    This maps a Gaussian one-step contrast to the simplex, not a fresh EC
    maximum-likelihood or posterior usage fit. A common scalar offset in log
    anchor cancels. A/B log-ratio difference equals contrast exactly before
    machine underflow. No strength or direction is selected with LR outcomes.
    """
    anchor, contrast = np.asarray(anchor, dtype=float), np.asarray(contrast, dtype=float)
    if anchor.ndim != 1 or len(anchor) < 2 or contrast.shape != (len(anchor) - 1,) or not np.isfinite(anchor).all() or not np.isfinite(contrast).all() or (anchor < 0).any() or anchor.sum() <= 0:
        raise ValueError("finite nonnegative path anchor and aligned ILR contrast required")
    logits = np.log(np.maximum(anchor / anchor.sum(), 1e-12))
    displacement = differential.helmert_basis(len(anchor)) @ contrast / 2
    return np.stack([special.softmax(logits - displacement), special.softmax(logits + displacement)])


@dataclass
class PathScoreComponents:
    """M scores of length D and M D-by-D information/heterogeneity matrices."""

    scores: np.ndarray
    information: np.ndarray
    biological_shapes: np.ndarray
    subject_ids: np.ndarray
    levels: tuple
    null_fits: list
    reporting_proportions: list


def efficient_shared_path_score(counts, designs, theta, path_index, observed_levels, n_levels):
    """Profile common path usage and independent type-specific nuisance.

    theta is C_observed by T. Reference-coded effects have dimension
    D=(C_global-1)*(S-1), with each contrast equal to type minus reference ILR.
    Expected-count whitening and orthogonal projection avoid subtracting nearly
    equal Fisher matrices. Nuisance directions use the linear transcript
    simplex, so nearly absent transcripts do not vanish solely because of a
    log-coordinate derivative. Zero-mass unobserved EC rows carry no score.
    """
    theta = np.asarray(theta, dtype=float)
    path_index, observed_levels = np.asarray(path_index, dtype=int), np.asarray(observed_levels, dtype=int)
    paths = np.unique(path_index[path_index >= 0])
    size, types, transcripts = len(paths), *theta.shape
    if size < 2 or not np.array_equal(paths, np.arange(size)) or path_index.shape != (transcripts,) or observed_levels.shape != (types,) or (observed_levels < 0).any() or (observed_levels >= n_levels).any():
        raise ValueError("aligned consecutive paths and global type indices required")
    basis = differential.helmert_basis(size)
    dimension = (n_levels - 1) * (size - 1)
    nuisance_dimension = transcripts - size
    psi = differential.path_proportions(theta[0], path_index)
    if (psi <= 0).any():
        raise ValueError("positive shared null path proportions required")
    if any(not np.allclose(differential.path_proportions(row, path_index), psi, atol=1e-8) for row in theta):
        raise ValueError("score must be evaluated at a shared-path null")
    # Analytic linear-simplex null space. Numerical null_space of the ILR
    # Jacobian loses nominal rank when one fitted path is nearly absent.
    groups = [np.flatnonzero(path_index == path) for path in paths]
    outside = np.flatnonzero(path_index < 0)
    nuisance_columns = []
    for group in groups + ([outside] if len(outside) else []):
        if len(group) > 1:
            embedded = np.zeros((transcripts, len(group) - 1))
            embedded[group] = differential.helmert_basis(len(group))
            nuisance_columns.extend(embedded.T)
    if len(outside):
        mass_direction = np.zeros(transcripts)
        for path, group in enumerate(groups):
            mass_direction[group] = psi[path] / len(group)
        mass_direction[outside] = -1 / len(outside)
        nuisance_columns.append(mass_direction / linalg.norm(mass_direction))
    nuisance = np.asarray(nuisance_columns).T if nuisance_columns else np.zeros((transcripts, 0))
    target_jacobian = np.zeros((transcripts, size - 1))
    target_jacobian[path_index >= 0] = basis[path_index[path_index >= 0]] - psi @ basis
    target_blocks, nuisance_blocks, residuals = [], [], []
    for local, level in enumerate(observed_levels):
        target = theta[local, :, None] * target_jacobian
        contrast = np.zeros((transcripts, dimension))
        if level > 0:
            contrast[:, (level - 1) * (size - 1):level * (size - 1)] = target
        null = np.zeros((transcripts, size - 1 + types * nuisance_dimension))
        null[:, :size - 1] = target
        null[:, size - 1 + local * nuisance_dimension:size - 1 + (local + 1) * nuisance_dimension] = nuisance
        for observed, mapping in zip(counts, designs):
            observed = np.asarray(observed, dtype=float)[local]
            total = float(observed.sum())
            if total <= 0:
                continue
            mapping = np.asarray(sparse.csr_matrix(mapping).todense())
            mass = mapping @ theta[local]
            normalizer = float(mass.sum())
            if normalizer <= 0 or ((mass <= 0) & (observed > 0)).any():
                raise ValueError("positive EC counts require positive null probability")
            active = mass > 0
            probability = mass[active] / normalizer
            joint = np.column_stack([contrast, null])
            derivative = (mapping[active] @ joint - probability[:, None] * (mapping.sum(axis=0) @ joint)) / normalizer
            whitened = np.sqrt(total / probability)[:, None] * derivative
            target_blocks.append(whitened[:, :dimension])
            nuisance_blocks.append(whitened[:, dimension:])
            residuals.append((observed[active] - total * probability) / np.sqrt(total * probability))
    if not target_blocks:
        raise ValueError("positive primer totals required")
    target, null, residual = np.vstack(target_blocks), np.vstack(nuisance_blocks), np.concatenate(residuals)
    norms = linalg.norm(null, axis=0)
    nonzero = norms > np.finfo(float).tiny
    q = linalg.orth(null[:, nonzero] / norms[nonzero], rcond=1e-10)
    efficient = target - q @ (q.T @ target)
    score = efficient.T @ residual
    information = efficient.T @ efficient
    score[np.abs(score) < 64 * np.finfo(float).eps * max(sum(np.asarray(x).sum() for x in counts), 1.)] = 0.
    path_shape = (basis.T / psi) @ basis
    biological_shape = np.kron(np.eye(n_levels - 1) + np.ones((n_levels - 1, n_levels - 1)), path_shape)
    return score, (information + information.T) / 2, biological_shape


def shared_path_score_components(data, path_index, labels, subjects, *, baseline=None, max_iter=300, reporting_concentration=None):
    """Fit each subject null once; failures invalidate the complete hypothesis.

    Missing/zero-total types are structural exclusions, not optimizer-based
    subject selection. Independent weak-prior usage is optional reporting;
    its failures produce unavailable reports without removing inference rows.
    """
    labels, subjects = np.asarray(labels), np.asarray(subjects)
    if labels.shape != subjects.shape or labels.shape != (len(data.counts[0]),):
        raise ValueError("labels and subjects must align with observations")
    levels, encoded = np.unique(labels, return_inverse=True)
    if len(levels) < 2:
        raise ValueError("at least two cell types required")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    if reporting_concentration is not None and (not np.isfinite(reporting_concentration) or reporting_concentration < 0):
        raise ValueError("finite nonnegative reporting concentration required")
    size = len(np.unique(np.asarray(path_index)[np.asarray(path_index) >= 0]))
    dimension = (len(levels) - 1) * (size - 1)
    scores, information, shapes, retained, fits, reports = [], [], [], [], [], []
    for subject in np.unique(subjects):
        local_levels = np.unique(encoded[subjects == subject])
        counts = [np.array([np.asarray(matrix[(subjects == subject) & (encoded == level)], dtype=float).sum(axis=0) for level in local_levels]) for matrix in data.counts]
        positive = sum(value.sum(axis=1) for value in counts) > 0
        local_levels, counts = local_levels[positive], tuple(value[positive] for value in counts)
        if len(local_levels) < 2:
            continue
        shared = SharedPathNullProblem(counts, data.compatibility, baseline, path_index).fit(max_iter=max_iter)
        if not shared.converged:
            raise ValueError(f"shared-path null failed for subject {subject}, iterations={shared.iterations}, gradient={shared.gradient_norm:.6g}, termination={shared.termination_message}")
        score, info, shape = efficient_shared_path_score(counts, data.compatibility, shared.theta, path_index, local_levels, len(levels))
        scores.append(score)
        information.append(info)
        shapes.append(shape)
        retained.append(subject)
        fits.append(shared)
        local_reports = []
        if reporting_concentration is not None:
            for local, level in enumerate(local_levels):
                try:
                    fit = differential.fit_free_isoform_paths(tuple(value[local] for value in counts), data.compatibility, baseline, path_index, path_pseudocount=reporting_concentration, path_pseudocount_scaling="total", max_iter=max_iter, tolerance=1e-12)
                    proportions = fit.path_proportions if fit.converged else np.full(size, np.nan)
                except (ValueError, np.linalg.LinAlgError):
                    proportions = np.full(size, np.nan)
                local_reports.append((int(level), proportions))
        reports.append(local_reports)
    return PathScoreComponents(np.asarray(scores).reshape(-1, dimension), np.asarray(information).reshape(-1, dimension, dimension), np.asarray(shapes).reshape(-1, dimension, dimension), np.asarray(retained), tuple(levels), fits, reports)


def mixed_score_test(scores, information, biological_shapes=None, *, biological_variance=None):
    """REML aggregation of g_u ~ N(I_u beta, I_u + tau^2 I_u B_u I_u).

    PSD/singular information is represented only in identifiable directions;
    missing contrasts are not assigned zero measurement variance. beta is D,
    information is M by D by D. The primary reference is an experimental
    modified-Knapp-Hartung-style F(D,M_informative-1), with residual scale at
    least one. Only the scalar full-rank case is the usual meta-analysis mKH;
    the multivariate/estimated-count-information extension is approximate.
    """
    scores, information = np.asarray(scores, dtype=float), np.asarray(information, dtype=float)
    if scores.ndim != 2 or scores.shape[1] < 1 or information.shape != (*scores.shape, scores.shape[1]) or not np.isfinite(scores).all() or not np.isfinite(information).all():
        raise ValueError("finite M-by-D scores and M-by-D-by-D information required")
    dimension = scores.shape[1]
    if biological_shapes is None:
        biological_shapes = np.broadcast_to(np.eye(dimension), information.shape)
    biological_shapes = np.asarray(biological_shapes, dtype=float)
    if biological_shapes.shape != information.shape or not np.isfinite(biological_shapes).all():
        raise ValueError("finite aligned biological covariance shapes required")
    if biological_variance is not None and (not np.isfinite(biological_variance) or biological_variance < 0):
        raise ValueError("finite nonnegative biological variance required")
    observations = []
    for score, info, shape in zip(scores, information, biological_shapes):
        eigenvalues, vectors = linalg.eigh((info + info.T) / 2)
        tolerance = max(float(eigenvalues[-1]), 1.) * 1e-10
        if eigenvalues[0] < -tolerance or linalg.eigvalsh((shape + shape.T) / 2)[0] < -1e-10:
            raise ValueError("information and biological shapes must be positive semidefinite")
        keep = eigenvalues > tolerance
        if not keep.any():
            continue
        design = vectors[:, keep].T
        values = design @ score / eigenvalues[keep]
        measurement = np.diag(1 / eigenvalues[keep])
        biological = design @ shape @ design.T
        observations.append((values, design, measurement, biological))
    clusters = len(observations)
    residual_df = sum(len(values) for values, _, _, _ in observations) - dimension
    if clusters < 4 or residual_df <= 0:
        raise ValueError("four informative subject clusters and positive residual df required")

    def fit(variance):
        precision, rhs, quadratic, logdet = np.zeros((dimension, dimension)), np.zeros(dimension), 0., 0.
        for values, design, measurement, biological in observations:
            covariance = measurement + variance * biological
            factor = linalg.cho_factor(covariance, lower=True)
            weighted_design = linalg.cho_solve(factor, design)
            weighted_values = linalg.cho_solve(factor, values)
            precision += design.T @ weighted_design
            rhs += design.T @ weighted_values
            quadratic += float(values @ weighted_values)
            logdet += 2 * np.log(np.diag(factor[0])).sum()
        factor = linalg.cho_factor(precision, lower=True)
        mean = linalg.cho_solve(factor, rhs)
        covariance = linalg.cho_solve(factor, np.eye(dimension))
        residual = max(quadratic - float(rhs @ mean), 0.)
        objective = logdet + 2 * np.log(np.diag(factor[0])).sum() + residual
        return objective, mean, covariance, residual, precision

    # Structural aggregate rank is checked once; individual subjects may lack
    # some contrasts. No small ridge fabricates information in missing axes.
    if np.linalg.matrix_rank(sum(design.T @ design for _, design, _, _ in observations)) < dimension:
        raise ValueError("aggregate score design does not identify every requested contrast")
    if biological_variance is None:
        scale = np.median([np.trace(measurement) / max(np.trace(biological), 1e-12) for _, _, measurement, biological in observations])
        search = optimize.minimize_scalar(lambda value: fit(np.exp(value))[0], bounds=(np.log(scale) - 24, np.log(scale) + 16), method="bounded", options={"xatol": 1e-6})
        if not search.success:
            raise ValueError("biological-variance REML optimization failed")
        biological_variance = min((0., float(np.exp(search.x))), key=lambda value: fit(value)[0])
    objective, mean, covariance, residual, precision = fit(biological_variance)
    inflation = max(1., residual / residual_df)
    wald = max(float(mean @ precision @ mean), 0.)
    statistic = wald / (dimension * inflation)
    return {"p_value": float(stats.f.sf(statistic, dimension, clusters - 1)), "statistic": statistic, "degrees_of_freedom": dimension, "denominator_degrees_of_freedom": clusters - 1, "n_subjects": clusters, "converged": True, "mean_difference": mean, "mean_covariance": inflation * covariance, "biological_variance": float(biological_variance), "residual_inflation": inflation, "residual_degrees_of_freedom": residual_df, "restricted_objective": objective, "chi_square_p_value": float(stats.chi2.sf(wald, dimension)), "residual_F_p_value": float(stats.f.sf(statistic, dimension, residual_df))}


def mixed_path_score_test(data, path_index, labels, subjects, **kwargs):
    """General C-type experimental EC score test, with optional usage reports."""
    components = shared_path_score_components(data, path_index, labels, subjects, **kwargs)
    result = mixed_score_test(components.scores, components.information, components.biological_shapes)
    # Unweighted one-step responses are diagnostic only, NOT the fitted mean.
    differences = np.array([linalg.pinvh(info, rtol=1e-10) @ score for score, info in zip(components.scores, components.information)])
    return {**result, "differences": differences, "components": components, "subject_ids": components.subject_ids, "levels": components.levels, "n_fitted_subjects": len(components.subject_ids)}
