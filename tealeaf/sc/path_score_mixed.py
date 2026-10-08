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
    score_coordinate: str = "ilr"
    reference_information: np.ndarray | None = None


def path_proportion_covariance(proportions):
    """Stable S by S categorical covariance, including nearly absent paths.

    Pair products avoid subtracting nearly equal diagonal terms at a simplex
    boundary. This is observation/Dirichlet covariance shape, not a prior.
    """
    proportions = np.asarray(proportions, dtype=float)
    if proportions.ndim != 1 or len(proportions) < 2 or not np.isfinite(proportions).all() or (proportions < 0).any() or proportions.sum() <= 0:
        raise ValueError("finite nonnegative path proportions required")
    proportions = proportions / proportions.sum()
    products = np.outer(proportions, proportions)
    np.fill_diagonal(products, 0.)
    covariance = -products
    np.fill_diagonal(covariance, products.sum(axis=1))
    return covariance


def efficient_shared_path_score(counts, designs, theta, path_index, observed_levels, n_levels, *, score_coordinate="ilr", return_reference=False):
    """Profile common path usage and independent type-specific nuisance.

    theta is C_observed by T. Reference-coded effects have dimension
    D=(C_global-1)*(S-1), with each contrast equal to type minus reference ILR.
    Expected-count whitening and orthogonal projection avoid subtracting nearly
    equal Fisher matrices. Nuisance directions use the linear transcript
    simplex, so nearly absent transcripts do not vanish solely because of a
    log-coordinate derivative. Zero-mass unobserved EC rows carry no score.
    Optional proportion coordinates test H.T @ (psi_type - psi_reference)
    rather than ILR differences. Subject-dependent transformations change
    the common random-effects estimand; they are not a numerical-only fix.
    """
    theta = np.asarray(theta, dtype=float)
    path_index, observed_levels = np.asarray(path_index, dtype=int), np.asarray(observed_levels, dtype=int)
    paths = np.unique(path_index[path_index >= 0])
    size, types, transcripts = len(paths), *theta.shape
    if size < 2 or not np.array_equal(paths, np.arange(size)) or path_index.shape != (transcripts,) or observed_levels.shape != (types,) or (observed_levels < 0).any() or (observed_levels >= n_levels).any():
        raise ValueError("aligned consecutive paths and global type indices required")
    if score_coordinate not in ("ilr", "proportion"):
        raise ValueError("score coordinate must be ilr or proportion")
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
    if score_coordinate == "ilr":
        target_jacobian[path_index >= 0] = basis[path_index[path_index >= 0]] - psi @ basis
    else:
        target_jacobian[path_index >= 0] = basis[path_index[path_index >= 0]] / psi[path_index[path_index >= 0], None]
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
    path_shape = (basis.T / psi) @ basis if score_coordinate == "ilr" else basis.T @ path_proportion_covariance(psi) @ basis
    biological_shape = np.kron(np.eye(n_levels - 1) + np.ones((n_levels - 1, n_levels - 1)), path_shape)
    result = score, (information + information.T) / 2, biological_shape
    return (*result, target.T @ target) if return_reference else result


def shared_path_score_components(data, path_index, labels, subjects, *, baseline=None, max_iter=300, reporting_concentration=None, null_multistart=False, score_coordinate="ilr", count_likelihood="multinomial"):
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
    if score_coordinate not in ("ilr", "proportion"):
        raise ValueError("score coordinate must be ilr or proportion")
    if count_likelihood not in ("multinomial", "conditional") or (count_likelihood == "conditional" and (score_coordinate != "ilr" or len(levels) != 2)):
        raise ValueError("conditional count prototype requires two types and ILR coordinates")
    if baseline is None:
        baseline = pooled_isoform_weights(data)
    if reporting_concentration is not None and (not np.isfinite(reporting_concentration) or reporting_concentration < 0):
        raise ValueError("finite nonnegative reporting concentration required")
    size = len(np.unique(np.asarray(path_index)[np.asarray(path_index) >= 0]))
    dimension = (len(levels) - 1) * (size - 1)
    scores, information, shapes, retained, fits, reports, reference = [], [], [], [], [], [], []
    for subject in np.unique(subjects):
        local_levels = np.unique(encoded[subjects == subject])
        counts = [np.array([np.asarray(matrix[(subjects == subject) & (encoded == level)], dtype=float).sum(axis=0) for level in local_levels]) for matrix in data.counts]
        positive = sum(value.sum(axis=1) for value in counts) > 0
        local_levels, counts = local_levels[positive], tuple(value[positive] for value in counts)
        if len(local_levels) < 2:
            continue
        if count_likelihood == "conditional":
            from .conditional_path_score import ConditionalPathNullProblem, efficient_conditional_path_score
            conditional = ConditionalPathNullProblem(counts, data.compatibility, baseline, path_index)
            shared, offsets = conditional.fit(max_iter=max_iter, multistart=null_multistart)
        else:
            shared = SharedPathNullProblem(counts, data.compatibility, baseline, path_index).fit(max_iter=max_iter, multistart=null_multistart)
        if not shared.converged:
            raise ValueError(f"shared-path null failed for subject {subject}, iterations={shared.iterations}, gradient={shared.gradient_norm:.6g}, termination={shared.termination_message}")
        if count_likelihood == "conditional":
            score, info, shape, target_info = efficient_conditional_path_score(conditional, shared.theta, offsets)
        else:
            score, info, shape, target_info = efficient_shared_path_score(counts, data.compatibility, shared.theta, path_index, local_levels, len(levels), score_coordinate=score_coordinate, return_reference=True)
        scores.append(score)
        information.append(info)
        shapes.append(shape)
        reference.append(target_info)
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
    return PathScoreComponents(np.asarray(scores).reshape(-1, dimension), np.asarray(information).reshape(-1, dimension, dimension), np.asarray(shapes).reshape(-1, dimension, dimension), np.asarray(retained), tuple(levels), fits, reports, score_coordinate, np.asarray(reference).reshape(-1, dimension, dimension))


def mixed_score_test(scores, information, biological_shapes=None, *, biological_variance=None, reference_information=None, scalar_fast=False):
    """REML aggregation of g_u ~ N(I_u beta, I_u + tau^2 I_u B_u I_u).

    PSD/singular information is represented only in identifiable directions;
    missing contrasts are not assigned zero measurement variance. beta is D,
    information is M by D by D. The primary reference is an experimental
    modified-Knapp-Hartung-style F(D,M_informative-1), with residual scale at
    least one. Only the scalar full-rank case is the usual meta-analysis mKH;
    the multivariate/estimated-count-information extension is approximate.
    Optional unprofiled target information defines a unit-invariant numerical
    rank metric. The mean, measurement covariance and biological model remain
    in the original target coordinates. No reference information is added to
    the profiled information, and exact aliasing remains uninformative.
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
    if reference_information is not None:
        reference_information = np.asarray(reference_information, dtype=float)
        if reference_information.shape != information.shape or not np.isfinite(reference_information).all():
            raise ValueError("finite aligned reference information required")
    if scalar_fast and dimension == 1:
        return _scalar_mixed_score_test(scores[:, 0], information[:, 0, 0], biological_shapes[:, 0, 0], biological_variance, None if reference_information is None else reference_information[:, 0, 0])
    observations = []
    for index, (score, info, shape) in enumerate(zip(scores, information, biological_shapes)):
        if reference_information is not None:
            info_values = linalg.eigvalsh((info + info.T) / 2)
            if info_values[0] < -max(float(info_values[-1]), 1.) * 1e-10:
                raise ValueError("information must be positive semidefinite")
            reference = reference_information[index]
            ref_values, ref_vectors = linalg.eigh((reference + reference.T) / 2)
            ref_tolerance = max(float(ref_values[-1]), np.finfo(float).tiny) * 1e-10
            if ref_values[0] < -ref_tolerance:
                raise ValueError("reference information must be positive semidefinite")
            supported = ref_values > ref_tolerance
            if not supported.any():
                continue
            root = ref_vectors[:, supported] * np.sqrt(ref_values[supported])
            whitening = (ref_vectors[:, supported] / np.sqrt(ref_values[supported])).T
            scaled = whitening @ ((info + info.T) / 2) @ whitening.T
            eigenvalues, vectors = linalg.eigh((scaled + scaled.T) / 2)
            tolerance = max(float(eigenvalues[-1]), 1.) * 1e-10
            if eigenvalues[0] < -tolerance or linalg.eigvalsh((shape + shape.T) / 2)[0] < -1e-10:
                raise ValueError("information and biological shapes must be positive semidefinite")
            if eigenvalues[-1] > 1 + 1e-8:
                raise ValueError("profiled information exceeds its target reference")
            keep = eigenvalues > tolerance
            if not keep.any():
                continue
            design = vectors[:, keep].T @ root.T
            values = vectors[:, keep].T @ (whitening @ score) / eigenvalues[keep]
            measurement = np.diag(1 / eigenvalues[keep])
            observations.append((values, design, measurement, design @ shape @ design.T))
            continue
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
    mean_scales = np.ones(dimension)
    if reference_information is not None:
        # A computational reparameterization of the SAME common mean, not
        # per-subject rescaling of its biological estimand. Undo on output.
        mean_scales = linalg.norm(np.vstack([design for _, design, _, _ in observations]), axis=0)
        if (mean_scales <= np.finfo(float).tiny).any():
            raise ValueError("aggregate score design does not identify every requested contrast")
        observations = [(values, design / mean_scales, measurement, biological) for values, design, measurement, biological in observations]

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
    return {"p_value": float(stats.f.sf(statistic, dimension, clusters - 1)), "statistic": statistic, "degrees_of_freedom": dimension, "denominator_degrees_of_freedom": clusters - 1, "n_subjects": clusters, "converged": True, "mean_difference": mean / mean_scales, "mean_covariance": inflation * covariance / mean_scales[:, None] / mean_scales[None, :], "biological_variance": float(biological_variance), "residual_inflation": inflation, "residual_degrees_of_freedom": residual_df, "restricted_objective": objective, "chi_square_p_value": float(stats.chi2.sf(wald, dimension)), "residual_F_p_value": float(stats.f.sf(statistic, dimension, residual_df))}


def _scalar_mixed_score_test(scores, information, shapes, biological_variance, reference):
    """Exact scalar specialization, same rank rule, REML search and F tail."""
    if (information < -np.maximum(information, 1.) * 1e-10).any() or (shapes < -1e-10).any():
        raise ValueError("information and biological shapes must be positive semidefinite")
    if reference is None:
        keep = information > np.maximum(information, 1.) * 1e-10
        values = scores[keep] / information[keep]
        design = np.ones(keep.sum())
        measurement, biological = 1 / information[keep], shapes[keep]
        mean_scale = 1.
    else:
        tolerance = np.maximum(reference, np.finfo(float).tiny) * 1e-10
        if (reference < -tolerance).any():
            raise ValueError("reference information must be positive semidefinite")
        supported = reference > tolerance
        eigenvalues = information[supported] / reference[supported]
        threshold = np.maximum(eigenvalues, 1.) * 1e-10
        if (eigenvalues < -threshold).any():
            raise ValueError("information and biological shapes must be positive semidefinite")
        if (eigenvalues > 1 + 1e-8).any():
            raise ValueError("profiled information exceeds its target reference")
        keep = supported.copy()
        keep[supported] = eigenvalues > threshold
        design = np.sqrt(reference[keep])
        values = (scores[keep] / design) / (information[keep] / reference[keep])
        measurement = reference[keep] / information[keep]
        biological = reference[keep] * shapes[keep]
        mean_scale = linalg.norm(design)
        if len(design) and mean_scale <= np.finfo(float).tiny:
            raise ValueError("aggregate score design does not identify every requested contrast")
        design = design / mean_scale if len(design) else design
    clusters = len(values)
    if clusters < 4:
        raise ValueError("four informative subject clusters and positive residual df required")

    def fit(variance):
        variances = measurement + variance * biological
        if (variances <= 0).any():
            raise np.linalg.LinAlgError("scalar covariance is not positive definite")
        weights = 1 / variances
        precision = float(np.sum(design**2 * weights))
        if precision <= 0:
            raise np.linalg.LinAlgError("scalar mean is not identified")
        rhs = float(np.sum(design * values * weights))
        mean, covariance = rhs / precision, 1 / precision
        residual = max(float(np.sum(values**2 * weights)) - rhs * mean, 0.)
        objective = float(np.log(variances).sum()) + np.log(precision) + residual
        return objective, mean, covariance, residual, precision

    if biological_variance is None:
        scale = np.median(measurement / np.maximum(biological, 1e-12))
        search = optimize.minimize_scalar(lambda value: fit(np.exp(value))[0], bounds=(np.log(scale) - 24, np.log(scale) + 16), method="bounded", options={"xatol": 1e-6})
        if not search.success:
            raise ValueError("biological-variance REML optimization failed")
        biological_variance = min((0., float(np.exp(search.x))), key=lambda value: fit(value)[0])
    objective, mean, covariance, residual, precision = fit(biological_variance)
    residual_df = clusters - 1
    inflation = max(1., residual / residual_df)
    wald = max(mean * precision * mean, 0.)
    statistic = wald / inflation
    return dict(p_value=float(stats.f.sf(statistic, 1, clusters - 1)), statistic=statistic, degrees_of_freedom=1, denominator_degrees_of_freedom=clusters - 1, n_subjects=clusters, converged=True, mean_difference=np.array([mean / mean_scale]), mean_covariance=np.array([[inflation * covariance / mean_scale**2]]), biological_variance=float(biological_variance), residual_inflation=inflation, residual_degrees_of_freedom=residual_df, restricted_objective=float(objective), chi_square_p_value=float(stats.chi2.sf(wald, 1)), residual_F_p_value=float(stats.f.sf(statistic, 1, residual_df)))


def aggregate_path_scores(components, *, information_metric="absolute", scalar_fast=False):
    """Aggregate already-fitted subject scores without repeating null fits."""
    if information_metric not in ("absolute", "reference") or (information_metric == "reference" and components.reference_information is None):
        raise ValueError("reference rank metric requires unprofiled target information")
    result = mixed_score_test(components.scores, components.information, components.biological_shapes, reference_information=components.reference_information if information_metric == "reference" else None, scalar_fast=scalar_fast)
    # Unweighted one-step responses are diagnostic only, NOT the fitted mean.
    differences = np.array([linalg.pinvh(info, rtol=1e-10) @ score for score, info in zip(components.scores, components.information)])
    return {**result, "differences": differences, "components": components, "subject_ids": components.subject_ids, "levels": components.levels, "n_fitted_subjects": len(components.subject_ids), "score_coordinate": components.score_coordinate, "information_metric": information_metric}


def signed_path_score_p_value(components, rng, *, information_metric="absolute", scalar_fast=False):
    """Paired subject-sign diagnostic, refitting REML for every realization.

    This is not an independent count-null calibration. The binary type
    contrast changes sign, its subject information and heterogeneity shape
    do not. Components are never modified in place.
    """
    if len(components.levels) != 2:
        raise ValueError("subject sign null requires exactly two types")
    if information_metric not in ("absolute", "reference") or (information_metric == "reference" and components.reference_information is None):
        raise ValueError("reference rank metric requires unprofiled target information")
    signs = rng.choice((-1., 1.), size=len(components.subject_ids))
    return mixed_score_test(components.scores * signs[:, None], components.information, components.biological_shapes, reference_information=components.reference_information if information_metric == "reference" else None, scalar_fast=scalar_fast)["p_value"]


def binary_score_subject_influence(components, *, information_metric="reference"):
    """Describe conditional subject precision, without changing inference.

    For M binary scores, d_u=S_u/I_u and w_u=1/(1/I_u+tau^2 B_u).
    Normalized w_u are precision shares of the common mean, not EC depth
    weights or an alternative degrees-of-freedom rule. Effective weighted
    subjects is 1/sum(share_u**2). Stable residuals are evaluated at the
    SAME fitted biological variance, not an independently optimized model.
    Returns M diagnostic records and a scalar summary; excludes no subjects
    beyond the existing numerical information-rank rule.
    """
    size = len(components.subject_ids)
    if len(components.levels) != 2 or components.scores.shape != (size, 1) or len(set(map(str, components.subject_ids))) != size:
        raise ValueError("unique aligned binary subject components required")
    fitted = aggregate_path_scores(components, information_metric=information_metric, scalar_fast=True)
    score, info, shape = (value[:, 0] if value.ndim == 2 else value[:, 0, 0] for value in (components.scores, components.information, components.biological_shapes))
    reference = None if components.reference_information is None else components.reference_information[:, 0, 0]
    if information_metric == "reference":
        supported = reference > np.maximum(reference, np.finfo(float).tiny) * 1e-10
        keep = supported.copy()
        relative = info[supported] / reference[supported]
        keep[supported] = relative > np.maximum(relative, 1.) * 1e-10
    else:
        keep = info > np.maximum(info, 1.) * 1e-10
    if int(keep.sum()) != fitted["n_subjects"]:
        raise ValueError("diagnostic information rank differs from fitted test")
    precision, effect = np.zeros(size), np.full(size, np.nan)
    effect[keep] = score[keep] / info[keep]
    precision[keep] = 1 / (1 / info[keep] + fitted["biological_variance"] * shape[keep])
    if not np.isfinite(precision).all() or not np.isfinite(effect[keep]).all() or (precision[keep] <= 0).any():
        raise ValueError("finite positive conditional subject precision required")
    shares = precision / precision.sum()
    mean = float(fitted["mean_difference"][0])
    reconstructed = float(shares[keep] @ effect[keep])
    if abs(reconstructed - mean) > 1e-7 * max(abs(mean), float(shares[keep] @ np.abs(effect[keep])), 1.):
        raise ValueError("conditional precision does not reconstruct the fitted mean")
    residual = float(np.sum(precision[keep] * (effect[keep] - mean)**2))
    inflation = max(1., residual / (int(keep.sum()) - 1))
    stable_statistic = mean**2 * precision.sum() / inflation
    dominant = int(np.argmax(shares))
    records = [dict(subject=str(subject), retained=bool(keep[index]), pseudo_effect=float(effect[index]), precision_share=float(shares[index]), weighted_mean_contribution=float(shares[index] * effect[index])) for index, subject in enumerate(components.subject_ids)]
    summary = dict(n_subjects=size, n_informative_subjects=int(keep.sum()), fitted_p_value=float(fitted["p_value"]), fitted_mean=mean, biological_variance=float(fitted["biological_variance"]), maximum_precision_share=float(shares[dominant]), effective_weighted_subjects=float(1 / np.sum(shares**2)), dominant_subject=str(components.subject_ids[dominant]), informative_subject_sign_agreement=float((np.sign(effect[keep]) == np.sign(mean)).mean()) if mean != 0 else np.nan, stable_residual_sum=residual, stable_residual_inflation=inflation, fitted_residual_inflation=float(fitted["residual_inflation"]), stable_fixed_variance_p_value=float(stats.f.sf(stable_statistic, 1, int(keep.sum()) - 1)))
    return records, summary


def binary_score_components_to_proportions(components):
    """Reexpress binary EC scores in the existing proportion coordinate.

    For null inclusion p_u, the archived ILR biological shape is
    B_u=1/(p_u*(1-p_u)). Local coordinate derivative a_u=2/B_u gives
    S'_u=S_u/a_u, I'_u=I_u/a_u**2, B'_u=4/B_u; reference information
    transforms like I_u. Subject-dependent a_u changes the common mean
    estimand from ILR to simplex-tangent difference, not numerical units.
    This cannot restore count fits, support multivariate events or supply
    bounded usages. Original component arrays and reports are unchanged.
    """
    size = len(components.subject_ids)
    if components.score_coordinate != "ilr" or len(components.levels) != 2 or components.scores.shape != (size, 1) or any(value.shape != (size, 1, 1) for value in (components.information, components.biological_shapes)):
        raise ValueError("aligned binary ILR EC components required")
    shape = components.biological_shapes[:, 0, 0]
    if not np.isfinite(shape).all() or (shape < 4 * (1 - 1e-8)).any():
        raise ValueError("binary EC Dirichlet geometry must have shape at least four")
    scale = (2 / shape)[:, None]
    scores = components.scores / scale
    information = (components.information / scale[:, :, None]) / scale[:, :, None]
    reference = None if components.reference_information is None else (components.reference_information / scale[:, :, None]) / scale[:, :, None]
    biological = (4 / shape)[:, None, None]
    if any(not np.isfinite(value).all() for value in (scores, information, biological)) or (reference is not None and not np.isfinite(reference).all()):
        raise ValueError("proportion-coordinate components exceed finite range")
    return PathScoreComponents(scores, information, biological, components.subject_ids.copy(), tuple(components.levels), list(components.null_fits), list(components.reporting_proportions), "proportion", reference)


def binary_information_geometry(information, biological_shape, reference_information):
    """Sign-invariant binary precision geometry, not fitted-score leverage.

    Inputs are aligned M-vectors of profiled information, positive biological
    shape and unprofiled reference information. Apply the existing relative
    information-rank rule, then evaluate precision shares at variance zero,
    infinity and an 81-point fixed log-variance grid around median 1/(I*B).
    The grid maximum is descriptive, not a certified continuous supremum.
    No score, fitted heterogeneity, p-value, effect or LR result is inspected.
    Global coordinate-unit transformations leave the geometry unchanged.
    """
    information, shape, reference = (np.asarray(value, dtype=float) for value in (information, biological_shape, reference_information))
    if information.ndim != 1 or shape.shape != information.shape or reference.shape != information.shape or any(not np.isfinite(value).all() for value in (information, shape, reference)) or (information < 0).any() or (reference < 0).any() or (shape <= 0).any():
        raise ValueError("aligned finite nonnegative information/reference and positive biological shape required")
    supported = reference > np.maximum(reference, np.finfo(float).tiny) * 1e-10
    keep = supported.copy()
    relative = information[supported] / reference[supported]
    if (relative > 1 + 1e-8).any():
        raise ValueError("profiled information exceeds its reference")
    keep[supported] = relative > np.maximum(relative, 1.) * 1e-10
    if keep.sum() < 4:
        raise ValueError("four informative subjects required for binary geometry")
    log_info, log_shape = np.log(information[keep]), np.log(shape[keep])
    log_scale = np.median(-log_info - log_shape)
    log_variance = log_scale + np.linspace(-24., 16., 81)
    log_precision = log_info[None, :] - np.logaddexp(0., log_variance[:, None] + log_info[None, :] + log_shape[None, :])
    shares = special.softmax(np.vstack([log_info, log_precision, -log_shape]), axis=1)
    maximum = float(shares.max())
    return dict(n_informative_subjects=int(keep.sum()), measurement_maximum_share=float(shares[0].max()), biological_limit_maximum_share=float(shares[-1].max()), grid_maximum_precision_share=maximum, grid_minimum_effective_subjects=float(np.min(1 / np.sum(shares**2, axis=1))), geometry_class="balanced" if maximum <= .5 else "intermediate" if maximum <= .9 else "dominated")


def paired_score_reporting(components):
    """Complete-subject arithmetic PSI reporting, separate from score testing.

    For two types and S paths, returns an S-vector of type1-minus-type0 means.
    Missing or failed reports make the complete estimand unavailable rather
    than silently defining a successful-subject reporting subset.
    """
    if len(components.levels) != 2 or len(components.reporting_proportions) != len(components.subject_ids):
        raise ValueError("paired reporting requires aligned two-type subject reports")
    size = components.scores.shape[1] + 1
    differences = []
    for reports in components.reporting_proportions:
        local = dict(reports)
        if len(local) != len(reports):
            raise ValueError("duplicate type in subject reporting")
        if set(local) != {0, 1}:
            continue
        values = np.asarray([local[level] for level in (0, 1)], dtype=float)
        if values.shape != (2, size):
            raise ValueError("paired reporting path dimension differs")
        if np.isfinite(values).all() and (values >= 0).all() and np.allclose(values.sum(axis=1), 1.):
            differences.append(values[1] - values[0])
    complete = len(differences) == len(components.subject_ids) and len(differences) > 0
    return {"effect": np.mean(differences, axis=0) if complete else np.full(size, np.nan), "n_reported_subjects": len(differences), "complete": complete}


def binary_subject_score_records(components, test_id):
    """Serializable binary subject scores, including unsuccessful reports.

    These retain actual fitted scores and both information matrices so rank
    diagnostics do not require another count fit. Missing reports remain NaN.
    This is not a successful-subject subset or a replacement inferential test.
    """
    size = len(components.subject_ids)
    if len(components.levels) != 2 or components.scores.shape != (size, 1) or any(value.shape != (size, 1, 1) for value in (components.information, components.biological_shapes)) or components.reference_information is None or components.reference_information.shape != (size, 1, 1) or len(components.reporting_proportions) != size:
        raise ValueError("aligned binary subject components and reference information required")
    output = []
    for index, subject in enumerate(components.subject_ids):
        reports = dict(components.reporting_proportions[index])
        if len(reports) != len(components.reporting_proportions[index]):
            raise ValueError("duplicate subject reporting level")
        proportions = []
        for level in (0, 1):
            value = np.asarray(reports.get(level, [np.nan, np.nan]), dtype=float)
            if value.shape != (2,):
                raise ValueError("binary report dimension differs")
            valid = np.isfinite(value).all() and (value >= 0).all() and np.isclose(value.sum(), 1.)
            proportions.append(float(value[0]) if valid else np.nan)
        output.append(dict(test_id=str(test_id), subject=str(subject), score=float(components.scores[index, 0]), information=float(components.information[index, 0, 0]), biological_shape=float(components.biological_shapes[index, 0, 0]), reference_information=float(components.reference_information[index, 0, 0]), report_inclusion_a=proportions[0], report_inclusion_b=proportions[1]))
    return output


def mixed_path_score_test(data, path_index, labels, subjects, *, information_metric="absolute", scalar_fast=False, **kwargs):
    """General C-type experimental EC score test, with optional usage reports."""
    return aggregate_path_scores(shared_path_score_components(data, path_index, labels, subjects, **kwargs), information_metric=information_metric, scalar_fast=scalar_fast)


def binary_score_components_from_records(records, *, levels=(0, 1), score_coordinate="ilr"):
    """Restore an entire binary subject score archive without fitting counts.

    Failed usage reports stay missing. Null optimizer objects are not restored,
    and this archive cannot support refitting a different count model.
    """
    records = list(records)
    subjects = [str(row["subject"]) for row in records]
    if not records or len(set(subjects)) != len(subjects) or len(levels) != 2 or score_coordinate not in ("ilr", "proportion"):
        raise ValueError("unique nonempty binary subject archive required")
    arrays = [np.asarray([row[name] for row in records], dtype=float) for name in ("score", "information", "biological_shape", "reference_information")]
    if any(not np.isfinite(value).all() for value in arrays) or any((value < 0).any() for value in arrays[1:]):
        raise ValueError("finite scores and nonnegative information/shapes required")
    reports = []
    for row in records:
        local = []
        for level, name in enumerate(("report_inclusion_a", "report_inclusion_b")):
            value = float(row[name])
            if np.isnan(value):
                local.append((level, np.array([np.nan, np.nan])))
            elif np.isfinite(value) and 0 <= value <= 1:
                local.append((level, np.array([value, 1 - value])))
            else:
                raise ValueError("reported inclusion must be missing or within the simplex")
        reports.append(local)
    size = len(records)
    return PathScoreComponents(arrays[0].reshape(size, 1), arrays[1].reshape(size, 1, 1), arrays[2].reshape(size, 1, 1), np.asarray(subjects), tuple(levels), [], reports, score_coordinate, arrays[3].reshape(size, 1, 1))
