"""Residual-based multivariate regression diagnostics for blocked contrasts."""

import numpy as np
import scipy.linalg
import scipy.stats


def regression_omnibus(values, design, tested_columns):
    """Return trace F, Pillai, and maximum-coordinate F diagnostics.

    Values have shape N by D, design N by P, tested columns length Q.
    These analytic reference distributions assume independent homoscedastic
    residuals. Repeated biological observations require design-matched
    permutation calibration; these are not replacements for that calibration.
    Trace F assumes spherical response covariance; Pillai accounts for full
    response covariance. Maximum-coordinate F is an unadjusted ordering
    statistic and its analytic minimum p-value is not a global p-value.
    """
    values = np.asarray(values, dtype=float)
    design = np.asarray(design, dtype=float)
    tested_columns = np.asarray(tested_columns, dtype=int)
    if values.ndim != 2 or design.ndim != 2 or design.shape[0] != len(values) or not np.isfinite(values).all() or not np.isfinite(design).all():
        raise ValueError("finite aligned regression matrices required")
    if not len(tested_columns) or len(np.unique(tested_columns)) != len(tested_columns) or np.any(tested_columns < 0) or np.any(tested_columns >= design.shape[1]):
        raise ValueError("valid nonempty tested columns required")
    null_design = np.delete(design, tested_columns, axis=1)
    rank = np.linalg.matrix_rank(design)
    q = rank - np.linalg.matrix_rank(null_design)
    residual_df = len(values) - rank
    if rank != design.shape[1] or q <= 0 or residual_df <= 0:
        raise ValueError("full-rank regression and residual degrees required")
    full_residual = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
    null_residual = values - null_design @ np.linalg.lstsq(null_design, values, rcond=None)[0]
    error = full_residual.T @ full_residual
    hypothesis = null_residual.T @ null_residual - error
    dimension = values.shape[1]
    improvement = max(float(np.trace(hypothesis)), 0.)
    trace_stat = (improvement / (q * dimension)) / max(float(np.trace(error)) / (residual_df * dimension), np.finfo(float).tiny)
    coordinate_stat = np.maximum(np.diag(hypothesis), 0) / q / np.maximum(np.diag(error) / residual_df, np.finfo(float).tiny)
    total = error + hypothesis
    if np.linalg.matrix_rank(total) != dimension:
        raise ValueError("response covariance is rank deficient")
    pillai = float(np.trace(scipy.linalg.solve(total, hypothesis, assume_a="sym")))
    s = min(dimension, q)
    numerator_df = s * (abs(dimension - q) + s)
    denominator_df = s * (residual_df - dimension + s)
    if denominator_df <= 0:
        raise ValueError("insufficient residual degrees for Pillai reference")
    pillai = np.clip(pillai, 0., s - 1e-12)
    pillai_stat = denominator_df / numerator_df * pillai / (s - pillai)
    return {
        "trace F": {"statistic": trace_stat, "p_value": float(scipy.stats.f.sf(trace_stat, q * dimension, residual_df * dimension)), "degrees_of_freedom": q * dimension},
        "Pillai": {"statistic": pillai_stat, "p_value": float(scipy.stats.f.sf(pillai_stat, numerator_df, denominator_df)), "degrees_of_freedom": numerator_df},
        "maximum-coordinate F": {"statistic": float(coordinate_stat.max()), "p_value": float(scipy.stats.f.sf(coordinate_stat.max(), q, residual_df)), "degrees_of_freedom": q},
    }


def cluster_max_f(values, design, tested_columns, subjects):
    """CR2-adjusted maximum-coordinate cluster Wald ordering statistic.

    Arrays have dimensions N by D, N by P, Q, and N respectively. Subject
    fixed effects give singular cluster residual-maker blocks; CR2 uses the
    inverse square root only on their positive-eigenvalue residual subspace.
    The minimum marginal analytic p-value is not a global p-value. Use a
    subject-level wild residual bootstrap for global calibration. The F
    reference with M-Q denominator degrees is a small-sample approximation.
    """
    values = np.asarray(values, dtype=float)
    design = np.asarray(design, dtype=float)
    subjects = np.asarray(subjects)
    tested_columns = np.asarray(tested_columns, dtype=int)
    if values.ndim != 2 or design.ndim != 2 or design.shape[0] != len(values) or subjects.shape != (len(values),):
        raise ValueError("cluster regression inputs must align")
    groups = np.unique(subjects)
    q = len(tested_columns)
    if len(groups) <= q or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("insufficient clusters or rank-deficient design")
    inverse = np.linalg.inv(design.T @ design)
    coefficients = inverse @ design.T @ values
    residuals = values - design @ coefficients
    scores = []
    for subject in groups:
        selected = subjects == subject
        local = design[selected]
        maker = np.eye(np.sum(selected)) - local @ inverse @ local.T
        eigenvalues, eigenvectors = np.linalg.eigh((maker + maker.T) / 2)
        root = np.zeros_like(eigenvalues)
        positive = eigenvalues > 1e-8
        root[positive] = 1 / np.sqrt(eigenvalues[positive])
        adjusted = (eigenvectors * root) @ eigenvectors.T @ residuals[selected]
        scores.append((inverse @ local.T @ adjusted)[tested_columns])
    scores = np.asarray(scores)
    statistics = []
    for coordinate in range(values.shape[1]):
        covariance = scores[:, :, coordinate].T @ scores[:, :, coordinate]
        if np.linalg.matrix_rank(covariance) != q:
            raise ValueError("unidentifiable cluster-tested coordinate")
        effect = coefficients[tested_columns, coordinate]
        statistics.append(float(effect @ np.linalg.solve(covariance, effect) / q))
    maximum = max(statistics)
    return {"statistic": maximum, "p_value": float(scipy.stats.f.sf(maximum, q, len(groups) - q)), "degrees_of_freedom": q}
