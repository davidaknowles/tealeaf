"""Experimental dispersion corrections for effective-path DM regression.

This module reuses the DM mean optimizer. Fractional effective counts remain
an approximation; neither a Cox-Reid correction nor a chi-square tail supplies
calibration by itself. Production defaults are not changed.
"""

import numpy as np
from scipy.special import digamma, polygamma, softmax
from scipy.stats import chi2, f

from . import differential


def dm_coefficient_information(counts, design, coefficients, concentration):
    """Observed negative Hessian for B, in row-major P by (S-1) order.

    Counts have shape N by S, design N by P, coefficients P by (S-1),
    and the returned information P(S-1) by P(S-1). Concentration is positive.
    The softmax second derivative is included, even away from the optimum.
    Positive determinant alone does not imply positive-definite information.
    """
    counts, design, coefficients = (np.asarray(value, dtype=float) for value in (counts, design, coefficients))
    concentration = float(concentration)
    if counts.ndim != 2 or counts.shape[1] < 2 or design.ndim != 2 or design.shape[0] != len(counts) or coefficients.shape != (design.shape[1], counts.shape[1] - 1):
        raise ValueError("incompatible counts, design, and coefficient dimensions")
    if not np.isfinite(concentration) or concentration <= 0 or not all(np.isfinite(value).all() for value in (counts, design, coefficients)) or np.any(counts < 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("finite positive concentration and valid counts/design required")
    means = softmax(np.column_stack([design @ coefficients, np.zeros(len(counts))]), axis=1)
    alpha = concentration * means
    score_alpha = concentration * (digamma(alpha + counts) - digamma(alpha))
    curvature_alpha = concentration ** 2 * (polygamma(1, alpha + counts) - polygamma(1, alpha))
    gradient = means * (score_alpha - (means * score_alpha).sum(axis=1, keepdims=True))
    identity = np.eye(counts.shape[1])
    jacobian = means[:, :, None] * identity - means[:, :, None] * means[:, None, :]
    hessian = gradient[:, :, None] * identity - gradient[:, :, None] * means[:, None, :] - means[:, :, None] * gradient[:, None, :]
    hessian += np.einsum("nab,nb,nbc->nac", jacobian, curvature_alpha, jacobian)
    dimensions = coefficients.size
    information = -np.einsum("ni,nab,nj->iajb", design, hessian[:, :-1, :-1], design).reshape(dimensions, dimensions)
    return (information + information.T) / 2


def _validate_designs(counts, null_design, alternative_design):
    counts, null_design, alternative_design = (np.asarray(value, dtype=float) for value in (counts, null_design, alternative_design))
    if counts.ndim != 2 or counts.shape[1] < 2 or not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts.sum(axis=1) <= 0):
        raise ValueError("valid positive-depth counts with at least two paths required")
    for design in (null_design, alternative_design):
        if design.ndim != 2 or design.shape[0] != len(counts) or design.shape[1] < 1 or not np.isfinite(design).all() or np.linalg.matrix_rank(design) != design.shape[1]:
            raise ValueError("finite full-rank designs aligned with counts required")
    if alternative_design.shape[1] <= null_design.shape[1] or not np.allclose(alternative_design[:, :null_design.shape[1]], null_design):
        raise ValueError("alternative design must begin with and extend the null design")
    if len(counts) <= alternative_design.shape[1]:
        raise ValueError("regression needs residual observations")
    return counts, null_design, alternative_design


def _fit_at_concentration(counts, design, concentration, initial=None, max_iter=250):
    fitted = differential._dirichlet_multinomial_fit(counts, design, initial=initial, max_iter=max_iter, fixed_concentration=concentration)
    if not fitted["converged"]:
        fitted = differential._dirichlet_multinomial_fit(counts, design, initial=fitted["parameters"], max_iter=4 * int(max_iter), fixed_concentration=concentration)
    return fitted


def cox_reid_profile(counts, design, *, concentrations=None, max_iter=250):
    """Select precision on a prespecified full-design adjusted-profile grid.

    At each kappa, refit B and maximize ell(Bhat, kappa)-log|I_BB|/2.
    The default 21-point grid is 100 times 2**linspace(-10,10,21).
    A profile with any failed/non-PD/boundary mean fit is rejected, not silently
    maximized over whichever grid points succeeded. No cross-block moderation
    or interpolation is performed. Grid-endpoint selections remain finite and
    are exported explicitly, not labeled exact multinomial limits.
    """
    counts, design = np.asarray(counts, dtype=float), np.asarray(design, dtype=float)
    if counts.ndim != 2 or counts.shape[1] < 2 or design.ndim != 2:
        raise ValueError("counts and design must be matrices with at least two paths")
    # Reuse the public information validation before starting optimization.
    dm_coefficient_information(counts, design, np.zeros((design.shape[1], counts.shape[1] - 1)), 1.)
    if np.linalg.matrix_rank(design) != design.shape[1] or len(counts) <= design.shape[1]:
        raise ValueError("full-rank design with residual observations required")
    grid = 100. * 2. ** np.linspace(-10, 10, 21) if concentrations is None else np.asarray(concentrations, dtype=float)
    if grid.ndim != 1 or len(grid) < 3 or not np.isfinite(grid).all() or np.any(grid <= 0) or np.any(np.diff(grid) <= 0):
        raise ValueError("at least three ordered positive concentrations required")
    initial = differential._multinomial_fit(counts, design, max_iter=max_iter)["parameters"]
    fits, adjusted, errors = [], [], []
    for concentration in grid:
        fitted = _fit_at_concentration(counts, design, concentration, initial, max_iter)
        fits.append(fitted)
        try:
            if not fitted["converged"] or not np.isfinite(fitted["objective"]) or np.any(np.abs(fitted["parameters"]) >= 19.99):
                raise ValueError("failed or coefficient-boundary mean fit")
            coefficients = fitted["parameters"].reshape(design.shape[1], counts.shape[1] - 1)
            information = dm_coefficient_information(counts, design, coefficients, concentration)
            factor = np.linalg.cholesky(information)
            adjusted.append(-fitted["objective"] - np.log(np.diag(factor)).sum())
            initial = fitted["parameters"]
            errors.append("")
        except (ValueError, np.linalg.LinAlgError) as exception:
            adjusted.append(np.nan)
            errors.append(str(exception))
    valid = np.isfinite(adjusted)
    if not valid.all():
        raise ValueError(f"Cox-Reid profile has {int((~valid).sum())}/{len(grid)} invalid points, {next(error for error in errors if error)}")
    selected = int(np.argmax(adjusted))
    return {"concentration": float(grid[selected]), "fit": fits[selected], "concentrations": grid, "adjusted_profile": np.asarray(adjusted), "profile_index": selected, "profile_boundary": selected in (0, len(grid) - 1)}


def corrected_dm_test(counts, null_design, alternative_design, *, concentration=None, concentrations=None, max_iter=250):
    """Fixed-precision DM LRT after full-design Cox-Reid selection.

    If concentration is provided, it is fixed instead of selected, intended
    for prespecified/oracle diagnostics. Otherwise precision is estimated
    under the full design and held identical under the null and alternative.
    The correction is NOT subtracted from the LR statistic. Both chi-square
    and a heuristic residual-df F tail are returned; neither is exact here.
    A label-permutation caller must reselect precision for each permutation.
    """
    counts, null_design, alternative_design = _validate_designs(counts, null_design, alternative_design)
    if concentration is None:
        profile = cox_reid_profile(counts, alternative_design, concentrations=concentrations, max_iter=max_iter)
        concentration, alternative = profile["concentration"], profile["fit"]
    else:
        concentration = float(concentration)
        if not np.isfinite(concentration) or concentration <= 0:
            raise ValueError("finite positive concentration required")
        profile = {"profile_index": -1, "profile_boundary": False}
        alternative = _fit_at_concentration(counts, alternative_design, concentration, max_iter=max_iter)
    coefficients = alternative["parameters"].reshape(alternative_design.shape[1], counts.shape[1] - 1)
    null = _fit_at_concentration(counts, null_design, concentration, initial=coefficients[:null_design.shape[1]].ravel(), max_iter=max_iter)
    difference = 2 * (null["objective"] - alternative["objective"])
    converged = bool(null["converged"] and alternative["converged"] and np.isfinite(difference) and difference >= -1e-6)
    statistic = max(float(difference), 0.) if converged else 0.
    degrees = (alternative_design.shape[1] - null_design.shape[1]) * (counts.shape[1] - 1)
    residual_degrees = (len(counts) - alternative_design.shape[1]) * (counts.shape[1] - 1)
    return {"model": "cox_reid_dm" if profile["profile_index"] >= 0 else "fixed_precision_dm", "statistic": statistic, "degrees_of_freedom": degrees, "residual_degrees_of_freedom": residual_degrees, "p_value": float(chi2.sf(statistic, degrees)) if converged else 1., "f_p_value": float(f.sf(statistic / degrees, degrees, residual_degrees)) if converged else 1., "null_concentration": concentration, "alternative_concentration": concentration, "null_converged": null["converged"], "alternative_converged": alternative["converged"], "null_objective": null["objective"], "null_coefficients": null["parameters"].reshape(null_design.shape[1], counts.shape[1] - 1), "alternative_coefficients": coefficients, "profile_index": profile["profile_index"], "profile_boundary": profile["profile_boundary"]}
