"""Evidence eligibility and numerical availability, never a test filter."""

import numpy as np

from .local_read_mixed import LocalReadMixed


def local_read_design_diagnostic(counts, min_subjects=4):
    """Classify original M x P x two-type x two-class count evidence.

    Apply the existing model's structural checks, without fitting, changing
    counts, adding pseudocounts or using p-values. A ready design is not a
    guarantee of a usable fit. Keep unavailable requests in every denominator.
    """
    if not isinstance(min_subjects, int) or isinstance(min_subjects, bool) or min_subjects < 1:
        raise ValueError('positive integer minimum subject count required')
    try:
        likelihood = LocalReadMixed(counts)
    except ValueError as exc:
        reason = str(exc)
        allowed = ('local read evidence with both marker classes required', 'primer and cell-type means are confounded')
        if reason not in allowed:
            raise
        return dict(design_status=reason, n_paired_marker_subjects=0, design_ready=False)
    ready = likelihood.n_paired_subjects >= min_subjects
    return dict(design_status='ready' if ready else 'too few paired marker subjects', n_paired_marker_subjects=likelihood.n_paired_subjects, design_ready=ready)


def local_read_failure_reason(row):
    """Explain unavailable fits separately from original design eligibility."""
    flag = str(row.converged).lower()
    if flag not in ('true', 'false'):
        raise ValueError('explicit boolean fitting availability required')
    if flag == 'true':
        return 'usable'
    error = getattr(row, 'error', '')
    if error is not None and str(error) and str(error).lower() not in ('nan', 'none'):
        return str(error)
    tolerance = 1e-4 if getattr(row, 'model', 'unconditional') == 'conditional' else 1e-3
    if getattr(row, 'quadrature_error', np.nan) > tolerance:
        return 'doubled-order quadrature mismatch'
    if str(getattr(row, 'parameter_boundary', False)).lower() == 'true':
        return 'nonvalidated parameter boundary'
    return 'other numerical unavailability'
