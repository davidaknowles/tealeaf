from types import SimpleNamespace

import numpy as np
import pytest

from tealeaf.sc.local_read_diagnostics import local_read_design_diagnostic, local_read_failure_reason


def test_design_uses_existing_marker_rules_without_fitting_or_changing_counts():
    counts = np.ones((4, 2, 2, 2), dtype=np.int64)
    original = counts.copy()
    assert local_read_design_diagnostic(counts) == dict(design_status='ready', n_paired_marker_subjects=4, design_ready=True)
    np.testing.assert_array_equal(counts, original)
    counts[0] = 0
    assert local_read_design_diagnostic(counts) == dict(design_status='too few paired marker subjects', n_paired_marker_subjects=3, design_ready=False)


def test_globally_constant_primers_do_not_invent_paired_marker_subjects():
    counts = np.ones((4, 2, 2, 2), dtype=np.int64)
    counts[:, 1, :, 1] = 0
    counts[0, 0] = 0
    assert not local_read_design_diagnostic(counts)['design_ready']
    counts[..., 1] = 0
    assert local_read_design_diagnostic(counts)['design_status'] == 'local read evidence with both marker classes required'
    assert not local_read_design_diagnostic(np.zeros_like(counts))['design_ready']


def test_primer_type_confounding_is_not_called_ready():
    counts = np.ones((4, 2, 2, 2), dtype=np.int64)
    counts[:, 0, 1] = 0
    counts[:, 1, 0] = 0
    assert local_read_design_diagnostic(counts)['design_status'] == 'primer and cell-type means are confounded'


@pytest.mark.parametrize('bad', [-1, .5, np.nan])
def test_invalid_count_inputs_are_not_mistaken_for_low_evidence(bad):
    counts = np.ones((4, 2, 2, 2))
    counts[0, 0, 0, 0] = bad
    with pytest.raises(ValueError):
        local_read_design_diagnostic(counts)


def test_failure_flags_are_explicit_and_censored_usable_effects_stay_usable():
    row = SimpleNamespace(converged='False', error=np.nan, quadrature_error=.1, parameter_boundary=False)
    assert local_read_failure_reason(row) == 'doubled-order quadrature mismatch'
    row.quadrature_error = .0005
    row.parameter_boundary = True
    assert local_read_failure_reason(row) == 'nonvalidated parameter boundary'
    row.converged = True
    assert local_read_failure_reason(row) == 'usable'
    row.converged = None
    with pytest.raises(ValueError):
        local_read_failure_reason(row)


def test_original_exception_and_model_specific_tolerances_are_preserved():
    row = SimpleNamespace(converged=False, error='too little evidence', quadrature_error=1.)
    assert local_read_failure_reason(row) == 'too little evidence'
    row.error = ''
    row.quadrature_error = .0005
    assert local_read_failure_reason(row) == 'other numerical unavailability'
    row.model = 'conditional'
    assert local_read_failure_reason(row) == 'doubled-order quadrature mismatch'
