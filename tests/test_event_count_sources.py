import numpy as np
import pytest
from scipy import sparse

from extra_scripts.audit_event_count_sources import binary_transcript_support, count_fingerprint


def test_binary_support_is_not_assignment_and_filters_ECs_before_transcripts():
    support, kept = binary_transcript_support(sparse.csr_matrix([[.3, .7, 0], [0, 0, 1]]), [100, 2])
    np.testing.assert_array_equal(support, [100, 100, 2])
    np.testing.assert_array_equal(kept, [100, 100, 0])


def test_binary_support_requires_matching_nonnegative_count_totals():
    membership = sparse.eye(2, format="csr")
    with pytest.raises(ValueError):
        binary_transcript_support(membership, [1])
    with pytest.raises(ValueError):
        binary_transcript_support(membership, [-1, 0])


def test_count_fingerprint_checks_values_and_row_order_not_just_totals():
    values = sparse.csr_matrix([[1., 2], [0, 3]])
    assert count_fingerprint(values) == count_fingerprint(values.copy())
    assert count_fingerprint(values) != count_fingerprint(values[[1, 0]])
    assert count_fingerprint(values) != count_fingerprint(sparse.csr_matrix([[2., 1], [0, 3]]))
