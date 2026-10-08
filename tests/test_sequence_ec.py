import numpy as np
import pytest

from tealeaf.sc.ec_diagnostics import pooled_ec_mixture_diagnostics
from tealeaf.sc.sequence_ec import exact_sequence_words, exact_sequence_ec_kernels, project_sequence_ec_kernels, sequence_ec_background_kernels, replace_ec_kernel_rows_inplace


def test_exact_words_match_naive_strings_without_hash_collisions():
    sequence = "ACGTACNNNNATCGTACGTA"
    complement = str.maketrans("ACGT", "TGCA")
    for length in (1, 4, 10):
        observed = list(exact_sequence_words(sequence, length))
        expected = []
        for start in range(len(sequence) - length + 1):
            word = sequence[start:start + length]
            if "N" not in word:
                expected.append((start, min(word, word.translate(complement)[::-1])))
        keys = [word for _, word in observed]
        assert [start for start, _ in observed] == [start for start, _ in expected]
        for left in range(len(keys)):
            for right in range(len(keys)):
                assert (keys[left] == keys[right]) == (expected[left][1] == expected[right][1])


def test_shared_read_opportunities_are_not_shared_class_counts():
    kernel = exact_sequence_ec_kernels(["A" * 100 + "CGCG", "A" * 100 + "TCTC"], 4, unstranded=False)
    row = kernel.masks.index(3)
    assert kernel.rh[row, 0] == kernel.rh[row, 1] == 97
    assert kernel.rh_starts.tolist() == [101, 101]
    np.testing.assert_allclose(kernel.dt.sum(axis=0), [1., 1.])
    # The retained categorical probabilities exactly match this physical
    # uniform-start working model, while a class-uniform map does not.
    counts = kernel.dt @ np.array([500., 500.])
    correct = pooled_ec_mixture_diagnostics(counts, kernel.dt)
    uniform = np.array([[(mask & (1 << t)) > 0 for t in range(2)] for mask in kernel.masks], dtype=float)
    wrong = pooled_ec_mixture_diagnostics(counts, uniform)
    assert correct["KL"] < 1e-10
    assert wrong["KL_lower_bound"] > .3


def test_terminal_start_window_changes_dt_only_and_preserves_repeated_positions():
    full = exact_sequence_ec_kernels(["AAAACCCC", "AAAAGGGG"], 4, unstranded=False)
    end = exact_sequence_ec_kernels(["AAAACCCC", "AAAAGGGG"], 4, end_window=2, unstranded=False)
    assert full.masks == end.masks
    np.testing.assert_array_equal(full.rh, end.rh)
    assert end.dt_starts.tolist() == [2, 2]
    assert end.dt[end.masks.index(3)].sum() == 0
    np.testing.assert_allclose(end.dt.sum(axis=0), [1., 1.])


def test_unknown_observed_classes_are_not_given_opportunities():
    kernel = exact_sequence_ec_kernels(["AAAA", "AAAA"], 4)
    maps, matched = project_sequence_ec_kernels(kernel, [[1, 0], [1, 1]])
    assert matched.tolist() == [False, True]
    assert maps[0][0].sum() == maps[1][0].sum() == 0
    with pytest.raises(ValueError, match="unique"):
        project_sequence_ec_kernels(kernel, [[1, 1], [1, 1]])


def test_missing_valid_sequence_and_invalid_configuration_are_explicit():
    with pytest.raises(ValueError, match="no valid"):
        exact_sequence_ec_kernels(["NNNN", "AC"], 4)
    with pytest.raises(ValueError, match="positive"):
        exact_sequence_ec_kernels(["AAAA"], 0)
    with pytest.raises(ValueError, match="positive"):
        exact_sequence_ec_kernels(["AAAA"], 4, end_window=0)


def test_class_memberships_do_not_overflow_at_sixty_four_transcripts():
    kernel = exact_sequence_ec_kernels(["AAAA"] * 70, 4)
    maps, matched = project_sequence_ec_kernels(kernel, np.ones((1, 70)))
    assert kernel.masks == ((1 << 70) - 1,)
    assert matched.all()
    np.testing.assert_array_equal(maps[0], np.ones((1, 70)))


def test_background_retains_unknown_counts_without_changing_capture_units():
    kernel = exact_sequence_ec_kernels(["AAAA", "AAAA"], 4)
    membership = np.array([[1, 0], [0, 1], [1, 1]])
    maps, matched = sequence_ec_background_kernels(kernel, membership, .01)
    assert matched.tolist() == [False, False, True]
    assert (maps[0][membership > 0] > 0).all()
    np.testing.assert_allclose(maps[0].sum(axis=0), [1, 1])
    np.testing.assert_allclose(maps[1].sum(axis=0), kernel.rh_starts)
    # Censoring correct-read classes must not be renormalized away.
    censored, _ = sequence_ec_background_kernels(kernel, [[1, 0], [0, 1]], .01)
    np.testing.assert_allclose(censored[0].sum(axis=0), [.01, .01])


def test_background_cannot_invent_capture_for_sequence_without_valid_starts():
    kernel = exact_sequence_ec_kernels(["AAAA", "NNNN"], 4)
    with pytest.raises(ValueError, match="capture exposure"):
        sequence_ec_background_kernels(kernel, [[1, 0], [0, 1]], .01)


def test_sparse_kernel_injection_preserves_global_support_and_arbitrary_order():
    from scipy.sparse import csr_matrix

    original = csr_matrix([[1., 0, 2], [0, 4, 0], [3, 0, 0]])
    designs = (original.copy(), original.copy())
    kernel = np.array([[0, 7.], [8., 9.]])
    replace_ec_kernel_rows_inplace(designs, [2, 0], [2, 0], (kernel, kernel * 10))
    np.testing.assert_array_equal(designs[0].toarray(), [[9, 0, 8], [0, 4, 0], [7, 0, 0]])
    np.testing.assert_array_equal(designs[1].toarray(), [[90, 0, 80], [0, 4, 0], [70, 0, 0]])
    for value in designs:
        np.testing.assert_array_equal(value.indptr, original.indptr)
        np.testing.assert_array_equal(value.indices, original.indices)
    np.testing.assert_array_equal(original.toarray(), [[1, 0, 2], [0, 4, 0], [3, 0, 0]])


@pytest.mark.parametrize("bad", [np.array([[1., 0.]]), np.array([[1., -1.]]), np.array([[1., np.nan]]), np.array([[1.]]), np.array([[1e-100, 2.]])])
def test_injection_validates_all_primers_before_mutation(bad):
    from scipy.sparse import csr_matrix

    designs = (csr_matrix([[1., 2.]], dtype=np.float32), csr_matrix([[3., 4.]], dtype=np.float32))
    before = [value.data.copy() for value in designs]
    with pytest.raises(ValueError):
        replace_ec_kernel_rows_inplace(designs, [0], [0, 1], (np.array([[7., 8.]]), bad))
    for value, expected in zip(designs, before):
        np.testing.assert_array_equal(value.data, expected)


def test_injection_rejects_foreign_transcripts_and_duplicate_rows():
    from scipy.sparse import csr_matrix

    design = csr_matrix([[1., 2.]])
    with pytest.raises(ValueError, match="outside declared gene"):
        replace_ec_kernel_rows_inplace((design,), [0], [0], (np.array([[3.]]),))
    with pytest.raises(ValueError, match="unique"):
        replace_ec_kernel_rows_inplace((design,), [0, 0], [0, 1], (np.ones((2, 2)),))
