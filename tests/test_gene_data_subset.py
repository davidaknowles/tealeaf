import numpy as np
import pytest
from scipy import sparse

from tealeaf.sc.ec_glmm import subset_gene_data


@pytest.mark.parametrize("drop_zero", [True, False])
def test_column_first_selection_matches_old_row_first_bit_for_bit(drop_zero):
    rng = np.random.default_rng(192)
    counts = [sparse.csr_matrix(rng.integers(0, 50, (20, 300))) for _ in range(2)]
    maps = [sparse.csr_matrix(rng.integers(0, 2, (300, 30)).astype(float)) for _ in range(2)]
    rows, ecs, transcripts = np.array([9, 4, 17, 0]), np.array([200, 9, 5, 41]), np.array([19, 2, 7])
    for value in counts:
        value[rows[1], ecs] = 0
    fixed, clusters = np.ones((4, 1)), np.array(["b", "a", "b", "a"])
    old = subset_gene_data(tuple(value[rows] for value in counts), maps, transcripts, ecs, fixed, clusters, drop_zero=drop_zero)
    new = subset_gene_data(tuple(value.tocsc() for value in counts), maps, transcripts, ecs, fixed, clusters, rows=rows, drop_zero=drop_zero)
    for first, second in zip(old[0].counts + old[0].compatibility, new[0].counts + new[0].compatibility):
        np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(old[1], new[1])
    np.testing.assert_array_equal(old[2], new[2])
    np.testing.assert_array_equal(old[0].design, new[0].design)
    np.testing.assert_array_equal(old[0].clusters, new[0].clusters)
