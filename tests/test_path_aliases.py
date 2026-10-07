import numpy as np
import pytest

from tealeaf.sc.path_aliases import exact_path_column_groups


def test_aliases_require_every_primer_and_the_same_path():
    first = np.array([[1., 1., 1., 2., 2.], [0., 0., 0., 0., 0.]])
    second = np.array([[3., 3., 3., 4., 4.1]])
    groups = exact_path_column_groups((first, second), [0, 0, 1, -1, -1])
    assert [list(group) for group in groups] == [[0, 1], [2], [3], [4]]
    assert len(exact_path_column_groups((first, second * 2), [0, 0, 1, -1, -1])) == 4


def test_near_duplicates_are_not_aliases():
    groups = exact_path_column_groups((np.array([[0., -0., 1e-15]]),), [0, 0, 0])
    assert [list(group) for group in groups] == [[0, 1], [2]]
    with pytest.raises(ValueError, match="aligned"):
        exact_path_column_groups((np.ones((3, 2)),), [0, 1, 2])
