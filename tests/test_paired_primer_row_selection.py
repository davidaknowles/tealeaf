import numpy as np
import pytest

from tealeaf.sc.glm_cv import paired_primer_row_selection


def test_half_cell_selection_preserves_pair_order_and_requires_both_thresholds():
    pairs = [("one", "p", "h"), ("low", "lowp", "lowh"), ("missing", "absent", "h2"), ("threshold", "p2", "h2b")]
    barcodes = ["h", "p", "lowp", "lowh", "h2", "h2b", "p2", "other"]
    complete, groups = paired_primer_row_selection(barcodes, [501, 500, 499, 900, 700, 500, 500, 800], pairs)
    assert complete == [("one", 1, 0), ("threshold", 6, 5)]
    np.testing.assert_array_equal(groups, [1, 0, 0, 1, 1, 1, 0, -1])


@pytest.mark.parametrize("barcodes,totals,pairs", [(["p", "p"], [500, 500], [("one", "p", "h")]), (["p", "h"], [500], [("one", "p", "h")]), (["p", "h"], [500, -1], [("one", "p", "h")]), (["p", "h"], [500, np.nan], [("one", "p", "h")]), (["p", "h"], [500, 500], [("one", "p", "h"), ("two", "p", "missing")]), (["p", "h"], [499, 500], [("one", "p", "h")])])
def test_invalid_or_unretained_half_cell_families_fail(barcodes, totals, pairs):
    with pytest.raises(ValueError):
        paired_primer_row_selection(barcodes, totals, pairs)
