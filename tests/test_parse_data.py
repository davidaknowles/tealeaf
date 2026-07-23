import pytest

from tealeaf.data.parse import parse_primer_pairs


def test_parse_primer_pairs_preserves_batch_and_ligation_prefix():
    rt = ["AAAA", "CCCC", "GGGG", "TTTT"]
    cells = ["b1:ACGTAAAA", "b1:ACGTGGGG", "b2:TGCAAAAA"]
    assert parse_primer_pairs(cells, rt) == [
        ("b1:ACGTAAAA", "b1:ACGTAAAA", "b1:ACGTGGGG"),
        ("b2:TGCAAAAA", "b2:TGCAAAAA", "b2:TGCAGGGG"),
    ]


def test_parse_primer_pairs_rejects_ambiguous_variable_suffixes():
    with pytest.raises(ValueError, match="ambiguous"):
        parse_primer_pairs(["prefixAAAA"], ["AAAA", "CCCC", "AA", "GG"])
