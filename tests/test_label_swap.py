from pathlib import Path

from tealeaf.sc.local_path_test import read_swapped_pairs, swap_entry_labels


def test_swap_entry_labels(tmp_path: Path):
    table = tmp_path / "swapped.tsv"
    table.write_text("level_a\tlevel_b\tsubject\nB\tA\ts1\n")
    swapped = read_swapped_pairs(table)
    entries = [("s1", "A", "poly(dT)", 1, 2.), ("s1", "B", "poly(dT)", 2, 3.), ("s2", "A", "poly(dT)", 1, 4.)]
    out = swap_entry_labels(entries, ("A", "B"), swapped)
    assert out == [("s1", "B", "poly(dT)", 1, 2.), ("s1", "A", "poly(dT)", 2, 3.), ("s2", "A", "poly(dT)", 1, 4.)]
    assert swap_entry_labels(entries, ("A", "C"), swapped) == entries
