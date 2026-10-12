import numpy as np

from tealeaf.sc.local_path_reads import block_structure, compatibility_mask, group_paths, path_read_opportunities, path_structures, same_junction_groups, variable_junctions

SKIP = [(100, 200), (500, 600)]
INCLUDE = [(100, 200), (300, 400), (500, 600)]


def classify(blocks, junctions=()):
    return compatibility_mask(blocks, list(junctions), path_structures([SKIP, INCLUDE]))


def test_read_classes():
    assert classify([(150, 200), (300, 350)], [(200, 300)]) == 2
    assert classify([(150, 200), (500, 550)], [(200, 500)]) == 1
    assert classify([(120, 180)]) == 3
    assert classify([(320, 380)]) == 2
    assert classify([(250, 280)]) == 0
    assert classify([(180, 220)]) == 0
    assert classify([(700, 750)]) == 0


def test_window_clipping_ignores_outside_junctions():
    assert classify([(50, 90), (120, 170)], [(90, 120)]) == 3
    assert classify([(550, 600), (650, 700)], [(600, 650)]) == 3


def test_opportunities_cover_each_path_once():
    length = 50
    opportunities = path_read_opportunities([SKIP, INCLUDE], length)
    totals = sum(opportunities.values())
    assert np.allclose(totals, [200 + length - 1, 300 + length - 1])
    assert opportunities[1][1] == 0 and opportunities[2][0] == 0
    # skip-junction reads: starts covering positions 99 and 100 of the skip chain
    assert opportunities[1][0] == length - 1
    assert opportunities[3][0] == opportunities[3][1]


def test_pooled_path_shares_recovers_effective_length_composition():
    from tealeaf.sc.local_path_reads import pooled_path_shares
    opportunities = path_read_opportunities([SKIP, INCLUDE], 50)
    lengths = sum(opportunities.values())
    truth = np.array([.3, .7])
    expected = {mask: 1e6 * (vector * truth / lengths).sum() / truth.sum() for mask, vector in opportunities.items()}
    assert np.allclose(pooled_path_shares(expected, opportunities, pseudocount=0.), truth, atol=1e-6)


def test_projected_opportunities_equal_recomputed():
    from tealeaf.sc.local_path_reads import project_opportunities
    ALT = [(100, 200), (300, 350), (500, 600)]
    PRE = [(100, 600)]
    full = path_read_opportunities([SKIP, INCLUDE, ALT, PRE], 50)
    for kept, precursor in (([0, 1], True), ([0, 2], False), ([1, 2], True)):
        direct = path_read_opportunities([[SKIP, INCLUDE, ALT][i] for i in kept] + ([PRE] if precursor else []), 50)
        projected = project_opportunities(full, kept, 3, precursor)
        assert set(direct) == set(projected)
        for mask in direct:
            assert np.allclose(direct[mask], projected[mask])


def test_block_structure():
    cassette = [[[0, 100], [200, 300], [400, 500]], [[0, 100], [400, 500]]]
    assert block_structure(cassette) == {"structure": "internal", "same_junctions": False}
    first_exons = [[[0, 100], [400, 500]], [[200, 300], [400, 500]]]
    assert block_structure(first_exons) == {"structure": "terminal", "same_junctions": False}
    ends = [[[0, 100], [400, 500]], [[0, 100], [400, 800]]]
    assert block_structure(ends) == {"structure": "terminal", "same_junctions": True}


def test_junction_only_opportunities_and_groups():
    paths = [[(0, 100), (400, 500)], [(200, 300), (400, 500)], [(250, 300), (400, 500)]]
    variable = variable_junctions(paths)
    assert variable == {(100, 400), (300, 400)}
    groups = same_junction_groups(paths)
    assert groups == [[0], [1, 2]]
    opportunities = path_read_opportunities(paths, 151, variable)
    full = path_read_opportunities(paths, 151)
    assert all(np.all(opportunities[mask] <= full[mask]) for mask in opportunities)
    assert opportunities[1][0] == 150
    entries = [("s", "a", "dT", 1, 3.), ("s", "a", "dT", 6, 2.), ("s", "a", "dT", 2, 1.)]
    grouped_entries, grouped, kept = group_paths(entries, opportunities, groups, 3)
    assert kept == groups
    assert [entry[3] for entry in grouped_entries] == [1, 2, 2]
    np.testing.assert_allclose(sum(grouped.values()), [150., 150.])
