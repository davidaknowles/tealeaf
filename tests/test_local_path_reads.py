import numpy as np

from tealeaf.sc.local_path_reads import compatibility_mask, path_read_opportunities, path_structures

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
