import numpy as np

from tealeaf.sc.junction_paths import block_path_exons, path_junction_map, junction_identifiable, junction_path_usage


def cassette_paths():
    left, right = (100, 200), (500, 600)
    return [block_path_exons(left, [], right), block_path_exons(left, [(300, 400)], right)]


def test_cassette_map_and_classical_psi():
    keys, weights = path_junction_map(cassette_paths(), "+")
    assert keys == [(200, 300), (200, 500), (400, 500)]
    assert np.allclose(weights, [[0, 1], [1, 0], [0, 1]])
    assert junction_identifiable(weights)
    psi = junction_path_usage([30., 20., 50.], weights, concentration=0.)
    inclusion = (30 + 50) / 2
    assert np.allclose(psi, [20 / (20 + inclusion), inclusion / (20 + inclusion)], atol=1e-6)


def test_minus_strand_keys_are_genomic():
    keys, _ = path_junction_map(cassette_paths(), "-")
    assert keys == [(200, 300), (200, 500), (400, 500)]


def test_retained_intron_path_is_not_identifiable():
    left, right = (100, 200), (500, 600)
    paths = [block_path_exons(left, [], right), [(100, 600)]]
    _, weights = path_junction_map(paths, "+")
    assert not junction_identifiable(weights)


def test_prior_keeps_unobserved_path_positive():
    _, weights = path_junction_map(cassette_paths(), "+")
    psi = junction_path_usage([0., 10., 0.], weights, concentration=1.)
    assert psi[1] > 0 and psi[0] > psi[1]
    assert np.isclose(psi.sum(), 1.)
