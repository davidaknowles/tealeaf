import numpy as np
from scipy.special import expit

from tealeaf.sc.local_path_reads import path_read_opportunities
from tealeaf.sc.path_marginal import BinaryECPathLikelihood
from tealeaf.sc.path_marginal_grid import estimate_primer_offset
from tealeaf.sc.position_bias import a_run_starts, primed_path_read_opportunities, priming_kernel, priming_weights


def test_a_run_starts():
    # the start is the first 10-base window with at least 8 A
    assert a_run_starts("CCCCAAAAAAAACCCC").tolist() == [2]
    assert a_run_starts("CCCCAAAAAAACCCCCAAAAAAAAAA").tolist() == [14]
    assert a_run_starts("ACGT" * 10).tolist() == []


def test_priming_kernel_window_and_floor():
    delta = np.arange(-1000, 1001)
    enrichment = np.where((delta >= 0) & (delta < 100), 4., .9)
    kernel_delta, kappa = priming_kernel(enrichment, delta, smooth=1, low=-300, high=400)
    assert kernel_delta.min() == -300 and kernel_delta.max() == 400
    assert np.all(kappa >= 0)
    assert np.allclose(kappa[(kernel_delta >= 0) & (kernel_delta < 100)], 3.)
    assert np.allclose(kappa[kernel_delta >= 100], 0.)


def test_priming_weights_place_excess_upstream_of_a_run():
    sequence = "C" * 200 + "A" * 10 + "C" * 200
    kernel = (np.arange(0, 50), np.ones(50))
    weights = priming_weights(sequence, np.arange(0, 410), kernel)
    start = int(a_run_starts(sequence)[0])
    assert start == 198
    assert np.allclose(weights[start - 49:start + 1], 2.)  # read ends 0..49 bases upstream of the A-run start
    assert np.allclose(weights[:start - 49], 1.) and np.allclose(weights[start + 1:], 1.)


def test_primed_opportunities_flat_primer_matches_uniform_and_primed_favours_a_run_path():
    paths = [[(0, 300), (1000, 1100), (2000, 2300)], [(0, 300), (2000, 2300)]]
    exon = "A" * 12 + "C" * 88
    sequences = ["G" * 300 + exon + "G" * 300, "G" * 600]
    kernel = (np.arange(0, 100), np.full(100, 2.))
    for strand in ("+", "-"):
        result = primed_path_read_opportunities(paths, 50, strand, sequences if strand == "+" else [s[::-1] for s in sequences], kernel)
        uniform = path_read_opportunities(paths, 50)
        assert set(result["random hexamer"]) == set(uniform)
        for mask, vector in uniform.items():
            assert np.allclose(result["random hexamer"][mask], vector)
        primed = {mask: np.asarray(vector) for mask, vector in result["poly(dT)"].items()}
        lengths = sum(primed.values())
        assert lengths[0] > lengths[1]  # the path with the A-run gains poly(dT) opportunities
        flat = sum(np.asarray(vector) for vector in result["random hexamer"].values())
        assert lengths[0] / flat[0] > lengths[1] / flat[1]


def test_primer_offset_prior_shrinks_unidentified_offset():
    rng = np.random.default_rng(3)
    rows = []
    for subject in range(10):
        # poly(dT) sees only class 1, hexamer only class 2: the offset is unbounded without a prior
        rows.append(np.r_[rng.integers(1, 3), 0, 0, rng.integers(1, 3)])
    components = np.array([[0., .5, 0.], [0., 0., .4]])
    counts = np.array(rows, dtype=float)
    like = BinaryECPathLikelihood((counts[:, :2], counts[:, 2:]), (components, components), np.array([f"s{i}" for i in range(10)]), np.array(["a"] * 10), np.ones((10, 2)), np.zeros(10))
    free = estimate_primer_offset(like)
    shrunk = estimate_primer_offset(like, prior_sd=1.)
    assert abs(free) > 3.5
    assert abs(shrunk) < abs(free)


def test_primer_offset_prior_keeps_identified_offset():
    rng = np.random.default_rng(5)
    rows = []
    for subject in range(10):
        psi = expit(rng.normal(0, 1))
        shifted = expit(np.log(psi / (1 - psi)) + 1.2)
        rows.append(np.r_[rng.multinomial(400, [.5 * psi / (.5 * psi + .4 * (1 - psi)), .4 * (1 - psi) / (.5 * psi + .4 * (1 - psi))]), rng.multinomial(400, [.5 * shifted / (.5 * shifted + .4 * (1 - shifted)), .4 * (1 - shifted) / (.5 * shifted + .4 * (1 - shifted))])])
    components = np.array([[0., .5, 0.], [0., 0., .4]])
    like = BinaryECPathLikelihood((np.array(rows)[:, :2], np.array(rows)[:, 2:]), (components, components), np.array([f"s{i}" for i in range(10)]), np.array(["a"] * 10), np.ones((10, 2)), np.zeros(10))
    assert abs(estimate_primer_offset(like, prior_sd=1.) - 1.2) < .2
