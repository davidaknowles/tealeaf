import numpy as np
import pytest
from scipy.special import expit

from tealeaf.sc.local_path_reads import path_read_opportunities
from tealeaf.sc.local_path_test import PRIMERS, ProfiledPrecursorLikelihood, fit_local_block, local_class_masks, local_read_likelihood, prepare_block_test

SKIP = [(100, 200), (500, 600)]
INCLUDE = [(100, 200), (300, 400), (500, 600)]
PRECURSOR = [(100, 600)]
OPPORTUNITIES = path_read_opportunities([SKIP, INCLUDE, PRECURSOR], 50)


def test_profiled_precursor_matches_brute_force():
    masks = local_class_masks(OPPORTUNITIES, 2, True)
    components = np.column_stack([[OPPORTUNITIES[m][2] for m in masks], [OPPORTUNITIES[m][0] for m in masks], [OPPORTUNITIES[m][1] for m in masks]])
    counts = np.random.default_rng(1).integers(0, 40, size=(1, len(masks))).astype(float)
    like = ProfiledPrecursorLikelihood((counts, counts), (components, components), np.array(["s"]), np.array(["a"]), np.ones((1, 2)), np.zeros(1))
    t = np.array([-1., .3, 2.])
    profiled = like.primer_log_likelihood(0, t, 0)
    for index, value in enumerate(t):
        psi = expit(value)
        brute = max(counts[0] @ np.log(x * components[:, 0] + psi * components[:, 1] + (1 - psi) * components[:, 2]) - counts[0].sum() * np.log(x * components[:, 0].sum() + psi * components[:, 1].sum() + (1 - psi) * components[:, 2].sum()) for x in np.exp(np.linspace(-15, 8, 20001)))
        assert profiled[index] >= brute - 1e-9 and profiled[index] - brute < 1e-4


def simulate(rng, precursor_fraction):
    """Same mature inclusion in both types; type b carries extra precursor."""
    lookup, subjects, labels = {}, [], []
    lengths = sum(OPPORTUNITIES.values())
    masks = sorted(OPPORTUNITIES)
    for subject in range(8):
        psi = expit(rng.normal(0, .3))
        for label, rho in enumerate(precursor_fraction):
            molecules = np.array([(1 - rho) * psi, (1 - rho) * (1 - psi), rho])
            probabilities = np.array([OPPORTUNITIES[m] @ (molecules / lengths) for m in masks])
            for primer in PRIMERS:
                for mask, value in zip(masks, rng.multinomial(1500, probabilities / probabilities.sum())):
                    lookup[("b", f"s{subject}", "AB"[label], primer, mask)] = value
            subjects.append(f"s{subject}")
            labels.append(label)
    return lookup, np.array(subjects), np.array(labels)


def test_precursor_model_ignores_unspliced_excess():
    rng = np.random.default_rng(7)
    lookup, subjects, labels = simulate(rng, (.05, .6))
    with_precursor, _ = fit_local_block(lookup, "b", OPPORTUNITIES, 2, subjects, labels, ("A", "B"), anchors=False, precursor=True, primer_offset=False)
    mature_only = {mask & 3: 0 for mask in OPPORTUNITIES if mask & 3}
    collapsed = {}
    for (block, subject, level, primer, mask), value in lookup.items():
        if mask & 3:
            key = (block, subject, level, primer, mask & 3)
            collapsed[key] = collapsed.get(key, 0) + value
    plain_opportunities = path_read_opportunities([SKIP, INCLUDE], 50)
    without, _ = fit_local_block(collapsed, "b", plain_opportunities, 2, subjects, labels, ("A", "B"), anchors=False, precursor=False, primer_offset=False)
    assert with_precursor["converged"] and without["converged"]
    assert with_precursor["p_value"] > .01
    assert without["p_value"] < 1e-4


def test_compiled_precursor_likelihood_matches_reference():
    rng = np.random.default_rng(5)
    components = np.column_stack([rng.uniform(0, 3, 6), rng.uniform(0, 2, 6), rng.uniform(0, 2, 6)])
    components[0, 0] = 0.
    counts = rng.poisson(8, size=(3, 6)).astype(float)
    counts[1, 2] = 0.
    like = ProfiledPrecursorLikelihood((counts, counts), (components, components), np.array(["s", "s", "u"]), np.array([0, 1, 0]), np.ones((3, 2)), np.zeros(3))
    t = np.linspace(-18, 18, 241)
    scales = np.exp(np.linspace(-25, 15, 40001))
    for row in range(3):
        compiled = like.primer_log_likelihood(row, t, 0)
        # At least as high as the grid-and-Newton reference, and the maximum of a dense scale grid.
        assert (compiled >= like.primer_log_likelihood_reference(row, t, 0) - 1e-9).all()
        observed = counts[row]
        for index in range(0, len(t), 20):
            psi = expit(t[index])
            mass = scales[:, None] * components[:, 0] + psi * components[:, 1] + (1 - psi) * components[:, 2]
            dense = (np.log(np.maximum(mass, 1e-300)) @ observed - observed.sum() * np.log(mass.sum(axis=1))).max()
            assert dense - 1e-9 <= compiled[index] <= dense + 1e-6


def test_eligibility_rechecked_after_dropping_precursor():
    paths = [[(0, 100), (300, 400)], [(0, 100), (200, 250), (300, 400)], [(0, 400)]]
    opportunities = path_read_opportunities(paths, 50)
    precursor_only = 1 << 2
    mature = 0b011
    levels = ("a", "b")
    entries = [(f"s{index}", level, PRIMERS[0], precursor_only, 10.) for index in range(6) for level in levels]
    entries += [(f"s{index}", level, PRIMERS[0], mask, 5.) for index in range(2) for level in levels for mask in (mature, 0b001, 0b010)]
    prepare_block_test("b", entries, opportunities, 2, levels, precursor=True)
    with pytest.raises(ValueError):
        prepare_block_test("b", entries, opportunities, 2, levels, precursor=False)


def test_identical_primer_opportunities_match_shared():
    rng = np.random.default_rng(1)
    lookup = {}
    subjects, labels = [], []
    for index in range(6):
        for label, level in enumerate(("a", "b")):
            subjects.append(f"s{index}")
            labels.append(label)
            for primer in PRIMERS:
                for mask in OPPORTUNITIES:
                    if mask & 3:
                        lookup[("b", f"s{index}", level, primer, mask)] = float(rng.poisson(20))
    mature = {mask: vector[:2] for mask, vector in OPPORTUNITIES.items() if mask & 3}
    subjects, labels = np.array(subjects), np.array(labels)
    shared = local_read_likelihood(lookup, "b", mature, subjects, labels, ("a", "b"), True, n_paths=2)
    split = local_read_likelihood(lookup, "b", mature, subjects, labels, ("a", "b"), True, n_paths=2, primer_opportunities={primer: mature for primer in PRIMERS})
    for row in range(len(shared.counts[0])):
        assert np.allclose(shared.row_log_likelihood(row, np.array([.3, .6])), split.row_log_likelihood(row, np.array([.3, .6])))


def test_rest_component_splits_molecules_as_pooled_shares():
    # paths of different effective lengths; the collapsed rest must give each
    # non-target path its pooled molecule share, so weights are share / length
    paths = [SKIP, INCLUDE, [(100, 200), (250, 450), (500, 600)]]
    opportunities = path_read_opportunities(paths, 50)
    lookup = {("b", "s0", level, primer, mask): 1. for level in ("a", "b") for primer in PRIMERS for mask in opportunities}
    shares = np.array([.2, .3, .5])
    likelihood = local_read_likelihood(lookup, "b", opportunities, np.array(["s0", "s0"]), np.array([0, 1]), ("a", "b"), True, target=0, shares=shares, n_paths=3)
    masks = local_class_masks(opportunities, 3, True)
    n = np.array([opportunities[mask] for mask in masks], dtype=float)
    weights = np.linalg.lstsq(n[:, 1:], likelihood.components[0][:, 2], rcond=None)[0]
    molecules = weights * n[:, 1:].sum(axis=0)
    assert np.allclose(molecules / molecules.sum(), shares[1:] / shares[1:].sum())
    assert not np.isclose(n[:, 1].sum(), n[:, 2].sum())


def test_transcript_path_assignments():
    from tealeaf.sc.local_path_reads import transcript_path_assignments

    transcripts = {"skip": [(0, 200), (500, 700)], "include": [(50, 200), (300, 400), (500, 900)], "outside": [(0, 200), (250, 280), (500, 700)], "partial": [(0, 150)]}
    assignments = transcript_path_assignments([[SKIP], [INCLUDE]], transcripts)
    # "partial" fits both paths and "outside" has an exon in no path
    assert assignments == {"skip": 0, "include": 1}
    assert transcript_path_assignments([[SKIP, INCLUDE]], transcripts) == {"skip": 0, "include": 0, "partial": 0}


def test_unpaired_eligibility_counts_subjects_per_level():
    mature = 0b011
    entries = [(f"s{index}", "a" if index < 4 else "b", PRIMERS[0], mask, 5.) for index in range(8) for mask in (mature, 0b001, 0b010)]
    paths = [[(0, 100), (300, 400)], [(0, 100), (200, 250), (300, 400)]]
    opportunities = path_read_opportunities(paths, 50)
    with pytest.raises(ValueError):
        prepare_block_test("b", entries, opportunities, 2, ("a", "b"))
    prepared = prepare_block_test("b", entries, opportunities, 2, ("a", "b"), unpaired=True)
    assert sorted(prepared["labels"].tolist()) == [0] * 4 + [1] * 4
    with pytest.raises(ValueError):
        prepare_block_test("b", entries, opportunities, 2, ("a", "b"), min_subjects=5, unpaired=True)
