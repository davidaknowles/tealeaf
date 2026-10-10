import numpy as np
from scipy.special import expit

from tealeaf.sc.local_path_reads import path_read_opportunities
from tealeaf.sc.local_path_test import PRIMERS, ProfiledPrecursorLikelihood, fit_local_block, local_class_masks, local_read_likelihood

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
