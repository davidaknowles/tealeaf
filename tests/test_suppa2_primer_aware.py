import numpy as np
from scipy import sparse
from scipy.stats import wilcoxon

from extra_scripts.run_suppa2_primer_aware import (
    event_matrices,
    fit_primer_aware_psi,
)
from tealeaf.sc import sc_utils
from extra_scripts.run_suppa2_full_data_comparison import (
    event_test_pvalues,
    fast_paired_wilcoxon,
    hybrid_exact_paired_wilcoxon,
)


def test_primer_aware_event_mapping_and_shared_logit():
    catalog = __import__("pandas").DataFrame(
        [
            {
                "feature_id": "e1",
                "event_id": "e1",
                "event_type": "SE",
                "gene_id": "g1",
                "gene_name": "",
                "included": "tx1.1",
                "excluded": "tx2.1",
            },
            {
                "feature_id": "e2",
                "event_id": "e2",
                "event_type": "SE",
                "gene_id": "g1",
                "gene_name": "",
                "included": "missing",
                "excluded": "tx2.1",
            },
        ]
    )
    kept, included, excluded = event_matrices(catalog, np.array(["tx1.1", "tx2.1"]))
    assert kept.event_id.tolist() == ["e1"]
    assert included.toarray().tolist() == [[1.0, 0.0]]
    assert excluded.toarray().tolist() == [[0.0, 1.0]]
    psi, offset = fit_primer_aware_psi(
        np.array([[80.0, 40.0], [60.0, 30.0]]),
        np.array([[20.0, 60.0], [40.0, 70.0]]),
    )
    assert np.all(np.isfinite(psi))
    assert offset > 0


def test_grouped_probability_raw_mode_preserves_ec_rows(tmp_path):
    membership = sparse.csr_matrix([[1.0, 1.0], [1.0, 0.0]])
    sidecar = tmp_path / "probs.tsv.gz"
    import gzip

    with gzip.open(sidecar, "wt") as handle:
        handle.write("cell_idx\teqid\tumi_rank\tprobs\n")
        handle.write("0\t0\t0\t0.25,0.75\n")
        handle.write("0\t1\t0\t1.0\n")
        handle.write("1\t0\t0\t0.5,0.5\n")
        handle.write("1\t1\t0\t1.0\n")
    raw = sc_utils.grouped_ec_probability_matrices(
        sidecar, membership, np.array([0, 1]), 2, normalize_columns=False
    )
    assert np.allclose(raw[0].toarray(), [[0.25, 0.75], [1.0, 0.0]])
    assert np.allclose(raw[1].toarray(), [[0.5, 0.5], [1.0, 0.0]])


def test_vectorized_wilcoxon_averages_tied_ranks():
    differences = np.array([[1.0, 1.0, -2.0, 0.0, np.nan]])
    valid = np.isfinite(differences)
    observed = fast_paired_wilcoxon(differences, valid)[0]
    expected = wilcoxon(differences[0, valid[0]], zero_method="wilcox", method="approx").pvalue
    assert np.isclose(observed, expected)


def test_exact_wilcoxon_ties_matches_exhaustive_signs():
    from itertools import product
    from scipy.stats import rankdata

    differences = np.array([[1., 1., -2., 3., 0., np.nan]])
    ranks = rankdata(np.abs(differences[0, :4]))
    sums = np.array([np.dot(signs, ranks) for signs in product((0, 1), repeat=4)])
    observed_sum = ranks[differences[0, :4] > 0].sum()
    expected = min(1, 2 * np.mean(sums <= min(observed_sum, ranks.sum() - observed_sum)))
    actual = hybrid_exact_paired_wilcoxon(differences, np.isfinite(differences))[0]
    assert actual == expected
    assert actual == hybrid_exact_paired_wilcoxon(differences[:, ::-1], np.isfinite(differences[:, ::-1]))[0]


def test_exact_wilcoxon_many_subjects_does_not_overflow():
    from tealeaf.sc.event_tests import signed_rank_cdf

    for n in (8, 24, 64, 83):
        cdf = signed_rank_cdf(n)
        assert np.all(np.diff(cdf) >= 0)
        assert np.isclose(cdf[-1], 1)
        assert cdf[0] == 2. ** (-n)
        differences = np.arange(1, n + 1, dtype=float)[None, :]
        assert hybrid_exact_paired_wilcoxon(differences, np.ones_like(differences, dtype=bool))[0] == 2. ** (1 - n)


def test_signed_rank_missing_zero_and_tied_extreme_tail():
    from tealeaf.sc.event_tests import paired_signed_rank

    differences = np.array([[0, 0, np.nan], [np.nan, np.nan, np.nan], [1, 1, 1]], dtype=float)
    valid = np.ones_like(differences, dtype=bool)
    for exact in (False, True):
        p = paired_signed_rank(differences, valid, exact=exact)
        assert p[0] == 1
        assert np.isnan(p[1])
    assert paired_signed_rank(differences, valid, exact=True)[2] == .25


def test_normal_wilcoxon_matches_scipy_with_random_ties():
    rng = np.random.default_rng(123)
    differences = rng.integers(-4, 5, size=(30, 24)).astype(float)
    valid = rng.random(differences.shape) > .1
    observed = fast_paired_wilcoxon(differences, valid)
    expected = [wilcoxon(row[mask], zero_method="wilcox", method="approx", correction=False).pvalue for row, mask in zip(differences, valid)]
    np.testing.assert_allclose(observed, expected)


def test_merged_tealeaf_reference_uses_calibrated_pvalues(tmp_path):
    import pandas as pd
    from extra_scripts.compare_suppa2_primer_aware_tealeaf import load_merged_tealeaf

    table = pd.DataFrame({"method": ["local_path"] * 3, "gene_id": ["g1", "g2", "g3"], "level_a": ["a"] * 3, "level_b": ["b"] * 3, "converged": [True, False, True], "n_subjects": [8, 8, 2], "p_value": [.04, .01, .001], "raw_p_value": [.00001, .00002, .00003]})
    path = tmp_path / "paired_path.tsv"
    table.to_csv(path, sep="\t", index=False)
    actual = load_merged_tealeaf(path)
    assert actual.gene_id.tolist() == ["g1"]
    assert actual.p_value.tolist() == [.04]


def test_event_paired_t_matches_scipy():
    from scipy.stats import ttest_rel

    first = np.array([
        [0.1, 0.2, 0.4, 0.3],
        [0.2, 0.2, np.nan, 0.2],
    ])
    second = np.array([
        [0.3, 0.5, 0.7, 0.6],
        [0.2, 0.2, np.nan, 0.2],
    ])
    valid = np.isfinite(first) & np.isfinite(second)
    observed = event_test_pvalues(first, second, valid, "paired_t")
    assert np.isclose(observed[0], ttest_rel(second[0], first[0]).pvalue)
    assert observed[1] == 1.0


def test_hybrid_exact_wilcoxon_matches_untied_exact_tail():
    differences = np.array([[1.0, 2.0, 3.0, 4.0, -5.0]])
    valid = np.ones_like(differences, dtype=bool)
    observed = hybrid_exact_paired_wilcoxon(differences, valid)[0]
    expected = wilcoxon(differences[0], method="exact").pvalue
    assert np.isclose(observed, expected)
