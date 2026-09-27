import numpy as np
from scipy import sparse
from scipy.stats import wilcoxon

from extra_scripts.run_suppa2_primer_aware import (
    event_matrices,
    fit_primer_aware_psi,
)
from tealeaf.sc import sc_utils
from extra_scripts.run_suppa2_full_data_comparison import fast_paired_wilcoxon


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
