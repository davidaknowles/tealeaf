import gzip

import numpy as np
import scipy.sparse as sp

from tealeaf.data.salmon import build_positional_ec_design


def test_build_positional_design_collapses_rows_and_falls_back(tmp_path):
    features = ["tx0", "tx1", "tx2"]
    membership = sp.csr_matrix([
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 0],
    ])
    quant = tmp_path / "quant"
    (quant / "aux_info").mkdir(parents=True)
    (quant / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "tx0\t100\t10\t0\t0\n"
        "tx1\t200\t20\t0\t0\n"
        "tx2\t300\t30\t0\t0\n"
    )
    (quant / "aux_info" / "meta_info.json").write_text(
        '{"num_mapped": 6}\n'
    )
    # Salmon target order differs from alevin. The first two rows are
    # range-factorized instances of the same transcript set.
    with gzip.open(quant / "aux_info" / "eq_classes.txt.gz", "wt") as handle:
        handle.write(
            "3\n3\n"
            "tx2\n"
            "tx0\n"
            "tx1\n"
            "2\t1\t2\t0.8\t0.2\t3\n"
            "2\t2\t1\t0.6\t0.4\t1\n"
            "1\t1\t1.0\t2\n"
        )
    quant2 = tmp_path / "quant2"
    (quant2 / "aux_info").mkdir(parents=True)
    (quant2 / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "tx0\t100\t20\t0\t0\n"
        "tx1\t200\t40\t0\t0\n"
        "tx2\t300\t60\t0\t0\n"
    )
    (quant2 / "aux_info" / "meta_info.json").write_text(
        '{"num_mapped": 4}\n'
    )
    with gzip.open(quant2 / "aux_info" / "eq_classes.txt.gz", "wt") as handle:
        handle.write(
            "3\n1\n"
            "tx0\n"
            "tx1\n"
            "tx2\n"
            "2\t0\t1\t0.2\t0.8\t2\n"
        )

    design, stats, effective_lengths = build_positional_ec_design(
        membership,
        features,
        [
            quant / "aux_info" / "eq_classes.txt.gz",
            quant2 / "aux_info" / "eq_classes.txt.gz",
        ],
        [quant, quant2],
        ec_counts=np.array([8, 2, 4]),
    )
    # Before column normalization, EC0 is the count-weighted average
    # across range-factorized rows and runs. EC1 uses the averaged effective
    # lengths as a fallback, and EC2 is exact.
    raw = np.array([
        [3.2 / 6, 2.8 / 6, 0.0],
        [0.0, 0.6, 0.4],
        [1.0, 0.0, 0.0],
    ])
    expected = raw / raw.sum(axis=0)
    np.testing.assert_allclose(design.toarray(), expected, rtol=1e-6)
    assert stats["salmon_duplicate_rows"] == 2
    assert stats["salmon_quantifications"] == 2
    assert stats["matched_alevin_ecs"] == 2
    assert stats["fallback_alevin_ecs"] == 1
    assert stats["matched_molecule_fraction"] == 12 / 14
    np.testing.assert_allclose(effective_lengths, [14, 28, 42])
