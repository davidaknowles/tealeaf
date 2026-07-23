import gzip

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

from tealeaf.data.alevin import (
    load_alevin_counts,
    load_alevin_structure,
    merge_alevin_quantifications,
    validate_alevin_quantification,
)


def _write_quant(path, features, barcodes, counts, membership, probabilities):
    path.mkdir()
    sp.save_npz(path / "geqc_counts.npz", sp.csr_matrix(counts))
    sp.save_npz(path / "gene_eqclass.npz", sp.csr_matrix(membership))
    (path / "quants_mat_cols.txt").write_text("\n".join(features) + "\n")
    (path / "quants_mat_rows.txt").write_text("\n".join(barcodes) + "\n")
    with gzip.open(path / "gene_eqclass_probs.tsv.gz", "wt") as handle:
        handle.write("cell_idx\teqid\tumi_rank\tprobs\n")
        for row in probabilities:
            handle.write("\t".join(map(str, row)) + "\n")


def test_merge_aligns_features_ecs_and_probabilities(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _write_quant(
        first,
        ["tx1", "tx2"],
        ["cell"],
        [[2, 3]],
        [[1, 1], [1, 0]],
        [(0, 0, 0, "0.25000000000000001,0.74999999999999999")],
    )
    _write_quant(
        second,
        ["tx2", "tx1"],
        ["cell"],
        [[5, 7]],
        [[1, 1], [1, 0]],
        [(0, 0, 0, "0.6,0.4")],
    )
    output = tmp_path / "merged"
    report = merge_alevin_quantifications(
        [("b1", first), ("b2", second)], output
    )
    counts = sp.load_npz(output / "geqc_counts.npz").toarray()
    assert report["cells"] == 2
    assert report["equivalence_classes"] == 3
    assert (output / "quants_mat_rows.txt").read_text().splitlines() == [
        "b1:cell",
        "b2:cell",
    ]
    np.testing.assert_array_equal(counts[0], [2, 3, 0])
    np.testing.assert_array_equal(counts[1], [5, 0, 7])
    with gzip.open(output / "gene_eqclass_probs.tsv.gz", "rt") as handle:
        rows = handle.read().splitlines()
    assert rows[1].split("\t")[-1] == "0.25000000000000001,0.74999999999999999"
    np.testing.assert_allclose(
        np.fromstring(rows[2].split("\t")[-1], sep=","),
        [0.4, 0.6],
    )


def test_merge_accepts_alevin_fry_output_root(tmp_path):
    first = tmp_path / "first" / "alevin"
    second = tmp_path / "second" / "alevin"
    first.parent.mkdir()
    second.parent.mkdir()
    _write_quant(first, ["tx1"], ["cell1"], [[2]], [[1]], [])
    _write_quant(second, ["tx1"], ["cell2"], [[3]], [[1]], [])

    output = tmp_path / "merged"
    report = merge_alevin_quantifications(
        [("run1", first.parent), ("run2", second.parent)], output
    )

    assert report["cells"] == 2
    assert (output / "quants_mat_rows.txt").read_text().splitlines() == [
        "run1:cell1",
        "run2:cell2",
    ]


def test_loaders_accept_uncached_alevin_fry_output_root(tmp_path):
    root = tmp_path / "quant"
    quant = root / "alevin"
    quant.mkdir(parents=True)
    (quant / "quants_mat_cols.txt").write_text("tx1\ntx2\n")
    (quant / "quants_mat_rows.txt").write_text("cell1\ncell2\n")
    scipy.io.mmwrite(
        quant / "geqc_counts.mtx",
        sp.csr_matrix([[2, 0], [0, 3]], dtype=np.int64),
    )
    with gzip.open(quant / "gene_eqclass.txt.gz", "wt") as handle:
        handle.write("2\n2\n0\t0\n1\t1\n")

    features, membership = load_alevin_structure(root)
    barcodes, counts = load_alevin_counts(root)

    assert features == ["tx1", "tx2"]
    assert barcodes == ["cell1", "cell2"]
    np.testing.assert_array_equal(membership.toarray(), np.eye(2))
    np.testing.assert_array_equal(counts.toarray(), [[2, 0], [0, 3]])


def test_validate_merged_quantification(tmp_path):
    quant = tmp_path / "quant"
    _write_quant(
        quant,
        ["tx1"],
        ["run1:cell1", "run1:cell2"],
        [[500], [20]],
        [[1]],
        [(0, 0, 0, "1")],
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text(
        "cell_id\tpolydt_barcode\tranhex_barcode\n"
        "run1:cell1\trun1:cell1\trun1:cell2\n"
    )

    report = validate_alevin_quantification(
        quant,
        expected_prefixes={"run1"},
        reference_ids={"run1:cell1"},
        primer_pair_file=pairs,
        min_cell_umis=500,
    )

    assert report["eligible_cells"] == 1
    assert report["reference_overlap"] == 1
    assert report["complete_primer_pairs"] == 1
    assert report["max_total_molecules"] is None

    with pytest.raises(ValueError, match="exceeding the limit"):
        validate_alevin_quantification(
            quant,
            min_cell_umis=1,
            max_total_molecules=519,
        )
