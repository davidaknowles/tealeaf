import gzip

import numpy as np
import scipy.sparse as sp

from tealeaf.data.alevin import merge_alevin_quantifications


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
        [(0, 0, 0, "0.25,0.75")],
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
    assert rows[1].split("\t")[-1] == "0.25,0.75"
    assert rows[2].split("\t")[-1] == "0.4,0.6"
