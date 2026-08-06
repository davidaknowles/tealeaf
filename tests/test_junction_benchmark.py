from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from tealeaf.sc.junction_benchmark import (
    JunctionBundle,
    aggregate_starsolo_runs,
    annotate_scquint_groups,
    normalize_starsolo_barcode,
    prepare_splitseq_whitelists,
    stratified_subject_folds,
    write_leafcutter_junctions,
)


def test_normalize_complex_starsolo_barcode():
    assert normalize_starsolo_barcode("AAACATCG_AACAACCA_AACCGAGA") == (
        "AAACATCGAACAACCAAACCGAGA"
    )


def write_solo(root: Path, barcodes, features, counts):
    sj = root / "Solo.out" / "SJ" / "raw"
    sj.mkdir(parents=True)
    mmwrite(sj / "matrix.mtx", sparse.csr_matrix(counts).T)
    (sj / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
    pd.DataFrame(features).to_csv(sj / "features.tsv", sep="\t", header=False, index=False)


def test_prepare_splitseq_whitelists(tmp_path):
    paths = prepare_splitseq_whitelists(
        ["AAAAAAAACCCCCCCCGGGGGGGG", "AAAAAAAATTTTTTTTGGGGGGGG"], tmp_path
    )
    assert paths[0].read_text() == "AAAAAAAA\n"
    assert paths[1].read_text() == "CCCCCCCC\nTTTTTTTT\n"
    assert paths[2].read_text() == "GGGGGGGG\n"


def test_aggregate_and_leafcutter_export(tmp_path):
    features = [
        ["chr1", 101, 200, 1, 1, 1, 0, 0, 0],
        ["chr1", 101, 250, 1, 1, 1, 0, 0, 0],
    ]
    write_solo(tmp_path / "run1", ["A", "B"], features, [[2, 0], [1, 3]])
    write_solo(tmp_path / "run2", ["C"], features, [[4, 5]])
    metadata = pd.DataFrame(
        {
            "run": ["run1", "run1", "run2"],
            "barcode": ["A", "B", "C"],
            "cell_type": ["X", "X", "Y"],
            "condition": ["case", "case", "control"],
            "subject": ["m1", "m1", "m2"],
        }
    )
    bundle = aggregate_starsolo_runs(
        {"run1": tmp_path / "run1", "run2": tmp_path / "run2"}, metadata
    )
    np.testing.assert_array_equal(bundle.counts.toarray(), [[3, 3], [4, 5]])
    prefix = tmp_path / "bundle"
    bundle.save(prefix)
    loaded = JunctionBundle.load(prefix)
    np.testing.assert_array_equal(loaded.counts.toarray(), bundle.counts.toarray())
    paths = write_leafcutter_junctions(bundle, tmp_path / "leaf")
    first = Path(paths.read_text().splitlines()[0]).read_text().splitlines()[0].split("\t")
    assert first[:3] == ["chr1", "100", "200"]
    assert first[4:6] == ["3", "+"]


def test_annotate_scquint_groups_uses_exon_boundaries(tmp_path):
    bundle = JunctionBundle(
        sparse.csr_matrix([[2, 3]]),
        pd.DataFrame(
            [
                ["chr1", 101, 200, "+"],
                ["chr1", 151, 200, "+"],
            ],
            columns=["chromosome", "start", "end", "strand"],
        ),
        pd.DataFrame(
            [["m1__X", "X", "A", "m1", 1]],
            columns=["sample_id", "cell_type", "condition", "subject", "n_cells"],
        ),
    )
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\tt\texon\t1\t100\t.\t+\t.\tgene_id "g1"; gene_name "G";\n'
        '1\tt\texon\t201\t300\t.\t+\t.\tgene_id "g1"; gene_name "G";\n'
        '1\tt\texon\t1\t150\t.\t+\t.\tgene_id "g1"; gene_name "G";\n'
    )
    result = annotate_scquint_groups(bundle, gtf)
    assert result.junctions.intron_group.nunique() == 1
    assert result.junctions.intron_group_size.tolist() == [2, 2]


def test_stratified_subject_folds_balances_conditions_and_subjects():
    samples = pd.DataFrame(
        {
            "subject": np.repeat([f"s{i}" for i in range(8)], 2),
            "condition": np.repeat(["A"] * 4 + ["B"] * 4, 2),
            "cell_type": ["X", "Y"] * 8,
        }
    )
    folds = stratified_subject_folds(samples, n_folds=2, seed=3)
    assert folds.subject.is_unique
    assert folds.groupby(["fold", "condition"]).size().eq(2).all()
    repeated = stratified_subject_folds(samples, n_folds=2, seed=3)
    pd.testing.assert_frame_equal(folds, repeated)
