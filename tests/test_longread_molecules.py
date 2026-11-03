import pandas as pd

from tealeaf.sc import longread_molecules as lm


def test_exon_chain_splits_on_introns_only():
    # 5S 10M 2D 5M 100N 20M 1I 5M
    cigar = [(4, 5), (0, 10), (2, 2), (0, 5), (3, 100), (0, 20), (1, 1), (0, 5)]
    assert lm.exon_chain(1000, cigar) == [(1000, 1017), (1117, 1142)]
    assert lm.parse_chain(lm.format_chain([(1, 2), (5, 9)])) == [(1, 2), (5, 9)]


def test_gene_index_strand_and_overlap():
    genes = {
        "g1": {"chromosome": "chr1", "strand": "+", "transcripts": {"t1": [(100, 200), (300, 400)]}},
        "g2": {"chromosome": "chr1", "strand": "+", "transcripts": {"t2": [(380, 600)]}},
        "g3": {"chromosome": "chr1", "strand": "-", "transcripts": {"t3": [(100, 400)]}},
    }
    index = lm.GeneIndex(genes)
    assert index.assign("chr1", "+", [(150, 200), (300, 390)]) == "g1"
    assert index.assign("chr1", "+", [(390, 590)]) == "g2"
    assert index.assign("chr1", "-", [(150, 200)]) == "g3"
    assert index.assign("chr1", "+", [(210, 290)]) is None
    assert index.assign("chr2", "+", [(150, 200)]) is None


def test_collapse_molecules_merges_one_mismatch_umis():
    frame = pd.DataFrame({
        "cell": ["c1"] * 4 + ["c2"],
        "gene": ["g"] * 5,
        "umi": ["AAAAAAAAAA", "AAAAAAAAAC", "AAAAAAAAAC", "TTTTTTTTTT", "AAAAAAAAAA"],
        "chain": ["1-5", "1-5;8-9", "1-5;8-9", "1-5", "1-5"],
        "span": [4, 8, 8, 4, 4],
    })
    molecules = lm.collapse_molecules(frame)
    assert len(molecules) == 3
    first = molecules.loc[molecules.cell.eq("c1") & molecules.n_reads.eq(3)]
    assert len(first) == 1 and first.chain.iloc[0] == "1-5;8-9"


def test_molecule_source_counts_chains_by_type_and_replicate():
    molecules = pd.DataFrame({"gene": ["ENSG1.2"] * 4 + ["ENSG3.1"], "chain": ["1-5;8-9", "1-5;8-9", "1-9", "1-5;8-9", "20-30"], "cell": ["a", "b", "a", "c", "d"]})
    labels = pd.DataFrame({"cell_type": ["X", "X", "Y", "Y"], "replicate": [1, 2, 1, 1]}, index=["a", "b", "c", "d"])
    matrix, features, columns, chains = lm.molecule_source(molecules, labels)
    assert matrix.shape == (3, 3) and matrix.sum() == 5
    row = features.set_index("transcript_id").row["ENSG1|1-5;8-9"]
    column = columns.set_index(["tealeaf_cell_type", "replicate"]).column
    assert matrix[row, column[("X", 1)]] == 1 and matrix[row, column[("X", 2)]] == 1 and matrix[row, column[("Y", 1)]] == 1
    assert chains["ENSG1|1-5;8-9"] == [(1, 5), (8, 9)]


def test_cell_labels_clean_and_stranded_chains(tmp_path):
    table = pd.DataFrame({"cell": ["A1_B2_C3", "A2_B2_C3"], "low": ["Tcm/Naive helper T cells", "CD16+ NK cells"]})
    table.to_csv(tmp_path / "labels.tsv", sep="\t", index=False)
    labels = lm.cell_labels(tmp_path / "labels.tsv", "low", clean=True)
    assert labels.cell_type.tolist() == ["Tcm_Naive_helper_T_cells", "CD16__NK_cells"] and labels.replicate.tolist() == [2, 1]
    gtf = tmp_path / "a.gtf"
    gtf.write_text('chr1\tx\texon\t11\t20\t.\t-\t.\tgene_id "ENSG1.2"; transcript_id "ENST1.1";\n')
    chains = lm.stranded_chains({"ENSG1|10-20": [(10, 20)], "ENSG9|1-5": [(1, 5)]}, gtf)
    assert chains == {"ENSG1|10-20": ("chr1", "-", ((10, 20),))}
