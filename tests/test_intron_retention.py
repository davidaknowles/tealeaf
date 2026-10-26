import numpy as np

from tealeaf.sc.intron_retention import a_run_positions, contained_introns, downstream_run_distance, intron_catalogue


def test_catalogue_clean_and_retained_labels():
    transcripts = {
        "t1": {"gene_id": "g", "chromosome": "chr1", "strand": "+", "type": "protein_coding", "exons": [(0, 100), (200, 300), (400, 500)]},
        # skips the middle exon: its intron [100, 400) holds exon boundaries of t1, so it is not clean
        "t2": {"gene_id": "g", "chromosome": "chr1", "strand": "+", "type": "protein_coding", "exons": [(0, 100), (400, 500)]},
        # retains intron [300, 400)
        "t3": {"gene_id": "g", "chromosome": "chr1", "strand": "+", "type": "retained_intron", "exons": [(200, 500)]},
        # an exon of another gene on the other strand inside intron [100, 200)
        "t4": {"gene_id": "h", "chromosome": "chr1", "strand": "-", "type": "protein_coding", "exons": [(120, 150), (160, 180)]},
    }
    table = intron_catalogue(transcripts, genes={"g"}).set_index(["start", "end"])
    assert list(table.index) == [(100, 200), (100, 400), (300, 400)]
    assert table.clean.to_dict() == {(100, 200): False, (100, 400): False, (300, 400): True}
    assert table.retained_by_ri.to_dict() == {(100, 200): False, (100, 400): False, (300, 400): True}
    assert table.loc[(100, 400), "n_transcripts"] == 1


def test_a_runs_and_downstream_distance():
    codes = np.array([1] * 5 + [0] * 9 + [1] * 5 + [3] * 10, dtype=np.uint8)
    assert a_run_positions(codes, 0, window=10, min_count=8).tolist() == [3, 4, 5, 6]
    assert a_run_positions(codes, 3, window=10, min_count=8).tolist() == [17, 18, 19]
    runs = np.array([50, 200])
    assert downstream_run_distance([10, 50, 120, 300], runs, plus=True).tolist() == [40, 0, 80, np.inf]
    assert downstream_run_distance([10, 50, 120, 300], runs, plus=False).tolist() == [np.inf, 0, 70, 100]


def test_contained_introns_with_overhang():
    intron_starts, intron_ends = np.array([100, 300, 600]), np.array([200, 400, 700])
    block, intron = contained_introns([90, 50, 250, 0], [210, 450, 420, 800], intron_starts, intron_ends, overhang=10)
    assert list(zip(block.tolist(), intron.tolist())) == [(0, 0), (1, 0), (1, 1), (2, 1), (3, 0), (3, 1), (3, 2)]
