import numpy as np

from tealeaf.sc.intron_retention import a_run_positions, contained_introns, downstream_run_distance, calibrated_p_value, intron_catalogue, null_inflation, retention_test


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


def test_retention_test_offset_shifts_intercept_only():
    rng = np.random.default_rng(3)
    subjects, labels = np.repeat(np.arange(8), 2), np.tile([0, 1], 8)
    offset = rng.normal(-2, .5, 16)
    probability = 1 / (1 + np.exp(-(-1 + offset + 1.5 * labels)))
    total = rng.poisson(40, 16) + 5
    retained = rng.binomial(total, probability)
    base = retention_test(retained, total - retained, subjects, labels, offset=offset)
    shifted = retention_test(retained, total - retained, subjects, labels, offset=offset + 1.)
    assert base["converged"] and base["p_value"] < 1e-3 and base["coefficients"][1] > .5
    assert np.isclose(shifted["statistic"], base["statistic"], atol=1e-3)
    assert np.isclose(shifted["coefficients"][0], base["coefficients"][0] - 1., atol=1e-2)
    # a gene-wide change carried by the offset is not a retention effect
    gene_wide = retention_test(retained, total - retained, subjects, labels, offset=offset + 1.5 * labels)
    assert np.isclose(gene_wide["coefficients"][1], base["coefficients"][1] - 1.5, atol=1e-2)


def test_null_inflation_recovers_scale():
    rng = np.random.default_rng(0)
    statistics = 1.3 * rng.chisquare(1, 200000)
    inflation = null_inflation(statistics)
    assert abs(inflation - 1.3) < .05
    assert np.isclose(calibrated_p_value(1.3 * 3.841, 1.3), .05, atol=1e-3)
