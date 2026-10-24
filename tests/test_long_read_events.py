from tealeaf.sc.long_read_events import chain_introns, chain_junction_rows, event_transcript_classes, rmats_event, suppa2_event


def test_suppa2_skipped_exon_forms_from_exon_chains():
    event = suppa2_event("G;SE:chr1:100-201:250-301:+")
    include = ((49, 100), (200, 250), (300, 400))
    skip = ((49, 100), (300, 400))
    other = ((49, 120), (300, 400))
    assert chain_introns(include) == ((101, 200), (251, 300))
    assert event.classify(include) == 1 and event.classify(skip) == 0 and event.classify(other) is None
    chains = {"a": ("chr1", "+", include), "b": ("chr1", "+", skip), "c": ("chr1", "-", skip)}
    assert event_transcript_classes([event], chains) == [(["a"], ["b"])]


def test_retained_intron_and_terminal_forms():
    retained = suppa2_event("G;RI:chr1:50:100-201:300:+")
    assert retained.classify(((49, 300),)) == 1 and retained.classify(((49, 100), (200, 300))) == 0
    first = suppa2_event("G;AF:chr1:10:50-201:120:150-201:+")
    assert first.classify(((5, 50), (200, 300))) == 1 and first.classify(((100, 150), (200, 300))) == 0
    assert first.classify(((1, 20), (100, 150), (200, 300))) is None
    assert first.classify(((1, 20), (30, 50), (200, 300))) is None


def test_rmats_mxe_minus_strand_uses_transcript_5prime_exon():
    row = {"ID": "0", "GeneID": "g", "geneSymbol": "g", "chr": "chr1", "strand": "+", "1stExonStart_0base": "200", "1stExonEnd": "250", "2ndExonStart_0base": "300", "2ndExonEnd": "350", "upstreamES": "50", "upstreamEE": "100", "downstreamES": "400", "downstreamEE": "500"}
    first = ((50, 100), (200, 250), (400, 500))
    second = ((50, 100), (300, 350), (400, 500))
    assert rmats_event("MXE", row).classify(first) == 1
    assert rmats_event("MXE", {**row, "strand": "-"}).classify(first) == 0
    assert rmats_event("MXE", {**row, "strand": "-"}).classify(second) == 1


def test_junction_rows_cover_both_strands():
    chains = {"p": ("chr1", "+", ((0, 10), (20, 30))), "m": ("chr1", "-", ((0, 10), (20, 30)))}
    rows = chain_junction_rows(chains, {"p": [0], "m": [1, 2]}, {("chr1", 11, 20, "+"), ("chr1", 11, 20, "-")})
    assert rows[("chr1", 11, 20, "+")].tolist() == [0] and rows[("chr1", 11, 20, "-")].tolist() == [1, 2]
