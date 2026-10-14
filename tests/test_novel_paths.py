from tealeaf.sc.novel_paths import novel_block_paths, read_star_junctions

INCLUDED = [(0, 100), (200, 300), (400, 500)]
SKIPPED = [(0, 100), (400, 500)]
ANNOTATED = {(100, 200): (500, 8), (300, 400): (500, 8), (100, 400): (300, 8)}


def chains(paths, support, **options):
    return [[tuple(exon) for exon in chain] for chain, _ in novel_block_paths(paths, {**ANNOTATED, **support}, **options)]


def test_alternative_acceptor_uses_path_without_extension():
    assert chains([INCLUDED, SKIPPED], {(100, 220): (100, 8)}) == [[(0, 100), (220, 300), (400, 500)]]


def test_new_exon_from_junction_pair_without_extended_singles():
    found = chains([INCLUDED, SKIPPED], {(100, 330): (60, 8), (360, 400): (50, 8)})
    assert found == [[(0, 100), (330, 360), (400, 500)]]


def test_support_filters_and_window():
    assert chains([INCLUDED, SKIPPED], {(100, 220): (19, 8)}) == []
    assert chains([INCLUDED, SKIPPED], {(100, 220): (100, 3)}) == []
    assert chains([INCLUDED, SKIPPED], {(100, 220): (100, 8)}, min_fraction=.5) == []
    assert chains([INCLUDED, SKIPPED], {(100, 600): (100, 8)}) == []


def test_exitron_and_annotated_paths_not_repeated():
    found = chains([INCLUDED, SKIPPED], {(420, 460): (100, 8)})
    assert found == [[(0, 100), (200, 300), (400, 420), (460, 500)]]


def test_read_star_junctions(tmp_path):
    first, second = tmp_path / "a.tab", tmp_path / "b.tab"
    first.write_text("chr1\t101\t200\t1\t1\t0\t5\t0\t30\nchr1\t101\t300\t1\t0\t0\t9\t0\t30\n")
    second.write_text("chr1\t101\t200\t1\t1\t0\t7\t2\t30\n")
    assert read_star_junctions([[str(first)], [str(second)]]) == {("chr1", 100, 200): (12, 2, "+")}


def test_new_block_for_skipped_constitutive_exon():
    from tealeaf.sc.novel_paths import novel_junction_blocks
    transcripts = [[(0, 100), (200, 300), (400, 500), (600, 700)], [(50, 100), (200, 300), (400, 500), (600, 650)]]
    support = {(100, 200): (500, 8, "+"), (300, 400): (500, 8, "+"), (500, 600): (500, 8, "+"), (300, 600): (80, 8, "+")}
    blocks = novel_junction_blocks(transcripts, support, windows=[(0, 300)])
    assert blocks == [([[[200, 300], [400, 500], [600, 650]]], [([[200, 300], [600, 650]], 80)])]
    assert novel_junction_blocks(transcripts, support, windows=[(250, 650)]) == []


def test_novel_first_and_last_exons():
    included, skipped = [(0, 100), (1000, 1100), (2000, 2100)], [(0, 100), (2000, 2100)]
    annotated = {(100, 1000): (500, 8), (1100, 2000): (500, 8), (100, 2000): (300, 8)}

    def found(support, **options):
        return [[tuple(exon) for exon in chain] for chain, _ in novel_block_paths([included, skipped], {**annotated, **support}, **options)]

    # deep donor at 1600 spliced to the last exon: new first exon, clipped at the preceding exon's end
    assert found({(1600, 2000): (100, 8)}) == []
    assert found({(1600, 2000): (100, 8)}, terminal_length=30) == [[(1570, 1600), (2000, 2100)]]
    assert found({(1600, 2000): (100, 8)}, terminal_length=1000) == [[(1100, 1600), (2000, 2100)]]
    # deep acceptor: new last exon
    assert found({(1100, 1500): (100, 8)}, terminal_length=30) == [[(0, 100), (1000, 1100), (1500, 1530)]]
    # a junction pair still makes an internal new exon
    assert found({(100, 1400): (60, 8), (1500, 2000): (50, 8)}, terminal_length=30) == [[(0, 100), (1400, 1500), (2000, 2100)]]


def test_new_block_for_novel_first_exon():
    from tealeaf.sc.novel_paths import novel_junction_blocks
    transcripts = [[(0, 100), (2000, 2100), (3000, 3100)]]
    support = {(100, 2000): (500, 8, "+"), (2100, 3000): (500, 8, "+"), (1500, 2000): (80, 8, "+")}
    assert novel_junction_blocks(transcripts, support, windows=[]) == []
    blocks = novel_junction_blocks(transcripts, support, windows=[], terminal_length=151)
    assert blocks == [([[[0, 100], [2000, 2100]]], [([[1349, 1500], [2000, 2100]], 80)])]
