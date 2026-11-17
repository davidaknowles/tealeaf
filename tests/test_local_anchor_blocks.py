import gzip
from pathlib import Path

from tealeaf.sc.differential import _constitutive_intervals, _local_anchor_intervals, build_splice_blocks
from tealeaf.sc.junction_paths import block_path_exons


def write_gtf(path, transcripts, strand="+"):
    lines = []
    for name, exons in transcripts.items():
        for start, end in exons:
            lines.append(f'chr1\tt\texon\t{start + 1}\t{end}\t.\t{strand}\t.\tgene_id "G"; transcript_id "{name}";')
    path.write_text("\n".join(lines) + "\n")


# a CD45-like gene: exon 3, cassette exons 4-6 in combinations, exon 7; plus a transcript starting downstream at exon 9
EXONS = {3: (100, 150), 4: (300, 400), 5: (500, 600), 6: (700, 800), 7: (900, 950), 8: (1100, 1150), 9: (1300, 1350), 10: (1500, 1600)}
TRANSCRIPTS = {
    "RABC": [EXONS[i] for i in (3, 4, 5, 6, 7, 8, 9, 10)],
    "RA": [EXONS[i] for i in (3, 4, 7, 8, 9, 10)],
    "RB": [EXONS[i] for i in (3, 5, 7, 8, 9, 10)],
    "RO": [EXONS[i] for i in (3, 7, 8, 9, 10)],
    "short": [(1250, 1350), EXONS[10]],
}


def test_local_anchors_survive_a_transcript_starting_downstream():
    chains = [tuple(c) for c in TRANSCRIPTS.values()]
    assert (100, 150) not in _constitutive_intervals(chains)
    local = _local_anchor_intervals(chains)
    assert (100, 150) in local and (900, 950) in local


def test_local_block_tests_the_cassette_combinations(tmp_path: Path):
    gtf = tmp_path / "g.gtf"
    write_gtf(gtf, TRANSCRIPTS)
    blocks = build_splice_blocks(gtf, local_anchors=True)
    local = [b for b in blocks if ":L" in b.block_id]
    cassette = [b for b in local if any(len(s) > 2 for s in b.path_signatures) and min(a for s in b.path_signatures for a, _ in s) == 100]
    assert len(cassette) == 1
    paths = {tuple(block_path_exons(None, s, None)) for s in cassette[0].path_signatures}
    assert paths == {((100, 150), (300, 400), (500, 600), (700, 800), (900, 950)), ((100, 150), (300, 400), (900, 950)), ((100, 150), (500, 600), (900, 950)), ((100, 150), (900, 950))}
    assert cassette[0].path_index[list(TRANSCRIPTS).index("short")] == -1
    # the global blocks are unchanged and come first
    plain = build_splice_blocks(gtf)
    assert [b.block_id for b in blocks[:len(plain)]] == [b.block_id for b in plain]
    assert build_splice_blocks(gtf, local_anchors=True, local_only=True) == local


def test_minus_strand_local_blocks(tmp_path: Path):
    gtf = tmp_path / "g.gtf"
    write_gtf(gtf, TRANSCRIPTS, strand="-")
    local = [b for b in build_splice_blocks(gtf, local_anchors=True) if ":L" in b.block_id]
    assert any(len(b.path_signatures) == 4 for b in local)


def test_local_rule_is_one_consistent_block_set(tmp_path: Path):
    gtf = tmp_path / "g.gtf"
    write_gtf(gtf, TRANSCRIPTS)
    blocks = build_splice_blocks(gtf, anchor_rule="local")
    assert all(":B" in b.block_id for b in blocks)
    contents = [{tuple(block_path_exons(None, s, None)) for s in b.path_signatures} for b in blocks]
    # the cassette combinations form one block, as with the added local blocks
    cassette = {((100, 150), (300, 400), (500, 600), (700, 800), (900, 950)), ((100, 150), (300, 400), (900, 950)), ((100, 150), (500, 600), (900, 950)), ((100, 150), (900, 950))}
    assert cassette in contents
    # no block spans the whole gene any more
    assert all(max(e for p in c for _, e in p) - min(a for p in c for a, _ in p) < 1500 for c in contents)
