import numpy as np
import pytest

from tealeaf.sc import native_align as na
from tealeaf.sc.local_path_reads import compatibility_mask, path_structures


class Genome:
    def __init__(self, sequences):
        self.sequences = sequences
        self.references = list(sequences)

    def fetch(self, chromosome, start=None, end=None):
        return self.sequences[chromosome][start:end]


def revcomp(sequence):
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


EXONS = [(1000, 1200), (3000, 3150), (5000, 5300), (8000, 8400)]
NOVEL_EXON = (4000, 4100)


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(1)
    chr1 = list("".join(rng.choice(list("ACGT"), 60000)))
    chr2 = list("".join(rng.choice(list("ACGT"), 60000)))
    for a, b in EXONS + [NOVEL_EXON]:
        chr1[a - 2:a] = "AG"
        chr1[b:b + 2] = "GT"
    chr1 = "".join(chr1)
    chr2 = "".join(chr2)
    # an intronic stretch of chr1 repeated outside the target loci
    chr2 = chr2[:50000] + chr1[2000:2300] + chr2[50300:]
    # a copy of the third exon inside the second target gene
    chr2 = chr2[:25000] + chr1[5000:5300] + chr2[25300:]
    genome = Genome({"chr1": chr1, "chr2": chr2})
    genes = {
        "g1": {"chromosome": "chr1", "strand": "+", "transcripts": {"t1": EXONS, "t2": [EXONS[0], EXONS[2], EXONS[3]]}},
        "g2": {"chromosome": "chr2", "strand": "+", "transcripts": {"t3": [(20000, 20500), (30000, 30500)]}},
    }
    reference = na.build_reference(genome, genes, pad=500)
    index = na.build_index(reference)
    unique = na.genome_uniqueness(reference, index, genome)
    return genome, reference, index, unique


def align(setup, sequences):
    _, reference, index, unique = setup
    reads = np.stack([na.encode(s) for s in sequences])
    return na.align_reads(reads, reference["ref"], index["index_codes"], index["index_positions"], index["buckets"], unique, na.K, reference["locus_ref_start"], reference["locus_ref_end"], reference["tx_space_start"], reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], reference["exon_ref_start"], reference["exon_length"], reference["junctions"])


def genomic(setup, n_blocks, blocks, read):
    _, reference, _, _ = setup
    out = []
    for i in range(n_blocks[read]):
        start, end = blocks[read, i]
        locus = int(np.searchsorted(reference["locus_ref_start"], start, side="right") - 1)
        out.append((na.locus_to_genome(reference, locus, start)[1], na.locus_to_genome(reference, locus, end)[1]))
    return out


def spliced(chromosome, segments):
    return "".join(chromosome[a:b] for a, b in segments)


def test_unspliced_spliced_and_reverse(setup):
    chr1 = setup[0].sequences["chr1"]
    reads = [chr1[8100:8251],
             spliced(chr1, [(1130, 1200), (3000, 3081)]),
             revcomp(spliced(chr1, [(1130, 1200), (3000, 3081)])),
             spliced(chr1, [(1054, 1200), (5000, 5005)])]
    status, orientation, _, _, mismatches, n_blocks, blocks = align(setup, reads)
    assert list(status) == [na.UNIQUE] * 4
    assert genomic(setup, n_blocks, blocks, 0) == [(8100, 8251)]
    assert genomic(setup, n_blocks, blocks, 1) == [(1130, 1200), (3000, 3081)]
    assert genomic(setup, n_blocks, blocks, 2) == [(1130, 1200), (3000, 3081)] and orientation[2] == 1
    assert genomic(setup, n_blocks, blocks, 3) == [(1054, 1200), (5000, 5005)]
    assert list(mismatches) == [0, 0, 0, 0]


def test_minimum_junction_overhang(setup):
    """Terminal blocks shorter than 3 bases (annotated junction) or 5 bases (novel) are clipped, as in STAR."""
    chr1 = setup[0].sequences["chr1"]
    reads = [spliced(chr1, [(1051, 1200), (3000, 3002)]), spliced(chr1, [(1053, 1200), (3000, 3004)]),
             spliced(chr1, [(1196, 1200), (3000, 3147)]), spliced(chr1, [(1053, 1200), (4000, 4004)])]
    status, _, _, _, _, n_blocks, blocks = align(setup, reads)
    assert list(status) == [na.UNIQUE] * 4
    assert n_blocks[0] == 1 and genomic(setup, n_blocks, blocks, 0)[0][0] == 1051
    assert genomic(setup, n_blocks, blocks, 1) == [(1053, 1200), (3000, 3004)]
    assert genomic(setup, n_blocks, blocks, 2) == [(1196, 1200), (3000, 3147)]
    assert n_blocks[3] == 1


def test_novel_exon_and_mismatches(setup):
    chr1 = setup[0].sequences["chr1"]
    read = spliced(chr1, [(1170, 1200), (4000, 4100), (5000, 5021)])
    noisy = list(chr1[8100:8251])
    for position in (10, 70, 120):
        noisy[position] = "A" if noisy[position] != "A" else "C"
    status, _, _, _, mismatches, n_blocks, blocks = align(setup, [read, "".join(noisy)])
    assert list(status) == [na.UNIQUE, na.UNIQUE]
    assert genomic(setup, n_blocks, blocks, 0) == [(1170, 1200), (4000, 4100), (5000, 5021)]
    assert mismatches[1] == 3


def test_repeats_multi_and_unmapped(setup):
    chr1, chr2 = setup[0].sequences["chr1"], setup[0].sequences["chr2"]
    rng = np.random.default_rng(5)
    reads = [chr1[2050:2201], chr1[5100:5251], "".join(rng.choice(list("ACGT"), 151)), chr1[8100:8170] + "".join(rng.choice(list("ACGT"), 81))]
    status, *_ = align(setup, reads)
    assert list(status) == [na.NOT_GENOME_UNIQUE, na.MULTI, na.UNMAPPED, na.TOO_SHORT]


def test_block_masks_match_compatibility_mask(setup):
    genome, reference, _, _ = setup
    chr1 = genome.sequences["chr1"]
    window = (EXONS[0][0], EXONS[2][1])
    paths = [[EXONS[0], EXONS[1], EXONS[2]], [EXONS[0], EXONS[2]], [window]]
    blocks = {"k": {"block_id": "g1:B1", "gene_id": "g1", "chromosome": "chr1", "signatures": [0, 0], "paths": paths}}
    keys, tables = na.block_tables(blocks, reference)
    reads = [chr1[1100:1251], spliced(chr1, [(1130, 1200), (3000, 3081)]), spliced(chr1, [(1130, 1200), (5000, 5081)]), chr1[3050:3201], chr1[8100:8251]]
    status, _, _, _, _, n_blocks, aligned = align(setup, reads)
    locus = np.zeros(len(reads), dtype=np.int64)
    offsets, starts, ends = na.flatten_blocks(n_blocks, aligned)
    read, block, mask, informative = na.block_masks(locus, offsets, starts, ends, tables["locus"], tables["window_start"], tables["window_end"], tables["order"], tables["path_offsets"], tables["exon_offsets"], tables["exon_start"], tables["exon_end"], tables["junction_offsets"], tables["junction_keys"], tables["variable_offsets"], tables["variable_keys"])
    structures = path_structures(paths)
    expected = {}
    for r in range(len(reads)):
        segments = genomic(setup, n_blocks, aligned, r)
        junctions = [(segments[i][1], segments[i + 1][0]) for i in range(len(segments) - 1)]
        value = compatibility_mask(segments, junctions, structures)
        if any(a < window[1] and b > window[0] for a, b in segments):
            expected[r] = value
    assert dict(zip(read.tolist(), mask.tolist())) == expected
    assert dict(zip(read.tolist(), informative.tolist())) == {0: False, 1: True, 2: True, 3: False}


def test_barcodes_and_umis():
    table = na.segment_table(["AAAAAAAA", "CCCCCCCC", "AAAACCCC"], max_distance=2)
    code = lambda word: int(na.pack(na.encode(word)[None, :])[0][0])
    assert table[code("AAAAAAAA")] == 0 and table[code("AAAAAAAT")] == 0
    assert table[code("AAAACCCA")] == 2 and table[code("AAAAACCC")] == 2 and table[code("AAAAAACC")] == -1
    groups = np.array([0, 0, 0, 1])
    umis = np.array([0, 1, 5, 1], dtype=np.int64)
    out = na.collapse_umis(groups, umis, 5)
    # 0 (AAAAA) and 1 (AAAAC) are one mismatch apart; 5 (AAACC) is one from 1: one component
    assert out.tolist() == [0, 0, 0, 1]


RUST_LIBRARY = __import__("pathlib").Path(__file__).resolve().parents[1] / "rust/native_align/target/release/libnative_align.so"


@pytest.mark.skipif(not RUST_LIBRARY.exists(), reason="Rust aligner not built (cargo build --release in rust/native_align)")
def test_rust_kernel_matches_numba(setup):
    genome, reference, index, unique = setup
    chr1, chr2 = genome.sequences["chr1"], genome.sequences["chr2"]
    rng = np.random.default_rng(9)
    sequences = [chr1[8100:8251], spliced(chr1, [(1130, 1200), (3000, 3081)]), revcomp(spliced(chr1, [(1054, 1200), (5000, 5005)])),
                 spliced(chr1, [(1170, 1200), (4000, 4100), (5000, 5021)]), chr1[2050:2201], chr1[5100:5251], "".join(rng.choice(list("ACGT"), 151)), chr2[20100:20251]]
    for _ in range(200):
        start = int(rng.integers(500, 9000))
        read = list(chr1[start:start + 151])
        for position in rng.integers(0, 151, 2):
            read[position] = "ACGT"[int(rng.integers(4))]
        sequences.append("".join(read))
    reads = np.stack([na.encode(s) for s in sequences])
    arguments = (reads, reference["ref"], index["index_codes"], index["index_positions"], index["buckets"], unique, na.K, reference["locus_ref_start"], reference["locus_ref_end"], reference["tx_space_start"], reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], reference["exon_ref_start"], reference["exon_length"], reference["junctions"])
    expected = na.align_reads(*arguments)
    observed = na.rust_aligner(RUST_LIBRARY, threads=2, chunk=16)(*arguments)
    for a, b in zip(expected, observed):
        np.testing.assert_array_equal(a, b)


def test_simulated_opportunities_track_geometry(setup):
    from tealeaf.sc.local_path_reads import path_read_opportunities

    genome, reference, index, unique = setup
    window = (EXONS[0][0], EXONS[3][1])
    paths = [[EXONS[0], EXONS[1], EXONS[3]], [EXONS[0], EXONS[3]], [window]]
    blocks = {"k": {"block_id": "g1:B1", "gene_id": "g1", "chromosome": "chr1", "strand": "+", "signatures": [0, 0], "paths": paths}}
    tables = (reference["ref"], index["index_codes"], index["index_positions"], index["buckets"], unique, na.K, reference["locus_ref_start"], reference["locus_ref_end"], reference["tx_space_start"], reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], reference["exon_ref_start"], reference["exon_length"], reference["junctions"])
    simulated, junction = na.simulated_opportunities(blocks, reference, lambda reads: na.align_reads(reads, *tables), 151, batch=5000)["k"]
    geometric = path_read_opportunities(paths, 151)
    masks = set(simulated) | set(geometric)
    difference = sum(np.abs(simulated.get(m, np.zeros(3)) - geometric.get(m, np.zeros(3))).sum() for m in masks)
    total = sum(v.sum() for v in geometric.values())
    assert difference / total < .02
    # the skip junction is variable; reads across it are junction-informative for path 2 only
    assert junction[2][1] > 0 and 2 in junction


def test_known_junctions_second_pass(setup):
    """A first-pass novel junction added with add_junctions allows annotated-length (3-base) overhangs."""
    genome, reference, index, unique = setup
    chr1 = genome.sequences["chr1"]
    read = spliced(chr1, [(1053, 1200), (4000, 4004)])
    junctions = na.add_junctions(reference, ["chr1", "chr2"], [1200, 1200], [4000, 4000])
    assert len(junctions) == len(reference["junctions"]) + 1
    second = dict(reference, junctions=junctions)
    status, _, _, _, _, n_blocks, blocks = na.align_reads(np.stack([na.encode(read)]), second["ref"], index["index_codes"], index["index_positions"], index["buckets"], unique, na.K, second["locus_ref_start"], second["locus_ref_end"], second["tx_space_start"], second["tx_ref_start"], second["tx_length"], second["exon_offsets"], second["exon_tx_start"], second["exon_ref_start"], second["exon_length"], second["junctions"])
    assert status[0] == na.UNIQUE and genomic(setup, n_blocks, blocks, 0) == [(1053, 1200), (4000, 4004)]


@pytest.mark.skipif(not RUST_LIBRARY.exists(), reason="Rust aligner not built (cargo build --release in rust/native_align)")
def test_rust_fastq_reader_matches_numpy(tmp_path):
    import gzip
    from extra_scripts.native_collect import SEGMENTS, UMI, numpy_batches

    rng = np.random.default_rng(3)
    words = ["".join(rng.choice(list("ACGT"), 8)) for _ in range(6)]
    tables = [na.segment_table(words) for _ in range(3)]
    size = len(words)
    lookup = np.full(size ** 3, -1, dtype=np.int64)
    lookup[rng.choice(size ** 3, 50, replace=False)] = np.arange(50)
    records = [[], []]
    for i in range(1000):
        cdna = "".join(rng.choice(list("ACGTN"), 151, p=[.24, .24, .24, .24, .04]))
        umi = "A" * 10 if i % 50 == 0 else "".join(rng.choice(list("ACGTN"), 10, p=[.245, .245, .245, .245, .02]))
        segments = [words[rng.integers(size)] for _ in range(3)]
        segments = [s[:3] + "T" + s[4:] if rng.random() < .2 else s for s in segments]
        barcode = umi + segments[0] + "G" * 30 + segments[1] + "C" * 22 + segments[2] + "A" * 65
        for mate, sequence in enumerate((cdna, barcode)):
            records[mate].append(f"@r{i}\n{sequence}\n+\n{'F' * len(sequence)}\n")
    paths = [tmp_path / f"x_{mate}.fastq.gz" for mate in (1, 2)]
    for path, lines in zip(paths, records):
        with gzip.open(path, "wt") as handle:
            handle.write("".join(lines))
    expected = list(numpy_batches(paths[0], paths[1], 300, tables, lookup, size))
    observed = list(na.rust_fastq_pairs(RUST_LIBRARY, paths[0], paths[1], tables, lookup, size, UMI, SEGMENTS, batch=300))
    assert [len(b[0]) for b in observed] == [300, 300, 300, 100]
    for a, b in zip(expected, observed):
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)
    assert sum((b[1] >= 0).sum() for b in observed) > 0 and sum((~b[3]).sum() for b in observed) > 20
    with gzip.open(paths[1], "wt") as handle:
        handle.write("".join(records[1][:999]))
    with pytest.raises(ValueError, match="different record counts"):
        list(na.rust_fastq_pairs(RUST_LIBRARY, paths[0], paths[1], tables, lookup, size, UMI, SEGMENTS, batch=300))


@pytest.mark.skipif(not RUST_LIBRARY.exists(), reason="Rust aligner not built (cargo build --release in rust/native_align)")
def test_v3_barcode_layout(tmp_path):
    import gzip
    from extra_scripts.native_collect import BARCODE_LAYOUTS, UMI, numpy_batches

    rng = np.random.default_rng(4)
    words = ["".join(rng.choice(list("ACGT"), 8)) for _ in range(5)]
    tables = [na.segment_table(words) for _ in range(3)]
    size = len(words)
    lookup = np.arange(size ** 3, dtype=np.int64)
    triples, records = [], [[], []]
    for i in range(200):
        triple = rng.integers(size, size=3)
        triples.append((triple[0] * size + triple[1]) * size + triple[2])
        barcode = "".join(rng.choice(list("ACGT"), 10)) + words[triple[0]] + "ATGAGGGGTCAG" + words[triple[1]] + "TCCAACCACCTC" + words[triple[2]]
        assert len(barcode) == 58
        for mate, sequence in enumerate(("".join(rng.choice(list("ACGT"), 64)), barcode)):
            records[mate].append(f"@r{i}\n{sequence}\n+\n{'F' * len(sequence)}\n")
    paths = [tmp_path / f"x_{mate}.fastq.gz" for mate in (1, 2)]
    for path, lines in zip(paths, records):
        with gzip.open(path, "wt") as handle:
            handle.write("".join(lines))
    layout = BARCODE_LAYOUTS["v3"]
    expected = list(numpy_batches(paths[0], paths[1], 1000, tables, lookup, size, layout))
    observed = list(na.rust_fastq_pairs(RUST_LIBRARY, paths[0], paths[1], tables, lookup, size, UMI, layout, batch=1000))
    np.testing.assert_array_equal(expected[0][1], triples)
    for x, y in zip(expected[0], observed[0]):
        np.testing.assert_array_equal(x, y)
