import numpy as np

from tealeaf.sc.local_pseudoalign import canonical_kmers, kmer_read_opportunities

rng = np.random.default_rng(3)
ANCHOR_LEFT, CASSETTE, ANCHOR_RIGHT = (rng.choice(list(b"ACGT"), size=n).astype(np.uint8).tobytes() for n in (200, 100, 200))
INTRON = rng.choice(list(b"ACGT"), size=150).astype(np.uint8).tobytes()
SKIP = ANCHOR_LEFT + ANCHOR_RIGHT
INCLUDE = ANCHOR_LEFT + CASSETTE + ANCHOR_RIGHT
PRECURSOR = ANCHOR_LEFT + INTRON + CASSETTE + INTRON[::-1] + ANCHOR_RIGHT


def test_canonical_kmers_strand_invariant():
    sequence = rng.choice(list(b"ACGT"), size=60).astype(np.uint8).tobytes()
    complement = sequence.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]
    forward, _ = canonical_kmers(sequence)
    reverse, _ = canonical_kmers(complement)
    assert np.array_equal(np.sort(forward), np.sort(reverse))


def test_kmer_opportunities_classes():
    length, k = 60, 31
    opportunities = kmer_read_opportunities([SKIP, INCLUDE, PRECURSOR], length, k)
    starts = lambda size: size - k + 1 + (length - k)
    totals = sum(opportunities.values())
    assert np.allclose(totals, [starts(len(SKIP)), starts(len(INCLUDE)), starts(len(PRECURSOR))])
    # reads wholly in an anchor are compatible with all three targets
    assert opportunities[7][0] > 0 and opportunities[7][2] > 0
    # skip-junction reads are skip-only; cassette-body reads include + precursor
    assert opportunities[1][0] > 0 and opportunities[6][1] > 0
    # intronic reads are precursor-only
    assert opportunities[4][2] > 0 and opportunities[4][0] == 0
