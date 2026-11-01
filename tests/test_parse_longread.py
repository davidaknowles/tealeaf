import numpy as np

from tealeaf.sc import parse_longread as pl

TSO = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
ADAPTER = "CTACACGACGCTCTTCCGATCT"


def molecule(rng, words, umi, cdna, primer_tail):
    """Sense-strand S-read: TSO, cDNA, primer tail, then the read-2 barcode block reverse-complemented, adapter."""
    b3, b2, b1 = words
    read2 = ADAPTER + umi + b3 + pl.LINKER_3 + b2 + pl.LINKER_2 + b1
    return TSO + cdna + primer_tail + pl.reverse_complement(read2)


def test_call_reads_both_strands_and_errors():
    rng = np.random.default_rng(0)
    random = lambda n: "".join(rng.choice(list("ACGT"), n))
    whitelists = [[random(8) for _ in range(20)] for _ in range(3)]
    reads, truth = [], []
    for i in range(60):
        idx = rng.integers(20, size=3)
        words = [whitelists[k][idx[k]] for k in range(3)]
        umi, cdna = random(10), "G" + random(300) + "C"
        tail = "A" * 25 if i % 2 else ""
        observed = list(words)
        if i % 3 == 1:
            observed[0] = observed[0][:3] + observed[0][4:]  # deletion in bc3
        if i % 3 == 2:
            observed[2] = observed[2][:5] + "A" + observed[2][5:]  # insertion in bc1
        read = molecule(rng, observed, umi, cdna, tail)
        if i % 5 == 0:
            j = len(read) - len(ADAPTER) - 10 - 8 - 4
            read = read[:j] + {"A": "C", "C": "G", "G": "T", "T": "A"}[read[j]] + read[j + 1:]  # linker substitution
        if i % 4 == 0:
            read = pl.reverse_complement(read)
        reads.append(read)
        truth.append((idx, umi, cdna, len(tail), i % 3 == 2))
    calls = pl.call_reads(reads, *whitelists)
    assert (calls[:, 0] == 0).all(), [pl.STATUS[s] for s in calls[:, 0]]
    for read, row, (idx, umi, cdna, tail, bc1_indel) in zip(reads, calls, truth):
        assert tuple(row[2:5]) == tuple(idx) and row[12] == 1
        u, insert, quality = pl.extract(read, "I" * len(read), row)
        assert u == umi and len(quality) == len(insert)
        if bc1_indel:  # an indel at the bc1 end can equally be a substitution: insert end within 1 nt
            assert insert[:len(cdna) - 1] == cdna[:-1] and abs(len(insert) - len(cdna)) <= 1 and abs(row[11] - tail) <= 1
        else:
            assert insert == cdna and row[11] == tail


def test_call_reads_rejects_reads_without_barcodes():
    rng = np.random.default_rng(1)
    random = lambda n: "".join(rng.choice(list("ACGT"), n))
    whitelists = [[random(8) for _ in range(20)] for _ in range(3)]
    calls = pl.call_reads([TSO + random(500) for _ in range(50)], *whitelists)
    assert (calls[:, 0] == 1).mean() > .8 and not (calls[:, 0] == 0).any()
