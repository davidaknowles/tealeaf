"""Native spliced alignment of single-cell reads to target gene loci.

A dependency-free replacement for the genome alignments (STARsolo) that feed
block-local path classes. Only the loci of the target genes are indexed, so
the aligner reports, for each read, a genomic alignment (aligned blocks and
splice junctions) at most at one target position, together with the same
uniqueness and quality filters that the genome alignments provide.

Reference. Each target gene span (optionally padded by `pad` bases; padding
by 1 kb changed under 0.2% of block molecules) is merged with
overlapping spans into a locus; loci and the spliced sequences of the genes'
annotated transcripts are concatenated, separated by runs of N, into one
2-bit coded reference. Every k-mer of the loci and every transcript k-mer
spanning an exon junction is indexed (sorted codes with a bucket table on the
leading BUCKET_BITS bits). Transcript k-mers let reads with short exon
overhangs or short annotated exons seed across junctions; locus k-mers let
unspliced (exonic or intronic) reads and reads across unannotated junctions
seed in each aligned segment.

Uniqueness. The genome has no repeat annotation, so multi-mapping is judged
from k-mer multiplicity: every indexed k-mer is counted over the whole genome
(canonical, both strands; genome_uniqueness). A read is kept only if its
alignment contains, inside one aligned block, a k-mer that matches the read
exactly and occurs once in the genome (or, for a transcript junction k-mer,
not at all). A read from a paralogous or repeated copy elsewhere cannot
contain such a k-mer without a sequencing error. A read is also dropped when
two target alignments at different positions score within one of the best.

Alignment (align_reads). Seeds at a fixed stride (and the read end) are
looked up for the read and its reverse complement. Transcript hits are
extended ungapped along the transcript and projected to genomic blocks. Locus
hits are grouped by diagonal and chained (increasing read and reference
position, within one locus); consecutive segments are joined by a splice
junction (gap >= MIN_INTRON, canonical motif or annotated junction, split
point minimising mismatches), a deletion or an insertion, and the ends are
extended with soft clipping (match +1, mismatch -1). Scores follow STAR's:
matches minus mismatches, junction penalties (GC/AG -4, AT/AC -8, annotated
+2), indel open and per-base penalties of -2, and a genomic-length penalty
ceil(-0.25 log2(span) - 0.5). Alignments need score and matched bases at
least 0.66 (L - 1), at most 10 mismatches, and at most 0.3 mismatches per
aligned base.

Coordinates in the kernels are reference positions; locus_to_genome maps
locus positions back to (chromosome, 0-based position).
"""

from collections import defaultdict

import numpy as np
import pandas as pd
from numba import njit, prange, types
from numba.typed import List

K = 23
BUCKET_BITS = 24
SEPARATOR = 300
MIN_INTRON = 21
MAX_INTRON = 1_000_000
MAX_INDEL = 6
MIN_OVERHANG = 5  # STAR alignSJoverhangMin: minimum terminal block length after a novel junction
MIN_OVERHANG_ANNOTATED = 3  # STAR alignSJDBoverhangMin, for annotated junctions
MAX_OCCURRENCES = 32
MAX_HITS = 1024
MAX_SEGMENTS = 128
MAX_BLOCKS = 8
MAX_CANDIDATES = 24
STRIDE = 4
MAX_MISMATCHES = 10

CODE = np.full(256, 4, dtype=np.uint8)
for _base, _value in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    CODE[_base] = _value

# status codes
UNMAPPED, UNIQUE, MULTI, NOT_GENOME_UNIQUE, TOO_SHORT = 0, 1, 2, 3, 4
# junction motifs as in STAR SJ.out.tab: 0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT
MOTIF_STRAND = np.array([0, 1, 2, 1, 2, 1, 2], dtype=np.int8)


def encode(sequence):
    """2-bit codes (N and other symbols as 4) of a bytes or str sequence."""
    if isinstance(sequence, str):
        sequence = sequence.encode()
    return CODE[np.frombuffer(sequence, dtype=np.uint8)]


# ---------------------------------------------------------------- reference


def target_loci(genes, pad=0):
    """Merged padded gene spans: list of (chromosome, start, end, [gene ids])."""
    spans = defaultdict(list)
    for gene, record in genes.items():
        chains = list(record["transcripts"].values())
        start = min(chain[0][0] for chain in chains)
        end = max(chain[-1][1] for chain in chains)
        spans[record["chromosome"]].append((max(0, start - pad), end + pad, gene))
    loci = []
    for chromosome in sorted(spans):
        for start, end, gene in sorted(spans[chromosome]):
            if loci and loci[-1][0] == chromosome and start <= loci[-1][2]:
                loci[-1][2] = max(loci[-1][2], end)
                loci[-1][3].append(gene)
            else:
                loci.append([chromosome, start, end, [gene]])
    return [tuple(locus) for locus in loci]


def build_reference(genome, genes, pad=0):
    """Concatenated 2-bit reference of target loci and annotated transcripts.

    genome is a pysam.FastaFile (or any object with fetch(chromosome, start,
    end) returning a str), genes a read_gtf_exons-style dict of target genes.
    Returns a dict of arrays: ref (uint8 codes), locus tables, transcript
    tables (exon offsets in transcript and reference coordinates), the start
    of transcript space and sorted annotated junction keys
    (donor << 32 | acceptor, reference coordinates).
    """
    loci = target_loci(genes, pad)
    pieces, position = [], 0
    locus_ref_start, locus_ref_end, locus_genome_start, locus_chromosome = [], [], [], []
    for chromosome, start, end, _ in loci:
        pieces.append(np.full(SEPARATOR, 4, dtype=np.uint8))
        position += SEPARATOR
        sequence = encode(genome.fetch(chromosome, start, end))
        locus_ref_start.append(position)
        locus_genome_start.append(start)
        locus_chromosome.append(chromosome)
        pieces.append(sequence)
        position += len(sequence)
        locus_ref_end.append(position)
    locus_index = {}
    for index, (chromosome, start, end, members) in enumerate(loci):
        for gene in members:
            locus_index[gene] = index
    tx_space_start = position
    tx_ref_start, tx_length, tx_gene, exon_offsets = [], [], [], [0]
    exon_tx_start, exon_ref_start, exon_length, junctions = [], [], [], set()
    for gene in sorted(genes):
        record = genes[gene]
        locus = locus_index[gene]
        shift = locus_ref_start[locus] - locus_genome_start[locus]
        for transcript in sorted(record["transcripts"]):
            chain = sorted(record["transcripts"][transcript])
            if len(chain) < 2:
                continue
            pieces.append(np.full(SEPARATOR, 4, dtype=np.uint8))
            position += SEPARATOR
            tx_ref_start.append(position)
            offset = 0
            for a, b in chain:
                pieces.append(encode(genome.fetch(record["chromosome"], a, b)))
                exon_tx_start.append(offset)
                exon_ref_start.append(a + shift)
                exon_length.append(b - a)
                offset += b - a
            for (_, b), (a, _) in zip(chain, chain[1:]):
                junctions.add(((b + shift) << 32) | (a + shift))
            exon_offsets.append(len(exon_tx_start))
            tx_length.append(offset)
            tx_gene.append(gene)
            position += offset
    pieces.append(np.full(SEPARATOR, 4, dtype=np.uint8))
    return {
        "ref": np.concatenate(pieces),
        "locus_ref_start": np.asarray(locus_ref_start, dtype=np.int64),
        "locus_ref_end": np.asarray(locus_ref_end, dtype=np.int64),
        "locus_genome_start": np.asarray(locus_genome_start, dtype=np.int64),
        "locus_chromosome": np.asarray(locus_chromosome),
        "locus_genes": np.asarray([";".join(members) for *_, members in loci]),
        "tx_space_start": np.int64(tx_space_start),
        "tx_ref_start": np.asarray(tx_ref_start, dtype=np.int64),
        "tx_length": np.asarray(tx_length, dtype=np.int64),
        "tx_gene": np.asarray(tx_gene),
        "exon_offsets": np.asarray(exon_offsets, dtype=np.int64),
        "exon_tx_start": np.asarray(exon_tx_start, dtype=np.int64),
        "exon_ref_start": np.asarray(exon_ref_start, dtype=np.int64),
        "exon_length": np.asarray(exon_length, dtype=np.int64),
        "junctions": np.asarray(sorted(junctions), dtype=np.int64),
    }


@njit(cache=True)
def _kmer_codes(seq, k):
    """Forward k-mer codes and validity for every start of seq."""
    n = len(seq) - k + 1
    codes = np.zeros(max(n, 0), dtype=np.uint64)
    valid = np.zeros(max(n, 0), dtype=np.bool_)
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    code = np.uint64(0)
    run = 0
    for i in range(len(seq)):
        c = seq[i]
        if c > 3:
            run = 0
            code = np.uint64(0)
        else:
            code = ((code << np.uint64(2)) | np.uint64(c)) & mask
            run += 1
        if i >= k - 1:
            codes[i - k + 1] = code
            valid[i - k + 1] = run >= k
    return codes, valid


@njit(cache=True)
def _reverse_complement_code(code, k):
    out = np.uint64(0)
    for _ in range(k):
        out = (out << np.uint64(2)) | (np.uint64(3) - (code & np.uint64(3)))
        code >>= np.uint64(2)
    return out


@njit(cache=True)
def _junction_positions(tx_ref_start, tx_length, exon_offsets, exon_tx_start, k):
    """Reference starts of transcript k-mers that span an exon junction."""
    total = 0
    for t in range(len(tx_ref_start)):
        total += max(tx_length[t] - k + 1, 0)
    out = np.empty(total, dtype=np.int64)
    n = 0
    for t in range(len(tx_ref_start)):
        first, last = exon_offsets[t], exon_offsets[t + 1]
        for e in range(first + 1, last):
            boundary = exon_tx_start[e]
            low = max(boundary - k + 1, 0)
            if n > 0 and out[n - 1] >= tx_ref_start[t] + low:
                low = out[n - 1] - tx_ref_start[t] + 1
            for p in range(low, min(boundary, tx_length[t] - k + 1)):
                out[n] = tx_ref_start[t] + p
                n += 1
    return out[:n]


def _bucket_table(codes, k):
    shift = np.uint64(2 * k - BUCKET_BITS)
    buckets = (codes >> shift).astype(np.int64)
    return np.searchsorted(buckets, np.arange(2 ** BUCKET_BITS + 1)).astype(np.int64)


def build_index(reference, k=K):
    """k-mer index: sorted codes, reference positions and bucket offsets."""
    ref = reference["ref"]
    codes, valid = _kmer_codes(ref, k)
    locus_end = int(reference["tx_space_start"])
    locus_positions = np.flatnonzero(valid[:locus_end])
    junction_positions = _junction_positions(reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], k)
    junction_positions = junction_positions[valid[junction_positions]]
    positions = np.concatenate([locus_positions, junction_positions])
    selected = codes[positions]
    order = np.argsort(selected, kind="stable")
    index_codes = selected[order]
    return {"index_codes": index_codes, "index_positions": positions[order].astype(np.int64), "buckets": _bucket_table(index_codes, k)}


@njit(cache=True)
def _count_genome(seq, unique_codes, buckets, counts, k, shift):
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    forward = np.uint64(0)
    reverse = np.uint64(0)
    top = np.uint64(2 * (k - 1))
    run = 0
    for i in range(len(seq)):
        c = seq[i]
        if c > 3:
            run = 0
            forward = np.uint64(0)
            reverse = np.uint64(0)
            continue
        forward = ((forward << np.uint64(2)) | np.uint64(c)) & mask
        reverse = (reverse >> np.uint64(2)) | (np.uint64(3 - c) << top)
        run += 1
        if run < k:
            continue
        code = forward if forward < reverse else reverse
        b = np.int64(code >> shift)
        low, high = buckets[b], buckets[b + 1]
        while low < high:
            middle = (low + high) // 2
            if unique_codes[middle] < code:
                low = middle + 1
            else:
                high = middle
        if low < buckets[b + 1] and unique_codes[low] == code and counts[low] < 255:
            counts[low] += 1


def genome_uniqueness(reference, index, genome, k=K):
    """uint8 flag per reference position: 1 if the k-mer starting there is genome-unique.

    Locus k-mers are unique when their canonical code occurs once in the
    genome; transcript junction k-mers when it does not occur at all.
    """
    ref = reference["ref"]
    codes = index["index_codes"]
    canonical = np.minimum(codes, np.array([_reverse_complement_code(code, k) for code in codes], dtype=np.uint64)) if len(codes) < 1000 else _canonical_all(codes, k)
    unique_codes, inverse = np.unique(canonical, return_inverse=True)
    buckets = _bucket_table(unique_codes, k)
    counts = np.zeros(len(unique_codes), dtype=np.uint8)
    shift = np.uint64(2 * k - BUCKET_BITS)
    for chromosome in genome.references:
        _count_genome(encode(genome.fetch(chromosome)), unique_codes, buckets, counts, k, shift)
    positions = index["index_positions"]
    in_loci = positions < reference["tx_space_start"]
    observed = counts[inverse]
    flags = np.zeros(len(ref), dtype=np.uint8)
    flags[positions] = np.where(in_loci, observed == 1, observed == 0)
    return flags


@njit(cache=True, parallel=True)
def _canonical_all(codes, k):
    out = np.empty(len(codes), dtype=np.uint64)
    for i in prange(len(codes)):
        r = _reverse_complement_code(codes[i], k)
        out[i] = codes[i] if codes[i] < r else r
    return out


# ---------------------------------------------------------------- alignment


@njit(cache=True)
def _lookup(code, index_codes, buckets, shift):
    b = np.int64(code >> shift)
    low, high = buckets[b], buckets[b + 1]
    end = high
    while low < high:
        middle = (low + high) // 2
        if index_codes[middle] < code:
            low = middle + 1
        else:
            high = middle
    first = low
    high = end
    while low < high:
        middle = (low + high) // 2
        if index_codes[middle] <= code:
            low = middle + 1
        else:
            high = middle
    return first, low


@njit(cache=True)
def _find(sorted_values, value):
    """Index of the last element <= value, or -1."""
    low, high = 0, len(sorted_values)
    while low < high:
        middle = (low + high) // 2
        if sorted_values[middle] <= value:
            low = middle + 1
        else:
            high = middle
    return low - 1


@njit(cache=True)
def _contains(sorted_values, value):
    i = _find(sorted_values, value)
    return i >= 0 and sorted_values[i] == value


@njit(cache=True)
def _motif(ref, donor, acceptor):
    """STAR motif code of the intron [donor, acceptor) on the forward reference."""
    a, b, c, d = ref[donor], ref[donor + 1], ref[acceptor - 2], ref[acceptor - 1]
    if a == 2 and b == 3 and c == 0 and d == 2:
        return 1  # GT/AG
    if a == 1 and b == 3 and c == 0 and d == 1:
        return 2  # CT/AC
    if a == 2 and b == 1 and c == 0 and d == 2:
        return 3  # GC/AG
    if a == 1 and b == 3 and c == 2 and d == 1:
        return 4  # CT/GC
    if a == 0 and b == 3 and c == 0 and d == 1:
        return 5  # AT/AC
    if a == 2 and b == 3 and c == 0 and d == 3:
        return 6  # GT/AT
    return 0


@njit(cache=True)
def junction_motifs(ref, donors, acceptors):
    """STAR motif codes (_motif) of introns [donors[i], acceptors[i])."""
    out = np.zeros(len(donors), dtype=np.int8)
    for i in range(len(donors)):
        out[i] = _motif(ref, donors[i], acceptors[i])
    return out


@njit(cache=True)
def _junction_score(motif, annotated):
    if annotated:
        return 2
    if motif == 1 or motif == 2:
        return 0
    if motif == 3 or motif == 4:
        return -4
    return -8


@njit(cache=True)
def _mismatch(q, r, ref, diagonal):
    """1 if read base r mismatches the reference on this diagonal, 0 otherwise (read N is neutral)."""
    b = q[r]
    if b > 3:
        return 0
    return 1 if ref[diagonal + r] != b else 0


@njit(cache=True)
def _extend_left(q, ref, diagonal, start, floor, budget):
    """Leftmost read start of a soft-clipped extension from read position start (exclusive) down to floor, with at most budget mismatches."""
    best, score, best_start = 0, 0, start
    r = start - 1
    while r >= floor:
        if ref[diagonal + r] > 3:
            break
        if _mismatch(q, r, ref, diagonal):
            if budget == 0:
                break
            budget -= 1
            score -= 1
        elif q[r] <= 3:
            score += 1
        if score > best:
            best, best_start = score, r
        if score < best - 8:
            break
        r -= 1
    return best_start


@njit(cache=True)
def _extend_right(q, ref, diagonal, end, ceiling, budget):
    best, score, best_end = 0, 0, end
    r = end
    while r < ceiling:
        if ref[diagonal + r] > 3:
            break
        if _mismatch(q, r, ref, diagonal):
            if budget == 0:
                break
            budget -= 1
            score -= 1
        elif q[r] <= 3:
            score += 1
        if score > best:
            best, best_end = score, r + 1
        if score < best - 8:
            break
        r += 1
    return best_end


@njit(cache=True)
def _score_blocks(q, ref, read_starts, read_ends, diagonals, n):
    """(matches, mismatches) of n ungapped read segments."""
    matches, mismatches = 0, 0
    for i in range(n):
        for r in range(read_starts[i], read_ends[i]):
            if q[r] > 3:
                continue
            if ref[diagonals[i] + r] == q[r]:
                matches += 1
            else:
                mismatches += 1
    return matches, mismatches


@njit(cache=True)
def _genome_unique(q, ref, unique, read_starts, read_ends, diagonals, n, k):
    """True if some aligned block contains an exactly matching genome-unique k-mer."""
    for i in range(n):
        run = 0
        for r in range(read_starts[i], read_ends[i]):
            if q[r] <= 3 and ref[diagonals[i] + r] == q[r]:
                run += 1
            else:
                run = 0
            if run >= k and unique[diagonals[i] + r - k + 1]:
                return True
    return False


@njit(cache=True)
def _transcript_candidate(q, ref, diagonal, seed_start, seed_end, tx_start, tx_end, t, exon_offsets, exon_tx_start, exon_ref_start, exon_length, read_starts, read_ends, diagonals):
    """Ungapped transcript extension projected to genomic segments; returns segment count (0 if none)."""
    length = len(q)
    floor = max(0, tx_start - diagonal)
    ceiling = min(length, tx_end - diagonal)
    start = _extend_left(q, ref, diagonal, seed_start, floor, MAX_MISMATCHES)
    spent = 0
    for r in range(start, seed_end):
        spent += _mismatch(q, r, ref, diagonal)
    end = _extend_right(q, ref, diagonal, seed_end, ceiling, max(0, MAX_MISMATCHES - spent))
    n = 0
    for e in range(exon_offsets[t], exon_offsets[t + 1]):
        e_start = tx_start + exon_tx_start[e]
        e_end = e_start + exon_length[e]
        low = max(diagonal + start, e_start)
        high = min(diagonal + end, e_end)
        if low >= high:
            continue
        if n >= MAX_BLOCKS:
            return 0
        read_starts[n] = low - diagonal
        read_ends[n] = high - diagonal
        diagonals[n] = exon_ref_start[e] - e_start + diagonal
        n += 1
    return n


@njit(cache=True)
def _join(q, ref, a_diagonal, b_diagonal, low, high, junctions, locus_start, locus_end):
    """Best split s in [low, high] joining diagonals a and b; returns (s, mismatches, motif, annotated, ok)."""
    best_s, best_mm, best_motif, best_annotated, best_rank = -1, 1 << 30, 0, False, -1
    gap = b_diagonal - a_diagonal
    for s in range(low, high + 1):
        mm = 0
        for r in range(low, s):
            mm += _mismatch(q, r, ref, a_diagonal)
        for r in range(s, high):
            mm += _mismatch(q, r, ref, b_diagonal)
        motif, annotated, rank = 0, False, 0
        if gap >= MIN_INTRON:
            donor = a_diagonal + s
            acceptor = b_diagonal + s
            if donor < locus_start or acceptor > locus_end:
                continue
            annotated = _contains(junctions, (donor << 32) | acceptor)
            motif = _motif(ref, donor, acceptor)
            if not annotated and motif == 0:
                continue
            rank = 3 if annotated else (2 if motif <= 2 else 1)
        if mm < best_mm or (mm == best_mm and rank > best_rank):
            best_s, best_mm, best_motif, best_annotated, best_rank = s, mm, motif, annotated, rank
    return best_s, best_mm, best_motif, best_annotated, best_s >= 0


@njit(cache=True)
def _rescue_right(q, ref, diagonal, block_start, end, acceptors, junctions, locus_end):
    """Best junction from the clipped right end to an annotated acceptor.

    Returns (split, acceptor, delta score) with delta <= 0 when nothing improves
    on soft clipping. The tail [split, len(q)) aligns ungapped from the acceptor.
    """
    length = len(q)
    best_split, best_acceptor, best_delta = -1, -1, 0
    if length - end < MIN_OVERHANG_ANNOTATED:
        return best_split, best_acceptor, best_delta
    first = _find(acceptors, diagonal + end + MIN_INTRON - 1) + 1
    for i in range(first, len(acceptors)):
        acceptor = acceptors[i]
        if acceptor > diagonal + end + MAX_INTRON or acceptor + (length - end) > locus_end:
            break
        for split in range(max(block_start + 1, end - 6), end + 1):
            donor = diagonal + split
            annotated = _contains(junctions, (donor << 32) | acceptor)
            motif = _motif(ref, donor, acceptor)
            if not annotated and (motif == 0 or length - split < MIN_OVERHANG):
                continue
            gain, lost, tail_mm = 0, 0, 0
            for r in range(split, length):
                b = q[r]
                if b > 3:
                    continue
                if ref[acceptor + r - split] == b:
                    gain += 1
                else:
                    gain -= 1
                    tail_mm += 1
            for r in range(split, end):
                b = q[r]
                if b <= 3:
                    lost += 1 if ref[diagonal + r] == b else -1
            if tail_mm > max(1, (length - split) // 10):
                continue
            delta = gain - lost + _junction_score(motif, annotated)
            if delta > best_delta:
                best_split, best_acceptor, best_delta = split, acceptor, delta
    return best_split, best_acceptor, best_delta


@njit(cache=True)
def _rescue_left(q, ref, diagonal, start, block_end, donors, junctions, locus_start):
    """Best junction from an annotated donor into the clipped left end; returns (split, donor, delta)."""
    best_split, best_donor, best_delta = -1, -1, 0
    if start < MIN_OVERHANG_ANNOTATED:
        return best_split, best_donor, best_delta
    last = _find(donors, diagonal + start - MIN_INTRON + 1)
    i = last
    while i >= 0:
        donor = donors[i]
        if donor < diagonal + start - MAX_INTRON or donor - start < locus_start:
            break
        for split in range(start, min(start + 6, block_end - 1) + 1):
            acceptor = diagonal + split
            annotated = _contains(junctions, (donor << 32) | acceptor)
            motif = _motif(ref, donor, acceptor)
            if not annotated and (motif == 0 or split < MIN_OVERHANG):
                continue
            gain, lost, head_mm = 0, 0, 0
            for r in range(split):
                b = q[r]
                if b > 3:
                    continue
                if ref[donor - split + r] == b:
                    gain += 1
                else:
                    gain -= 1
                    head_mm += 1
            for r in range(start, split):
                b = q[r]
                if b <= 3:
                    lost += 1 if ref[diagonal + r] == b else -1
            if head_mm > max(1, split // 10):
                continue
            delta = gain - lost + _junction_score(motif, annotated)
            if delta > best_delta:
                best_split, best_donor, best_delta = split, donor, delta
        i -= 1
    return best_split, best_donor, best_delta


@njit(cache=True)
def _locus_of(locus_ref_start, locus_ref_end, position):
    i = _find(locus_ref_start, position)
    if i >= 0 and position < locus_ref_end[i]:
        return i
    return -1


@njit(cache=True)
def _align_one(q, ref, index_codes, index_positions, buckets, shift, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions, donors, acceptors, hits_diag, hits_read, seg_diag, seg_start, seg_end, seg_n, cand_score, cand_matches, cand_mm, cand_start, cand_end, cand_n, cand_rs, cand_re, cand_dg, cand_unique, n_cand, work_rs, work_re, work_dg, order, cand_orient, o):
    """Add alignment candidates of read orientation o (sequence q); returns the new candidate count."""
    length = len(q)
    n_hits = 0
    p = 0
    while p <= length - k:
        code = np.uint64(0)
        good = True
        for i in range(k):
            c = q[p + i]
            if c > 3:
                good = False
                break
            code = (code << np.uint64(2)) | np.uint64(c)
        if good:
            first, last = _lookup(code, index_codes, buckets, shift)
            if 0 < last - first <= MAX_OCCURRENCES:
                for j in range(first, last):
                    if n_hits < MAX_HITS:
                        hits_diag[n_hits] = index_positions[j] - p
                        hits_read[n_hits] = p
                        n_hits += 1
        if p == length - k:
            break
        p = min(p + STRIDE, length - k)
    if n_hits == 0:
        return n_cand
    keys = hits_diag[:n_hits] * 256 + hits_read[:n_hits]
    sorted_hits = np.argsort(keys, kind="mergesort")
    n_seg = 0
    for h in sorted_hits:
        d, r = hits_diag[h], hits_read[h]
        if n_seg > 0 and seg_diag[n_seg - 1] == d:
            seg_end[n_seg - 1] = max(seg_end[n_seg - 1], r + k)
            seg_n[n_seg - 1] += 1
        elif n_seg < MAX_SEGMENTS:
            seg_diag[n_seg], seg_start[n_seg], seg_end[n_seg], seg_n[n_seg] = d, r, r + k, 1
            n_seg += 1
    # transcript-space candidates
    for s in range(n_seg):
        d = seg_diag[s]
        position = d + seg_start[s]
        if position < tx_space_start:
            continue
        t = _find(tx_ref_start, position)
        if t < 0 or position >= tx_ref_start[t] + tx_length[t]:
            continue
        n = _transcript_candidate(q, ref, d, seg_start[s], seg_end[s], tx_ref_start[t], tx_ref_start[t] + tx_length[t], t, exon_offsets, exon_tx_start, exon_ref_start, exon_length, work_rs, work_re, work_dg)
        if n == 0:
            continue
        n_cand = _add_candidate(q, ref, unique, k, junctions, work_rs, work_re, work_dg, n, cand_score, cand_matches, cand_mm, cand_start, cand_end, cand_n, cand_rs, cand_re, cand_dg, cand_unique, n_cand, cand_orient, o)
    # locus-space chaining
    n_locus = 0
    for s in range(n_seg):
        if seg_diag[s] + seg_start[s] < tx_space_start:
            order[n_locus] = s
            n_locus += 1
    if n_locus == 0:
        return n_cand
    by_start = order[:n_locus][np.argsort(seg_start[order[:n_locus]], kind="mergesort")]
    dp = np.zeros(n_locus, dtype=np.int64)
    previous = np.full(n_locus, -1, dtype=np.int64)
    loci = np.empty(n_locus, dtype=np.int64)
    for j in range(n_locus):
        sj = by_start[j]
        loci[j] = _locus_of(locus_ref_start, locus_ref_end, seg_diag[sj] + seg_start[sj])
        dp[j] = seg_end[sj] - seg_start[sj]
        for i in range(j):
            si = by_start[i]
            if loci[i] != loci[j] or loci[j] < 0:
                continue
            gap = seg_diag[sj] - seg_diag[si]
            if gap == 0 or gap < -MAX_INDEL or gap > MAX_INTRON or (0 < gap < MIN_INTRON and gap > MAX_INDEL):
                continue
            if seg_start[si] >= seg_start[sj] or seg_end[si] > seg_end[sj]:
                continue
            value = dp[i] + seg_end[sj] - max(seg_start[sj], seg_end[si])
            if value > dp[j]:
                dp[j], previous[j] = value, i
    used = np.zeros(n_locus, dtype=np.bool_)
    ends = np.argsort(-dp, kind="mergesort")
    for e in ends:
        if used[e] or dp[e] < 2 * k // 2:
            continue
        chain_len = 0
        j = e
        while j >= 0:
            used[j] = True
            chain_len += 1
            j = previous[j]
        chain = np.empty(chain_len, dtype=np.int64)
        j = e
        for i in range(chain_len - 1, -1, -1):
            chain[i] = by_start[j]
            j = previous[j]
        locus = loci[e]
        l_start, l_end = locus_ref_start[locus], locus_ref_end[locus]
        n = 0
        first = chain[0]
        work_dg[0] = seg_diag[first]
        work_rs[0] = _extend_left(q, ref, seg_diag[first], seg_start[first], max(0, l_start - seg_diag[first]), MAX_MISMATCHES)
        ok = True
        for c in range(1, chain_len):
            a, b = chain[c - 1], chain[c]
            gap = seg_diag[b] - seg_diag[a]
            low = min(seg_end[a], seg_start[b])
            high = max(seg_end[a], seg_start[b])
            low = max(low - 4, seg_start[a] + 1, work_rs[n] + 1)
            high = min(high + 4, seg_end[b] - 1)
            if gap < 0:
                high = min(high, seg_end[b] + gap - 1)
            if low > high:
                ok = False
                break
            split, mm, motif, annotated, good = _join(q, ref, seg_diag[a], seg_diag[b], low, high, junctions, l_start, l_end)
            if not good or n + 1 >= MAX_BLOCKS:
                ok = False
                break
            work_re[n] = split
            n += 1
            work_rs[n] = split - gap if gap < 0 else split
            work_dg[n] = seg_diag[b]
        if not ok:
            continue
        last = chain[chain_len - 1]
        work_re[n] = seg_end[last]
        spent = _score_blocks(q, ref, work_rs, work_re, work_dg, n + 1)[1]
        work_re[n] = _extend_right(q, ref, seg_diag[last], seg_end[last], min(length, l_end - seg_diag[last]), max(0, MAX_MISMATCHES - spent))
        n += 1
        if n < MAX_BLOCKS:
            split, acceptor, delta = _rescue_right(q, ref, work_dg[n - 1], work_rs[n - 1], work_re[n - 1], acceptors, junctions, l_end)
            if delta > 0:
                work_re[n - 1] = split
                work_rs[n] = split
                work_re[n] = length
                work_dg[n] = acceptor - split
                n += 1
        if n < MAX_BLOCKS:
            split, donor, delta = _rescue_left(q, ref, work_dg[0], work_rs[0], work_re[0], donors, junctions, l_start)
            if delta > 0:
                for i in range(n, 0, -1):
                    work_rs[i], work_re[i], work_dg[i] = work_rs[i - 1], work_re[i - 1], work_dg[i - 1]
                work_rs[1] = split
                work_rs[0] = 0
                work_re[0] = split
                work_dg[0] = donor - split
                n += 1
        valid = True
        for i in range(n):
            if work_re[i] <= work_rs[i]:
                valid = False
        if not valid:
            continue
        n_cand = _add_candidate(q, ref, unique, k, junctions, work_rs, work_re, work_dg, n, cand_score, cand_matches, cand_mm, cand_start, cand_end, cand_n, cand_rs, cand_re, cand_dg, cand_unique, n_cand, cand_orient, o)
    return n_cand


@njit(cache=True)
def _trim_overhangs(ref, junctions, rs, re, dg, n):
    """First and last kept block after dropping terminal blocks shorter than the minimum junction overhang (STAR's 3 for annotated, 5 for novel junctions)."""
    first, last = 0, n - 1
    while last > first and dg[first + 1] - dg[first] >= MIN_INTRON:
        donor, acceptor = dg[first] + re[first], dg[first + 1] + rs[first + 1]
        minimum = MIN_OVERHANG_ANNOTATED if _contains(junctions, (donor << 32) | acceptor) else MIN_OVERHANG
        if re[first] - rs[first] >= minimum:
            break
        first += 1
    while last > first and dg[last] - dg[last - 1] >= MIN_INTRON:
        donor, acceptor = dg[last - 1] + re[last - 1], dg[last] + rs[last]
        minimum = MIN_OVERHANG_ANNOTATED if _contains(junctions, (donor << 32) | acceptor) else MIN_OVERHANG
        if re[last] - rs[last] >= minimum:
            break
        last -= 1
    return first, last


@njit(cache=True)
def _gap_score(ref, junctions, rs, re, dg, n):
    """Junction scores plus indel penalties between consecutive blocks."""
    score = 0
    for i in range(n - 1):
        gap = dg[i + 1] - dg[i]
        if gap >= MIN_INTRON:
            donor, acceptor = dg[i] + re[i], dg[i + 1] + rs[i + 1]
            score += _junction_score(_motif(ref, donor, acceptor), _contains(junctions, (donor << 32) | acceptor))
        else:
            score -= 2 + 2 * abs(gap)
    return score


@njit(cache=True)
def _add_candidate(q, ref, unique, k, junctions, rs, re, dg, n, cand_score, cand_matches, cand_mm, cand_start, cand_end, cand_n, cand_rs, cand_re, cand_dg, cand_unique, n_cand, cand_orient, o):
    """Store a candidate; a full table drops its lowest-scoring candidate if the new one scores higher."""
    first, last = _trim_overhangs(ref, junctions, rs, re, dg, n)
    rs, re, dg, n = rs[first:last + 1], re[first:last + 1], dg[first:last + 1], last + 1 - first
    matches, mismatches = _score_blocks(q, ref, rs, re, dg, n)
    span = dg[n - 1] + re[n - 1] - (dg[0] + rs[0])
    length_penalty = min(0, int(np.ceil(np.log2(max(span, 1)) * -0.25 - 0.5)))
    score = matches - mismatches + _gap_score(ref, junctions, rs, re, dg, n) + length_penalty
    start, end = dg[0] + rs[0], dg[n - 1] + re[n - 1]
    for c in range(n_cand):
        if cand_start[c] == start and cand_end[c] == end and cand_n[c] == n and cand_orient[c] == o:
            if score > cand_score[c]:
                cand_score[c] = score
            return n_cand
    j = n_cand
    if n_cand >= MAX_CANDIDATES:
        j = 0
        for c in range(1, n_cand):
            if cand_score[c] < cand_score[j]:
                j = c
        if cand_score[j] >= score:
            return n_cand
    cand_score[j], cand_matches[j], cand_mm[j] = score, matches, mismatches
    cand_start[j], cand_end[j], cand_n[j], cand_orient[j] = start, end, n, o
    for i in range(n):
        cand_rs[j, i], cand_re[j, i], cand_dg[j, i] = rs[i], re[i], dg[i]
    cand_unique[j] = _genome_unique(q, ref, unique, rs, re, dg, n, k)
    return max(n_cand, j + 1)


def align_reads(reads, ref, index_codes, index_positions, buckets, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions, chunk=4096):
    """Align reads (N x L uint8 codes) to the target reference.

    Returns status (N), orientation (N, 1 = reverse complement), score, matches,
    mismatches, block counts and blocks (N x MAX_BLOCKS x 2, reference
    coordinates of aligned bases, half-open). Status is UNMAPPED, UNIQUE,
    MULTI (another position within one of the best score), NOT_GENOME_UNIQUE
    or TOO_SHORT. Clipped read ends of at least 5 bases are tried across a
    junction to the annotated acceptors (right end) or from the annotated
    donors (left end) of the locus.
    """
    donors = np.unique(junctions >> 32)
    acceptors = np.unique(junctions & ((1 << 32) - 1))
    return _align_reads(reads, ref, index_codes, index_positions, buckets, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions, donors, acceptors, chunk)


@njit(cache=True, parallel=True)
def _align_reads(reads, ref, index_codes, index_positions, buckets, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions, donors, acceptors, chunk):
    n_reads, length = reads.shape
    shift = np.uint64(2 * k - BUCKET_BITS)
    status = np.zeros(n_reads, dtype=np.int8)
    orientation = np.zeros(n_reads, dtype=np.int8)
    scores = np.zeros(n_reads, dtype=np.int16)
    matched = np.zeros(n_reads, dtype=np.int16)
    mismatched = np.zeros(n_reads, dtype=np.int16)
    n_blocks = np.zeros(n_reads, dtype=np.int8)
    blocks = np.zeros((n_reads, MAX_BLOCKS, 2), dtype=np.int64)
    minimum = int(0.66 * (length - 1))
    n_chunks = (n_reads + chunk - 1) // chunk
    for c in prange(n_chunks):
        hits_diag = np.empty(MAX_HITS, dtype=np.int64)
        hits_read = np.empty(MAX_HITS, dtype=np.int64)
        seg_diag = np.empty(MAX_SEGMENTS, dtype=np.int64)
        seg_start = np.empty(MAX_SEGMENTS, dtype=np.int64)
        seg_end = np.empty(MAX_SEGMENTS, dtype=np.int64)
        seg_n = np.empty(MAX_SEGMENTS, dtype=np.int64)
        order = np.empty(MAX_SEGMENTS, dtype=np.int64)
        cand_score = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_matches = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_mm = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_start = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_end = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_n = np.empty(MAX_CANDIDATES, dtype=np.int64)
        cand_rs = np.empty((MAX_CANDIDATES, MAX_BLOCKS), dtype=np.int64)
        cand_re = np.empty((MAX_CANDIDATES, MAX_BLOCKS), dtype=np.int64)
        cand_dg = np.empty((MAX_CANDIDATES, MAX_BLOCKS), dtype=np.int64)
        cand_unique = np.empty(MAX_CANDIDATES, dtype=np.bool_)
        cand_orient = np.empty(MAX_CANDIDATES, dtype=np.int8)
        work_rs = np.empty(MAX_BLOCKS, dtype=np.int64)
        work_re = np.empty(MAX_BLOCKS, dtype=np.int64)
        work_dg = np.empty(MAX_BLOCKS, dtype=np.int64)
        reverse = np.empty(length, dtype=np.uint8)
        for read in range(c * chunk, min((c + 1) * chunk, n_reads)):
            q = reads[read]
            for i in range(length):
                b = q[length - 1 - i]
                reverse[i] = 3 - b if b <= 3 else 4
            n_cand = 0
            for o in range(2):
                query = q if o == 0 else reverse
                n_cand = _align_one(query, ref, index_codes, index_positions, buckets, shift, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions, donors, acceptors, hits_diag, hits_read, seg_diag, seg_start, seg_end, seg_n, cand_score, cand_matches, cand_mm, cand_start, cand_end, cand_n, cand_rs, cand_re, cand_dg, cand_unique, n_cand, work_rs, work_re, work_dg, order, cand_orient, o)
            if n_cand == 0:
                continue
            best = 0
            for j in range(1, n_cand):
                if cand_score[j] > cand_score[best]:
                    best = j
            if cand_score[best] < minimum or cand_matches[best] < minimum or cand_mm[best] > MAX_MISMATCHES or cand_mm[best] > 0.3 * (cand_matches[best] + cand_mm[best]):
                status[read] = TOO_SHORT
                continue
            multi = False
            for j in range(n_cand):
                if j != best and cand_score[j] >= cand_score[best] - 1 and (cand_end[j] <= cand_start[best] or cand_start[j] >= cand_end[best] or cand_orient[j] != cand_orient[best]):
                    multi = True
            orientation[read] = cand_orient[best]
            scores[read], matched[read], mismatched[read] = cand_score[best], cand_matches[best], cand_mm[best]
            n_blocks[read] = cand_n[best]
            for i in range(cand_n[best]):
                blocks[read, i, 0] = cand_dg[best, i] + cand_rs[best, i]
                blocks[read, i, 1] = cand_dg[best, i] + cand_re[best, i]
            if multi:
                status[read] = MULTI
            elif not cand_unique[best]:
                status[read] = NOT_GENOME_UNIQUE
            else:
                status[read] = UNIQUE
    return status, orientation, scores, matched, mismatched, n_blocks, blocks


def rust_aligner(library, threads=None, chunk=4096):
    """align_reads backed by the Rust port of the kernel (rust/native_align, a C-ABI cdylib).

    library is the path of the built shared object. The returned function has
    the signature and outputs of align_reads.
    """
    import ctypes
    import os

    handle = ctypes.CDLL(str(library))
    function = handle.align_reads
    function.restype = ctypes.c_int32
    threads = threads or int(os.environ.get("NUMBA_NUM_THREADS", os.cpu_count()))

    def pointer(array):
        return array.ctypes.data_as(ctypes.c_void_p)

    def run(reads, ref, index_codes, index_positions, buckets, unique, k, locus_ref_start, locus_ref_end, tx_space_start, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions):
        arrays = [np.ascontiguousarray(a, dtype=t) for a, t in ((reads, np.uint8), (ref, np.uint8), (index_codes, np.uint64), (index_positions, np.int64), (buckets, np.int64), (unique, np.uint8), (locus_ref_start, np.int64), (locus_ref_end, np.int64), (tx_ref_start, np.int64), (tx_length, np.int64), (exon_offsets, np.int64), (exon_tx_start, np.int64), (exon_ref_start, np.int64), (exon_length, np.int64), (junctions, np.int64))]
        reads, ref, index_codes, index_positions, buckets, unique, locus_ref_start, locus_ref_end, tx_ref_start, tx_length, exon_offsets, exon_tx_start, exon_ref_start, exon_length, junctions = arrays
        donors = np.unique(junctions >> 32)
        acceptors = np.unique(junctions & ((1 << 32) - 1))
        n_reads, length = reads.shape
        outputs = (np.zeros(n_reads, dtype=np.int8), np.zeros(n_reads, dtype=np.int8), np.zeros(n_reads, dtype=np.int16), np.zeros(n_reads, dtype=np.int16), np.zeros(n_reads, dtype=np.int16), np.zeros(n_reads, dtype=np.int8), np.zeros((n_reads, MAX_BLOCKS, 2), dtype=np.int64))
        size = ctypes.c_size_t
        function(size(n_reads), size(length), pointer(reads), pointer(ref), size(len(ref)), pointer(index_codes), pointer(index_positions), size(len(index_codes)), pointer(buckets), pointer(unique), size(k),
                 pointer(locus_ref_start), pointer(locus_ref_end), size(len(locus_ref_start)), ctypes.c_int64(int(tx_space_start)), pointer(tx_ref_start), pointer(tx_length), size(len(tx_ref_start)),
                 pointer(exon_offsets), pointer(exon_tx_start), pointer(exon_ref_start), pointer(exon_length), size(len(exon_tx_start)),
                 pointer(junctions), size(len(junctions)), pointer(donors), size(len(donors)), pointer(acceptors), size(len(acceptors)),
                 size(threads), size(chunk), *(pointer(a) for a in outputs))
        return outputs

    return run


def rust_fastq_pairs(library, read1, read2, tables, lookup, size, umi, segments, batch=1 << 22):
    """Yield (read-1 codes, cell, UMI code, UMI valid) per batch of a pair of gzipped FASTQ files.

    Reading is done by the Rust library (rust/native_align, fastq.rs): one
    background thread per mate decompresses and parses ahead of the caller.
    Read-1 codes are n x L uint8 as CODE; read 2 is reduced to the cell index
    lookup[corrected segment triple] (-1 when a segment does not correct or
    the triple is not in lookup), the UMI code (pack) and whether the UMI has
    no N and is not a homopolymer. tables are the three segment correction
    tables (segment_table), size the whitelist size, umi and segments the
    (start, end) of the UMI and the three barcode segments in read 2.
    """
    import ctypes

    handle = ctypes.CDLL(str(library))
    handle.fastq_pair_open.restype = ctypes.c_void_p
    handle.fastq_pair_peek.restype = ctypes.c_int64
    handle.fastq_pair_error.restype = ctypes.c_char_p
    for name in ("fastq_pair_peek", "fastq_pair_take", "fastq_pair_error", "fastq_pair_close"):
        getattr(handle, name).argtypes = [ctypes.c_void_p] + {"fastq_pair_peek": [ctypes.c_void_p], "fastq_pair_take": [ctypes.c_void_p] * 4}.get(name, [])
    segment_length = segments[0][1] - segments[0][0]
    stacked = np.ascontiguousarray(np.concatenate(tables), dtype=np.int64)
    lookup = np.ascontiguousarray(lookup, dtype=np.int64)
    starts = np.array([a for a, _ in segments], dtype=np.uintp)
    ends = np.array([b for _, b in segments], dtype=np.uintp)
    if len(stacked) != 3 * 4 ** segment_length or len(lookup) != size ** 3 or any(b - a != segment_length for a, b in segments):
        raise ValueError("barcode tables do not match the segment layout")

    def pointer(array):
        return array.ctypes.data_as(ctypes.c_void_p)

    pair = handle.fastq_pair_open(str(read1).encode(), str(read2).encode(), ctypes.c_size_t(batch), pointer(stacked), ctypes.c_size_t(segment_length), pointer(starts), pointer(ends), pointer(lookup), ctypes.c_int64(size), ctypes.c_size_t(umi[0]), ctypes.c_size_t(umi[1]))
    try:
        length = ctypes.c_size_t(0)
        while True:
            n = handle.fastq_pair_peek(pair, ctypes.byref(length))
            if n < 0:
                raise ValueError(handle.fastq_pair_error(pair).decode())
            if n == 0:
                return
            out = (np.empty((n, length.value), dtype=np.uint8), np.empty(n, dtype=np.int64), np.empty(n, dtype=np.int64), np.empty(n, dtype=np.uint8))
            handle.fastq_pair_take(pair, *(pointer(a) for a in out))
            yield out[0], out[1], out[2], out[3].astype(bool)
    finally:
        handle.fastq_pair_close(pair)


# ---------------------------------------------------------------- barcodes and UMIs


def segment_table(whitelist, max_distance=2):
    """Lookup over all 8-mers (4^L): index of the unique nearest whitelist entry within max_distance, else -1."""
    length = len(whitelist[0])
    codes = np.array([[int(c) for c in encode(word)] for word in whitelist], dtype=np.int64)
    every = np.arange(4 ** length, dtype=np.int64)
    digits = (every[:, None] >> (2 * np.arange(length - 1, -1, -1))[None, :]) & 3
    distance = (digits[:, None, :] != codes[None, :, :]).sum(axis=2)
    nearest = distance.argmin(axis=1)
    best = distance[every, nearest]
    distance[every, nearest] = length + 1
    second = distance.min(axis=1)
    return np.where((best <= max_distance) & (second > best), nearest, -1).astype(np.int64)


def pack(codes):
    """Integer code of each row of 2-bit codes (N mapped to A) and a flag for rows with N."""
    codes = np.asarray(codes)
    has_n = (codes > 3).any(axis=1)
    values = np.minimum(codes, 3).astype(np.int64)
    weights = 4 ** np.arange(codes.shape[1] - 1, -1, -1, dtype=np.int64)
    return values @ weights, has_n


@njit(cache=True)
def collapse_umis(groups, umis, umi_length):
    """1MM_All UMI collapsing: within each group (sorted by group, then UMI), UMIs one
    mismatch apart are merged transitively; returns the component representative
    (smallest UMI) of every row."""
    n = len(umis)
    out = umis.copy()
    start = 0
    while start < n:
        end = start
        while end < n and groups[end] == groups[start]:
            end += 1
        size = end - start
        if size > 1:
            parent = np.arange(size)
            values = umis[start:end]
            for i in range(size):
                for position in range(umi_length):
                    shift = 2 * (umi_length - 1 - position)
                    base = (values[i] >> shift) & 3
                    for alternative in range(4):
                        if alternative == base:
                            continue
                        neighbour = values[i] + ((alternative - base) << shift)
                        low, high = 0, size
                        while low < high:
                            middle = (low + high) // 2
                            if values[middle] < neighbour:
                                low = middle + 1
                            else:
                                high = middle
                        if low < size and values[low] == neighbour:
                            a, b = i, low
                            while parent[a] != a:
                                a = parent[a]
                            while parent[b] != b:
                                b = parent[b]
                            if a != b:
                                if a < b:
                                    parent[b] = a
                                else:
                                    parent[a] = b
            for i in range(size):
                a = i
                while parent[a] != a:
                    a = parent[a]
                out[start + i] = values[a]
        start = end
    return out


# ---------------------------------------------------------------- block classes


def block_tables(blocks, reference):
    """Flattened block structures in reference coordinates for block_masks.

    blocks is the collect_local_path_reads.py block dict (path key -> block).
    Returns (keys, tables) where tables holds per block its locus, gene strand
    (+1 or -1), window, path exon intervals and path junction keys (flattened
    with offsets) and variable junction keys.
    """
    locus_index = {}
    for index, genes in enumerate(reference["locus_genes"]):
        for gene in str(genes).split(";"):
            locus_index[gene] = index
    keys, locus, window_start, window_end, strand = [], [], [], [], []
    path_offsets, exon_offsets, exon_start, exon_end = [0], [0], [], []
    junction_offsets, junction_keys, variable_offsets, variable_keys = [0], [], [0], []
    for key, block in blocks.items():
        if block["gene_id"] not in locus_index:
            continue
        index = locus_index[block["gene_id"]]
        shift = int(reference["locus_ref_start"][index] - reference["locus_genome_start"][index])
        chains = [sorted((int(a) + shift, int(b) + shift) for a, b in path) for path in block["paths"]]
        mature = chains[:len(block["signatures"])]
        keys.append(key)
        locus.append(index)
        strand.append(1 if block.get("strand", "+") == "+" else -1)
        window_start.append(min(chain[0][0] for chain in chains))
        window_end.append(max(chain[-1][1] for chain in chains))
        path_junctions = []
        for chain in chains:
            for a, b in chain:
                exon_start.append(a)
                exon_end.append(b)
            exon_offsets.append(len(exon_start))
            js = sorted(((x[1]) << 32) | y[0] for x, y in zip(chain, chain[1:]) if x[1] < y[0])
            path_junctions.append(set(js))
            junction_keys.extend(js)
            junction_offsets.append(len(junction_keys))
        path_offsets.append(len(exon_offsets) - 1)
        mature_sets = path_junctions[:len(mature)]
        variable = sorted(set().union(*mature_sets) - set.intersection(*mature_sets)) if mature_sets else []
        variable_keys.extend(variable)
        variable_offsets.append(len(variable_keys))
    order = np.lexsort((np.asarray(window_start), np.asarray(locus)))
    as_array = lambda values, dtype=np.int64: np.asarray(values, dtype=dtype)
    tables = {"locus": as_array(locus), "strand": as_array(strand), "window_start": as_array(window_start), "window_end": as_array(window_end), "path_offsets": as_array(path_offsets), "exon_offsets": as_array(exon_offsets), "exon_start": as_array(exon_start), "exon_end": as_array(exon_end), "junction_offsets": as_array(junction_offsets), "junction_keys": as_array(junction_keys), "variable_offsets": as_array(variable_offsets), "variable_keys": as_array(variable_keys), "order": as_array(order)}
    return keys, tables


@njit(cache=True)
def _read_mask(starts, ends, n, b, window_start, window_end, path_offsets, exon_offsets, exon_start, exon_end, junction_offsets, junction_keys, variable_offsets, variable_keys, min_intron):
    """(mask, informative) of one read in block b, as compatibility_mask: clipped
    aligned blocks within one exon of a path and junctions inside the window
    among the path's junctions. mask is -1 when no aligned block overlaps the
    window (the read does not enter the block)."""
    w0, w1 = window_start[b], window_end[b]
    overlap = False
    for i in range(n):
        if starts[i] < w1 and ends[i] > w0:
            overlap = True
    if not overlap:
        return -1, False
    informative = False
    for i in range(n - 1):
        if starts[i + 1] - ends[i] >= min_intron:
            key = (ends[i] << 32) | starts[i + 1]
            for v in range(variable_offsets[b], variable_offsets[b + 1]):
                if variable_keys[v] == key:
                    informative = True
    mask = 0
    first_path = path_offsets[b]
    for s in range(path_offsets[b + 1] - first_path):
        p = first_path + s
        ok = True
        for i in range(n):
            a, e = max(starts[i], w0), min(ends[i], w1)
            if a >= e:
                continue
            inside = False
            for x in range(exon_offsets[p], exon_offsets[p + 1]):
                if exon_start[x] <= a and e <= exon_end[x]:
                    inside = True
                    break
            if not inside:
                ok = False
                break
        if ok:
            for i in range(n - 1):
                if starts[i + 1] - ends[i] < min_intron:
                    continue
                d, acc = ends[i], starts[i + 1]
                if d >= w0 and acc <= w1:
                    key = (d << 32) | acc
                    found = False
                    for j in range(junction_offsets[p], junction_offsets[p + 1]):
                        if junction_keys[j] == key:
                            found = True
                            break
                    if not found:
                        ok = False
                        break
        if ok:
            mask |= 1 << s
    return mask, informative


def flatten_blocks(n_blocks, blocks):
    """(offsets, starts, ends) flat arrays from per-read block counts and an N x MAX_BLOCKS x 2 array."""
    n_blocks = np.asarray(n_blocks, dtype=np.int64)
    offsets = np.r_[0, np.cumsum(n_blocks)].astype(np.int64)
    keep = np.arange(blocks.shape[1])[None, :] < n_blocks[:, None]
    return offsets, blocks[:, :, 0][keep].astype(np.int64), blocks[:, :, 1][keep].astype(np.int64)


@njit(cache=True)
def block_masks(read_locus, offsets, flat_starts, flat_ends, locus, window_start, window_end, order, path_offsets, exon_offsets, exon_start, exon_end, junction_offsets, junction_keys, variable_offsets, variable_keys, min_intron=MIN_INTRON):
    """Per (read, block) records: read index, block index, mask, informative.

    Read r's aligned blocks are flat_starts/flat_ends[offsets[r]:offsets[r + 1]]
    (reference coordinates). Blocks are visited in (locus, window start)
    order; a read enters every block of its locus whose window one of its
    aligned blocks overlaps.
    """
    n_reads = len(read_locus)
    sorted_locus = locus[order]
    sorted_start = window_start[order]
    out_read, out_block, out_mask = List.empty_list(types.int64), List.empty_list(types.int64), List.empty_list(types.int64)
    out_informative = List.empty_list(types.boolean)
    for r in range(n_reads):
        first = offsets[r]
        n = offsets[r + 1] - first
        if n == 0:
            continue
        starts = flat_starts[first:first + n]
        ends = flat_ends[first:first + n]
        low, high = 0, len(order)
        while low < high:
            middle = (low + high) // 2
            if sorted_locus[middle] < read_locus[r]:
                low = middle + 1
            else:
                high = middle
        j = low
        while j < len(order) and sorted_locus[j] == read_locus[r]:
            b = order[j]
            if sorted_start[j] >= ends[n - 1]:
                break
            mask, informative = _read_mask(starts, ends, n, b, window_start, window_end, path_offsets, exon_offsets, exon_start, exon_end, junction_offsets, junction_keys, variable_offsets, variable_keys, min_intron)
            if mask >= 0:
                out_read.append(r)
                out_block.append(b)
                out_mask.append(mask)
                out_informative.append(informative)
            j += 1
    n_out = len(out_read)
    read_out = np.empty(n_out, dtype=np.int64)
    block_out = np.empty(n_out, dtype=np.int64)
    mask_out = np.empty(n_out, dtype=np.int64)
    informative_out = np.empty(n_out, dtype=np.bool_)
    for i in range(n_out):
        read_out[i], block_out[i], mask_out[i], informative_out[i] = out_read[i], out_block[i], out_mask[i], out_informative[i]
    return read_out, block_out, mask_out, informative_out


def simulated_opportunities(blocks, reference, align, read_length, batch=4_000_000):
    """Aligner-consistent read opportunities {path key: (all-read, junction-only)}.

    Each value maps mask -> length-S vector (S = number of block paths,
    precursor included), as path_read_opportunities. Every start position of
    a read overlapping path s by at least one base (path sequence with genomic
    flanks, as in the geometric opportunities) gives an error-free read in the
    gene's sense orientation; it is aligned with align (align_reads or
    rust_aligner), kept if UNIQUE, and classified into its own block with
    block_masks, so its class follows the same alignment, clipping and
    uniqueness rules as the data. The junction-only opportunities count the
    reads spanning a variable junction. Precursor reads that overlap no mature
    exon can only be precursor-only and are counted without alignment.
    """
    keys, tables = block_tables(blocks, reference)
    ref = reference["ref"]
    index_of = {key: i for i, key in enumerate(keys)}
    n_paths = {key: len(blocks[key]["paths"]) for key in keys}
    totals = {key: (defaultdict(lambda n=n_paths[key]: np.zeros(n)), defaultdict(lambda n=n_paths[key]: np.zeros(n))) for key in keys}
    flank = read_length - 1
    pending, sources = [], []
    held = 0
    tables_args = (tables["locus"], tables["window_start"], tables["window_end"], tables["order"], tables["path_offsets"], tables["exon_offsets"], tables["exon_start"], tables["exon_end"], tables["junction_offsets"], tables["junction_keys"], tables["variable_offsets"], tables["variable_keys"])

    def flush():
        reads = np.concatenate(pending)
        source_block = np.concatenate([np.full(n, b) for b, _, n in sources])
        source_path = np.concatenate([np.full(n, p) for _, p, n in sources])
        status, _, _, _, _, n_blocks, aligned = align(reads)
        keep = status == UNIQUE
        offsets, starts, ends = flatten_blocks(n_blocks[keep], aligned[keep])
        locus = np.searchsorted(reference["locus_ref_start"], aligned[keep, 0, 0], side="right") - 1
        read, block, mask, informative = block_masks(locus, offsets, starts, ends, *tables_args)
        kept = np.flatnonzero(keep)[read]
        own = (block == source_block[kept]) & (mask > 0)
        frame = pd.DataFrame({"block": block[own], "path": source_path[kept[own]], "mask": mask[own], "informative": informative[own]})
        for (b, m, p, informative_value), count in frame.groupby(["block", "mask", "path", "informative"]).size().items():
            key = keys[b]
            totals[key][0][m][p] += count
            if informative_value:
                totals[key][1][m][p] += count
        pending.clear()
        sources.clear()

    for key in keys:
        b = index_of[key]
        block = blocks[key]
        shift = int(reference["locus_ref_start"][tables["locus"][b]] - reference["locus_genome_start"][tables["locus"][b]])
        chains = [sorted((int(a) + shift, int(e) + shift) for a, e in path) for path in block["paths"]]
        n_mature = len(block["signatures"])
        mature_exons = sorted({exon for chain in chains[:n_mature] for exon in chain})
        for p, chain in enumerate(chains):
            sequence = np.concatenate([ref[chain[0][0] - flank:chain[0][0]]] + [ref[a:e] for a, e in chain] + [ref[chain[-1][1]:chain[-1][1] + flank]])
            windows = np.lib.stride_tricks.sliding_window_view(sequence, read_length)
            if p >= n_mature and len(chain) == 1:
                genomic_start = chain[0][0] - flank + np.arange(len(windows))
                touches = np.zeros(len(windows), dtype=bool)
                for a, e in mature_exons:
                    touches |= (genomic_start < e) & (genomic_start + read_length > a)
                totals[key][0][1 << p][p] += int((~touches).sum())
                windows = windows[touches]
            if not len(windows):
                continue
            reads = windows if block.get("strand", "+") == "+" else np.where(windows[:, ::-1] <= 3, 3 - windows[:, ::-1], 4).astype(np.uint8)
            pending.append(np.ascontiguousarray(reads))
            sources.append((b, p, len(reads)))
            held += len(reads)
            if held >= batch:
                flush()
                held = 0
    if pending:
        flush()
    return {key: tuple({mask: vector for mask, vector in table.items()} for table in totals[key]) for key in keys}


def add_junctions(reference, chromosomes, donors, acceptors):
    """Sorted junction keys of the reference plus genomic junctions [donor, acceptor) that fall inside one target locus."""
    keys = [reference["junctions"]]
    genome_start = reference["locus_genome_start"]
    genome_end = genome_start + reference["locus_ref_end"] - reference["locus_ref_start"]
    chromosomes, donors, acceptors = np.asarray(chromosomes).astype(str), np.asarray(donors, dtype=np.int64), np.asarray(acceptors, dtype=np.int64)
    for locus, chromosome in enumerate(reference["locus_chromosome"].astype(str)):
        inside = (chromosomes == chromosome) & (donors >= genome_start[locus]) & (acceptors <= genome_end[locus])
        shift = reference["locus_ref_start"][locus] - genome_start[locus]
        keys.append(((donors[inside] + shift) << 32) | (acceptors[inside] + shift))
    return np.unique(np.concatenate(keys)).astype(np.int64)


def locus_to_genome(reference, locus, position):
    """(chromosome, genomic position) of a reference position in a locus."""
    return str(reference["locus_chromosome"][locus]), int(position - reference["locus_ref_start"][locus] + reference["locus_genome_start"][locus])
