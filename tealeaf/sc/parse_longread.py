"""Cell barcode, UMI and RT-primer calling for Parse (Split-seq) cDNA sequenced as full-length long reads.

A full-length Parse WT v3 molecule, written in read-2 orientation (the reverse complement of the cDNA, as Illumina read 2
reads it), is

    adapter | UMI (10) | bc3 (8) | LINKER_3 (12) | bc2 (8) | LINKER_2 (12) | bc1 (8) | poly(T) or hexamer | cDNA | TSO

with the template-switch oligo (TSO) at the far end. Long reads (e.g. PacBio Kinnex S-reads) come in either strand and
carry indels, so the two linkers are placed by semi-global edit distance in both strands, and each barcode is corrected
to its whitelist by edit distance with the linker-side end fixed. bc1 encodes the RT well and primer (poly(dT) or random
hexamer). The insert returned is the cDNA in sense orientation, with the TSO and the poly(A) tail removed.
"""
import numpy as np
from numba import njit

from tealeaf.sc.native_align import CODE, TSO_CORE, TSO_TAIL, encode

LINKER_3 = "ATGAGGGGTCAG"  # read-2 orientation, between bc3 and bc2
LINKER_2 = "TCCAACCACCTC"  # between bc2 and bc1
UMI_LENGTH = 10
BARCODE_LENGTH = 8
STATUS = ("ok", "no linker", "no second linker", "bc3", "bc2", "bc1", "short")
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@njit(cache=True)
def _place(text, lo, hi, pattern):
    """Best placement of pattern (all bases) inside text[lo:hi] by edit distance: (start, end, distance).

    Ties keep the leftmost end; among equal paths a substitution is preferred over an indel.
    """
    m = len(pattern)
    prev, cur = np.arange(m + 1), np.zeros(m + 1, dtype=np.int64)
    prev_start, cur_start = np.full(m + 1, lo), np.zeros(m + 1, dtype=np.int64)
    best, best_start, best_end = m + 1, -1, -1
    for j in range(lo, hi):
        cur[0], cur_start[0] = 0, j + 1
        for i in range(1, m + 1):
            cost, start = prev[i - 1] + (text[j] != pattern[i - 1]), prev_start[i - 1]
            if prev[i] + 1 < cost:
                cost, start = prev[i] + 1, prev_start[i]
            if cur[i - 1] + 1 < cost:
                cost, start = cur[i - 1] + 1, cur_start[i - 1]
            cur[i], cur_start[i] = cost, start
        if cur[m] < best:
            best, best_start, best_end = cur[m], cur_start[m], j + 1
        prev, cur, prev_start, cur_start = cur, prev, cur_start, prev_start
    return best_start, best_end, best


@njit(cache=True)
def _edit(text, lo, hi, word, row):
    """Global edit distance between text[lo:hi] and word (row is scratch of length len(word) + 1)."""
    m = len(word)
    for i in range(m + 1):
        row[i] = i
    for j in range(lo, hi):
        diagonal, row[0] = row[0], j - lo + 1
        for i in range(1, m + 1):
            value = min(diagonal + (text[j] != word[i - 1]), row[i] + 1, row[i - 1] + 1)
            diagonal, row[i] = row[i], value
    return row[m]


@njit(cache=True)
def _correct(text, anchor, side, lengths, words, max_edits, row):
    """Whitelist index of the barcode next to anchor, and its distance and length; index -1 if none or ambiguous.

    side 1 reads text[anchor:anchor + L], side -1 reads text[anchor - L:anchor], side 0 reads text[anchor:lengths[0]]
    (a barcode between two placed linkers). Each word takes its best length; ties between words are ambiguous.
    """
    best, second, index, length = max_edits + 1, max_edits + 1, -1, 0
    for w in range(words.shape[0]):
        word_best, word_length = max_edits + 1, 0
        for k in range(1 if side == 0 else len(lengths)):
            if side == 0:
                lo, hi = anchor, lengths[0]
            elif side > 0:
                lo, hi = anchor, anchor + lengths[k]
            else:
                lo, hi = anchor - lengths[k], anchor
            if lo < 0 or hi > len(text) or hi < lo:
                continue
            d = _edit(text, lo, hi, words[w], row)
            if d < word_best or (d == word_best and hi - lo == words.shape[1]):
                word_best, word_length = d, hi - lo
        if word_best < best:
            second, best, index, length = best, word_best, w, word_length
        elif word_best == best:
            second = word_best
    if best > max_edits or second == best:
        return -1, best, 0
    return index, best, length


@njit(cache=True)
def _poly_t(text, start):
    """Length of the poly(T) run starting at start."""
    end = start
    while end < len(text) and text[end] == 3:
        end += 1
    return end - start


@njit(cache=True)
def _call(codes, offsets, linker3, linker2, tso, bc3, bc2, bc1, window, max_linker_edits, max_edits, out):
    lengths = np.array([8, 7, 9])
    row = np.zeros(bc1.shape[1] + 1, dtype=np.int64)
    for r in range(len(offsets) - 1):
        forward = codes[offsets[r]:offsets[r + 1]]
        n = len(forward)
        reverse = np.empty(n, dtype=np.uint8)
        for j in range(n):
            reverse[j] = COMPLEMENT[forward[n - 1 - j]]
        out[r, :] = -1
        hi = min(window, n)
        s_f, e_f, d_f = _place(forward, 0, hi, linker3)
        s_r, e_r, d_r = _place(reverse, 0, hi, linker3)
        flip = d_r <= d_f  # forward is the cDNA strand (TSO first) for most reads
        text = reverse if flip else forward
        s3, e3, d3 = (s_r, e_r, d_r) if flip else (s_f, e_f, d_f)
        out[r, 0], out[r, 1] = 1, flip
        if d3 > max_linker_edits:
            continue
        s2, e2, d2 = _place(text, e3 + 4, min(e3 + BARCODE_LENGTH + len(linker2) + 4, n), linker2)
        if d2 > max_linker_edits:
            out[r, 0] = 2
            continue
        i3, k3, l3 = _correct(text, s3, -1, lengths, bc3, max_edits, row)
        if i3 < 0:
            out[r, 0] = 3
            continue
        i2, k2, l2 = _correct(text, e3, 0, np.array([s2]), bc2, max_edits, row)
        if i2 < 0:
            out[r, 0] = 4
            continue
        i1, k1, l1 = _correct(text, e2, 1, lengths, bc1, max_edits, row)
        if i1 < 0:
            out[r, 0] = 5
            continue
        umi = s3 - l3 - UMI_LENGTH
        start = e2 + l1
        poly = _poly_t(text, start)
        if poly == 0 and _poly_t(text, start + 1) >= 5:  # one base left over by an ambiguous bc1 indel
            poly = 1 + _poly_t(text, start + 1)
        tail = max(start + poly, n - 60)
        t_start, t_end, t_d = _place(text, tail, n, tso)
        end = t_start - TSO_TAIL if t_d <= 3 else n
        if umi < 0 or end - (start + poly) < 1:
            out[r, 0] = 6
            continue
        out[r, 0] = 0
        out[r, 2], out[r, 3], out[r, 4] = i3, i2, i1
        out[r, 5], out[r, 6], out[r, 7] = k3 + d3, k2 + d2, k1
        out[r, 8], out[r, 9], out[r, 10], out[r, 11] = umi, start + poly, end, poly
        out[r, 12] = t_d <= 3


COLUMNS = ("status", "flip", "bc3", "bc2", "bc1", "edits3", "edits2", "edits1", "umi", "insert_start", "insert_end", "poly_t", "tso")


def word_codes(words):
    return np.stack([encode(word) for word in words]).astype(np.uint8)


def call_reads(sequences, bc3, bc2, bc1, window=160, max_linker_edits=2, max_edits=1):
    """Call barcodes for read sequences (strings, either strand).

    bc3, bc2, bc1 are whitelists (lists of 8-nt words). Returns an int64 array R x len(COLUMNS): status (index into
    STATUS), flip (1 when read-2 orientation is the reverse complement of the read), whitelist indices, edits (barcode
    plus its linker), and UMI start, insert start/end, poly(T) length and TSO flag, all in read-2-orientation
    coordinates. window bounds the linker search from the read-2 start.
    """
    joined = np.frombuffer("".join(sequences).encode(), dtype=np.uint8)
    codes = CODE[joined]
    offsets = np.zeros(len(sequences) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum([len(s) for s in sequences])
    out = np.zeros((len(sequences), len(COLUMNS)), dtype=np.int64)
    _call(codes, offsets, encode(LINKER_3), encode(LINKER_2), encode(reverse_complement(TSO_CORE)), word_codes(bc3), word_codes(bc2), word_codes(bc1), window, max_linker_edits, max_edits, out)
    return out


def extract(sequence, quality, row):
    """(UMI, insert, insert quality) of a called read, the insert in cDNA (sense) orientation."""
    status, flip, umi, start, end = row[0], row[1], row[8], row[9], row[10]
    text = reverse_complement(sequence) if flip else sequence
    qual = quality[::-1] if flip else quality
    return text[umi:umi + UMI_LENGTH], reverse_complement(text[start:end]), qual[start:end][::-1]
