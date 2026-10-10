"""Block-local pseudoalignment classes and k-mer read opportunities.

Each block contributes one target per mature path (its local spliced exon
chain) and one precursor target (its genomic window). Reads are pseudoaligned
to all targets with kallisto; a read's class in a block is the bitmask of that
block's targets in its equivalence class, and a UMI molecule's class is the
AND over its reads that hit the block. Opportunities use the same rule on
error-free reads: a read's class is the AND over its canonical k-mers present
in the block of their target masks, so model weights match the pseudoaligner.
"""

from collections import defaultdict
import gzip
from pathlib import Path

import numpy as np

CODE = np.full(256, 4, dtype=np.uint8)
for _base, _value in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    CODE[_base] = _value


def path_sequences(fasta, chromosome, paths):
    """Spliced sequence of each exon chain (0-based half-open, genomic order).

    fasta is a pysam.FastaFile.
    """
    return [b"".join(fasta.fetch(chromosome, start, end).upper().encode() for start, end in sorted(path)) for path in paths]


def canonical_kmers(sequence, k=31):
    """Canonical 2-bit k-mer codes and validity (no N) for every start."""
    codes = CODE[np.frombuffer(sequence, dtype=np.uint8)]
    n = len(codes) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=bool)
    windows = np.lib.stride_tricks.sliding_window_view(codes, k)
    valid = (windows < 4).all(axis=1)
    values = np.minimum(windows, 3).astype(np.uint64)
    shifts = (2 * np.arange(k - 1, -1, -1)).astype(np.uint64)
    forward = (values << shifts).sum(axis=1, dtype=np.uint64)
    reverse = ((3 - values[:, ::-1]) << shifts).sum(axis=1, dtype=np.uint64)
    return np.minimum(forward, reverse), valid


def kmer_read_opportunities(sequences, read_length, k=31):
    """n[mask] -> per-target start counts under the k-mer intersection rule.

    sequences are the block's target sequences (mature paths, then precursor).
    Reads start at every position overlapping a target by at least k bases;
    their k-mers outside the target are absent from the block and ignored.
    """
    kmers = [canonical_kmers(sequence, k) for sequence in sequences]
    pooled = np.concatenate([codes[valid] for codes, valid in kmers])
    owners = np.concatenate([np.full(valid.sum(), index) for index, (_, valid) in enumerate(kmers)])
    unique, inverse = np.unique(pooled, return_inverse=True)
    membership = np.zeros(len(unique), dtype=np.int64)
    np.bitwise_or.at(membership, inverse, (1 << owners).astype(np.int64))
    opportunities = defaultdict(lambda: np.zeros(len(sequences)))
    bits, span, offset = len(sequences), read_length - k + 1, 0
    for index, (codes, valid) in enumerate(kmers):
        masks = np.full(len(codes), -1, dtype=np.int64)
        masks[valid] = membership[inverse[offset:offset + valid.sum()]]
        offset += valid.sum()
        if not len(codes):
            continue
        present = masks >= 0
        starts = np.arange(1 - span, len(masks))
        low, high = np.maximum(starts, 0), np.minimum(starts + span, len(masks))
        found = np.r_[0, np.cumsum(present)]
        observed = found[high] - found[low] > 0
        read_masks = np.zeros(len(starts), dtype=np.int64)
        for bit in range(bits):
            missing = np.r_[0, np.cumsum(present & ((masks >> bit) & 1 == 0))]
            read_masks |= ((missing[high] - missing[low]) == 0).astype(np.int64) << bit
        values, frequencies = np.unique(read_masks[observed], return_counts=True)
        for value, frequency in zip(values.tolist(), frequencies.tolist()):
            opportunities[value][index] += frequency
    return dict(opportunities)


def write_block_reference(blocks, fasta_path, output):
    """FASTA of all block targets, named '<path key>|<target index>'."""
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    with open(output, "w") as handle:
        for key, block in blocks.items():
            for index, sequence in enumerate(path_sequences(fasta, block["chromosome"], block["paths"])):
                handle.write(f">{key}|{index}\n{sequence.decode()}\n")


def read_equivalence_classes(output_dir):
    """kallisto bus output -> EC id -> {path key: target bitmask}."""
    output_dir = Path(output_dir)
    targets = [line.strip().rsplit("|", 1) for line in open(output_dir / "transcripts.txt")]
    classes = {}
    for line in open(output_dir / "matrix.ec"):
        ec, members = line.split()
        masks = defaultdict(int)
        for member in members.split(","):
            key, index = targets[int(member)]
            masks[key] |= 1 << int(index)
        classes[int(ec)] = dict(masks)
    return classes


def molecule_classes(bus_text, classes, barcode_groups):
    """Sorted bustools text (barcode, UMI, EC, count) -> class counts.

    Per (barcode, UMI, block) the masks of all records are ANDed. Returns
    counts[(path key, *group, mask)], where group is the barcode's
    (subject, cell type, primer); unrequested barcodes are skipped.
    """
    counts = defaultdict(int)
    current, masks = None, {}
    opener = gzip.open if str(bus_text).endswith(".gz") else open

    def flush():
        group = barcode_groups.get(current[0]) if current else None
        if group is not None:
            for key, mask in masks.items():
                if mask:
                    counts[(key, *group, mask)] += 1

    with opener(bus_text, "rt") as handle:
        for line in handle:
            barcode, umi, ec, _ = line.split("\t")
            if (barcode, umi) != current:
                flush()
                current, masks = (barcode, umi), {}
            for key, mask in classes[int(ec)].items():
                masks[key] = masks[key] & mask if key in masks else mask
        flush()
    return dict(counts)
