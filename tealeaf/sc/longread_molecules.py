"""Molecules from aligned long reads that carry cell barcode and UMI tags: exon chains, gene assignment, UMI collapsing.

Exon chains are zero-based half-open genomic intervals, as in `differential.read_gtf_exons`, so that a molecule's chain can
be assigned to tested paths with `local_path_reads.transcript_path_assignments`.
"""
from collections import Counter, defaultdict

import numpy as np

from tealeaf.sc.native_align import collapse_umis, encode, pack

REFERENCE_CONSUMING = {0, 2, 7, 8}  # M, D, =, X: extend the current exon
INTRON = 3  # N


def exon_chain(reference_start, cigartuples):
    """Exons (start, end) of an alignment: N operations split exons, deletions do not."""
    exons, start, position = [], reference_start, reference_start
    for op, length in cigartuples:
        if op in REFERENCE_CONSUMING:
            position += length
        elif op == INTRON:
            if position > start:
                exons.append((start, position))
            position += length
            start = position
    if position > start:
        exons.append((start, position))
    return exons


def format_chain(exons):
    return ";".join(f"{a}-{b}" for a, b in exons)


def parse_chain(text):
    return [tuple(int(x) for x in exon.split("-")) for exon in text.split(";")]


class GeneIndex:
    """Exonic overlap of exon chains with annotated genes on the same chromosome and strand.

    genes is the output of `differential.read_gtf_exons`; a gene's exons are the union of its transcripts' exons.
    """

    def __init__(self, genes):
        intervals = defaultdict(list)
        for gene, record in genes.items():
            exons = sorted({exon for chain in record["transcripts"].values() for exon in chain})
            merged = []
            for a, b in exons:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            for a, b in merged:
                intervals[(record["chromosome"], record["strand"])].append((a, b, gene))
        self.index = {}
        for key, rows in intervals.items():
            rows.sort()
            starts = np.array([r[0] for r in rows])
            ends = np.array([r[1] for r in rows])
            self.index[key] = (starts, ends, np.maximum.accumulate(ends), [r[2] for r in rows])

    def assign(self, chromosome, strand, exons):
        """Gene with the largest exonic overlap, or None when no gene overlaps or two genes tie."""
        if (chromosome, strand) not in self.index:
            return None
        starts, ends, running_end, names = self.index[(chromosome, strand)]
        overlap = Counter()
        for a, b in exons:
            j = np.searchsorted(starts, b) - 1
            while j >= 0 and running_end[j] > a:
                shared = min(b, ends[j]) - max(a, starts[j])
                if shared > 0:
                    overlap[names[j]] += shared
                j -= 1
        if not overlap:
            return None
        (best, value), *rest = overlap.most_common(2)
        return None if rest and rest[0][1] == value else best


def collapse_molecules(frame, umi_length=10):
    """One row per molecule from a read table with cell, gene, umi, chain and span columns.

    UMIs one mismatch apart within a (cell, gene) are merged transitively (`native_align.collapse_umis`); UMIs with N are
    dropped. A molecule takes the most frequent exon chain among its reads (ties: the longest span) and counts its reads.
    """
    frame = frame.loc[~frame.umi.str.contains("N")].copy()
    group = frame.cell.astype(str) + "\t" + frame.gene.astype(str)
    codes, _ = pack(np.stack([encode(u) for u in frame.umi]) if len(frame) else np.zeros((0, umi_length), dtype=np.uint8))
    frame["group"] = group.astype("category").cat.codes.to_numpy(np.int64)
    frame["umi_code"] = codes
    pairs = frame[["group", "umi_code"]].drop_duplicates().sort_values(["group", "umi_code"])
    pairs["molecule_umi"] = collapse_umis(pairs.group.to_numpy(), pairs.umi_code.to_numpy(), umi_length)
    frame = frame.merge(pairs, on=["group", "umi_code"])
    frame["n_chain"] = frame.groupby(["group", "molecule_umi", "chain"]).chain.transform("size")
    frame = frame.sort_values(["group", "molecule_umi", "n_chain", "span"], ascending=[True, True, False, False])
    reads = frame.groupby(["group", "molecule_umi"]).size().rename("n_reads")
    molecules = frame.drop_duplicates(["group", "molecule_umi"]).merge(reads, left_on=["group", "molecule_umi"], right_index=True)
    return molecules.drop(columns=["group", "umi_code", "molecule_umi", "n_chain"]).reset_index(drop=True)
