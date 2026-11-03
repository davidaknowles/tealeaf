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


def molecule_source(molecules, labels):
    """Long-read source for `long_read_assessment`: molecules per distinct exon chain and (cell type, replicate).

    molecules has gene, chain and cell columns; labels maps cell -> (cell type, replicate). Returns a features x columns
    CSR matrix of molecule counts, the feature table (stable_gene_id, transcript_id = feature key, row), the column table
    (tealeaf_cell_type, replicate, column) and the exon chain of each feature key.
    """
    from scipy import sparse

    frame = molecules.loc[molecules.cell.isin(labels.index)]
    frame = frame.assign(stable_gene_id=frame.gene.str.split(".").str[0], cell_type=frame.cell.map(labels.cell_type), replicate=frame.cell.map(labels.replicate))
    frame["transcript_id"] = frame.stable_gene_id + "|" + frame.chain
    counts = frame.groupby(["stable_gene_id", "transcript_id", "cell_type", "replicate"]).size().rename("molecules").reset_index()
    features = counts[["stable_gene_id", "transcript_id"]].drop_duplicates().reset_index(drop=True)
    features["row"] = range(len(features))
    columns = counts[["cell_type", "replicate"]].drop_duplicates().sort_values(["cell_type", "replicate"]).reset_index(drop=True).rename(columns={"cell_type": "tealeaf_cell_type"})
    columns["column"] = range(len(columns))
    row = counts.transcript_id.map(features.set_index("transcript_id").row)
    column = counts.set_index(["cell_type", "replicate"]).index.map(columns.set_index(["tealeaf_cell_type", "replicate"]).column)
    matrix = sparse.csr_matrix((counts.molecules.to_numpy(float), (row.to_numpy(), np.asarray(column))), shape=(len(features), len(columns)))
    chains = {key: parse_chain(key.split("|", 1)[1]) for key in features.transcript_id}
    return matrix, features, columns, chains


def clean_label(label):
    """Cell-type name as written into pseudobulk read groups: characters other than letters, digits and _.- become _."""
    return "".join(character if character.isalnum() or character in "_.-" else "_" for character in label)


def cell_labels(path, column, merge=(), suffix="", clean=False):
    """cell -> (cell_type, replicate) from a label table (e.g. CellTypist), replicate = 1 + RT-well column parity.

    Spaces in labels become underscores; merge holds old=new renames applied after that; with clean, labels are also
    passed through clean_label (the naming of junction-comparator pseudobulks).
    """
    import pandas as pd

    table = pd.read_csv(path, sep="\t", index_col=0)
    types = table[column].str.replace(" ", "_")
    for pair in merge:
        old, new = pair.split("=")
        types = types.replace(old, new)
    if clean:
        types = types.map(clean_label)
    types = types + suffix
    well = pd.Series(table.index.str.split("_").str[0], index=table.index)
    replicate = 1 + well.map(lambda w: int(w[1:]) % 2)
    return pd.DataFrame({"cell_type": types, "replicate": replicate})


def stranded_chains(chains, gtf_path):
    """{feature key: (chromosome, strand, exons)} for molecule_source chains (keys stable_gene_id|chain), taking each
    gene's chromosome and strand from the annotation; features of genes missing from it are dropped."""
    from tealeaf.sc.differential import read_gtf_exons

    genes = {gene.split(".")[0]: (record["chromosome"], record["strand"]) for gene, record in read_gtf_exons(gtf_path).items()}
    output = {}
    for key, exons in chains.items():
        location = genes.get(key.split("|", 1)[0])
        if location is not None:
            output[key] = (*location, tuple(exons))
    return output
