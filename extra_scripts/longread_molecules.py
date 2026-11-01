"""Long-read molecules from aligned, barcode-tagged reads (e.g. Kinnex Parse reads from parse_longread_barcodes.py).

reads: one BAM (primary alignments, tags CB/UB/pr) -> <output>.parquet with cell, UMI, primer, chromosome, strand,
       span, exon chain, gene (largest exonic overlap on the read strand) and MAPQ.
collapse: all read tables -> <output-dir>/molecules.parquet (UMIs collapsed within cell and gene, chromosome by
       chromosome, for cells in --cells with at least --min-cell-reads reads and reads with MAPQ >= --min-mapq) and
       <output-dir>/genes/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz} (molecules per gene and cell).
"""
import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pysam
from scipy import io as sio
from scipy import sparse

from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.longread_molecules import GeneIndex, collapse_molecules, exon_chain, format_chain


def reads(args):
    index = GeneIndex(read_gtf_exons(args.gtf))
    rows = []
    with pysam.AlignmentFile(args.bam, threads=2) as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            exons = exon_chain(read.reference_start, read.cigartuples)
            strand = "-" if read.is_reverse else "+"
            rows.append((read.get_tag("CB"), read.get_tag("UB"), read.get_tag("pr"), read.reference_name, strand, exons[0][0], exons[-1][1], format_chain(exons), index.assign(read.reference_name, strand, exons) or "", read.mapping_quality))
    frame = pd.DataFrame(rows, columns=["cell", "umi", "primer", "chromosome", "strand", "start", "end", "chain", "gene", "mapq"])
    frame["span"] = frame.end - frame.start
    frame.to_parquet(args.output, index=False)
    print(f"{len(frame):,} primary alignments, {frame.gene.ne('').mean():.3f} in a gene, MAPQ>=5 {frame.mapq.ge(5).mean():.3f}")


def collapse(args):
    cells = pd.read_csv(args.cells, sep="\t")
    keep = set(cells.loc[cells.reads >= args.min_cell_reads, "cell"])
    dataset = ds.dataset(sorted(map(str, args.reads)), format="parquet")
    chromosomes = sorted(set(pq.read_table(args.reads[0], columns=["chromosome"]).column(0).to_pylist()) | {f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]})
    writer, counts = None, []
    for chromosome in chromosomes:
        frame = dataset.to_table(columns=["cell", "umi", "primer", "chromosome", "strand", "start", "end", "chain", "gene", "mapq", "span"], filter=ds.field("chromosome") == chromosome).to_pandas()
        frame = frame.loc[frame.cell.isin(keep) & frame.gene.ne("") & frame.mapq.ge(args.min_mapq)]
        if not len(frame):
            continue
        molecules = collapse_molecules(frame)
        print(chromosome, f"{len(frame):,} reads -> {len(molecules):,} molecules", flush=True)
        table = pa.Table.from_pandas(molecules, preserve_index=False)
        writer = writer or pq.ParquetWriter(args.output_dir / "molecules.parquet", table.schema)
        writer.write_table(table)
        counts.append(molecules.groupby(["gene", "cell"]).size().rename("molecules").reset_index())
    writer.close()
    counts = pd.concat(counts, ignore_index=True)
    genes, barcodes = sorted(set(counts.gene)), sorted(keep)
    gene_row, cell_column = {g: i for i, g in enumerate(genes)}, {c: j for j, c in enumerate(barcodes)}
    matrix = sparse.coo_matrix((counts.molecules.to_numpy(), (counts.gene.map(gene_row), counts.cell.map(cell_column))), shape=(len(genes), len(barcodes)))
    out = args.output_dir / "genes"
    out.mkdir(parents=True, exist_ok=True)
    sio.mmwrite(out / "matrix.mtx", matrix)
    pd.Series(genes).to_csv(out / "features.tsv.gz", index=False, header=False)
    pd.Series(barcodes).to_csv(out / "barcodes.tsv.gz", index=False, header=False)
    with open(out / "matrix.mtx", "rb") as source, gzip.open(out / "matrix.mtx.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    (out / "matrix.mtx").unlink()
    per_cell = np.asarray(matrix.tocsr().sum(axis=0)).ravel()
    print(f"{len(barcodes):,} cells, {len(genes):,} genes, median molecules per cell {np.median(per_cell):,.0f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="stage", required=True)
    one = sub.add_parser("reads")
    one.add_argument("--bam", required=True)
    one.add_argument("--gtf", required=True)
    one.add_argument("--output", required=True, type=Path)
    two = sub.add_parser("collapse")
    two.add_argument("--reads", nargs="+", required=True, type=Path)
    two.add_argument("--cells", required=True, type=Path, help="cells.tsv.gz from parse_longread_barcodes.py")
    two.add_argument("--min-cell-reads", type=int, default=1000)
    two.add_argument("--min-mapq", type=int, default=5)
    two.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.stage == "reads":
        reads(args)
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        collapse(args)


if __name__ == "__main__":
    main()
