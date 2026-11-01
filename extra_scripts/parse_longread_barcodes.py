"""Call Parse cell barcodes, UMIs and RT primers in long reads (PacBio Kinnex S-reads of a Parse WT library).

Writes <out>/reads.fastq.gz (cDNA inserts in sense orientation; the comment carries SAM tags for `minimap2 -y`: CB cell
as bc1-well_bc2-well_bc3-well, UB UMI, pr primer T/R, pa poly(T) length, ts TSO found), <out>/cells.tsv.gz (reads per
cell and primer), and <out>/summary.tsv (reads by status, edits, orientation).
Usage: parse_longread_barcodes.py --bam segmented.bam --output-dir out [--tables dir] [--kit n102_R1_v3_3,v1,R3_v3]
"""
import argparse
import csv
import gzip
import multiprocessing as mp
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from tealeaf.sc import parse_longread as pl

TABLES = Path("/gpfs/commons/home/daknowles/knowles_lab/data/sc/ParsePipeline/newvolume/ParseBiosciences-Pipeline.1.2.1/splitpipe/barcodes")
PRIMERS = ("T", "R")


def kit_tables(tables, kit):
    """bc3, bc2, bc1 tables (lists of csv rows); bc1 keeps the poly(dT) and random-hexamer entries."""
    bc1, bc2, bc3 = (list(csv.DictReader(open(Path(tables) / f"bc_data_{name}.csv"))) for name in kit.split(","))
    return bc3, bc2, [r for r in bc1 if r["stype"] in PRIMERS]


def _init(tables):
    global TABLE, WORDS, WELLS
    TABLE = tables
    WORDS = [[r["sequence"] for r in t] for t in tables]
    WELLS = sorted({r["well"] for r in tables[2]}, key=lambda w: (w[0], int(w[1:])))


def _work(batch):
    bc3, bc2, bc1 = TABLE
    calls = pl.call_reads([s for _, s, _ in batch], *WORDS)
    lines, cells, stats = [], Counter(), Counter()
    for (name, sequence, quality), row in zip(batch, calls):
        status = pl.STATUS[row[0]]
        stats[("status", status)] += 1
        if status != "ok":
            continue
        stats[("flip", int(row[1]))] += 1
        stats[("edits", int(row[5] + row[6] + row[7]))] += 1
        stats[("tso", int(row[12]))] += 1
        r1, r2, r3 = bc1[row[4]], bc2[row[3]], bc3[row[2]]
        stats[("primer", r1["stype"])] += 1
        stats[(f"poly_t_ge10_{r1['stype']}", int(row[11] >= 10))] += 1
        cell = f"{r1['well']}_{r2['well']}_{r3['well']}"
        cells[(cell, r1["stype"])] += 1
        umi, insert, qual = pl.extract(sequence, quality, row)
        lines.append(f"@{name}\tCB:Z:{cell}\tUB:Z:{umi}\tpr:Z:{r1['stype']}\tpa:i:{row[11]}\tts:i:{row[12]}\n{insert}\n+\n{qual}\n")
    return gzip.compress("".join(lines).encode(), compresslevel=1), cells, stats


def batches(bam, size, limit):
    batch, count = [], 0
    with pysam.AlignmentFile(bam, check_sq=False, threads=4) as handle:
        for read in handle.fetch(until_eof=True):
            values = read.query_qualities  # some S-reads have qualities above 93, which FASTQ cannot hold
            quality = "I" * read.query_length if values is None else (np.minimum(np.frombuffer(values, dtype=np.uint8), 93) + 33).tobytes().decode()
            batch.append((read.query_name, read.query_sequence, quality))
            count += 1
            if len(batch) == size or count == limit:
                yield batch
                batch = []
            if count == limit:
                return
    if batch:
        yield batch


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tables", default=TABLES, type=Path, help="split-pipe barcode tables (bc_data_<name>.csv)")
    parser.add_argument("--kit", default="n102_R1_v3_3,v1,R3_v3", help="bc1,bc2,bc3 table names")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--batch", type=int, default=20000)
    parser.add_argument("--limit", type=int, default=0, help="stop after this many reads (0: all)")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables = kit_tables(args.tables, args.kit)
    cells, stats = Counter(), Counter()
    with open(args.output_dir / "reads.fastq.gz", "wb") as fastq, mp.Pool(args.workers, _init, (tables,)) as pool:
        for k, (chunk, c, s) in enumerate(pool.imap(_work, batches(args.bam, args.batch, args.limit or -1), chunksize=1)):
            fastq.write(chunk)
            cells.update(c)
            stats.update(s)
            if k % 500 == 0:
                print(f"{(k + 1) * args.batch:,} reads, ok {stats[('status', 'ok')] / max(1, sum(v for (a, _), v in stats.items() if a == 'status')):.3f}", flush=True)
    frame = pd.Series(cells, name="reads").rename_axis(["cell", "primer"]).reset_index()
    frame = frame.pivot(index="cell", columns="primer", values="reads").fillna(0).astype(int).reindex(columns=list(PRIMERS), fill_value=0)
    frame["reads"] = frame.sum(axis=1)
    frame.sort_values("reads", ascending=False).to_csv(args.output_dir / "cells.tsv.gz", sep="\t")
    summary = pd.DataFrame([(a, b, v) for (a, b), v in stats.items()], columns=["field", "value", "reads"]).sort_values(["field", "value"])
    summary["fraction"] = summary.reads / summary.groupby("field").reads.transform("sum")
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
