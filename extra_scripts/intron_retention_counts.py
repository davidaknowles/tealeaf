#!/usr/bin/env python3
"""Intron-body, boundary and spliced molecules per annotated intron from native alignments of one library.

Introns are the junctions of the --gtf transcripts (use the comprehensive annotation, so that clean introns
overlap no exon of any transcript) of the genes of the native reference. For each clean intron, a single-block
read inside the interior (MARGIN bases from each end) is an intron-body read, binned by the distance from its
sense 3' end to the next downstream sense A-run (internal poly(dT) priming); an unspliced read across the intron
start or end with OVERHANG bases on each side is a boundary read; a read with a junction equal to the intron
splices it. Molecules are (cell, UMI) pairs, deduplicated within an alignment chunk.

Writes introns.tsv.gz (catalogue with clean and retained labels, intron = row index), opportunities.tsv.gz
(read 3'-end positions per clean intron and distance bin for reads of --read-length) and
counts_<library>.tsv.gz (kind, intron, subject, cell_type, primer, sense, bin, count).
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from tealeaf.sc.intron_retention import a_run_positions, downstream_run_distance, intron_catalogue, read_gtf_transcripts

MARGIN, OVERHANG, WINDOW, MIN_A = 5, 8, 10, 8
BINS = np.array([0, 100, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 5000, np.inf])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", type=Path, required=True, help="native reference.npz")
    parser.add_argument("--gtf", type=Path, required=True, help="comprehensive annotation GTF")
    parser.add_argument("--recipe", type=Path, required=True, help="native recipe.json (libraries, runs, barcode groups)")
    parser.add_argument("--runs-dir", type=Path, required=True, help="directory of per-run alignment chunks (e.g. runs_pass2)")
    parser.add_argument("--library-index", type=int, required=True)
    parser.add_argument("--read-length", type=int, default=150, help="read length for intron-body opportunities")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    reference = np.load(args.reference)
    ref = reference["ref"]
    locus_ref_start, locus_ref_end = reference["locus_ref_start"], reference["locus_ref_end"]
    locus_genome_start, locus_chromosome = reference["locus_genome_start"], reference["locus_chromosome"]
    targets = {gene for genes in reference["locus_genes"] for gene in genes.split(";")}

    table = intron_catalogue(read_gtf_transcripts(args.gtf), targets)
    table["intron"] = np.arange(len(table))
    clean = table[table.clean].drop_duplicates(["chromosome", "start", "end"])
    index = {c: (g.start.to_numpy(), g.end.to_numpy(), g.intron.to_numpy()) for c, g in clean.groupby("chromosome")}
    plus_of = table.strand.eq("+").to_numpy()

    # sense A-run centres per chromosome: A-rich windows serve plus genes, T-rich windows minus genes
    runs = {"+": defaultdict(list), "-": defaultdict(list)}
    for locus in range(len(locus_ref_start)):
        codes = ref[locus_ref_start[locus]:locus_ref_end[locus]]
        shift = locus_genome_start[locus] + WINDOW // 2
        runs["+"][locus_chromosome[locus]].append(a_run_positions(codes, 0, WINDOW, MIN_A) + shift)
        runs["-"][locus_chromosome[locus]].append(a_run_positions(codes, 3, WINDOW, MIN_A) + shift)
    runs = {strand: {c: np.unique(np.concatenate(v)) for c, v in by.items()} for strand, by in runs.items()}
    del ref
    print("A-run windows", {s: sum(len(v) for v in by.values()) for s, by in runs.items()}, flush=True)

    def distance_bin(chromosome, positions, plus):
        distance = np.empty(len(positions))
        for strand, sel in (("+", plus), ("-", ~plus)):
            distance[sel] = downstream_run_distance(positions[sel], runs[strand].get(chromosome, np.array([], dtype=np.int64)), strand == "+")
        return np.searchsorted(BINS, distance, side="right") - 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.library_index == 0:
        rows = []
        for chromosome, local in clean.groupby("chromosome"):
            for row in local.itertuples():
                plus = row.strand == "+"
                lo, hi = (row.start + MARGIN + args.read_length - 1, row.end - MARGIN) if plus else (row.start + MARGIN, row.end - MARGIN - args.read_length + 1)
                if hi <= lo:
                    continue
                positions = np.arange(lo, hi)
                bins = distance_bin(chromosome, positions, np.full(len(positions), plus))
                for b, n in zip(*np.unique(bins, return_counts=True)):
                    rows.append((row.intron, b, n))
        table.to_csv(args.output_dir / "introns.tsv.gz", sep="\t", index=False)
        pd.DataFrame(rows, columns=["intron", "bin", "positions"]).to_csv(args.output_dir / "opportunities.tsv.gz", sep="\t", index=False)

    library = json.load(open(args.recipe))["bams"][args.library_index]
    groups = [tuple(library["barcode_groups"][b]) for b in sorted(library["barcode_groups"])]
    subject = np.array([g[0] for g in groups], dtype=object)
    cell_type = np.array([g[1] for g in groups], dtype=object)
    primer = np.array([g[2] for g in groups], dtype=object)
    counts = defaultdict(int)

    def add(kind, introns, reads, bins, cell, umi, orient):
        key = pd.DataFrame({"intron": introns, "subject": subject[cell[reads]], "ct": cell_type[cell[reads]], "primer": primer[cell[reads]], "sense": (orient[reads] == 0) == plus_of[introns], "bin": bins, "mol": (cell[reads].astype(np.int64) << 32) | umi[reads]}).drop_duplicates(["intron", "mol"])
        for k, v in key.groupby(["intron", "subject", "ct", "primer", "sense", "bin"]).size().items():
            counts[(kind,) + k] += v

    shift_of = locus_genome_start - locus_ref_start
    for run in library["runs"]:
        for path in sorted((args.runs_dir / Path(run["path"]).parent.name).glob("chunk_*.npz")):
            z = np.load(path)
            cell, umi, locus, orient, offs, starts, ends = (z[k] for k in ("cell", "umi", "locus", "orientation", "offsets", "starts", "ends"))
            nblocks = np.diff(offs)
            shift = shift_of[locus]
            chrom = locus_chromosome[locus]
            single = np.flatnonzero(nblocks == 1)
            s, e = starts[offs[single]] + shift[single], ends[offs[single]] + shift[single]
            for c in np.unique(chrom[single]):
                if c not in index:
                    continue
                m = chrom[single] == c
                r, rs, re_ = single[m], s[m], e[m]
                cs, ce, ci = index[c]
                pos = np.searchsorted(cs, rs - MARGIN, side="right") - 1
                ok = pos >= 0
                ok[ok] &= re_[ok] <= ce[pos[ok]] - MARGIN
                ii = ci[pos[ok]]
                three_prime = np.where(plus_of[ii], re_[ok] - 1, rs[ok])
                add("body", ii, r[ok], distance_bin(c, three_prime, plus_of[ii]), cell, umi, orient)
                for kind, edge in (("left", cs), ("right", ce)):
                    pos = np.searchsorted(edge, rs + OVERHANG, side="left")
                    ok = pos < len(edge)
                    ok[ok] &= edge[pos[ok]] <= re_[ok] - OVERHANG
                    ii = ci[pos[ok]]
                    add(kind, ii, r[ok], np.zeros(ok.sum(), dtype=int), cell, umi, orient)
            multi = np.flatnonzero(nblocks > 1)
            k = np.concatenate([np.arange(offs[i], offs[i + 1] - 1) for i in multi]) if len(multi) else np.array([], dtype=np.int64)
            owner = np.repeat(multi, nblocks[multi] - 1)
            junctions = pd.DataFrame({"chromosome": chrom[owner], "start": ends[k] + shift[owner], "end": starts[k + 1] + shift[owner], "read": owner})
            junctions = junctions.merge(table[["chromosome", "start", "end", "intron"]], on=["chromosome", "start", "end"])
            add("spliced", junctions.intron.to_numpy(), junctions.read.to_numpy(), np.zeros(len(junctions), dtype=int), cell, umi, orient)
            print(path, len(counts), flush=True)

    pd.DataFrame([k + (v,) for k, v in counts.items()], columns=["kind", "intron", "subject", "cell_type", "primer", "sense", "bin", "count"]).to_csv(args.output_dir / f"counts_{args.library_index}.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
