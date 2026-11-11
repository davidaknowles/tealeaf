#!/usr/bin/env python3
"""Primer-specific read-end profiles from single-transcript genes (label-blind).

For a random set of protein-coding genes with one annotated transcript of at
least --min-length spliced bases, the 3'-most aligned base of every primary,
uniquely mapped, gene-assigned read in a QC cell (one per barcode, UMI and
position) is located on the transcript. Per library and primer, read ends are
binned by spliced distance from the 3' and 5' ends; --collect pools libraries
and writes tealeaf.sc.position_bias.end_profiles per primer to profiles.json.
--by-cell-type keys histograms and profiles by "primer|cell type" instead
(a diagnostic: the opportunities use the label-blind primer profiles).
With --genome, each gene histogram also records the reads whose end lies near an
A-run of the spliced sense sequence (an A-run start, tealeaf.sc.position_bias.
a_run_starts, at offset A-run start - read end in [--a-run-low, --a-run-high)),
as "three_arun", "five_arun" and "middle_arun" (a subset of the totals).
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path
import re

import numpy as np

from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode
from tealeaf.sc.position_bias import a_run_starts, end_profiles, transcript_end_distance

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def single_transcript_genes(gtf, min_length):
    exons, transcripts = defaultdict(list), defaultdict(set)
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "exon" or 'gene_type "protein_coding"' not in fields[8]:
                continue
            gene = re.search('gene_id "([^"]+)"', fields[8]).group(1)
            transcript = re.search('transcript_id "([^"]+)"', fields[8]).group(1)
            transcripts[gene].add(transcript)
            exons[gene].append((fields[0], fields[6], int(fields[3]) - 1, int(fields[4])))
    genes = {}
    for gene, names in transcripts.items():
        if len(names) == 1 and sum(b - a for _, _, a, b in exons[gene]) >= min_length:
            genes[gene] = (exons[gene][0][0], exons[gene][0][1], sorted((a, b) for _, _, a, b in exons[gene]))
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--library-index", type=int, default=0)
    parser.add_argument("--genes", type=int, default=3000)
    parser.add_argument("--min-length", type=int, default=4500)
    parser.add_argument("--bin-width", type=int, default=25)
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--by-cell-type", action="store_true")
    parser.add_argument("--genome", type=Path, help="genome FASTA; records reads near A-runs")
    parser.add_argument("--a-run-low", type=int, default=-200)
    parser.add_argument("--a-run-high", type=int, default=200)
    args = parser.parse_args()
    middle = args.bin_width * args.bins
    if args.collect:
        parts = sorted(args.output_dir.glob("library_*.json"))
        rows = defaultdict(list)
        for path in parts:
            for primer, genes in json.loads(path.read_text()).items():
                rows[primer].extend((value["length"], np.array(value["three"]), np.array(value["five"]), value["middle"]) for value in genes)
        profiles = {}
        for primer, values in rows.items():
            f, g = end_profiles(values, args.bin_width, args.bins, middle)
            profiles[primer] = {"three_prime": f.tolist(), "five_prime": g.tolist()}
        result = {"bin_width": args.bin_width, "libraries": len(parts), "genes": {primer: len(values) for primer, values in rows.items()}, "profiles": profiles}
        (args.output_dir / "profiles.json").write_text(json.dumps(result, indent=1) + "\n")
        for primer, value in sorted(profiles.items()):
            print(primer, result["genes"][primer], "3':", np.round(value["three_prime"][:12], 2).tolist(), "5':", np.round(value["five_prime"][:12], 2).tolist(), flush=True)
        return
    import pysam

    genes = single_transcript_genes(args.gtf, args.min_length)
    names = sorted(genes)
    chosen = [names[index] for index in np.random.default_rng(20261012).permutation(len(names))[:args.genes]]
    library = json.loads(args.recipe.read_text())["bams"][args.library_index]
    groups = library["barcode_groups"]
    sources = [pysam.AlignmentFile(run["path"], "rb", index_filename=run["index"]) for run in library["runs"]]
    fasta = pysam.FastaFile(str(args.genome)) if args.genome else None
    label = (lambda group: f"{group[2]}|{group[1]}") if args.by_cell_type else (lambda group: group[2])
    result = defaultdict(list)
    for gene in chosen:
        chromosome, strand, exons = genes[gene]
        length = sum(b - a for a, b in exons)
        seen = defaultdict(set)
        if fasta is not None:
            sequence = "".join(fasta.fetch(chromosome, a, b).upper() for a, b in exons)
            runs = a_run_starts(sequence if strand == "+" else sequence.translate(COMPLEMENT)[::-1])
        for source in sources:
            if chromosome not in source.references:
                continue
            for alignment in source.fetch(chromosome, exons[0][0], exons[-1][1]):
                if alignment.is_secondary or alignment.is_supplementary or not all(alignment.has_tag(tag) for tag in ("GX", "NH", "CB", "UB")) or alignment.get_tag("NH") != 1:
                    continue
                if alignment.get_tag("GX").split(".")[0] != gene.split(".")[0]:
                    continue
                barcode = normalize_starsolo_barcode(alignment.get_tag("CB"))
                if barcode not in groups:
                    continue
                end = alignment.reference_end - 1 if strand == "+" else alignment.reference_start
                distance = transcript_end_distance(exons, strand, end)
                if distance is not None:
                    seen[label(groups[barcode])].add((barcode, alignment.get_tag("UB"), distance))
        for primer, reads in seen.items():
            hist = {key: [np.zeros(args.bins), np.zeros(args.bins), 0] for key in ("", "_arun")}
            for _, _, (d3, d5) in reads:
                keys = [""]
                if fasta is not None and len(runs):
                    first = np.searchsorted(runs, d5 + args.a_run_low)
                    if first < len(runs) and runs[first] < d5 + args.a_run_high:
                        keys.append("_arun")
                for key in keys:
                    if d3 < middle:
                        hist[key][0][d3 // args.bin_width] += 1
                    if d5 < middle:
                        hist[key][1][d5 // args.bin_width] += 1
                    if d3 >= middle and d5 >= middle:
                        hist[key][2] += 1
            if reads:
                row = {"gene": gene, "length": length}
                for key in hist if fasta is not None else ("",):
                    row.update({f"three{key}": hist[key][0].tolist(), f"five{key}": hist[key][1].tolist(), f"middle{key}": hist[key][2]})
                result[primer].append(row)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"library_{args.library_index}.json").write_text(json.dumps(result) + "\n")
    print(f"library {library['library']}: {len(chosen)} genes, {sum(len(value) for value in result.values())} gene-primer histograms", flush=True)


if __name__ == "__main__":
    main()
