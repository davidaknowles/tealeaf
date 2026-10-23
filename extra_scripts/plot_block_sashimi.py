#!/usr/bin/env python3
"""Sashimi plots of example blocks: Parse coverage by primer and long-read coverage, under the gene's isoforms.

For each example of export_block_examples.py, the plotted region is the
block window. Parse tracks are UMI-deduplicated native alignments
(native_collect.py align chunks) of the two cell types, by primer, per 1,000
half-cells. Long-read tracks are the coverage and junctions implied by the
full-length annotated-transcript UMI counts of the long-read data, summed
over replicates, per --long-read-scale long-read UMIs of the cell type (default one million). The annotation
track shows the gene's annotated transcripts and the block's tested paths.
"""

import argparse
import gzip
import json
import re
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_tilgner_long_read_replication import read_tilgner_matrix, stable_identifier
from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.sashimi import SashimiEvent, collect_native_event_support, long_read_event_support, merge_support, write_ggsashimi_inputs

PRIMERS = ("poly(dT)", "random hexamer")
LONG = "long reads"
COLORS = ("#2C7FB8", "#D95F0E")


def library_support(library, runs_dir, locus_chromosome, locus_shift, events):
    barcodes = sorted(library["barcode_groups"])
    groups = {index: tuple(library["barcode_groups"][barcode]) for index, barcode in enumerate(barcodes)}
    chunks = [path for run in library["runs"] for path in sorted((runs_dir / Path(run["path"]).parent.name).glob("chunk_*.npz"))]
    return collect_native_event_support(chunks, groups, locus_chromosome, locus_shift, events)


def gene_gtf_lines(gtf, genes):
    """GTF transcript and exon lines of the given genes, with transcript_id replaced by the transcript name."""
    lines = {gene: [] for gene in genes}
    gene_pattern, name_pattern, id_pattern = re.compile(r'gene_id "([^".]+)'), re.compile(r'transcript_name "([^"]+)"'), re.compile(r'transcript_id "[^"]+"')
    with open(gtf) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon"):
                continue
            gene = gene_pattern.search(fields[8])
            if gene and gene.group(1) in lines:
                name = name_pattern.search(fields[8])
                fields[8] = id_pattern.sub(f'transcript_id "{name.group(1) if name else "transcript"}"', fields[8])
                lines[gene.group(1)].append("\t".join(fields))
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--examples", type=Path, required=True)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--native-dir", type=Path, required=True)
    parser.add_argument("--runs", default="runs_pass2")
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--long-read-matrix", type=Path, required=True)
    parser.add_argument("--ggsashimi", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--long-read-scale", type=float, default=1e6)
    parser.add_argument("--min-junction", type=int, default=1, help="ggsashimi minimum normalized junction coverage")
    parser.add_argument("--plot-setup", default="", help="shell commands run before ggsashimi (making Rscript with ggplot2, data.table and gridExtra, and a python with pandas, available)")
    args = parser.parse_args()
    examples = pd.read_csv(args.examples, sep="\t")
    blocks = {block["block_id"]: block for block in json.load(gzip.open(args.blocks, "rt")).values()}
    events, info = [], {}
    for example in examples.itertuples():
        block = blocks[example.block_id]
        mature = block["paths"][:len(block["signatures"])]
        event_id = f"{example.gene_name}_{example.block_id.split(':')[-1]}"
        events.append(SashimiEvent(event_id, block["chromosome"], min(a for path in mature for a, _ in path), max(b for path in mature for _, b in path), block["strand"]))
        info[event_id] = (example, mature)
    reference = np.load(args.native_dir / "reference.npz")
    locus_chromosome = reference["locus_chromosome"].astype(str)
    locus_shift = (reference["locus_genome_start"] - reference["locus_ref_start"]).astype(np.int64)
    recipe = json.loads(args.recipe.read_text())
    with ProcessPoolExecutor(max_workers=len(recipe["bams"])) as executor:
        results = list(executor.map(library_support, recipe["bams"], repeat(args.native_dir / args.runs), repeat(locus_chromosome), repeat(locus_shift), repeat(tuple(events))))
    exon_blocks, junctions, _ = merge_support(results)
    cells = pd.DataFrame([group for library in recipe["bams"] for group in library["barcode_groups"].values()], columns=["subject", "cell_type", "primer"])
    cell_sizes = cells.groupby(["subject", "cell_type", "primer"]).size().rename("cells").reset_index()
    matrix, features, columns = read_tilgner_matrix(args.long_read_matrix, args.gtf)
    chains = {}
    for gene, record in read_gtf_exons(args.gtf).items():
        if gene.split(".")[0] in set(examples.gene_id):
            chains.update({stable_identifier(transcript): sorted(chain) for transcript, chain in record["transcripts"].items()})
    long_sizes = []
    lines = gene_gtf_lines(args.gtf, set(examples.gene_id))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for event in events:
        example, mature = info[event.event_id]
        levels = [example.level_a, example.level_b]
        counts = {}
        for cell_type in levels:
            selected = columns.index[columns.tealeaf_cell_type.eq(cell_type)].to_numpy()
            totals = np.asarray(matrix[:, selected].sum(axis=1)).ravel()
            long_sizes.append({"subject": LONG, "cell_type": cell_type, "primer": LONG, "cells": float(totals.sum())})
            rows = features.index[features.stable_gene_id.eq(example.gene_id) & features.transcript_id.notna()]
            counts.update({(cell_type, stable_identifier(features.transcript_id[row])): float(totals[row]) for row in rows})
        long_blocks, long_junctions = long_read_event_support(event, counts, chains)
        sizes = pd.concat([cell_sizes, pd.DataFrame(long_sizes).drop_duplicates(["cell_type"])], ignore_index=True)
        event_dir = args.output_dir / event.event_id
        intron, exon, strand = write_ggsashimi_inputs(event_dir, event, exon_blocks + long_blocks, junctions + long_junctions, sizes, levels, primers=(*PRIMERS, LONG), scales={LONG: args.long_read_scale})
        gtf = event_dir / f"{event.event_id}.gtf"
        path_lines = []
        for index, chain in enumerate(mature, start=1):
            attributes = f'gene_id "{example.gene_id}"; transcript_id "Path {index}";'
            path_lines.append("\t".join((event.chromosome, "Tealeaf", "transcript", str(chain[0][0] + 1), str(chain[-1][1]), ".", event.strand, ".", attributes)))
            path_lines += ["\t".join((event.chromosome, "Tealeaf", "exon", str(start + 1), str(end), ".", event.strand, ".", attributes)) for start, end in chain]
        gtf.write_text("\n".join(path_lines) + "\n" + "".join(lines[example.gene_id]))
        palette = event_dir / "palette.txt"
        palette.write_text("".join(f"{color}\n" for color in COLORS for _ in range(3)))
        prefix = event_dir / f"{event.event_id}_sashimi"
        coordinates = f"{event.chromosome}:{event.start + 1}-{event.end}"
        command = ["python3", str(args.ggsashimi), "--intron", str(intron), "--exon", str(exon), "--strand_info", str(strand), "--data_type", "sc", "--aggregation", "False", "-c", coordinates, "-g", str(gtf), "-o", str(prefix), "-M", str(args.min_junction), "--shrink", "--height", "1.0", "--ann-height", "2.5", "--width", "9", "--base-size", "10", "-P", str(palette), "-F", "pdf"]
        subprocess.run(["bash", "-lc", f"{args.plot_setup}; {shlex.join(command)}" if args.plot_setup else shlex.join(command)], check=True)
        print(event.event_id, prefix.with_suffix(".pdf"), flush=True)


if __name__ == "__main__":
    main()
