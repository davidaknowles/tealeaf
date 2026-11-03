#!/usr/bin/env python3
"""Retained and spliced long-read molecules per clean intron and cell type, from barcode-tagged molecules.

Molecules come from longread_molecules.py collapse (one exon chain per molecule) and cells take their type from a
label table, as in assess_longread_molecules.py. Clean introns come from --gtf (comprehensive annotation). A
molecule retains an intron when one of its exons contains the intron with --overhang bases on each side, and
splices it when two consecutive exons are joined by exactly that intron (intron_retention.chain_intron_calls).

Writes counts.tsv.gz (intron, chromosome, start, end, gene_id, cell_type, replicate, retained, spliced).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from tealeaf.sc.intron_retention import chain_intron_calls, intron_catalogue, read_gtf_transcripts
from tealeaf.sc.longread_molecules import cell_labels


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--label-column", default="Immune_All_Low_majority")
    parser.add_argument("--merge", nargs="*", default=[])
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--overhang", type=int, default=20)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalogue = intron_catalogue(read_gtf_transcripts(args.gtf))
    catalogue["intron"] = np.arange(len(catalogue))
    clean = catalogue[catalogue.clean].drop_duplicates(["chromosome", "start", "end"])
    by_chromosome = {c: (g.start.to_numpy(), g.end.to_numpy(), g.intron.to_numpy()) for c, g in clean.groupby("chromosome")}
    labels = cell_labels(args.labels, args.label_column, args.merge)
    source = pq.ParquetFile(args.molecules)
    tables = []
    for group in range(source.num_row_groups):
        frame = source.read_row_group(group, columns=["cell", "chromosome", "chain"]).to_pandas()
        frame = frame.join(labels, on="cell", how="inner").reset_index(drop=True)
        exons = frame.chain.str.split(";").explode().rename("exon").reset_index().rename(columns={"index": "molecule"})
        bounds = exons.exon.str.split("-", expand=True).astype(np.int64)
        exons = exons.assign(chromosome=frame.chromosome.to_numpy()[exons.molecule], start=bounds[0], end=bounds[1]).sort_values(["molecule", "start"])
        calls = []
        for chromosome, local in exons.groupby("chromosome"):
            if chromosome in by_chromosome:
                calls.append(chain_intron_calls(local.molecule.to_numpy(), local.start.to_numpy(), local.end.to_numpy(), *by_chromosome[chromosome], overhang=args.overhang))
        if not calls:
            continue
        calls = pd.concat(calls, ignore_index=True).join(frame[["cell_type", "replicate"]], on="molecule")
        tables.append(calls.groupby(["intron", "cell_type", "replicate", "kind"]).size().rename("count").reset_index())
        print(f"row group {group}: {len(frame):,} molecules, {len(calls):,} intron calls", flush=True)
    counts = pd.concat(tables).groupby(["intron", "cell_type", "replicate", "kind"])["count"].sum().unstack("kind", fill_value=0).reset_index()
    for kind in ("retained", "spliced"):
        if kind not in counts:
            counts[kind] = 0
    counts = counts.merge(catalogue[["intron", "chromosome", "start", "end", "gene_id"]], on="intron")
    counts.to_csv(args.output_dir / "counts.tsv.gz", sep="\t", index=False)
    pooled = counts.groupby("intron")[["retained", "spliced"]].sum()
    print(f"{len(pooled):,} introns with calls; {int((pooled.retained >= 5).sum()):,} with >= 5 retaining molecules", flush=True)


if __name__ == "__main__":
    main()
