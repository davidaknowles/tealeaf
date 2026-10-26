#!/usr/bin/env python3
"""Retained and spliced long-read molecules per clean intron, by cell type.

Reads come from ScISOr-Seq2 AllInfo tables (read, gene, cell label, barcode,
UMI, intron chain, TSS, polyA, exon chain; exons as chr_start_end_strand,
1-based inclusive, on mm10). Molecules are deduplicated by sample, barcode,
UMI and gene. Exons are lifted to the annotation's build with UCSC liftOver
(--chain); a read with an exon that fails to lift, or that lifts to another
chromosome or strand, is dropped. Clean introns (tealeaf.sc.intron_retention)
come from --gtf. A molecule retains an intron when one of its exons contains
the intron with --overhang bases of exon on each side, and splices it when
two consecutive exons are joined by exactly that intron, so retention is
called from the reads and needs no retained-intron transcript. Cell labels
are mapped to Tealeaf cell types as in assess_tilgner_long_read_replication;
P56 cells without a counterpart are kept as "other".

Writes counts.tsv.gz (intron, chromosome, start, end, gene_id, sample,
replicate, cell_type, retained, spliced) and introns.tsv.gz (the catalogue
with pooled retained and spliced molecules over P56 cells).
"""

import argparse
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from assess_tilgner_long_read_replication import biological_replicate, tilgner_cell_type
from tealeaf.sc.intron_retention import contained_introns, intron_catalogue, read_gtf_transcripts

COLUMNS = ["read", "gene", "label", "barcode", "umi", "introns", "tss", "polya", "exons"]


def lift_exons(exons, chain, work):
    """mm39 0-based half-open (chromosome, start, end) for unique (chromosome, start, end, strand) exons; NaN where lifting fails."""
    bed = work / "exons.bed"
    exons.assign(name=np.arange(len(exons)), score=0)[["chromosome", "start", "end", "name", "score", "strand"]].to_csv(bed, sep="\t", header=False, index=False)
    subprocess.run(["liftOver", str(bed), str(chain), str(work / "lifted.bed"), str(work / "unmapped.bed")], check=True, capture_output=True)
    lifted = pd.read_csv(work / "lifted.bed", sep="\t", header=None, names=["chromosome", "start", "end", "name", "score", "strand"])
    lifted = lifted.drop_duplicates("name", keep=False).set_index("name")
    out = exons[["chromosome", "strand"]].copy()
    out[["new_chromosome", "new_start", "new_end", "new_strand"]] = lifted.reindex(np.arange(len(exons)))[["chromosome", "start", "end", "strand"]].to_numpy()
    good = out.new_chromosome.eq(out.chromosome) & out.new_strand.eq(out.strand)
    return out.new_start.where(good), out.new_end.where(good)


def read_allinfo(path, chunksize):
    with tarfile.open(path) as archive:
        member = next(m for m in archive.getmembers() if m.name.endswith(".tsv.gz"))
        handle = archive.extractfile(member)
        yield from pd.read_csv(handle, sep="\t", header=None, names=COLUMNS, usecols=["gene", "label", "barcode", "umi", "exons"], compression="gzip", chunksize=chunksize)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allinfo", type=Path, nargs="+", required=True, help="<sample>_AllInfo.tsv.tar files")
    parser.add_argument("--chain", type=Path, required=True, help="liftOver chain from the reads' build to the annotation's")
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--overhang", type=int, default=20)
    parser.add_argument("--chunksize", type=int, default=2_000_000)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalogue = intron_catalogue(read_gtf_transcripts(args.gtf))
    catalogue["intron"] = np.arange(len(catalogue))
    clean = catalogue[catalogue.clean].drop_duplicates(["chromosome", "start", "end"])
    by_chromosome = {c: (g.start.to_numpy(), g.end.to_numpy(), g.intron.to_numpy()) for c, g in clean.groupby("chromosome")}
    tables = []
    for path in args.allinfo:
        seen = set()
        for chunk in read_allinfo(path, args.chunksize):
            parts = chunk.label.str.split("::", expand=True)
            chunk = chunk.loc[parts[0].eq("P56")].copy()
            parts = parts.loc[chunk.index]
            chunk["cell_type"] = [tilgner_cell_type(*row) or "other" for row in parts[[0, 2, 4, 5]].itertuples(index=False)]
            chunk["sample"] = parts[1] + "_" + parts[2]
            chunk["replicate"] = [biological_replicate(region, sample) for region, sample in zip(parts[2], parts[1])]
            key = chunk["sample"] + "|" + chunk.barcode + "|" + chunk.umi + "|" + chunk.gene
            first = ~key.duplicated() & ~key.isin(seen)
            seen.update(key[first])
            chunk = chunk.loc[first].reset_index(drop=True)
            exons = chunk.exons.str.strip(";%").str.split(";%;").explode().rename("exon").reset_index().rename(columns={"index": "molecule"})
            fields = exons.exon.str.rsplit("_", n=3, expand=True)
            exons = exons.assign(chromosome=fields[0], start=fields[1].astype(np.int64) - 1, end=fields[2].astype(np.int64), strand=fields[3])
            unique = exons[["chromosome", "start", "end", "strand"]].drop_duplicates().reset_index(drop=True)
            with tempfile.TemporaryDirectory(dir=args.output_dir) as work:
                unique["new_start"], unique["new_end"] = lift_exons(unique, args.chain, Path(work))
            exons = exons.merge(unique, on=["chromosome", "start", "end", "strand"])
            failed = exons.loc[exons.new_start.isna(), "molecule"].unique()
            exons = exons.loc[~exons.molecule.isin(failed)].astype({"new_start": np.int64, "new_end": np.int64}).sort_values(["molecule", "new_start"])
            calls = []
            for chromosome, local in exons.groupby("chromosome"):
                if chromosome not in by_chromosome:
                    continue
                starts, ends, ids = by_chromosome[chromosome]
                block, intron = contained_introns(local.new_start.to_numpy(), local.new_end.to_numpy(), starts, ends, args.overhang)
                calls.append(pd.DataFrame({"molecule": local.molecule.to_numpy()[block], "intron": ids[intron], "kind": "retained"}))
                same = local.molecule.to_numpy()[1:] == local.molecule.to_numpy()[:-1]
                junctions = pd.DataFrame({"molecule": local.molecule.to_numpy()[1:][same], "start": local.new_end.to_numpy()[:-1][same], "end": local.new_start.to_numpy()[1:][same]})
                junctions = junctions.merge(pd.DataFrame({"start": starts, "end": ends, "intron": ids}), on=["start", "end"])
                calls.append(junctions[["molecule", "intron"]].assign(kind="spliced"))
            calls = pd.concat(calls, ignore_index=True).drop_duplicates()
            calls = calls.join(chunk[["sample", "replicate", "cell_type"]], on="molecule")
            tables.append(calls.groupby(["intron", "sample", "replicate", "cell_type", "kind"], dropna=False).size().rename("count").reset_index())
            print(f"{path.name}: {len(chunk)} molecules, {len(failed)} dropped by liftOver, {len(calls)} intron calls", flush=True)
    counts = pd.concat(tables).groupby(["intron", "sample", "replicate", "cell_type", "kind"], dropna=False)["count"].sum().unstack("kind", fill_value=0).reset_index()
    for kind in ("retained", "spliced"):
        if kind not in counts:
            counts[kind] = 0
    counts = counts.merge(catalogue[["intron", "chromosome", "start", "end", "gene_id"]], on="intron")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(args.output_dir / "counts.tsv.gz", sep="\t", index=False)
    pooled = counts.groupby("intron")[["retained", "spliced"]].sum()
    catalogue.join(pooled, on="intron").fillna({"retained": 0, "spliced": 0}).to_csv(args.output_dir / "introns.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
