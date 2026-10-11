#!/usr/bin/env python3
"""Block-local pseudoalignment of Parse reads with kallisto | bustools.

Stages: --stage reference writes the block-target FASTA and kallisto index
for the blocks of an existing collect_local_path_reads.py --annotated output;
--stage library pseudoaligns both runs of one physical library together,
corrects barcodes against that library's production-QC barcodes, and writes
molecule class counts; --stage collate sums libraries and stores k-mer read
opportunities for every block. Output matches collect_local_path_reads.py.
Without a D-list, reads from repeat copies elsewhere in the genome that share
one k-mer with a precursor window are assigned to it; --d-list with the genome
removes them.
"""

import argparse
import gzip
import json
from pathlib import Path
import subprocess

import pandas as pd

from tealeaf.sc.local_pseudoalign import kmer_read_opportunities, molecule_classes, path_sequences, read_equivalence_classes, write_block_reference

TECHNOLOGY = "1,10,18,1,48,56,1,78,86:1,0,10:0,0,0"
READ_LENGTH = 151


def run(command):
    print(" ".join(map(str, command)), flush=True)
    subprocess.run([str(part) for part in command], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("reference", "library", "collate"), required=True)
    parser.add_argument("--blocks", type=Path, required=True, help="blocks.json.gz from collect_local_path_reads.py --annotated")
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--fastq-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--library-index", type=int, default=0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--d-list", type=Path, help="FASTA of sequences to mask (the genome), so reads with genomic k-mers flanking a target are discarded")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    blocks = json.load(gzip.open(args.blocks, "rt"))
    index = args.output_dir / "blocks.idx"
    if args.stage == "reference":
        write_block_reference(blocks, args.genome, args.output_dir / "blocks.fa")
        masking = ["-d", args.d_list] if args.d_list else []
        run(["kallisto", "index", "-k", 31, "-t", args.threads, *masking, "-T", args.output_dir / "tmp", "-i", index, args.output_dir / "blocks.fa"])
        return
    recipe = json.loads(args.recipe.read_text())
    if args.stage == "library":
        library = recipe["bams"][args.library_index]
        out = args.output_dir / f"library_{library['library']}"
        out.mkdir(exist_ok=True)
        fastqs = []
        for entry in library["runs"]:
            accession = Path(entry["path"]).parent.name
            fastqs += [args.fastq_dir / f"{accession}_1.fastq.gz", args.fastq_dir / f"{accession}_2.fastq.gz"]
        run(["kallisto", "bus", "-i", index, "-o", out, "-x", TECHNOLOGY, "-t", args.threads, *fastqs])
        whitelist = out / "whitelist.txt"
        whitelist.write_text("\n".join(sorted(library["barcode_groups"])) + "\n")
        run(["bustools", "correct", "-w", whitelist, "-o", out / "corrected.bus", out / "output.bus"])
        run(["bustools", "sort", "-t", args.threads, "-o", out / "sorted.bus", out / "corrected.bus"])
        run(["bustools", "text", "-o", out / "sorted.txt", out / "sorted.bus"])
        groups = {barcode: tuple(group) for barcode, group in library["barcode_groups"].items()}
        counts = molecule_classes(out / "sorted.txt", read_equivalence_classes(out), groups)
        rows = [{"path_key": key[0], "subject": key[1], "cell_type": key[2], "primer": key[3], "mask": key[4], "count": value} for key, value in counts.items()]
        pd.DataFrame(rows, columns=["path_key", "subject", "cell_type", "primer", "mask", "count"]).to_csv(out / "counts.tsv.gz", sep="\t", index=False)
        for name in ("output.bus", "corrected.bus", "sorted.txt"):
            (out / name).unlink()
        print(f"library {library['library']}: {sum(counts.values())} molecules", flush=True)
        return
    parts = sorted(args.output_dir.glob("library_*/counts.tsv.gz"))
    if len(parts) != len(recipe["bams"]):
        raise ValueError(f"expected {len(recipe['bams'])} libraries, found {len(parts)}")
    table = pd.concat([pd.read_csv(path, sep="\t") for path in parts], ignore_index=True)
    table = table.groupby(["path_key", "subject", "cell_type", "primer", "mask"], as_index=False)["count"].sum()
    table.to_csv(args.output_dir / "counts.tsv.gz", sep="\t", index=False)
    with gzip.open(args.output_dir / "blocks.json.gz", "wt") as handle:
        json.dump(blocks, handle)
    import pysam
    fasta = pysam.FastaFile(str(args.genome))
    opportunities = {key: {str(mask): vector.tolist() for mask, vector in kmer_read_opportunities(path_sequences(fasta, block["chromosome"], block["paths"]), READ_LENGTH).items()} for key, block in blocks.items()}
    with gzip.open(args.output_dir / "opportunities.json.gz", "wt") as handle:
        json.dump(opportunities, handle)
    print(f"{len(blocks)} blocks, {len(table)} count rows, {int(table['count'].sum())} molecules", flush=True)


if __name__ == "__main__":
    main()
