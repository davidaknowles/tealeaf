"""Read-end position along transcripts by cell type and primer in GSE320507 (does 3' bias or internal priming vary
between cell types?).

For 2,000 random protein-coding genes with one annotated transcript of at least 4.5 kb spliced, the 3'-most aligned
base of every unique, gene-assigned read in the donor x cell type pseudobulk BAMs (one per barcode, UMI and position)
is placed on the transcript; the primer comes from the cell barcode (recipe barcode groups). Per BAM and primer,
per-gene histograms of read ends from the 3' and 5' ends (25-base bins over 2 kb), the reads in between, and the subset
of reads whose end lies within 200 bases of an A-run start (tealeaf.sc.position_bias.a_run_starts).
Usage: celltype_read_ends.py <task> <tasks> <output dir>  |  celltype_read_ends.py collect <output dir>
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "."); sys.path.insert(0, "extra_scripts")

B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark/majiq_bams"
REF = Path.home() / "tealeaf_runs/pbmc/ref"
BIN, BINS = 25, 80
MIDDLE = BIN * BINS
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

if sys.argv[1] == "collect":
    out = Path(sys.argv[2]); rows = []
    for path in sorted(out.glob("task_*.json")):
        for sample, primers in json.loads(path.read_text()).items():
            donor, cell_type = sample.split("__", 1)
            for primer, h in primers.items():
                three, five = np.array(h["three"]), np.array(h["five"]); total = three.sum() + five.sum() + h["middle"]
                body = total - three[:8].sum() - five[:8].sum()  # outside the first and last 200 bases
                rows.append(dict(donor=donor, cell_type=cell_type, primer=primer, reads=total, genes=h["genes"],
                                 last_200=three[:8].sum() / total, first_2kb=five.sum() / total, last_2kb=three.sum() / total,
                                 body_near_arun=(h["body_arun"] / body) if body else np.nan))
    d = pd.DataFrame(rows); d.to_csv(out / "per_sample.tsv", sep="\t", index=False)
    d = d[d.reads >= 2000]
    agg = d.groupby(["primer", "cell_type"]).agg(samples=("donor", "size"), reads=("reads", "sum"), last_200=("last_200", "median"), last_2kb=("last_2kb", "median"), first_2kb=("first_2kb", "median"), body_near_arun=("body_near_arun", "median"))
    agg.to_csv(out / "by_cell_type.tsv", sep="\t")
    print(agg.round(3).to_string())
    sys.exit(0)

import pysam
from estimate_read_end_profiles import single_transcript_genes
from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode
from tealeaf.sc.position_bias import a_run_starts, transcript_end_distance
task, tasks, out = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
recipe = json.loads(Path(".cache/pbmc/GSE320507/recipe.json").read_text())
primer_of = {barcode: group[2] for library in recipe["bams"] for barcode, group in library["barcode_groups"].items()}
genes = single_transcript_genes(REF / "gencode.v49.basic.annotation.gtf", 4500)
names = sorted(genes); chosen = [names[i] for i in np.random.default_rng(20261115).permutation(len(names))[:2000]]
fasta = pysam.FastaFile(str(REF / "GRCh38.primary_assembly.genome.fa"))
runs = {}
for g in chosen:
    chromosome, strand, exons = genes[g]
    seq = "".join(fasta.fetch(chromosome, a, b).upper() for a, b in exons)
    runs[g] = a_run_starts(seq if strand == "+" else seq.translate(COMPLEMENT)[::-1])
bams = sorted(B.glob("*.bam"))[task::tasks]
result = {}
for bam in bams:
    source = pysam.AlignmentFile(str(bam), "rb")
    hist = defaultdict(lambda: {"three": np.zeros(BINS), "five": np.zeros(BINS), "middle": 0, "body_arun": 0, "genes": 0})
    for g in chosen:
        chromosome, strand, exons = genes[g]
        if chromosome not in source.references:
            continue
        seen = set(); found = set()
        for a in source.fetch(chromosome, exons[0][0], exons[-1][1]):
            if a.is_secondary or a.is_supplementary or not a.has_tag("GX") or not a.has_tag("UB") or a.get_tag("NH") != 1 or a.get_tag("GX").split(".")[0] != g.split(".")[0]:
                continue
            primer = primer_of.get(normalize_starsolo_barcode(a.get_tag("CB")))
            if primer is None:
                continue
            end = a.reference_end - 1 if strand == "+" else a.reference_start
            key = (a.get_tag("CB"), a.get_tag("UB"), end)
            if key in seen:
                continue
            seen.add(key)
            dist = transcript_end_distance(exons, strand, end)
            if dist is None:
                continue
            d3, d5 = dist; h = hist[primer]; found.add(primer)
            if d3 < MIDDLE:
                h["three"][d3 // BIN] += 1
            if d5 < MIDDLE:
                h["five"][d5 // BIN] += 1
            if d3 >= MIDDLE and d5 >= MIDDLE:
                h["middle"] += 1
            if d3 >= 200 and d5 >= 200:
                r = runs[g]; i = np.searchsorted(r, d5 - 200)
                if i < len(r) and r[i] < d5 + 200:
                    h["body_arun"] += 1
        for primer in found:
            hist[primer]["genes"] += 1
    result[bam.stem] = {p: {"three": h["three"].tolist(), "five": h["five"].tolist(), "middle": h["middle"], "body_arun": h["body_arun"], "genes": h["genes"]} for p, h in hist.items()}
    print(bam.stem, {p: int(h["three"].sum() + h["five"].sum() + h["middle"]) for p, h in hist.items()}, flush=True)
out.mkdir(parents=True, exist_ok=True)
(out / f"task_{task}.json").write_text(json.dumps(result))
