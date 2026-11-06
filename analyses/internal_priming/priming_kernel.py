"""Diagnostic 3: read-end enrichment around A-runs on single-transcript genes, per primer (label-blind).

For each read (unique, gene-assigned, QC cell, one per barcode/UMI/position) the 3'-most aligned base is placed on the
spliced transcript; the cross-correlation of read ends with A-run starts (>= 8 A in 10 bases, sense) at offset
delta = A-run position - read end is compared with the same quantity for reads spread uniformly over the gene.
Usage: priming_kernel.py <library index> <output dir>   |   priming_kernel.py collect <output dir>
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
sys.path.insert(0, "extra_scripts"); sys.path.insert(0, ".")
MAX = 1500

if sys.argv[1] == "collect":
    out = Path(sys.argv[2]); obs, exp = defaultdict(lambda: np.zeros(2 * MAX + 1)), defaultdict(lambda: np.zeros(2 * MAX + 1))
    for p in sorted(out.glob("library_*.json")):
        for primer, v in json.loads(p.read_text()).items():
            obs[primer] += np.asarray(v["observed"]); exp[primer] += np.asarray(v["expected"])
    delta = np.arange(-MAX, MAX + 1)
    res = {}
    for primer in obs:
        ratio = obs[primer] / np.maximum(exp[primer], 1e-12)
        res[primer] = ratio.tolist()
        bins = [(-1500, -500), (-500, -200), (-200, 0), (0, 100), (100, 200), (200, 300), (300, 500), (500, 800), (800, 1500)]
        print(primer, {f"{a}..{b}": round(float(obs[primer][(delta >= a) & (delta < b)].sum() / exp[primer][(delta >= a) & (delta < b)].sum()), 3) for a, b in bins})
    (out / "kernel.json").write_text(json.dumps({"delta": delta.tolist(), "enrichment": res}) + "\n")
    sys.exit(0)

import pysam
from estimate_read_end_profiles import single_transcript_genes
from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode
from tealeaf.sc.position_bias import transcript_end_distance
lib_index, out = int(sys.argv[1]), Path(sys.argv[2])
GTF = "/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/gencode.vM32.basic.annotation.gtf"
FA = pysam.FastaFile("/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/GRCm39.primary_assembly.genome.fa")
COMP = str.maketrans("ACGTN", "TGCAN")
genes = single_transcript_genes(GTF, 4500)
names = sorted(genes)
chosen = [names[i] for i in np.random.default_rng(20261106).permutation(len(names))[:3000]]
library = json.loads(Path(".cache/native_cells/recipe.json").read_text())["bams"][lib_index]
groups = library["barcode_groups"]
sources = [pysam.AlignmentFile(run["path"], "rb", index_filename=run["index"]) for run in library["runs"]]
obs, exp = defaultdict(lambda: np.zeros(2 * MAX + 1)), defaultdict(lambda: np.zeros(2 * MAX + 1))
for gene in chosen:
    chrom, strand, exons = genes[gene]
    seq = "".join(FA.fetch(chrom, a, b).upper() for a, b in exons)
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
    L = len(seq)
    a = np.frombuffer(seq.encode(), dtype=np.uint8) == ord("A")
    win = np.convolve(a, np.ones(10, int), "valid") >= 8
    starts = np.flatnonzero(win & ~np.r_[False, win[:-1]])  # A-run start positions (5'->3' transcript coordinates)
    arun = np.zeros(L); arun[starts] = 1
    seen = defaultdict(set)
    for src in sources:
        if chrom not in src.references:
            continue
        for al in src.fetch(chrom, exons[0][0], exons[-1][1]):
            if al.is_secondary or al.is_supplementary or not all(al.has_tag(t) for t in ("GX", "NH", "CB", "UB")) or al.get_tag("NH") != 1:
                continue
            if al.get_tag("GX").split(".")[0] != gene.split(".")[0]:
                continue
            bc = normalize_starsolo_barcode(al.get_tag("CB"))
            if bc not in groups:
                continue
            end = al.reference_end - 1 if strand == "+" else al.reference_start
            d = transcript_end_distance(exons, strand, end)
            if d is not None:
                seen[groups[bc][2]].add((bc, al.get_tag("UB"), d[1]))  # d[1]: distance from the 5' end
    if not starts.size:
        continue
    for primer, reads in seen.items():
        r = np.zeros(L)
        for _, _, x in reads:
            if 0 <= x < L:
                r[x] += 1
        n = r.sum()
        if n == 0:
            continue
        # c[delta] = sum_x r[x] * arun[x + delta]
        full = np.correlate(arun, r, mode="full")  # index k corresponds to delta = k - (L - 1)
        u = np.full(L, n / L)
        fullu = np.correlate(arun, u, mode="full")
        lo, hi = L - 1 - MAX, L - 1 + MAX + 1
        idx = np.arange(-MAX, MAX + 1) + (L - 1)
        ok = (idx >= 0) & (idx < len(full))
        obs[primer][ok] += full[idx[ok]]; exp[primer][ok] += fullu[idx[ok]]
out.mkdir(parents=True, exist_ok=True)
(out / f"library_{lib_index}.json").write_text(json.dumps({p: {"observed": obs[p].tolist(), "expected": exp[p].tolist()} for p in obs}) + "\n")
print("library", lib_index, "done", flush=True)
