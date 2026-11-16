import sys, json
sys.path.insert(0, "."); sys.path.insert(0, "extra_scripts")
from pathlib import Path
import numpy as np, pysam
from estimate_read_end_profiles import single_transcript_genes
from tealeaf.sc.junction_benchmark import normalize_starsolo_barcode
from tealeaf.sc.position_bias import transcript_end_distance
REF = Path.home() / "tealeaf_runs/pbmc/ref"; B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark/majiq_bams"
recipe = json.loads(Path(".cache/pbmc/GSE320507/recipe.json").read_text())
primer_of = {b: g[2] for lib in recipe["bams"] for b, g in lib["barcode_groups"].items()}
genes = single_transcript_genes(REF / "gencode.v49.basic.annotation.gtf", 4500)
names = sorted(genes); chosen = [names[i] for i in np.random.default_rng(1).permutation(len(names))[:600]]
stats = {}
for bam in ["CO-109__Tcm_Naive_helper_T_cells", "CO-109__Classical_monocytes", "CO-112__Tcm_Naive_helper_T_cells"]:
    src = pysam.AlignmentFile(str(B / f"{bam}.bam"))
    for g in chosen:
        ch, strand, exons = genes[g]
        if ch not in src.references: continue
        for a in src.fetch(ch, exons[0][0], exons[-1][1]):
            if a.is_secondary or a.is_supplementary or not a.has_tag("GX") or a.get_tag("GX").split(".")[0] != g.split(".")[0]: continue
            p = primer_of.get(normalize_starsolo_barcode(a.get_tag("CB")))
            if p is None: continue
            sense = (a.is_reverse == (strand == "-"))
            for label, pos in (("3prime_of_read_in_tx", a.reference_end - 1 if strand == "+" else a.reference_start), ("5prime_of_read_in_tx", a.reference_start if strand == "+" else a.reference_end - 1)):
                d = transcript_end_distance(exons, strand, pos)
                if d is None: continue
                k = (p, label); s = stats.setdefault(k, [0, 0, 0, 0])
                s[0] += 1; s[1] += d[0] < 200; s[2] += d[0] < 2000; s[3] += d[1] < 2000
            stats.setdefault((p, "sense_reads"), [0, 0])[0] += sense; stats[(p, "sense_reads")][1] += 1
for k, v in sorted(stats.items()):
    if k[1] == "sense_reads": print(k, "sense fraction %.2f" % (v[0] / v[1]))
    else: print(k, "reads", v[0], "last200 %.3f last2kb %.3f first2kb %.3f" % (v[1] / v[0], v[2] / v[0], v[3] / v[0]))
