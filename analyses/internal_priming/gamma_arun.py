"""Diagnostic 1: primer offset of the inclusion path of two-path internal cassette tests vs A-runs in the alternative exon."""
import gzip, json, re, sys
import numpy as np, pandas as pd, pysam
sys.path.insert(0, ".")
from tealeaf.sc.local_path_reads import block_structure, structure_paths

FA = pysam.FastaFile("/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor/GRCm39.primary_assembly.genome.fa")
blocks = {b["block_id"]: b for b in json.load(gzip.open(".cache/native_v4/collection_parity_p2_v4/all/blocks.json.gz", "rt")).values()}
fits = pd.read_csv(".cache/binary_grid/internal_native_v4_full/merged/paired_path.tsv", sep="\t", usecols=["test_id", "block_id", "n_paths", "kept_paths", "primer_offsets", "converged", "p_value"])
fits = fits.loc[fits.n_paths.eq(2) & fits.converged.astype(bool)]
COMP = str.maketrans("ACGTN", "TGCAN")

def sense(chrom, start, end, strand):
    s = FA.fetch(chrom, start, end).upper()
    return s if strand == "+" else s.translate(COMP)[::-1]

def aruns(seq, k=8, w=10):
    """number of positions starting a window of w bases with at least k A (merged runs)."""
    a = np.frombuffer(seq.encode(), dtype=np.uint8) == ord("A")
    if len(a) < w:
        return 0
    c = np.convolve(a, np.ones(w, int), "valid") >= k
    return int(np.sum(c[1:] & ~c[:-1]) + c[0])

rows, cache = [], {}
for f in fits.itertuples(index=False):
    b = blocks.get(f.block_id)
    if b is None or block_structure(structure_paths(b))["structure"] != "internal":
        continue
    kept = json.loads(f.kept_paths)
    offs = json.loads(f.primer_offsets)
    if len(kept) != 2:
        continue
    chains = [set(map(tuple, b["paths"][k])) for k in kept]
    only = [chains[0] - chains[1], chains[1] - chains[0]]
    # cassette: one path has exactly one extra exon and the other none
    if not ((len(only[0]) == 1 and not only[1]) or (len(only[1]) == 1 and not only[0])):
        continue
    inc = 0 if only[0] else 1
    (start, end), = only[inc]
    key = (f.block_id, start, end)
    if key not in cache:
        seq = sense(b["chromosome"], start, end, b["strand"])
        down = sorted(chains[inc], key=lambda e: e[0])
        cache[key] = (len(seq), aruns(seq), seq.count("A") / max(len(seq), 1))
    length, n_arun, frac_a = cache[key]
    rows.append({"block_id": f.block_id, "gamma_inclusion": (offs[0] if inc == 0 else -offs[0]), "exon_length": length, "aruns": n_arun, "frac_a": frac_a, "p_value": f.p_value})
t = pd.DataFrame(rows)
t.to_csv(".cache/internal_priming/gamma_arun_tests.tsv.gz", sep="\t", index=False)
g = t.groupby("block_id").agg(gamma=("gamma_inclusion", "median"), aruns=("aruns", "first"), frac_a=("frac_a", "first"), length=("exon_length", "first"), tests=("gamma_inclusion", "size")).reset_index()
print("tests", len(t), "blocks", len(g), "fraction of test offsets at |gamma|>=3.99:", round((t.gamma_inclusion.abs() >= 3.99).mean(), 3))
print("gamma quantiles (tests):", t.gamma_inclusion.quantile([.05, .25, .5, .75, .95]).round(2).tolist())
g["arun_class"] = pd.cut(g.aruns, [-1, 0, 1, 100], labels=["0", "1", "2+"])
print(g.groupby("arun_class", observed=True).agg(blocks=("gamma", "size"), median_gamma=("gamma", "median"), mean_gamma=("gamma", "mean"), frac_neg=("gamma", lambda x: (x < 0).mean())).round(3).to_string())
from scipy.stats import spearmanr
print("spearman gamma vs A fraction:", spearmanr(g.gamma, g.frac_a), " vs A-runs per kb:", spearmanr(g.gamma, g.aruns / g.length * 1000), " vs length:", spearmanr(g.gamma, g.length))
