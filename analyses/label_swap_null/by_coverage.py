"""Swap-null rejection at p < .01 by feature coverage (molecules or reads in the two cell types of the pair):
LeafCutter clusters (junction reads summed over the pair's pseudobulks) and Tealeaf blocks (median molecules per row)."""
import gzip, json
from pathlib import Path
import numpy as np, pandas as pd
O = Path(".cache/label_swap"); PB = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
contrasts = json.loads((O / "contrasts_swapped.json").read_text())
counts = pd.read_csv(PB / "leafcutter/clustering/benchmark_perind_numers.counts.gz", sep=" ", index_col=0)
counts["cluster"] = [f"{name.split(':')[0]}:{name.split(':')[3]}" for name in counts.index]
by_cluster = counts.groupby("cluster").sum()
rows = []
for i, c in enumerate(contrasts):
    t = pd.read_csv(O / f"leafcutter/normalized/{i}.tsv", sep="\t")
    t = t[t.p_value < 1]
    samples = [s for s in c["samples_a"] + c["samples_b"] if s in by_cluster.columns]
    t["coverage"] = by_cluster.loc[t.feature_id, samples].sum(axis=1).to_numpy()
    rows.append(t[["feature_id", "p_value", "coverage"]].assign(pair=f"{c['level_a']}||{c['level_b']}"))
lc = pd.concat(rows)
bins = [0, 100, 300, 1000, 3000, 10000, np.inf]
lc["bin"] = pd.cut(lc.coverage, bins)
print("LeafCutter (tested clusters), junction reads in the pair")
print(lc.groupby("bin", observed=True).p_value.agg(tests="size", p01=lambda p: (p < .01).mean(), p001=lambda p: (p < .001).mean()).round(4).to_string())
tl = []
for s in ("internal", "terminal"):
    t = pd.read_csv(O / f"tealeaf_{s}/merged/paired_path.tsv", sep="\t", usecols=["p_value", "converged", "median_gene_umis"])
    t["p"] = np.where(t.converged.astype(str).str.lower().eq("true"), t.p_value, 1.)
    tl.append(t)
tl = pd.concat(tl)
tl["bin"] = pd.cut(tl.median_gene_umis, [0, 5, 10, 20, 50, 100, np.inf])
print("Tealeaf (all screened tests), median molecules per row")
print(tl.groupby("bin", observed=True).p.agg(tests="size", p01=lambda p: (p < .01).mean(), p001=lambda p: (p < .001).mean()).round(4).to_string())
