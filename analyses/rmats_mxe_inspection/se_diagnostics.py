"""Why do rMATS's long-read-evaluable skipped-exon (SE) tests in GSE320507 disagree with the long reads?

For each evaluable full-data SE test (rmats_full_data_lr.tsv.gz), from the STARsolo junction counts of the same
pseudobulks: inclusion junctions UE (upstream exon to the exon) and ED (exon to downstream exon) and skipping junction UD.
Directions (level b minus level a, the convention of the long-read table):
- pooled junction PSI, ((UE + ED) / 2) / ((UE + ED) / 2 + UD), with counts summed over the pseudobulks of each cell type;
- mean per-pseudobulk junction PSI (unweighted, as rMATS's IncLevelDifference averages per-sample inclusion levels);
- rMATS (-IncLevelDifference).
Agreement of each with the long-read direction, overall, for genes rMATS replicates (permutation FDR), and by long-read
depth and long-read effect size; and, on the same genes, agreement of Tealeaf's best evaluable test.
Usage: se_diagnostics.py <output.tsv>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
sys.path.insert(0, ".")
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table

B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
lr = evaluable(read_source_table(Path("analyses/pbmc/GSE320507_comparators/lr_event/rmats_full_data_lr.tsv.gz")))
lr = lr.loc[lr.feature_id.str.startswith("rMATS:SE:")].copy()
lr["id"] = lr.feature_id.str.split(":").str[2].astype(int)
ev = pd.read_csv(B / "rmats/count_pass/fromGTF.SE.txt", sep="\t").rename(columns={"exonStart_0base": "es", "exonEnd": "ee"}).set_index("ID")
lr = lr.join(ev[["chr", "strand", "es", "ee", "upstreamEE", "downstreamES"]], on="id")
lr["gene"] = lr.gene_id.astype(str).str.split(".").str[0]
lr["lr_effect"] = lr.long_read_effect.astype(float)
lr["short_effect"] = lr.short_read_effect.astype(float)

junctions = pd.read_csv(B / "pseudobulk_junctions_junctions.tsv.gz", sep="\t", usecols=["chromosome", "start", "end"])
index = {(c, s, e): i for i, (c, s, e) in enumerate(zip(junctions.chromosome, junctions.start, junctions.end))}
counts = sparse.load_npz(B / "pseudobulk_junctions_counts.npz").tocsc()
samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t")
rows_of = {ct: np.flatnonzero(samples.cell_type.eq(ct).to_numpy()) for ct in samples.cell_type.unique()}


def column(key):
    i = index.get(key)
    return np.zeros(counts.shape[0]) if i is None else counts[:, i].toarray().ravel()


out = []
for e in lr.itertuples(index=False):
    ue, ed, ud = column((e.chr, e.upstreamEE + 1, e.es)), column((e.chr, e.ee + 1, e.downstreamES)), column((e.chr, e.upstreamEE + 1, e.downstreamES))
    stats = {}
    for side, level in (("a", e.level_a), ("b", e.level_b)):
        r = rows_of.get(level, np.array([], int))
        inc, skip = (ue[r] + ed[r]) / 2, ud[r]
        stats[f"pooled_{side}"] = inc.sum() / (inc.sum() + skip.sum()) if inc.sum() + skip.sum() else np.nan
        ok = inc + skip > 0
        stats[f"mean_{side}"] = np.mean(inc[ok] / (inc[ok] + skip[ok])) if ok.any() else np.nan
        stats[f"reads_{side}"] = inc.sum() + skip.sum()
    out.append(dict(gene=e.gene, contrast=e.contrast_id, p=e.p_value, short_effect=e.short_effect, lr_effect=e.lr_effect, lr_depth=e.minimum_pooled_depth,
                    pooled_effect=stats["pooled_b"] - stats["pooled_a"], mean_effect=stats["mean_b"] - stats["mean_a"], reads=min(stats["reads_a"], stats["reads_b"])))
d = pd.DataFrame(out)
nz = d.lr_effect.ne(0)
for name in ("short_effect", "pooled_effect", "mean_effect"):
    d[f"{name}_agrees"] = np.where(nz & d[name].ne(0) & d[name].notna(), np.sign(d[name]) == np.sign(d.lr_effect), np.nan)
d["rmats_matches_pooled"] = np.sign(d.short_effect) == np.sign(d.pooled_effect)
genes = pd.read_csv("analyses/label_swap_null/perm171/replicated_genes_filtered.tsv", sep="\t")
replicated = set(genes.loc[genes.method.eq("rMATS (paired)") & genes.universe.eq("comparator") & genes.procedure.eq("permutation FDR"), "gene_id"])
d["replicated"] = d.gene.isin(replicated)
d.to_csv(sys.argv[1], sep="\t", index=False)
cols = ["short_effect_agrees", "pooled_effect_agrees", "mean_effect_agrees"]
print("evaluable SE tests with nonzero long-read effect:", int(nz.sum()))
print("agreement with long reads (rMATS / pooled junction PSI / mean per-sample junction PSI):")
print(d[nz].groupby("replicated")[cols].mean().round(3).to_string())
print("rMATS direction equals pooled junction direction: %.3f" % d.rmats_matches_pooled[nz].mean())
d["lr_size"] = pd.cut(d.lr_effect.abs(), [0, .02, .05, .1, .2, 1])
d["depth"] = pd.cut(d.lr_depth, [0, 50, 100, 300, np.inf])
print(d[nz].groupby("lr_size", observed=True)[cols].agg(["size", "mean"]).round(2).to_string())
print(d[nz].groupby("depth", observed=True)[cols].mean().round(2).to_string())
# Tealeaf on the same genes
t = pd.read_csv("analyses/pbmc/GSE320507_lr_sh/lr_mapping.tsv.gz", sep="\t", low_memory=False)
t = t.loc[t.eligible.astype(str).str.lower().eq("true")]
t["gene"] = t.gene_id.astype(str).str.split(".").str[0]
t["agree"] = t.pooled_replicated.astype(str).str.lower().eq("true")
best = t.sort_values("p_value").drop_duplicates("gene").set_index("gene").agree
for label, gs in (("rMATS-replicated genes", replicated), ("all genes with an evaluable rMATS SE test", set(d.gene))):
    hits = best.reindex(sorted(gs)).dropna()
    print(f"Tealeaf best evaluable test, {label}: {len(hits)} genes, agreement {hits.mean():.3f}")
