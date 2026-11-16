"""Long-read effect size and direction agreement of each method's replicated genes in GSE320507 (permutation FDR,
comparator universes, coverage-filtered): per gene the smallest-p long-read-evaluable full-data test; effect size is the
largest absolute change in long-read share over the test's forms (paths, junctions or event forms) between the two cell
types. Agreement within long-read effect-size bands shows whether a method's overall agreement reflects the size of the
differences it replicates. Usage: agreement_by_effect.py <output.tsv>"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table

O = Path("analyses/pbmc/GSE320507_comparators")
gene = lambda s: s.astype(str).str.split(".").str[0]


def vec(v):
    return np.asarray(json.loads(v), float) if isinstance(v, str) and v.startswith("[") else np.atleast_1d(np.asarray(v, float))


def share_change(row):
    """Largest absolute change in long-read share between the levels, from pooled replicate counts."""
    a = vec(row.counts_a_rep1) + vec(row.counts_a_rep2); b = vec(row.counts_b_rep1) + vec(row.counts_b_rep2)
    return float(np.max(np.abs(b / b.sum() - a / a.sum()))) if a.sum() and b.sum() else np.nan


t = pd.read_csv("analyses/pbmc/GSE320507_lr_sh/lr_mapping.tsv.gz", sep="\t", low_memory=False)
t = t.loc[t.eligible.astype(str).str.lower().eq("true")].copy()
t["size"] = [share_change(r) for r in t.itertuples(index=False)]
t = t.assign(method="Tealeaf")
j = evaluable(read_source_table(O / "lr_junction_leaf5/junction_replication_merged.tsv.gz"))
j["size"] = j.long_read_effect_norm.astype(float)  # junction-share change norm
r = evaluable(read_source_table(O / "lr_event/rmats_full_data_lr.tsv.gz")).assign(method="rMATS (paired)")
r["size"] = r.long_read_effect.astype(float).abs()
lr = pd.concat([x[["method", "gene_id", "p_value", "pooled_replicated", "size"]] for x in (t, j, r)], ignore_index=True)
lr["gene_id"] = gene(lr.gene_id); lr["agree"] = lr.pooled_replicated.astype(str).str.lower().eq("true")
best = lr.sort_values("p_value").drop_duplicates(["method", "gene_id"])
genes = pd.read_csv("analyses/label_swap_null/perm171/replicated_genes_filtered.tsv", sep="\t", dtype={"gene_id": str})
genes = genes.loc[genes.universe.eq("comparator") & genes.procedure.eq("permutation FDR")]
rows = []
for (comparison, method), g in genes.groupby(["comparison", "method"]):
    b = best.loc[best.method.eq(method) & best.gene_id.isin(set(g.gene_id))].copy()
    b["band"] = pd.cut(b["size"], [-1e-9, .02, .05, .1, .2, 1.01])
    row = dict(comparison=comparison, method=method, replicated=g.gene_id.nunique(), evaluable=len(b), median_lr_size=b["size"].median(), agreement=b.agree.mean())
    for band, x in b.groupby("band", observed=False):
        row[f"n {band}"] = len(x); row[f"agree {band}"] = x.agree.mean() if len(x) else np.nan
    rows.append(row)
out = pd.DataFrame(rows); out.to_csv(sys.argv[1], sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.round(2).to_string(index=False))
