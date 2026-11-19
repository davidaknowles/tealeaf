"""GSE320507 genes replicated (dataset split, analyses/pbmc/GSE320507_<tag>) under the current blocks but not under the
consistent local rule: for each, whether the block with its smallest full-data p under the current rule was kept by the
local rule (matched by content) or replaced, and its structure and size."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from extra_scripts.assess_paired_inference_audit import fit_table

def genes(tag):
    t = pd.read_csv(f"analyses/pbmc/GSE320507_{tag}/split_by_pair.tsv", sep="\t")
    return t
m = pd.read_csv(".cache/pbmc/GSE320507/localrule_match.tsv", sep="\t", dtype=str, keep_default_na=False)
kept = set(m.loc[m.source.eq("current"), "old_block_id"])
cur = pd.concat([fit_table(f".cache/binary_grid/pbmc_GSE320507_{s}_sh_full/merged/paired_path.tsv").assign(structure=s) for s in ("internal", "terminal")])
new = pd.concat([fit_table(f".cache/binary_grid/pbmc_GSE320507_{s}_lr_full/merged/paired_path.tsv").assign(structure=s) for s in ("internal", "terminal")])
for t in (cur, new):
    t["gene"] = t.gene_id.str.split(".").str[0]
best_cur = cur.sort_values("p_value").drop_duplicates("gene").set_index("gene")
best_new = new.sort_values("p_value").drop_duplicates("gene").set_index("gene")
import json, gzip
blocks = json.load(gzip.open(".cache/pbmc/GSE320507/collection/all/blocks.json.gz", "rt"))
npaths = {b["block_id"]: len(b["signatures"]) for b in blocks.values()}
sizes = {b["block_id"]: max(e for p in b["paths"] for _, e in p) - min(a for p in b["paths"] for a, _ in p) for b in blocks.values()}
# replicated genes: rerun the split metric gene lists is not stored; use genes with BH-significant full-data discoveries in both
print("genes with a full-data discovery: current", (cur.groupby("gene").fdr.min() <= .05).sum(), "local rule", (new.groupby("gene").fdr.min() <= .05).sum())
lost = sorted(set(cur.loc[cur.fdr <= .05, "gene"]) - set(new.loc[new.fdr <= .05, "gene"]))
rows = []
for g in lost:
    b = best_cur.loc[g]
    rows.append(dict(gene=g, best_block=b.block_id, structure=b.structure, kept=b.block_id in kept, paths=npaths.get(b.block_id), span=sizes.get(b.block_id), p_current=b.p_value, p_local=best_new.p_value.get(g, np.nan)))
d = pd.DataFrame(rows)
print("genes with a discovery under the current blocks only:", len(d))
print(d.groupby(["kept", "structure"]).size().to_string())
print("replaced best blocks: median paths", d.loc[~d.kept, "paths"].median(), "median span", d.loc[~d.kept, "span"].median(), "| kept:", d.loc[d.kept, "paths"].median(), d.loc[d.kept, "span"].median())
print("median best p current %.1e, local rule %.1e" % (d.p_current.median(), d.p_local.median()))
d.to_csv("analyses/local_anchors/lost_genes_GSE320507.tsv", sep="\t", index=False)
