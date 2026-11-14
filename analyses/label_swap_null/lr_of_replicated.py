"""Long-read direction agreement of the genes each method replicates in GSE320507 (perm_fdr_split.py gene lists).

For each method, replication procedure and gene set (BH, permutation FDR, and the genes added by the permutation FDR),
a gene is evaluable when the method has a long-read-evaluable full-data test or feature in it; its direction is that of
its smallest-p evaluable test. Reports evaluable genes and the fraction agreeing with the long reads, and the agreement
of all evaluable tests in the genes. rMATS and SUPPA2 have long-read mappings only for their 500 top-ranked events per
contrast. Usage: lr_of_replicated.py <replicated_genes.tsv> <output.tsv>"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table

O = Path("analyses/pbmc/GSE320507_comparators")
gene = lambda s: s.astype(str).str.split(".").str[0]
tealeaf = pd.read_csv("analyses/pbmc/GSE320507_lr_sh/lr_mapping.tsv.gz", sep="\t", low_memory=False)
tealeaf = tealeaf.loc[tealeaf.eligible.astype(str).str.lower().eq("true")].assign(method="Tealeaf")
junction = evaluable(read_source_table(O / "lr_junction_leaf5/junction_replication_merged.tsv.gz"))
rmats = evaluable(read_source_table(O / "lr_event/rmats_full_data_lr.tsv.gz")).assign(method="rMATS (paired)")
lr = pd.concat([t[["method", "gene_id", "p_value", "pooled_replicated"]] for t in (tealeaf, junction, rmats)], ignore_index=True)
lr["gene_id"] = gene(lr.gene_id)
lr["agree"] = lr.pooled_replicated.astype(str).str.lower().eq("true")
best = lr.sort_values("p_value").drop_duplicates(["method", "gene_id"]).set_index(["method", "gene_id"]).agree

genes = pd.read_csv(sys.argv[1], sep="\t", dtype={"gene_id": str})
rows = []
for (comparison, method, universe), g in genes.groupby(["comparison", "method", "universe"]):
    sets = {p: set(x.gene_id) for p, x in g.groupby("procedure")}
    sets["added by permutation FDR"] = sets.get("permutation FDR", set()) - sets.get("BH", set())
    for procedure, chosen in sets.items():
        m = method if method in set(lr.method) else None
        if m is None:
            continue
        hits = [best.get((m, x)) for x in chosen]
        hits = [h for h in hits if h is not None]
        tests = lr.loc[lr.method.eq(m) & lr.gene_id.isin(chosen)]
        rows.append(dict(comparison=comparison, method=method, universe=universe, procedure=procedure, genes=len(chosen), evaluable_genes=len(hits), gene_agreement=sum(hits) / len(hits) if hits else float("nan"), evaluable_tests=len(tests), test_agreement=tests.agree.mean() if len(tests) else float("nan")))
out = pd.DataFrame(rows)
out.to_csv(sys.argv[2], sep="\t", index=False)
print(out.loc[out.universe.eq("comparator")].round(3).to_string(index=False))
