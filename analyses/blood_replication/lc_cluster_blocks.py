"""For genes LeafCutter replicates on its matched blood universe: is the junction of LeafCutter's best cluster (per fold)
in a Tealeaf block, is that block testable, and how strong is Tealeaf's test of it on the same cell-type pair?"""
import sys, gzip, json, collections; sys.path.insert(0, ".")
import numpy as np, pandas as pd
from extra_scripts.assess_paired_inference_audit import fit_table
from tealeaf.sc.local_path_reads import path_structures, variable_junctions, block_structure, structure_paths, same_junction_groups
F = "analyses/coverage_null_filter_sh/pbmc"; LC = "/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark/leafcutter/clustering/benchmark_perind_numers.counts.gz"
g = pd.read_csv(f"{F}/split_filtered/split_gene_tests.tsv.gz", sep="\t"); g = g[g.universe.eq("matched") & g.comparison.eq("LeafCutter")]
t = g[g.method.eq("Tealeaf")].set_index("gene_id"); c = g[g.method.eq("LeafCutter")].set_index("gene_id")
rt = set(t.index[t.conjunction_q_value <= .05]); rc = set(c.index[c.conjunction_q_value <= .05])
lc = pd.concat(ch[ch.method.eq("LeafCutter") & ch.fold.astype(str).isin(["0", "1"]) & ch.gene_id.isin(rc)] for ch in pd.read_csv(f"{F}/comparator_tests_filtered.tsv.gz", sep="\t", chunksize=2_000_000, usecols=["gene_id", "pair_id", "feature_id", "p_value", "method", "fold"]))
lc["fold"] = lc.fold.astype(int)
best = lc.loc[lc.groupby(["gene_id", "fold"]).p_value.idxmin()]
clusters = set(best.feature_id)
junctions = collections.defaultdict(list)
with gzip.open(LC, "rt") as handle:
    next(handle)
    for line in handle:
        name = line.split(" ", 1)[0]; chrom, a, b, clu = name.split(":")
        key = f"{chrom}:{clu}"
        if key in clusters:
            junctions[key].append((chrom, int(a), int(b)))
blocks = {}
for coll in ("all", "junction"):
    for k, v in json.load(gzip.open(f".cache/pbmc/GSE320507/collection/{coll}/blocks.json.gz", "rt")).items():
        blocks.setdefault(v["block_id"], v)
by_gene = collections.defaultdict(list)
for bid, v in blocks.items():
    by_gene[v["block_id"].split(".")[0]].append(v)
fits = {}
for k in (0, 1):
    fits[k] = pd.concat([fit_table(f".cache/binary_grid/pbmc_GSE320507_{s}_sh_{k}/merged/paired_path.tsv", 4)[["gene_id", "pair_id", "block_id", "p_value", "converged"]] for s in ("internal", "terminal")])
    fits[k]["gene_id"] = fits[k].gene_id.str.split(".").str[0]
rows = []
for r in best.itertuples():
    js = junctions.get(r.feature_id, [])
    hit, status = [], "no block contains the junctions"
    for v in by_gene.get(r.gene_id, []):
        n = len(v["signatures"]); mp = v["paths"][:n]
        _, J, _ = path_structures(mp); var = variable_junctions(mp)
        if any(abs(x - a) <= 1 and abs(y - b) <= 1 for (x, y) in var for (_, a, b) in js):
            st = block_structure(structure_paths(v))["structure"]
            groups = [gr for gr in same_junction_groups(mp) if J[gr[0]] & var]
            testable = st == "internal" or len(groups) >= 2
            hit.append((v["block_id"], st, testable))
    f = fits[r.fold]; ft = f[f.gene_id.eq(r.gene_id) & f.pair_id.eq(r.pair_id)]
    tp = ft[ft.block_id.isin([h[0] for h in hit])]
    if hit:
        status = "block untestable (terminal, one junction group)" if not any(h[2] for h in hit) else ("block not tested on this pair" if tp.empty else "block tested")
    rows.append(dict(gene=r.gene_id, fold=r.fold, pair=r.pair_id, lc_p=r.p_value, status=status, structure=",".join(sorted({h[1] for h in hit})), tealeaf_block_p=tp.p_value.min() if len(tp) else np.nan, tealeaf_pair_best_p=ft.p_value.min() if len(ft) else np.nan, cls="LeafCutter only" if r.gene_id not in rt else "both"))
d = pd.DataFrame(rows); d.to_csv(".cache/adhoc/blood_rep/lc_cluster_blocks.tsv", sep="\t", index=False)
for cls, x in d.groupby("cls"):
    print(cls, len(x)); print(x.status.value_counts().to_string()); print(x.groupby("status").tealeaf_block_p.median().to_string())
    y = x[x.status.eq("block tested")]; print("  tested: block p < 1e-4 %.2f, structure %s" % ((y.tealeaf_block_p < 1e-4).mean(), y.structure.value_counts().to_dict()))
