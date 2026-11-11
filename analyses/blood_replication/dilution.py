"""Tealeaf vs LeafCutter per gene on the LeafCutter matched universe: tests per gene and best test p, by replication class."""
import sys; sys.path.insert(0, ".")
import pandas as pd, numpy as np
from extra_scripts.assess_paired_inference_audit import fit_table
cfg = {"pbmc": ("analyses/coverage_null_filter_sh/pbmc", ".cache/binary_grid/pbmc_GSE320507_{s}_sh_{k}/merged/paired_path.tsv", ("internal", "terminal")),
       "brain": ("analyses/coverage_null_filter_sh/brain", ".cache/binary_grid/{s}_native_v4sh_{k}/merged/paired_path.tsv", ("internal", "tj"))}
for tissue, (F, pattern, structures) in cfg.items():
    g = pd.read_csv(f"{F}/split_filtered/split_gene_tests.tsv.gz", sep="\t"); g = g[g.universe.eq("matched") & g.comparison.eq("LeafCutter")]
    t = g[g.method.eq("Tealeaf")].set_index("gene_id"); c = g[g.method.eq("LeafCutter")].set_index("gene_id")
    rt = set(t.index[t.conjunction_q_value <= .05]); rc = set(c.index[c.conjunction_q_value <= .05])
    lc = []
    for ch in pd.read_csv(f"{F}/comparator_tests_filtered.tsv.gz", sep="\t", chunksize=2_000_000, usecols=["gene_id", "pair_id", "feature_id", "p_value", "method", "fold"]):
        lc.append(ch[ch.method.eq("LeafCutter") & ch.fold.astype(str).isin(["0", "1"])])
    lc = pd.concat(lc); lc["fold"] = lc.fold.astype(int)
    tl = []
    for k in (0, 1):
        for s in structures:
            x = fit_table(pattern.format(s=s, k=k), 4)[["gene_id", "pair_id", "block_id", "p_value"]]
            x["gene_id"] = x.gene_id.str.split(".").str[0]; x["fold"] = k; tl.append(x)
    tl = pd.concat(tl)
    for name, genes in (("both", rc & rt), ("LeafCutter only", rc - rt), ("Tealeaf only", rt - rc)):
        rows = []
        for k in (0, 1):
            a = tl[tl.fold.eq(k) & tl.gene_id.isin(genes)].groupby("gene_id").p_value.agg(["size", "min"])
            b = lc[lc.fold.eq(k) & lc.gene_id.isin(genes)].groupby("gene_id").p_value.agg(["size", "min"])
            rows.append(a.join(b, lsuffix="_T", rsuffix="_L"))
        d = pd.concat(rows)
        print(tissue, name, len(genes), "median tests/gene T %.0f L %.0f; median best p T %.1e L %.1e; median Bonferroni-like n*min T %.1e L %.1e" % (
            d.size_T.median(), d.size_L.median(), d.min_T.median(), d.min_L.median(), (d.size_T * d.min_T).median(), (d.size_L * d.min_L).median()), flush=True)
