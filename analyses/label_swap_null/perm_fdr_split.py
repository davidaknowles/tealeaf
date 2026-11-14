"""Split replication in GSE320507 with a permutation FDR from the within-subject label-swap null, in place of BH.

Gene statistic, as in split_universe_assessment: in half k, test p-values are combined by Simes within (gene, cell-type
pair) and then within gene over the universe keys (missing keys p = 1), giving s_k; a gene's replication statistic is
r = max(s_0, s_1). The null statistic r* is computed the same way from the swapped-label fits of both halves (the same
swap per pair and donor in both halves), over the same keys. For threshold t the estimated FDR is
#{genes with r* <= t} / #{genes with r <= t} (the gene sets are equal), averaged over swap draws; a method replicates the
genes with r <= t for the largest real t with estimated FDR <= .05. The statistics only rank genes.

Usage: perm_fdr_split.py <output dir> <pair set: level_pairs.tsv> <null dir> [<null dir> ...]
Each null dir holds fold{0,1}/{leafcutter/normalized/*.tsv, scquint/*.tsv, rmats/*.tsv.gz, tealeaf_{internal,terminal}/merged/paired_path.tsv}.
Real tests: unfiltered comparator universe (split_universe_leaf5) and the pbmc_GSE320507_*_sh half fits; with
--coverage-filtered, the real tests are the coverage-null-filtered tables of analyses/coverage_null_filter_sh/pbmc, and each
swapped-label fit is filtered by the same rule from its own coverage-only null (submit_coverage_null.sh: fold{k}/tnull_*),
so the filter is part of the procedure that the permutation FDR calibrates."""
import glob, json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from extra_scripts.assess_paired_inference_audit import fit_table
from extra_scripts.evaluate_suppa2_statistics import grouped_pvalues
from tealeaf.sc.coverage_null import coverage_flagged, drop_flagged
from tealeaf.sc.junction_benchmark import benjamini_hochberg

args = [a for a in sys.argv[1:] if not a.startswith("--")]
FILTERED = "--coverage-filtered" in sys.argv
OUT, PAIRS, NULLS = Path(args[0]), Path(args[1]), [Path(a) for a in args[2:]]
OUT.mkdir(parents=True, exist_ok=True)
pairs = {"||".join(sorted(p)) for p in pd.read_csv(PAIRS, sep="\t")[["level_a", "level_b"]].itertuples(index=False)}
METHODS = {"LeafCutter": "leafcutter/normalized/*.tsv", "scQuint": "scquint/*.tsv", "rMATS (paired)": "rmats/*.tsv.gz"}


def gene(values):
    return pd.Series(values).astype(str).str.split(".").str[0].to_numpy()


def read_all(pattern):
    frames = []
    for path in sorted(glob.glob(pattern)):
        try:
            frames.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
    return pd.concat(frames, ignore_index=True)


def tealeaf_half(root, k, null_root=None):
    """Tealeaf tests of half k (real or swapped labels), without the tests flagged by the coverage-only null in
    null_root(structure, k) (null.tsv.gz of assess_targeted_count_null.py) when filtering."""
    t = pd.concat([fit_table(root(s, k)) for s in ("internal", "terminal")], ignore_index=True)
    if FILTERED:
        null = pd.concat([pd.read_csv(null_root(s, k), sep="\t") for s in ("internal", "terminal")], ignore_index=True)
        t = drop_flagged(t, coverage_flagged(null, ["test_id"]), ["test_id"])
    return pd.DataFrame({"gene_id": gene(t.gene_id), "pair_id": t.pair_id, "p_value": t.p_value})


def comparator_flags(root, k, method):
    """(pair_id, feature_id) of half k of a swapped-label draw that its coverage-only null calls."""
    contrasts = json.loads((root / f"fold{k}_{'rmats' if method.startswith('rMATS') else 'junction'}_contrasts.json").read_text())
    pair_of = {c["contrast_id"]: "||".join(sorted((c["level_a"], c["level_b"]))) for c in contrasts}
    name = {"LeafCutter": "leafcutter", "scQuint": "scquint", "rMATS (paired)": "rmats"}[method]
    table = read_all(str(root / f"fold{k}/tnull_{name}/*.tsv"))
    if not len(table):
        return set()
    if method.startswith("rMATS"):
        table["feature_id"] = "rMATS:" + table.event_type + ":" + table.event_id.astype(str)
    flagged = coverage_flagged(table, ["contrast_id", "feature_id"])
    return set(zip(flagged.contrast_id.map(pair_of), flagged.feature_id))


# real tests
universe = pd.read_csv("analyses/coverage_null_filter_sh/pbmc/comparator_tests_filtered.tsv.gz" if FILTERED else "analyses/pbmc/GSE320507_comparators/split_universe_leaf5/comparator_tests.tsv.gz", sep="\t", dtype={"gene_id": str}, low_memory=False)
universe["gene_id"] = gene(universe.gene_id)
universe = universe.loc[universe.pair_id.isin(pairs)]
real = {"Tealeaf": [tealeaf_half(lambda s, k: f".cache/binary_grid/pbmc_GSE320507_{s}_sh_{k}/merged/paired_path.tsv", k, lambda s, k: f".cache/targeted_null/pbmc_GSE320507_{s}_sh_{k}/null.tsv.gz") for k in (0, 1)]}
real["Tealeaf"] = [t.loc[t.pair_id.isin(pairs)] for t in real["Tealeaf"]]
genes_of = universe.drop_duplicates(["comparison", "feature_id"]).set_index(["comparison", "feature_id"]).gene_id

# null tests, per draw
nulls = []
for root in NULLS:
    draw = {"Tealeaf": [tealeaf_half(lambda s, k: root / f"fold{k}/tealeaf_{s}/merged/paired_path.tsv", k, lambda s, k: root / f"fold{k}/tnull_tealeaf_{s}/null.tsv.gz") for k in (0, 1)]}
    for method, pattern in METHODS.items():
        halves = []
        for k in (0, 1):
            t = read_all(str(root / f"fold{k}" / pattern))
            pair = ["||".join(sorted(x)) for x in zip(t.level_a, t.level_b)]
            g = [genes_of.get((method, f), np.nan) for f in t.feature_id]
            half = pd.DataFrame({"gene_id": g, "pair_id": pair, "feature_id": t.feature_id, "p_value": t.p_value.fillna(1.)}).dropna(subset=["gene_id"])
            if FILTERED:
                flagged = comparator_flags(root, k, method)
                half = half.loc[[(p, f) not in flagged for p, f in zip(half.pair_id, half.feature_id)]]
            halves.append(half[["gene_id", "pair_id", "p_value"]])
        draw[method] = halves
    nulls.append(draw)


def gene_stat(halves, keys):
    """max over halves of the gene-level Simes statistic over keys (gene_id, pair_id); missing keys p = 1."""
    stats = []
    for t in halves:
        t = t.merge(keys, on=["gene_id", "pair_id"])
        pair = keys.merge(grouped_pvalues(t, ["gene_id", "pair_id"], "simes"), on=["gene_id", "pair_id"], how="left").fillna({"p_value": 1.})
        stats.append(grouped_pvalues(pair, ["gene_id"], "simes").set_index("gene_id").p_value)
    s = pd.concat(stats, axis=1).fillna(1.)
    return s.max(axis=1)


def permutation_fdr(r, null_rs, level=.05):
    """Number of genes replicated at estimated FDR <= level, the threshold, and the null count at it."""
    order = np.sort(r.to_numpy())
    best = (0, np.nan, np.nan)
    for t in np.unique(order):
        n_real = np.searchsorted(order, t, side="right")
        n_null = np.mean([np.searchsorted(np.sort(x.to_numpy()), t, side="right") for x in null_rs])
        if n_null / n_real <= level:
            best = (n_real, t, n_null)
    return best


rows = []
tealeaf_keys = [set(zip(t.gene_id, t.pair_id)) for t in real["Tealeaf"]]
for comparison, comp in universe.groupby("comparison"):
    if comparison not in METHODS:
        continue
    comp_halves = [comp.loc[comp.fold.eq(k), ["gene_id", "pair_id", "p_value"]] for k in (0, 1)]
    comp_keys = set(zip(comp_halves[0].gene_id, comp_halves[0].pair_id)) & set(zip(comp_halves[1].gene_id, comp_halves[1].pair_id))
    for scope, keys in (("matched", comp_keys & tealeaf_keys[0] & tealeaf_keys[1]), ("comparator", comp_keys)):
        frame = pd.DataFrame(sorted(keys), columns=["gene_id", "pair_id"])
        for method, halves in (("Tealeaf", real["Tealeaf"]), (comparison, comp_halves)):
            r = gene_stat(halves, frame)
            null_rs = [gene_stat(draw[method], frame).reindex(r.index).fillna(1.) for draw in nulls]
            n, t, n_null = permutation_fdr(r, null_rs)
            bh = int((benjamini_hochberg(r.to_numpy()) <= .05).sum())
            null_bh = [int((benjamini_hochberg(x.to_numpy()) <= .05).sum()) for x in null_rs]
            rows.append(dict(comparison=comparison, method=method, universe=scope, genes=len(r), replicated_perm_fdr=n, threshold=t, null_genes_at_threshold=n_null, replicated_bh=bh, null_replicated_bh=np.mean(null_bh)))
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / f"perm_fdr_split{'_filtered' if FILTERED else ''}.tsv", sep="\t", index=False)
