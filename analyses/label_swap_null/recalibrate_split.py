"""Split replication in GSE320507 with each method's half-data p-values recalibrated by its within-subject label-swap null.

The null (make_contrasts_spread.py, submit_spread.sh) refits 21 pairs of seven cell types spread over abundance in each
half with each donor's two labels exchanged with probability 1/2. A test of pair P in half k of method m gets the
calibrated p-value min(1, c p^a), with F(t) = c t^a fitted to the lower tail (1e-4 to .05) of the null p-values of m in
half k whose pairs share P's abundance stratum (cells of the rarer cell type: >= 8,000, 1,500-8,000, < 1,500). Tealeaf, LeafCutter, scQuint and rMATS are
recalibrated; MAJIQ and the SUPPA2 forms keep their p-values. Replication is then scored as in coverage_filter_eval.py
(split_universe_assessment on the coverage-filtered tables). Also writes false replications: genes replicated across
the halves by each method with the swapped labels (keys tested in both halves of the 21 pairs).
Usage: recalibrate_split.py <output dir>; env TAIL_MIN (default 1e-4) sets the lowest fitted threshold, PENALIZE_ONLY=1
keeps p-values from decreasing."""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from extra_scripts.assess_paired_inference_audit import fit_table, split_universe_assessment
from extra_scripts.evaluate_suppa2_statistics import grouped_pvalues, metric_row
from tealeaf.sc.coverage_null import coverage_flagged, drop_flagged

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
N = Path(".cache/label_swap2"); F = Path("analyses/coverage_null_filter_sh/pbmc"); T = Path(".cache/targeted_null")
B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
cells = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t").groupby("cell_type").n_cells.sum()


def stratum(pair_ids):
    low = np.array([min(cells[a], cells[b]) for a, b in (p.split("||") for p in pair_ids)])
    return np.where(low >= 8000, "high", np.where(low >= 1500, "mid", "low"))


def read_all(pattern):
    frames = []
    for path in sorted(glob.glob(pattern)):
        try:
            frames.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:  # a contrast without testable features in this half
            print("empty:", path, flush=True)
    return pd.concat(frames, ignore_index=True)


def pair_id(a, b):
    return ["||".join(sorted(x)) for x in zip(a, b)]


def null_tables(k):
    """Null p-values of half k per method: pair_id, feature_id, gene_id (Tealeaf), p_value."""
    out = {}
    lc = read_all(str(N / f"fold{k}/leafcutter/normalized/*.tsv"))
    out["LeafCutter"] = pd.DataFrame({"pair_id": pair_id(lc.level_a, lc.level_b), "feature_id": lc.feature_id, "p_value": lc.p_value.fillna(1.)})
    sq = read_all(str(N / f"fold{k}/scquint/*.tsv"))
    out["scQuint"] = pd.DataFrame({"pair_id": pair_id(sq.level_a, sq.level_b), "feature_id": sq.feature_id, "p_value": sq.p_value.fillna(1.)})
    rm = read_all(str(N / f"fold{k}/rmats/*.tsv.gz"))
    out["rMATS (paired)"] = pd.DataFrame({"pair_id": pair_id(rm.level_a, rm.level_b), "feature_id": rm.feature_id, "p_value": rm.p_value.fillna(1.)})
    tl = pd.concat([fit_table(N / f"fold{k}/tealeaf_{s}/merged/paired_path.tsv") for s in ("internal", "terminal")], ignore_index=True)
    out["Tealeaf"] = tl[["pair_id", "feature_id", "gene_id", "p_value"]].copy() if "feature_id" in tl else tl.assign(feature_id=tl.test_id)[["pair_id", "feature_id", "gene_id", "p_value"]]
    for table in out.values():
        table["stratum"] = stratum(table.pair_id)
    return out


import os
TAIL = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2])
TAIL = TAIL[TAIL >= float(os.environ.get("TAIL_MIN", "1e-4"))]  # lowest threshold of the fit
PENALIZE_ONLY = os.environ.get("PENALIZE_ONLY") == "1"  # never make a p-value smaller (no credit for conservative tests)


def tail_fit(null_p):
    """Power-law fit of the null CDF in its lower tail, F(t) = c t^a, by least squares on log F(t) at TAIL thresholds
    with at least 5 null p-values below them. Returns (c, a)."""
    null_p = np.asarray(null_p, float)
    counts = np.array([(null_p < t).sum() for t in TAIL])
    keep = counts >= 5
    a, log_c = np.polyfit(np.log(TAIL[keep]), np.log(counts[keep] / len(null_p)), 1)
    return float(np.exp(log_c)), float(a)


FITS = {}


def recalibrate(table, null, label):
    """Calibrated p-values min(1, c p^a) of table (pair_id, p_value), with (c, a) fitted to the null p-values of the same
    abundance stratum (tail_fit); the power law extrapolates below the smallest null p-value."""
    p = table.p_value.to_numpy(float).copy()
    strata = stratum(table.pair_id)
    for s in np.unique(strata):
        ref = null.loc[null.stratum.eq(s), "p_value"].to_numpy(float)
        c, a = tail_fit(ref)
        FITS[(*label, s)] = (c, a, len(ref))
        mask = strata == s
        p[mask] = np.minimum(1., c * np.power(p[mask], a))
        if PENALIZE_ONLY:
            p[mask] = np.maximum(p[mask], table.p_value.to_numpy(float)[mask])
    return p


nulls = [null_tables(k) for k in (0, 1)]
summary = []
for k in (0, 1):
    for method, t in nulls[k].items():
        for s, g in t.groupby("stratum"):
            summary.append(dict(fold=k, method=method, stratum=s, tests=len(g), p01=(g.p_value < .01).mean(), p001=(g.p_value < .001).mean()))
summary = pd.DataFrame(summary); summary.to_csv(OUT / "null_rates.tsv", sep="\t", index=False)
print(summary.round(4).to_string(index=False), flush=True)

# coverage-filtered Tealeaf half tables and comparator universe, as in coverage_filter_eval.py (pbmc, SUFFIX=_sh)
folds = []
for k in (0, 1):
    fold = pd.concat([fit_table(Path(".cache/binary_grid") / f"pbmc_GSE320507_{s}_sh_{k}/merged/paired_path.tsv") for s in ("internal", "terminal")], ignore_index=True)
    null = pd.concat([pd.read_csv(T / f"pbmc_GSE320507_{s}_sh_{k}/null.tsv.gz", sep="\t") for s in ("internal", "terminal")], ignore_index=True)
    folds.append(drop_flagged(fold, coverage_flagged(null, ["test_id"]), ["test_id"]))
universe = pd.read_csv(F / "comparator_tests_filtered.tsv.gz", sep="\t", dtype={"gene_id": str}, low_memory=False)
calibrated_folds = [fold.assign(p_value=recalibrate(fold, nulls[k]["Tealeaf"], ("Tealeaf", k))) for k, fold in enumerate(folds)]
parts = []
for (comparison, k), local in universe.groupby(["comparison", "fold"]):
    local = local.copy()
    if comparison in ("LeafCutter", "scQuint", "rMATS (paired)"):
        local["p_value"] = recalibrate(local, nulls[int(k)][comparison], (comparison, int(k)))
    parts.append(local)
fits = pd.DataFrame([dict(method=m, fold=k, stratum=st, c=c, a=a, null_tests=n) for (m, k, st), (c, a, n) in FITS.items()]).drop_duplicates(["method", "fold", "stratum"])
fits.to_csv(OUT / "tail_fits.tsv", sep="\t", index=False)
print(fits.round(3).to_string(index=False), flush=True)
pd.concat(parts, ignore_index=True).to_csv(OUT / "comparator_tests_recalibrated.tsv.gz", sep="\t", index=False)
(OUT / "split_recalibrated").mkdir(exist_ok=True)
split_universe_assessment(calibrated_folds, OUT / "comparator_tests_recalibrated.tsv.gz", OUT / "split_recalibrated", "Tealeaf")

# false replications with the swapped labels
genes = universe.drop_duplicates(["comparison", "feature_id"]).set_index(["comparison", "feature_id"]).gene_id
rows = []
for method in ("Tealeaf", "LeafCutter", "scQuint", "rMATS (paired)"):
    tables = []
    for k in (0, 1):
        t = nulls[k][method].copy()
        if method != "Tealeaf":
            t["gene_id"] = [genes.get((method, f), np.nan) for f in t.feature_id]
        t["gene_id"] = t.gene_id.astype(str).str.split(".").str[0]
        tables.append(t.loc[t.gene_id.ne("nan")])
    keys = set(zip(tables[0].gene_id, tables[0].pair_id)) & set(zip(tables[1].gene_id, tables[1].pair_id))
    frame = pd.DataFrame(sorted(keys), columns=["gene_id", "pair_id"])
    pairs = [frame.merge(grouped_pvalues(t.merge(frame, on=["gene_id", "pair_id"]), ["gene_id", "pair_id"], "simes"), on=["gene_id", "pair_id"], how="left").fillna({"p_value": 1.}) for t in tables]
    metric, _ = metric_row(method, method, "simes", "simes", pairs)
    rows.append(dict(method=method, genes=metric["shared_genes"], half1_bh=metric["fold0_bh"], half2_bh=metric["fold1_bh"], replicated=metric["replicated_bh"]))
false = pd.DataFrame(rows); false.to_csv(OUT / "false_replications.tsv", sep="\t", index=False)
print(false.to_string(index=False))
