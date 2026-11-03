"""Split and long-read endpoints with the targeted coverage-only null filter applied to every method.

A test of any method is removed when it was simulated (observed p < .01) and the median of its three null p-values is
below .05 (tealeaf.sc.coverage_null). SUPPA2 and primer-aware SUPPA2 are not simulated and stay unchanged.
Usage: coverage_filter_eval.py brain|pbmc <output dir> [--skip-lr]
"""
import glob, gzip, json, sys
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from extra_scripts.assess_paired_inference_audit import fit_table, split_universe_assessment  # noqa: E402
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table  # noqa: E402
from tealeaf.sc.coverage_null import coverage_flagged, drop_flagged  # noqa: E402
from tealeaf.sc.local_path_reads import block_structure, structure_paths  # noqa: E402
from tealeaf.sc.majiq_coverage import lsv_pvalues_without_coverage_calls  # noqa: E402
from tealeaf.sc.replication_audit import rank_direction_table, ranked_direction_summary  # noqa: E402

DATA, OUT = sys.argv[1], Path(sys.argv[2])
SKIP_LR = "--skip-lr" in sys.argv
OUT.mkdir(parents=True, exist_ok=True)
T = REPO / ".cache/targeted_null"
CACHE = REPO / ".cache/binary_grid"
if DATA == "brain":
    ROOT = Path.home() / "tealeaf_runs/microglia_less/run/junction_benchmark"
    MODELS = ("internal_native_v4", "tj_native_v4")
    MAJIQ = {"full": REPO / "analyses/brain_majiq_coverage_null/full/majiq_coverage_null.parquet",
             "0": REPO / ".cache/brain/majiq_coverage_null_fold0/majiq_coverage_null.parquet", "1": REPO / ".cache/brain/majiq_coverage_null_fold1/majiq_coverage_null.parquet"}
    UNIVERSE = REPO / "analyses/split_universe/comparator_tests.tsv.gz"
else:
    ROOT = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
    MODELS = ("pbmc_GSE320507_internal", "pbmc_GSE320507_terminal")
    MAJIQ = {run: REPO / f".cache/pbmc/GSE320507/majiq_coverage_null_{'full' if run == 'full' else 'fold' + run}/majiq_coverage_null.parquet" for run in ("full", "0", "1")}
    UNIVERSE = REPO / "analyses/pbmc/GSE320507_comparators/split_universe/comparator_tests.tsv.gz"

pair_of = {}
for path in [ROOT / "contrasts.json", ROOT / "reproducibility/fold0/contrasts.json", ROOT / "reproducibility/fold1/contrasts.json"]:
    for c in json.loads(path.read_text()):
        pair_of[c["contrast_id"]] = "||".join(sorted([c["level_a"], c["level_b"]]))


def read_dir(pattern):
    frames = []
    for path in sorted(glob.glob(pattern)):
        try:
            frames.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def flags(run):
    """Flagged keys per method for one run (full, 0, 1): tests simulated and called by the coverage-only null."""
    out, counts = {}, []
    tealeaf = pd.concat([pd.read_csv(T / f"{model}_{run}/null.tsv.gz", sep="\t") for model in MODELS], ignore_index=True)
    out["Tealeaf"] = coverage_flagged(tealeaf, ["test_id"])
    counts.append(("Tealeaf", len(tealeaf), len(out["Tealeaf"])))
    for method, name in (("LeafCutter", "leafcutter"), ("scQuint", "scquint")):
        table = read_dir(str(T / f"{DATA}_{name}_{run}/*.tsv"))
        out[method] = coverage_flagged(table, ["contrast_id", "feature_id"]) if len(table) else pd.DataFrame(columns=["contrast_id", "feature_id"])
        counts.append((method, len(table), len(out[method])))
    table = read_dir(str(T / f"rmats_{DATA}_{run}/*.tsv"))
    if len(table):
        table["feature_id"] = "rMATS:" + table.event_type + ":" + table.event_id.astype(str)
        out["rMATS (paired)"] = coverage_flagged(table, ["contrast_id", "feature_id"])
    else:
        out["rMATS (paired)"] = pd.DataFrame(columns=["contrast_id", "feature_id"])
    counts.append(("rMATS (paired)", len(table), len(out["rMATS (paired)"])))
    majiq = pd.read_parquet(MAJIQ[run]).dropna(subset=["feature_id"])
    out["MAJIQ Heterogen"] = coverage_flagged(majiq, ["contrast_id", "feature_id"], p_column="raw_pvalue", threshold=.01)
    counts.append(("MAJIQ Heterogen", int((majiq.raw_pvalue < .01).sum()), len(out["MAJIQ Heterogen"])))
    for flagged in out.values():
        if "contrast_id" in flagged:
            flagged["pair_id"] = flagged.contrast_id.map(pair_of)
    return out, pd.DataFrame(counts, columns=["method", "simulated", "removed"]).assign(run=run), majiq


# Tealeaf fold tables (as in the reported split)
def tealeaf_folds():
    if DATA == "brain":
        blocks = json.load(gzip.open(REPO / ".cache/native_v4/collection_parity_p2_v4/all/blocks.json.gz", "rt"))
        internal = {b["block_id"] for b in blocks.values() if block_structure(structure_paths(b))["structure"] == "internal"}
        folds = []
        for k in (0, 1):
            a = fit_table(CACHE / f"internal_native_v4_{k}/merged/paired_path.tsv")
            t = fit_table(CACHE / f"tj_native_v4_{k}/merged/paired_path.tsv")
            folds.append(pd.concat([a.loc[a.block_id.isin(internal)], t.loc[~t.block_id.isin(internal)]], ignore_index=True))
        return folds
    return [pd.concat([fit_table(CACHE / f"{m}_{k}/merged/paired_path.tsv") for m in MODELS], ignore_index=True) for k in (0, 1)]


all_flags, removal = {}, []
for run in ("full", "0", "1"):
    all_flags[run], counts, _ = flags(run)
    removal.append(counts)
removal = pd.concat(removal, ignore_index=True)
removal.to_csv(OUT / "removed_tests.tsv", sep="\t", index=False)
print(removal.to_string(index=False), flush=True)

# split
folds = tealeaf_folds()
filtered_folds = [drop_flagged(fold, all_flags[str(k)]["Tealeaf"], ["test_id"]) for k, fold in enumerate(folds)]
universe = pd.read_csv(UNIVERSE, sep="\t", dtype={"gene_id": str}, low_memory=False)
parts = []
for comparison, local in universe.groupby("comparison"):
    if comparison in ("LeafCutter", "scQuint", "MAJIQ Heterogen", "rMATS (paired)"):
        kept = []
        for k in (0, 1):
            f = all_flags[str(k)][comparison][["pair_id", "feature_id"]]
            kept.append(drop_flagged(local.loc[local.fold.eq(k)], f, ["pair_id", "feature_id"]))
        local = pd.concat(kept, ignore_index=True)
        both = local.groupby(["gene_id", "pair_id"]).fold.nunique().pipe(lambda s: s[s == 2]).index
        local = local.set_index(["gene_id", "pair_id"]).loc[lambda d: d.index.isin(both)].reset_index()
    parts.append(local[universe.columns])
filtered_universe = pd.concat(parts, ignore_index=True)
filtered_universe.to_csv(OUT / "comparator_tests_filtered.tsv.gz", sep="\t", index=False)
(OUT / "split_unfiltered").mkdir(exist_ok=True)
(OUT / "split_filtered").mkdir(exist_ok=True)
split_universe_assessment(folds, UNIVERSE, OUT / "split_unfiltered", "Tealeaf")
split_universe_assessment(filtered_folds, OUT / "comparator_tests_filtered.tsv.gz", OUT / "split_filtered", "Tealeaf")

if SKIP_LR or DATA != "brain":
    sys.exit(0)

# long reads (brain)
lr = pd.read_csv(REPO / "analyses/native_align/internal_native_v4+tj_native_v4_lrchain_groupfix/lr_mapping.tsv.gz", sep="\t", low_memory=False)
lr = lr.loc[lr.eligible.astype(str).str.lower().eq("true")].assign(method="Tealeaf", feature_id=lambda d: d.test_id, contrast_id=lambda d: d.test_id)
lr["pooled_replicated"] = lr.pooled_replicated.astype(str).str.lower().eq("true")
tables = {"unfiltered": [lr], "filtered": [drop_flagged(lr, all_flags["full"]["Tealeaf"], ["test_id"])]}
junction = evaluable(read_source_table(REPO / "analyses/tilgner_junction_replication/all_tested_comprehensive/junction_replication.tsv.gz"))
rmats = evaluable(read_source_table(REPO / "analyses/comparator_suppa_rmats/rmats_paired_full_tilgner_comprehensive.tsv.gz"))
majiq_lr = evaluable(read_source_table(REPO / "analyses/tilgner_junction_replication/all_tested_comprehensive/majiq/junction_replication.tsv.gz"))
tables["unfiltered"] += [junction, rmats, majiq_lr]
f = all_flags["full"]
junction_f = pd.concat([drop_flagged(junction.loc[junction.method.eq(m)], f[m], ["contrast_id", "feature_id"]) for m in ("LeafCutter", "scQuint")], ignore_index=True)
rmats_f = drop_flagged(rmats, f["rMATS (paired)"], ["contrast_id", "feature_id"])
null = pd.read_parquet(MAJIQ["full"]).dropna(subset=["feature_id"])
null.loc[null.raw_pvalue >= .01, [c for c in null if c.startswith("p_null")]] = 1.
lsv = lsv_pvalues_without_coverage_calls(null)
majiq_f = majiq_lr.merge(lsv[["lsv", "p_value"]].rename(columns={"p_value": "p_kept"}), left_on="feature_id", right_on="lsv")
majiq_f["p_value"] = majiq_f["raw_p_value"] = majiq_f.p_kept
tables["filtered"] += [junction_f, rmats_f, majiq_f.drop(columns=["lsv", "p_kept"])]
summaries, ranks = [], []
for label, group in tables.items():
    ranked = pd.concat([rank_direction_table(t, 200) for t in group], ignore_index=True, sort=False)
    s = pd.DataFrame(ranked_direction_summary(ranked)).assign(table=label)
    s["n_evaluable"] = s.method.map(pd.concat([t.groupby("method").size() for t in group]))
    summaries.append(s)
    ranks.append(ranked.assign(table=label))
summary = pd.concat(summaries, ignore_index=True)
summary.to_csv(OUT / "lr_summary.tsv", sep="\t", index=False)
pd.concat(ranks, ignore_index=True).to_csv(OUT / "lr_ranks.tsv.gz", sep="\t", index=False)
print(summary.pivot_table(index=["method", "cutoff"], columns="table", values=["n_agree", "normalized_auc"]).to_string(), flush=True)
