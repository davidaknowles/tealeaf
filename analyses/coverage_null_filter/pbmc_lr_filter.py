"""GSE320507 long-read direction ranks with the targeted coverage-only null filter (as coverage_filter_eval.py for brain).
Tealeaf, LeafCutter, scQuint, rMATS and MAJIQ full-data tests the null calls are removed; SUPPA2 is not simulated.
Usage: pbmc_lr_filter.py <coverage_null_filter/pbmc dir>"""
import sys
from pathlib import Path
import pandas as pd
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from extra_scripts.summarize_tilgner_comparator_ranks import evaluable, read_source_table  # noqa: E402
from tealeaf.sc.coverage_null import drop_flagged  # noqa: E402
from tealeaf.sc.majiq_coverage import lsv_pvalues_without_coverage_calls  # noqa: E402
from tealeaf.sc.replication_audit import rank_direction_table, ranked_direction_summary  # noqa: E402

F = Path(sys.argv[1])
O = REPO / "analyses/pbmc/GSE320507_comparators"
flags = pd.read_csv(F / "flagged_full.tsv.gz", sep="\t")
keys = lambda method: flags.loc[flags.method.eq(method), "key"]
lr = pd.read_csv(REPO / __import__("os").environ.get("PBMC_LR_MAPPING", "analyses/pbmc/GSE320507_lr_groupfix/lr_mapping.tsv.gz"), sep="\t", low_memory=False)
lr = lr.loc[lr.eligible.astype(str).str.lower().eq("true")].assign(method="Tealeaf", feature_id=lambda d: d.test_id, contrast_id=lambda d: d.test_id)
lr["pooled_replicated"] = lr.pooled_replicated.astype(str).str.lower().eq("true")
junction = evaluable(read_source_table(O / ("lr_junction_leaf5/junction_replication_merged.tsv.gz" if __import__("os").environ.get("LEAF5") == "1" else "lr_junction/junction_replication.tsv.gz")))
events = [evaluable(read_source_table(O / f"lr_event/{name}_full_data_lr.tsv.gz")) for name in ("rmats", "suppa2", "suppa2_primer_aware")]
unfiltered = [lr, junction.loc[junction.method.ne("MAJIQ Heterogen")], junction.loc[junction.method.eq("MAJIQ Heterogen")], *events]
drop = lambda table, method: table.loc[~(table.contrast_id + "|" + table.feature_id.astype(str)).isin(set(keys(method)))]
tealeaf_f = lr.loc[~lr.test_id.isin(set(keys("Tealeaf")))]
junction_f = pd.concat([drop(junction.loc[junction.method.eq(m)], m) for m in ("LeafCutter", "scQuint")], ignore_index=True)
rmats_f = drop(events[0], "rMATS")
null = pd.read_parquet(F / "majiq_coverage_null_targeted.parquet").dropna(subset=["feature_id"])
lsv = lsv_pvalues_without_coverage_calls(null)
majiq = junction.loc[junction.method.eq("MAJIQ Heterogen")]
majiq_f = majiq.merge(lsv[["lsv", "p_value"]].rename(columns={"p_value": "p_kept"}), left_on="feature_id", right_on="lsv")
majiq_f["p_value"] = majiq_f["raw_p_value"] = majiq_f.p_kept
print("MAJIQ LSVs", len(majiq), "matched", len(majiq_f), flush=True)
tables = {"unfiltered": unfiltered, "filtered": [tealeaf_f, junction_f, majiq_f.drop(columns=["lsv", "p_kept"]), rmats_f, *events[1:]]}
DEPTH = int(__import__("os").environ.get("LR_DEPTH", "200"))  # ranks kept per method; outputs get a _depth<N> suffix unless 200
CUTOFFS = tuple(c for c in (40, 100, 200, 500, 1000, 2000) if c <= DEPTH)
SUFFIX = "" if DEPTH == 200 else f"_depth{DEPTH}"
summaries, ranks = [], []
for label, group in tables.items():
    ranked = pd.concat([rank_direction_table(t, DEPTH) for t in group], ignore_index=True, sort=False)
    s = pd.DataFrame(ranked_direction_summary(ranked, CUTOFFS)).assign(table=label)
    s["n_evaluable"] = s.method.map(pd.concat([t.groupby("method").size() for t in group]))
    summaries.append(s)
    ranks.append(ranked.assign(table=label))
summary = pd.concat(summaries, ignore_index=True)
summary.to_csv(F / f"lr_summary{SUFFIX}.tsv", sep="\t", index=False)
pd.concat(ranks, ignore_index=True).to_csv(F / f"lr_ranks{SUFFIX}.tsv.gz", sep="\t", index=False)
print(summary.pivot_table(index=["method", "cutoff"], columns="table", values=["n_agree", "normalized_auc", "n_evaluable"]).to_string(), flush=True)
