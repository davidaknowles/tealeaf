"""Full-data coverage-null flags (targeted rule) for plot_discoveries_by_event.py --coverage-flagged, and a MAJIQ null
parquet with p_null = 1 for tests with raw p >= .01 (not simulated under the targeted rule).
Usage: coverage_flags_full.py brain|pbmc <output dir>"""
import glob, sys
from pathlib import Path
import pandas as pd
from tealeaf.sc.coverage_null import coverage_flagged

DATA, OUT = sys.argv[1], Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
T = Path(".cache/targeted_null")
models = ("internal_native_v4", "tj_native_v4") if DATA == "brain" else ("pbmc_GSE320507_internal", "pbmc_GSE320507_terminal")
majiq = Path("analyses/brain_majiq_coverage_null/full/majiq_coverage_null.parquet") if DATA == "brain" else Path(".cache/pbmc/GSE320507/majiq_coverage_null_full/majiq_coverage_null.parquet")
read = lambda pattern: pd.concat([pd.read_csv(p, sep="\t") for p in sorted(glob.glob(pattern)) if Path(p).stat().st_size > 0], ignore_index=True)
parts = []
tealeaf = pd.concat([pd.read_csv(T / f"{m}_full/null.tsv.gz", sep="\t") for m in models], ignore_index=True)
parts.append(coverage_flagged(tealeaf, ["test_id"]).rename(columns={"test_id": "key"}).assign(method="Tealeaf"))
for method, name in (("LeafCutter", "leafcutter"), ("scQuint", "scquint")):
    f = coverage_flagged(read(str(T / f"{DATA}_{name}_full/*.tsv")), ["contrast_id", "feature_id"])
    parts.append(pd.DataFrame({"method": method, "key": f.contrast_id + "|" + f.feature_id.astype(str)}))
r = read(str(T / f"rmats_{DATA}_full/*.tsv"))
r["feature_id"] = "rMATS:" + r.event_type + ":" + r.event_id.astype(str)
f = coverage_flagged(r, ["contrast_id", "feature_id"])
parts.append(pd.DataFrame({"method": "rMATS", "key": f.contrast_id + "|" + f.feature_id}))
flags = pd.concat(parts, ignore_index=True)[["method", "key"]]
flags.to_csv(OUT / "flagged_full.tsv.gz", sep="\t", index=False)
print(flags.groupby("method").size().to_string())
null = pd.read_parquet(majiq)
null.loc[~(null.raw_pvalue < .01), [c for c in null if c.startswith("p_null")]] = 1.
null.to_parquet(OUT / "majiq_coverage_null_targeted.parquet")
print("MAJIQ flagged", int((null[[c for c in null if c.startswith("p_null")]].median(axis=1) < .05).sum()))
