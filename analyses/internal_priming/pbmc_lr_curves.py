"""GSE320507 long-read curve inputs from the filtered ranks of pbmc_lr_filter.py (shrunk-offset evaluation).
Usage: pbmc_lr_curves.py [depth] (default 200; other depths read and write the _depth<N> files)."""
import sys
import pandas as pd
F = "analyses/coverage_null_filter_sh/pbmc"
DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 else 200
SUFFIX = "" if DEPTH == 200 else f"_depth{DEPTH}"
r = pd.read_csv(f"{F}/lr_ranks{SUFFIX}.tsv.gz", sep="\t", low_memory=False)
r = r.loc[r.table.eq("filtered") & r["rank"].le(DEPTH)].copy()
r["method"] = r.method.replace({"SUPPA2 (full data); exact paired Wilcoxon": "SUPPA2 (full data)", "SUPPA2 (primer aware); exact paired Wilcoxon": "SUPPA2 (primer aware)"})
r.loc[r.method.ne("Tealeaf")].to_csv(f"{F}/sign_agreement_by_rank_filtered{SUFFIX}.tsv", sep="\t", index=False, na_rep="NA")
r.loc[r.method.eq("Tealeaf")].assign(scope="all tested", ranking="continuous raw").to_csv(f"{F}/tealeaf_lr_rank_filtered{SUFFIX}.tsv.gz", sep="\t", index=False)
print(r.groupby("method").size())
