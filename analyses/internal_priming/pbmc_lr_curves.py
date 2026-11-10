"""GSE320507 long-read curve inputs from the filtered ranks of pbmc_lr_filter.py (shrunk-offset evaluation)."""
import pandas as pd
F = "analyses/coverage_null_filter_sh/pbmc"
r = pd.read_csv(f"{F}/lr_ranks.tsv.gz", sep="\t", low_memory=False)
r = r.loc[r.table.eq("filtered") & r["rank"].le(200)].copy()
r["method"] = r.method.replace({"SUPPA2 (full data); exact paired Wilcoxon": "SUPPA2 (full data)", "SUPPA2 (primer aware); exact paired Wilcoxon": "SUPPA2 (primer aware)"})
r.loc[r.method.ne("Tealeaf")].to_csv(f"{F}/sign_agreement_by_rank_filtered.tsv", sep="\t", index=False, na_rep="NA")
r.loc[r.method.eq("Tealeaf")].assign(scope="all tested", ranking="continuous raw").to_csv(f"{F}/tealeaf_lr_rank_filtered.tsv.gz", sep="\t", index=False)
print(r.groupby("method").size())
