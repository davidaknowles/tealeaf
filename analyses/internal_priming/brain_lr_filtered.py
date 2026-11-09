"""Brain long-read curve inputs for the shrunk-offset evaluation: Tealeaf filtered ranks from coverage_filter_eval.py;
comparator curves (filtered, SUPPA2 exact rows) are unchanged and copied from the leaf5 evaluation."""
import shutil
import pandas as pd
F, OLD = "analyses/coverage_null_filter_sh/brain", "analyses/coverage_null_filter_leaf5/brain"
ranks = pd.read_csv(f"{F}/lr_ranks.tsv.gz", sep="\t", low_memory=False)
ranks = ranks.loc[ranks.table.eq("filtered") & ranks["rank"].le(200) & ranks.method.eq("Tealeaf")]
ranks.assign(scope="all tested", ranking="continuous raw").to_csv(f"{F}/tealeaf_lr_rank_filtered.tsv.gz", sep="\t", index=False)
shutil.copy(f"{OLD}/sign_agreement_by_rank_filtered.tsv", f"{F}/sign_agreement_by_rank_filtered.tsv")
print(len(ranks), ranks.agree.head(100).mean() if "agree" in ranks else "")
