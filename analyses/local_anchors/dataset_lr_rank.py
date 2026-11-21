"""Long-read direction ranking (top 200 evaluable tests by p) of a Tealeaf long-read mapping, as written by
assess_longread_molecules.py, for a combined mapping. Usage: dataset_lr_rank.py <lr dir with lr_mapping.tsv.gz>"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.replication_audit import rank_direction_table, ranked_direction_summary

d = Path(sys.argv[1])
m = pd.read_csv(d / "lr_mapping.tsv.gz", sep="\t", low_memory=False)
m = m.loc[m.eligible.astype(str).str.lower().eq("true")].assign(method="Tealeaf", feature_id=lambda t: t.test_id, contrast_id=lambda t: t.test_id)
m["pooled_replicated"] = m.pooled_replicated.astype(str).str.lower().eq("true")
ranked = rank_direction_table(m, 200)
ranked.assign(scope="all tested", ranking="continuous raw").to_csv(d / "lr_rank.tsv.gz", sep="\t", index=False)
s = pd.DataFrame(ranked_direction_summary(ranked)); s["n_evaluable"] = len(m)
s.to_csv(d / "lr_rank_summary.tsv", sep="\t", index=False)
print(d, s[["cutoff", "n_agree", "normalized_auc", "n_evaluable"]].to_string(index=False))
