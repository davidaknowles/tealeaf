"""Do called MXE events behave like a 5'/3' junction-coverage ratio? Per contrast, the fraction of calls whose
effect points toward the 3' form (exon 2 form on +, exon 1 form on -), by whether the two exons share a transcript."""
import sys
import numpy as np, pandas as pd
c = pd.read_csv(sys.argv[1], sep="\t")
c["pseudo"] = c.tx_both > 0
# effect = PSI(exon1 form) a - b; exon1 is 5' on + strand. toward3 = level a favours the 3' (downstream) form
c["a_toward_3prime"] = np.where(c.strand.eq("+"), c.effect < 0, c.effect > 0)
for k, d in c.groupby("fold"):
    print("fold", k, "calls", len(d), "pseudo-MXE (exons share a transcript)", round(d.pseudo.mean(), 3),
          "pseudo with MXE-specific junctions (exon1-down, up-exon2) unannotated", round((d.pseudo & (d.annotated_junctions <= 2)).mean(), 3))
    per = d.groupby(["contrast_id", "pseudo"]).a_toward_3prime.agg(["size", "mean"]).reset_index()
    per = per[per["size"] >= 20]
    per["consistency"] = (per["mean"] - .5).abs() * 2
    print(per.groupby("pseudo").consistency.describe().to_string())
    print(per.sort_values("size", ascending=False).head(8).to_string(index=False))
