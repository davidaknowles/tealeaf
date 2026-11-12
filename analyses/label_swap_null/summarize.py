"""Rejection rates and BH discoveries of LeafCutter and Tealeaf on the 21 pairs, with the real labels and with the
within-subject label swaps (same swaps for both methods). LeafCutter: clusters; Tealeaf: blocks (internal, terminal),
non-converged fits at p = 1. BH within method and table over the 21 pairs."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.junction_benchmark import benjamini_hochberg
O = Path(".cache/label_swap"); PB = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
contrasts = json.loads((O / "contrasts.json").read_text())
full = [c["contrast_id"] for c in json.loads((PB / "contrasts.json").read_text())]
pairs = {tuple(sorted((c["level_a"], c["level_b"]))) for c in contrasts}
tables = []
real = pd.concat(pd.read_csv(PB / f"leafcutter/tests_cov5/normalized/{full.index(c['contrast_id'])}.tsv", sep="\t") for c in contrasts)
swap = pd.concat(pd.read_csv(O / f"leafcutter/normalized/{i}.tsv", sep="\t") for i in range(len(contrasts)))
for label, t in (("real", real), ("swapped", swap)):
    tables.append(pd.DataFrame({"method": "LeafCutter", "labels": label, "pair": t.level_a + "||" + t.level_b, "p": t.p_value.fillna(1.).to_numpy()}))
for s in ("internal", "terminal"):
    for label, path in (("real", f".cache/binary_grid/pbmc_GSE320507_{s}_sh_full/merged/paired_path.tsv"), ("swapped", O / f"tealeaf_{s}/merged/paired_path.tsv")):
        t = pd.read_csv(path, sep="\t", usecols=["level_a", "level_b", "p_value", "converged"])
        t = t[[tuple(sorted(x)) in pairs for x in zip(t.level_a, t.level_b)]]
        p = np.where(t.converged.astype(str).str.lower().eq("true"), t.p_value.astype(float), 1.)
        tables.append(pd.DataFrame({"method": f"Tealeaf {s}", "labels": label, "pair": (t.level_a + "||" + t.level_b).to_numpy(), "p": p}))
d = pd.concat(tables, ignore_index=True)
d["method_all"] = d.method.str.replace(" internal| terminal", "", regex=True)
rows = []
for keys, g in list(d.groupby(["method", "labels"])) + [((m + " (both)", l), g) for (m, l), g in d[d.method.str.startswith("Tealeaf")].groupby(["method_all", "labels"])]:
    p = g.p.to_numpy()
    rows.append(dict(method=keys[0], labels=keys[1], tests=len(p), p05=(p < .05).mean(), p01=(p < .01).mean(), p001=(p < .001).mean(), bh=int((benjamini_hochberg(p) <= .05).sum())))
r = pd.DataFrame(rows); r.to_csv("analyses/label_swap_null/summary.tsv", sep="\t", index=False)
print(r.to_string(index=False, float_format="%.4f"))
per = d[d.labels.eq("swapped")].groupby(["method_all", "pair"]).p.apply(lambda p: (p < .01).mean()).unstack(0)
per.to_csv("analyses/label_swap_null/swapped_by_pair.tsv", sep="\t")
print(per.round(4).to_string())
