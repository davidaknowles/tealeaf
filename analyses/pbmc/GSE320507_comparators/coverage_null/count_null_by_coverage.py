"""Tealeaf count-null rejection by coverage ratio between the two cell types (total block molecules per level)."""
import gzip, json, sys
import numpy as np, pandas as pd
rows = []
for structure, coll in (("internal", "all"), ("terminal", "junction")):
    t = pd.read_csv(f".cache/pbmc/GSE320507/count_null/{structure}/trials.tsv.gz", sep="\t")
    t = t[t.model.eq("anchors")]
    blocks = json.load(gzip.open(f".cache/pbmc/GSE320507/collection/{coll}/blocks.json.gz", "rt"))
    key_of = {b["block_id"]: k for k, b in blocks.items()}
    parts = t.test_id.str.split("|", expand=True)
    t["block_id"], t["a"], t["b"] = parts[0], parts[2], parts[3]
    keys = {key_of[b] for b in t.block_id}
    c = pd.concat(ch[ch.path_key.isin(keys)] for ch in pd.read_csv(f".cache/pbmc/GSE320507/collection/{coll}/counts.tsv.gz", sep="\t", chunksize=2_000_000))
    tot = c.groupby(["path_key", "cell_type"])["count"].sum()
    t["cov_a"] = [tot.get((key_of[b], a), 0) for b, a in zip(t.block_id, t.a)]
    t["cov_b"] = [tot.get((key_of[b], x), 0) for b, x in zip(t.block_id, t.b)]
    t["ratio"] = np.abs(np.log2((t.cov_a + 1) / (t.cov_b + 1)))
    t["bin"] = pd.cut(t.ratio, [-1, 1, 2, 3, 100])
    for (sc, b), d in t.groupby(["scenario", "bin"], observed=True):
        cv = d[d.converged.astype(str).eq("True")]
        rows.append(dict(structure=structure, scenario=sc, abs_log2_coverage_ratio=str(b), trials=len(d), converged=len(cv),
                         reject05_all=(d.p_value < .05).mean(), reject05_converged=(cv.p_value < .05).mean() if len(cv) else np.nan))
r = pd.DataFrame(rows)
print(r.to_string(index=False))
r.to_csv(sys.argv[1], sep="\t", index=False)
