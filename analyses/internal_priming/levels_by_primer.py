"""Diagnostic 2: path usage levels by primer (poly(dT) only, random hexamer only, pooled) against long-read levels,
for long-read-evaluable two-path internal tests of the brain data."""
import gzip, json, sys
import numpy as np, pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.local_path_reads import block_structure, structure_paths, pooled_path_shares
from tealeaf.sc.local_path_test import prepare_block_test

C = ".cache/native_v4/collection_parity_p2_v4/all"
blocks = json.load(gzip.open(f"{C}/blocks.json.gz", "rt"))
key_of = {b["block_id"]: k for k, b in blocks.items()}
internal = {b["block_id"] for b in blocks.values() if block_structure(structure_paths(b))["structure"] == "internal"}
lr = pd.read_csv("analyses/native_align/internal_native_v4+tj_native_v4_lrchain_groupfix/lr_mapping.tsv.gz", sep="\t", low_memory=False)
lr = lr.loc[lr.eligible.astype(str).str.lower().eq("true") & lr.block_id.isin(internal) & lr.n_paths.eq(2) & lr.minimum_pooled_depth.ge(50)]
lr = lr.sample(n=min(3000, len(lr)), random_state=1)
print("tests", len(lr), flush=True)
keys = {key_of[b] for b in lr.block_id}
counts = pd.concat(c.loc[c.path_key.isin(keys)] for c in pd.read_csv(f"{C}/counts.tsv.gz", sep="\t", dtype={"subject": str}, chunksize=2_000_000))
grouped = {k: g for k, g in counts.groupby("path_key")}
opps = json.load(gzip.open(f"{C}/opportunities.json.gz", "rt"))
rows = []
for r in lr.itertuples(index=False):
    key = key_of[r.block_id]; b = blocks[key]
    levels = (r.level_a, r.level_b)
    g = grouped[key]; g = g.loc[g.cell_type.isin(levels)]
    entries = [(x.subject, x.cell_type, x.primer, int(x.mask), float(x.count)) for x in g.itertuples(index=False)]
    full = {int(m): np.asarray(v, float) for m, v in opps[key].items()}
    try:
        prep = prepare_block_test(key, entries, full, len(b["signatures"]), levels, .01, False, None)
    except ValueError:
        continue
    if len(prep["kept"]) != 2:
        continue
    o = prep["opportunities"]
    out = {"test_id": r.test_id, "block_id": r.block_id}
    for li, level in enumerate(levels):
        for primer in ("poly(dT)", "random hexamer", "both"):
            pooled = {}
            for (_, _, lev, p, mask), v in prep["lookup"].items():
                if lev == level and (primer == "both" or p == primer):
                    pooled[mask] = pooled.get(mask, 0.) + v
            try:
                shares = pooled_path_shares(pooled, o)
            except ValueError:
                shares = [np.nan, np.nan]
            out[f"{'ab'[li]}_{primer}"] = float(shares[0]) / max(float(shares[0] + shares[1]), 1e-12)
            out[f"{'ab'[li]}_{primer}_n"] = sum(pooled.values())
        vec = sum(np.asarray(json.loads(getattr(r, f"counts_{'ab'[li]}_rep{k}")), float) for k in (1, 2))
        out[f"{'ab'[li]}_lr"] = vec[0] / vec.sum() if vec.sum() > 0 else np.nan
    rows.append(out)
t = pd.DataFrame(rows)
t.to_csv(".cache/internal_priming/levels_by_primer.tsv.gz", sep="\t", index=False)
a = pd.read_csv(".cache/internal_priming/gamma_arun_tests.tsv.gz", sep="\t").groupby("block_id").aruns.first()
t["arun"] = t.block_id.map(a)
def summary(d, label):
    res = {"tests": len(d)}
    for p in ("poly(dT)", "random hexamer", "both"):
        lev = pd.concat([d[f"a_{p}"], d[f"b_{p}"]]); lrv = pd.concat([d.a_lr, d.b_lr]); ok = lev.notna() & lrv.notna()
        res[f"{p} level MAE"] = float((lev[ok] - lrv[ok]).abs().mean())
        res[f"{p} level r"] = float(np.corrcoef(lev[ok], lrv[ok])[0, 1])
        eff = d[f"b_{p}"] - d[f"a_{p}"]; effl = d.b_lr - d.a_lr; ok = eff.notna() & effl.notna()
        res[f"{p} effect r"] = float(np.corrcoef(eff[ok], effl[ok])[0, 1])
        res[f"{p} sign agree"] = float((np.sign(eff[ok]) == np.sign(effl[ok])).mean())
    print(label, pd.Series(res).round(3).to_dict(), flush=True)
summary(t, "all")
summary(t.loc[t.arun.eq(0)], "cassette, no A-run")
summary(t.loc[t.arun.ge(1)], "cassette, A-run in exon")
