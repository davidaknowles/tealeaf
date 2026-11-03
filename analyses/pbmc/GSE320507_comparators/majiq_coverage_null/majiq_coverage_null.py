"""MAJIQ heterogen coverage-only null on GSE320507 full-data cell-type contrasts.

For every LSV connection and contrast, reproduce heterogen's raw Mann-Whitney test on per-sample
posterior means (x + 1/J) / (n + 1), then simulate x_s ~ Binomial(n_s, p0) with each sample's own LSV
total n_s and one pooled PSI p0 (no difference between the levels) and rerun the same test. Calls that
the null also makes are driven by coverage, not by a difference in PSI.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import rna_majiq as nm
from scipy.special import ndtr
from scipy.stats import rankdata

B = Path(sys.argv[1])
OUT = Path(sys.argv[2])
N_SIM = 3
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)

p = nm.PsiCoverage.from_zarr(B / "majiq/all.psicov")
prefixes = list(p.prefixes)
total = np.asarray(p.df["raw_total"].values, dtype=np.float64)          # (ec, prefix)
psi = np.asarray(p.df["raw_psi"].values, dtype=np.float64)
passed = np.asarray(p.df["event_passed"].values)
alpha = np.asarray(p.alpha_prior.values, dtype=np.float64)              # (ec,) = 1/J
cover = np.rint(np.nan_to_num(psi) * total)
sg = nm.SpliceGraph.from_zarr(B / "majiq/combined.splicegraph")
ecdf = p.get_events(sg.introns, sg.junctions).ec_dataframe().reset_index()
print("connections", total.shape, flush=True)


def mwu(values, in_a):
    """Two-sided asymptotic Mann-Whitney with tie and continuity correction, NaN = not passed."""
    r = rankdata(values, axis=1, nan_policy="omit")
    ok = ~np.isnan(values)
    na, nb = (ok & in_a).sum(1), (ok & ~in_a).sum(1)
    n = na + nb
    u = np.nansum(np.where(in_a, r, np.nan), axis=1) - na * (na + 1) / 2
    s = np.sort(values, axis=1)
    m, k = s.shape
    start = np.ones((m, k), bool)
    start[:, 1:] = s[:, 1:] != s[:, :-1]
    run = np.cumsum(start.ravel()) - 1
    t = np.bincount(run, weights=~np.isnan(s).ravel())
    row = np.zeros(run.max() + 1, int)
    row[run] = np.repeat(np.arange(m), k)
    ties = np.bincount(row, weights=t ** 3 - t, minlength=m)
    with np.errstate(invalid="ignore", divide="ignore"):
        var = na * nb / 12 * ((n + 1) - ties / (n * (n - 1)))
        z = (np.abs(u - na * nb / 2) - 0.5) / np.sqrt(var)
    pv = 2 * (1 - ndtr(np.maximum(z, 0)))
    return np.minimum(pv, 1.0)


def bh(pv):
    q = np.full(pv.shape, np.nan)
    f = np.isfinite(pv)
    x = pv[f]
    o = np.argsort(x)
    adj = np.minimum.accumulate((x[o] * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    q[np.flatnonzero(f)[o]] = np.minimum(adj, 1)
    return q


contrasts = json.loads((B / "contrasts.json").read_text())[: int(sys.argv[3]) if len(sys.argv) > 3 else None]
rows = []
for ci, c in enumerate(contrasts):
    ia = [prefixes.index(s) for s in c["samples_a"] if s in prefixes]
    ib = [prefixes.index(s) for s in c["samples_b"] if s in prefixes]
    cols = np.array(ia + ib)
    in_a = np.zeros(len(cols), bool)
    in_a[: len(ia)] = True
    ok = passed[:, cols]
    keep = (ok[:, in_a].sum(1) >= 0.5 * len(ia)) & (ok[:, ~in_a].sum(1) >= 0.5 * len(ib))
    idx = np.flatnonzero(keep)
    if len(idx) == 0:
        continue
    n = np.where(ok[idx], total[idx][:, cols], np.nan)
    x = np.where(ok[idx], cover[idx][:, cols], np.nan)
    a = alpha[idx][:, None]
    post = (x + a) / (n + 1)
    p_obs = mwu(post, in_a)
    pa = np.nansum(x[:, in_a], 1) / np.nansum(n[:, in_a], 1)
    pb = np.nansum(x[:, ~in_a], 1) / np.nansum(n[:, ~in_a], 1)
    p0 = np.nansum(x, 1) / np.nansum(n, 1)
    d_obs = np.nanmedian(post[:, in_a], 1) - np.nanmedian(post[:, ~in_a], 1)
    rec = dict(contrast_index=ci, contrast_id=c["contrast_id"], ec_idx=idx, p_obs=p_obs, q_obs=bh(p_obs),
               pooled_effect=pa - pb, majiq_effect=d_obs, p0=p0, J=1 / alpha[idx],
               log2_cov_ratio=np.log2(np.nanmedian(n[:, in_a], 1) / np.nanmedian(n[:, ~in_a], 1)),
               median_cov_a=np.nanmedian(n[:, in_a], 1), median_cov_b=np.nanmedian(n[:, ~in_a], 1))
    nn = np.nan_to_num(n).astype(np.int64)
    for k in range(N_SIM):
        xs = rng.binomial(nn, np.clip(np.nan_to_num(p0)[:, None], 0, 1)).astype(np.float64)
        ps = np.where(ok[idx], (xs + a) / (n + 1), np.nan)
        rec[f"p_null{k}"] = mwu(ps, in_a)
        rec[f"effect_null{k}"] = np.nanmedian(ps[:, in_a], 1) - np.nanmedian(ps[:, ~in_a], 1)
    rows.append(pd.DataFrame(rec))
    if ci == 0 and len(rows):
        # validate the reproduction against MAJIQ's own raw_pvalue
        raw = pd.read_csv(B / "majiq/tests/raw" / (c["contrast_id"].replace("/", "_") + ".tsv"), sep="\t", comment="#")
        key = ["gene_id", "event_type", "ref_exon_start", "ref_exon_end", "start", "end"]
        mine = rows[-1][["ec_idx", "p_obs"]].merge(ecdf[["ec_idx"] + key], on="ec_idx")
        mine["event_type"] = mine["event_type"].astype(str).str.strip("b'")
        raw["event_type"] = raw["event_type"].astype(str)
        j = mine.merge(raw[key + ["raw_pvalue"]], on=key)
        print("validation rows", len(j), "of", len(raw), "spearman", j[["p_obs", "raw_pvalue"]].corr("spearman").iloc[0, 1],
              "median |log10 ratio|", np.nanmedian(np.abs(np.log10(j.p_obs / j.raw_pvalue))), flush=True)
    if ci % 20 == 0:
        print(ci, c["contrast_id"], len(idx), flush=True)

t = pd.concat(rows, ignore_index=True)
for k in range(N_SIM):
    t[f"q_null{k}"] = t.groupby("contrast_index")[f"p_null{k}"].transform(lambda s: pd.Series(bh(s.to_numpy()), index=s.index))
t.to_parquet(OUT / "majiq_coverage_null.parquet")

t["eff_bin"] = pd.cut(t["pooled_effect"].abs(), [-1e-9, .02, .05, .1, .2, 1])
t["cov_bin"] = pd.cut(t["log2_cov_ratio"].abs(), [-1e-9, .5, 1, 2, 3, 100])
call = t["q_obs"] < .05
summ = []
for k in range(N_SIM):
    summ.append(dict(sim=k, null_p05=(t[f"p_null{k}"] < .05).mean(), null_q05=int((t[f"q_null{k}"] < .05).sum())))
summ = pd.DataFrame(summ)
print("tests", len(t), "observed q<.05", int(call.sum()))
print(summ.to_string())
null_call = (t["q_null0"] < .05)
same = np.sign(t["majiq_effect"]) == np.sign(t[[f"effect_null{k}" for k in range(N_SIM)]].mean(1))
g = t.assign(call=call, null_call=null_call, null_p05=t["p_null0"] < .05, same=same).groupby("eff_bin", observed=True)
by_eff = pd.DataFrame({"tests": g.size(), "observed_q05": g["call"].sum(), "null_q05": g["null_call"].sum(),
                       "null_p05_rate": g["null_p05"].mean(),
                       "calls_also_null_q05": g.apply(lambda d: (d.call & d.null_call).sum()),
                       "calls_sign_as_null": g.apply(lambda d: d.loc[d.call, "same"].mean())})
print(by_eff.to_string())
g2 = t.assign(call=call, null_call=null_call).groupby("cov_bin", observed=True)
by_cov = pd.DataFrame({"tests": g2.size(), "observed_q05_rate": g2["call"].mean(), "null_q05_rate": g2["null_call"].mean()})
print(by_cov.to_string())
small = t[call & (t["pooled_effect"].abs() <= .05)]
print("small-effect calls", len(small), "with null q<.05", int((small["q_null0"] < .05).sum()),
      "null p<.05", float((small["p_null0"] < .05).mean()))
summ.to_csv(OUT / "null_calls.tsv", sep="\t", index=False)
by_eff.to_csv(OUT / "by_pooled_effect.tsv", sep="\t")
by_cov.to_csv(OUT / "by_coverage_ratio.tsv", sep="\t")
