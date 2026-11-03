"""MAJIQ heterogen coverage-only null.

For every LSV connection and contrast, reproduce heterogen's raw Mann-Whitney test on per-sample posterior means
(x + 1/J) / (n + 1) (x connection reads, n LSV total, J connections in the LSV), then simulate x_s ~ Binomial(n_s, p0)
with each sample's own LSV total n_s and one pooled PSI p0 for both levels (no difference) and rerun the same test.
Calls that the null also makes are driven by coverage, not by a difference in PSI. Each test is joined to MAJIQ's
own raw table, giving the normalized feature_id and MAJIQ's raw p-value.
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import rna_majiq as nm
from scipy.special import ndtr
from scipy.stats import rankdata

ap = argparse.ArgumentParser()
ap.add_argument("--build", type=Path, required=True, help="directory with all.psicov and combined.splicegraph")
ap.add_argument("--contrasts", type=Path, required=True)
ap.add_argument("--raw-dir", type=Path, required=True, help="heterogen raw tables <contrast stem>.tsv")
ap.add_argument("--min-experiments", type=float, required=True, help="per level; < 1 is a fraction of the level's samples")
ap.add_argument("--output", type=Path, required=True)
ap.add_argument("--n-sim", type=int, default=3)
ap.add_argument("--limit", type=int)
args = ap.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)
KEY = ["gene_id", "event_type", "ref_exon_start", "ref_exon_end", "start", "end"]
FEATURE = ["gene_id", "seqid", "strand", "event_type", "ref_exon_start", "ref_exon_end", "start", "end", "is_intron", "other_exon_start", "other_exon_end"]

p = nm.PsiCoverage.from_zarr(args.build / "all.psicov")
prefixes = {s: i for i, s in enumerate(p.prefixes)}
total = np.asarray(p.df["raw_total"].values, dtype=np.float64)          # (connection, sample)
psi = np.asarray(p.df["raw_psi"].values, dtype=np.float64)
passed = np.asarray(p.df["event_passed"].values)
alpha = np.asarray(p.alpha_prior.values, dtype=np.float64)              # (connection,) = 1/J
cover = np.rint(np.nan_to_num(psi) * total)
del psi
sg = nm.SpliceGraph.from_zarr(args.build / "combined.splicegraph")
ecdf = p.get_events(sg.introns, sg.junctions).ec_dataframe().reset_index()[["ec_idx"] + KEY]
ecdf["event_type"] = ecdf["event_type"].map(lambda v: v.decode() if isinstance(v, bytes) else str(v))
print("connections x samples", total.shape, flush=True)


def mwu(values, in_a):
    """Two-sided asymptotic Mann-Whitney with tie and continuity correction; NaN = sample not passed."""
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
    return np.minimum(2 * (1 - ndtr(np.maximum(z, 0))), 1.0)


def bh(pv):
    q = np.full(pv.shape, np.nan)
    f = np.isfinite(pv)
    x = pv[f]
    o = np.argsort(x)
    adj = np.minimum.accumulate((x[o] * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    q[np.flatnonzero(f)[o]] = np.minimum(adj, 1)
    return q


def need(n_level):
    return args.min_experiments * n_level if args.min_experiments < 1 else min(args.min_experiments, n_level)


contrasts = json.loads(args.contrasts.read_text())[: args.limit]
rows = []
for ci, c in enumerate(contrasts):
    raw_path = args.raw_dir / (c["contrast_id"].replace("/", "_") + ".tsv")
    if not raw_path.exists():
        continue
    ia = [prefixes[s] for s in c["samples_a"] if s in prefixes]
    ib = [prefixes[s] for s in c["samples_b"] if s in prefixes]
    cols = np.array(ia + ib)
    in_a = np.zeros(len(cols), bool)
    in_a[: len(ia)] = True
    ok = passed[:, cols]
    keep = (ok[:, in_a].sum(1) >= need(len(ia))) & (ok[:, ~in_a].sum(1) >= need(len(ib)))
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
    rec = dict(contrast_index=ci, contrast_id=c["contrast_id"], ec_idx=idx, p_obs=p_obs, q_obs=bh(p_obs),
               pooled_effect=pa - pb,
               majiq_effect=np.nanmedian(post[:, in_a], 1) - np.nanmedian(post[:, ~in_a], 1),
               log2_cov_ratio=np.log2(np.nanmedian(n[:, in_a], 1) / np.nanmedian(n[:, ~in_a], 1)),
               median_cov_a=np.nanmedian(n[:, in_a], 1), median_cov_b=np.nanmedian(n[:, ~in_a], 1))
    nn = np.nan_to_num(n).astype(np.int64)
    for k in range(args.n_sim):
        xs = rng.binomial(nn, np.clip(np.nan_to_num(p0)[:, None], 0, 1)).astype(np.float64)
        ps = np.where(ok[idx], (xs + a) / (n + 1), np.nan)
        rec[f"p_null{k}"] = mwu(ps, in_a)
        rec[f"q_null{k}"] = bh(rec[f"p_null{k}"])
        rec[f"effect_null{k}"] = np.nanmedian(ps[:, in_a], 1) - np.nanmedian(ps[:, ~in_a], 1)
    local = pd.DataFrame(rec).merge(ecdf, on="ec_idx", how="left")
    raw = pd.read_csv(raw_path, sep="\t", comment="#", low_memory=False)
    raw["event_type"] = raw["event_type"].astype(str)
    raw["feature_id"] = raw[FEATURE].astype(str).agg(":".join, axis=1)
    local = local.merge(raw[KEY + ["feature_id", "raw_pvalue"]], on=KEY, how="left")
    rows.append(local)
    if ci % 20 == 0:
        print(ci, c["contrast_id"], len(idx), "joined", int(local.feature_id.notna().sum()), "of raw", len(raw), flush=True)

t = pd.concat(rows, ignore_index=True)
t.drop(columns=KEY).to_parquet(args.output / "majiq_coverage_null.parquet")
j = t.dropna(subset=["raw_pvalue"])
null_p = [f"p_null{k}" for k in range(args.n_sim)]
null_q = [f"q_null{k}" for k in range(args.n_sim)]
validation = dict(tests=len(t), joined=len(j), spearman_p=float(j[["p_obs", "raw_pvalue"]].corr("spearman").iloc[0, 1]),
                  median_abs_log10_ratio=float(np.nanmedian(np.abs(np.log10(j.p_obs / j.raw_pvalue)))))
print(validation)
call = t["q_obs"] < .05
summ = pd.DataFrame([dict(sim=k, tests=len(t), observed_q05=int(call.sum()), null_p05=float((t[pk] < .05).mean()), null_q05=int((t[qk] < .05).sum()))
                     for k, (pk, qk) in enumerate(zip(null_p, null_q))])
print(summ.to_string())
t["eff_bin"] = pd.cut(t["pooled_effect"].abs(), [-1e-9, .02, .05, .1, .2, 1])
t["cov_bin"] = pd.cut(t["log2_cov_ratio"].abs(), [-1e-9, .5, 1, 2, 3, 100])
t["call"], t["null_call"] = call, t["q_null0"] < .05
t["same"] = np.sign(t["majiq_effect"]) == np.sign(t[[f"effect_null{k}" for k in range(args.n_sim)]].mean(1))
rows_e = []
for b, d in t.groupby("eff_bin", observed=True):
    rows_e.append(dict(pooled_effect=str(b), tests=len(d), observed_q05=int(d.call.sum()), null_q05=int(d.null_call.sum()),
                       null_p05_rate=float((d.p_null0 < .05).mean()), calls_also_null_q05=int((d.call & d.null_call).sum()),
                       calls_sign_as_null=float(d.loc[d.call, "same"].mean())))
by_eff = pd.DataFrame(rows_e)
by_cov = pd.DataFrame([dict(abs_log2_coverage_ratio=str(b), tests=len(d), observed_q05_rate=float(d.call.mean()), null_q05_rate=float(d.null_call.mean()))
                       for b, d in t.groupby("cov_bin", observed=True)])
print(by_eff.to_string())
print(by_cov.to_string())
pd.DataFrame([validation]).to_csv(args.output / "validation.tsv", sep="\t", index=False)
summ.to_csv(args.output / "null_calls.tsv", sep="\t", index=False)
by_eff.to_csv(args.output / "by_pooled_effect.tsv", sep="\t", index=False)
by_cov.to_csv(args.output / "by_coverage_ratio.tsv", sep="\t", index=False)
