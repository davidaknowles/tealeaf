"""Coverage-only null for MAJIQ Heterogen.

Heterogen's raw test compares, between two groups of samples, the per-sample posterior means
psi_s = (x_s + 1/J) / (n_s + 1) of an LSV connection (x_s connection reads, n_s LSV total, J connections of the LSV)
with a Mann-Whitney test. The prior pulls a low-coverage sample toward 1/J, so two groups with the same PSI but
different expression differ in the tested quantity. The coverage-only null keeps every sample's n_s, draws
x_s ~ Binomial(n_s, p0) with one PSI p0 pooled over both groups, and reruns the test. A MAJIQ test whose null
p-value is small is driven by coverage and is removed before evaluation (drop_coverage_null_calls).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import rankdata

from tealeaf.sc.junction_benchmark import benjamini_hochberg


def mann_whitney_rows(values: np.ndarray, in_a: np.ndarray) -> np.ndarray:
    """Two-sided asymptotic Mann-Whitney p-value per row with tie and continuity correction.

    values (tests x samples), NaN for samples that do not enter a test; in_a (samples,) bool, group membership.
    """
    values = np.asarray(values, dtype=float)
    ranks = rankdata(values, axis=1, nan_policy="omit")
    present = ~np.isnan(values)
    n_a, n_b = (present & in_a).sum(1), (present & ~in_a).sum(1)
    n = n_a + n_b
    u = np.nansum(np.where(in_a, ranks, np.nan), axis=1) - n_a * (n_a + 1) / 2
    ordered = np.sort(values, axis=1)
    rows, cols = ordered.shape
    start = np.ones((rows, cols), bool)
    start[:, 1:] = ordered[:, 1:] != ordered[:, :-1]
    run = np.cumsum(start.ravel()) - 1
    sizes = np.bincount(run, weights=~np.isnan(ordered).ravel())
    run_row = np.zeros(run.max() + 1 if run.size else 0, int)
    run_row[run] = np.repeat(np.arange(rows), cols)
    ties = np.bincount(run_row, weights=sizes ** 3 - sizes, minlength=rows)
    with np.errstate(invalid="ignore", divide="ignore"):
        variance = n_a * n_b / 12 * ((n + 1) - ties / (n * (n - 1)))
        z = (np.abs(u - n_a * n_b / 2) - 0.5) / np.sqrt(variance)
    return np.minimum(2 * (1 - ndtr(np.maximum(z, 0))), 1.0)


def bh_finite(pvalues: np.ndarray) -> np.ndarray:
    """BH over the finite entries; NaN elsewhere."""
    pvalues = np.asarray(pvalues, dtype=float)
    q = np.full(pvalues.shape, np.nan)
    finite = np.isfinite(pvalues)
    if finite.any():
        q[finite] = benjamini_hochberg(pvalues[finite])
    return q


def coverage_null_contrast(total: np.ndarray, cover: np.ndarray, passed: np.ndarray, alpha: np.ndarray, in_a: np.ndarray,
                           min_experiments: float, rng: np.random.Generator, n_sim: int = 3) -> pd.DataFrame:
    """Observed and coverage-only null Heterogen tests of one contrast.

    total, cover, passed (connections x samples of the contrast): LSV totals, connection reads, MAJIQ's per-sample
    quantifiability; alpha (connections,) = 1/J; in_a (samples,) group membership. min_experiments per group, a
    fraction of the group's samples when < 1 (MAJIQ's --min-experiments). Returns one row per tested connection
    (ec_idx = row index into the inputs) with the reproduced p-value, pooled and posterior-mean effects, the coverage
    ratio between the groups and n_sim null p-values.
    """
    def need(size):
        return min_experiments * size if min_experiments < 1 else min(min_experiments, size)

    keep = (passed[:, in_a].sum(1) >= need(in_a.sum())) & (passed[:, ~in_a].sum(1) >= need((~in_a).sum()))
    index = np.flatnonzero(keep)
    if index.size == 0:
        return pd.DataFrame()
    ok = passed[index]
    n = np.where(ok, total[index], np.nan)
    x = np.where(ok, cover[index], np.nan)
    prior = alpha[index][:, None]
    posterior = (x + prior) / (n + 1)
    p_observed = mann_whitney_rows(posterior, in_a)
    pooled_a = np.nansum(x[:, in_a], 1) / np.nansum(n[:, in_a], 1)
    pooled_b = np.nansum(x[:, ~in_a], 1) / np.nansum(n[:, ~in_a], 1)
    p0 = np.clip(np.nan_to_num(np.nansum(x, 1) / np.nansum(n, 1)), 0, 1)
    result = {"ec_idx": index, "p_obs": p_observed, "q_obs": bh_finite(p_observed), "pooled_effect": pooled_a - pooled_b,
              "majiq_effect": np.nanmedian(posterior[:, in_a], 1) - np.nanmedian(posterior[:, ~in_a], 1),
              "log2_cov_ratio": np.log2(np.nanmedian(n[:, in_a], 1) / np.nanmedian(n[:, ~in_a], 1)),
              "median_cov_a": np.nanmedian(n[:, in_a], 1), "median_cov_b": np.nanmedian(n[:, ~in_a], 1)}
    trials = np.nan_to_num(n).astype(np.int64)
    for draw in range(n_sim):
        simulated = rng.binomial(trials, p0[:, None]).astype(float)
        null_posterior = np.where(ok, (simulated + prior) / (n + 1), np.nan)
        result[f"p_null{draw}"] = mann_whitney_rows(null_posterior, in_a)
        result[f"q_null{draw}"] = bh_finite(result[f"p_null{draw}"])
        result[f"effect_null{draw}"] = np.nanmedian(null_posterior[:, in_a], 1) - np.nanmedian(null_posterior[:, ~in_a], 1)
    return pd.DataFrame(result)


def null_pvalue(table: pd.DataFrame) -> pd.Series:
    """Median coverage-only null p-value over the simulated draws (p_null0, p_null1, ...)."""
    return table[[column for column in table if column.startswith("p_null")]].median(axis=1)


def drop_coverage_null_calls(tests: pd.DataFrame, null: pd.DataFrame, keys=("contrast_id", "feature_id"), level: float = .05) -> pd.DataFrame:
    """MAJIQ tests whose coverage-only null p-value (median over draws) is at least `level`.

    tests and null share `keys`; tests without a null row are kept (they were not reproduced, so not shown to be
    coverage-driven).
    """
    keys = list(keys)
    flagged = null.assign(p_null=null_pvalue(null)).loc[lambda frame: frame.p_null < level, keys].drop_duplicates()
    merged = tests.merge(flagged.assign(_coverage_null=True), on=keys, how="left")
    return merged.loc[merged._coverage_null.isna()].drop(columns="_coverage_null").reset_index(drop=True)


def lsv_key(contrast_id: pd.Series, feature_id: pd.Series) -> pd.Series:
    """LSV identifier <contrast>|<gene>|<event type>|<reference exon start>|<end> of normalized MAJIQ edge features
    (gene:seqid:strand:event_type:ref_exon_start:ref_exon_end:...)."""
    parts = feature_id.str.split(":", expand=True)
    return contrast_id + "|" + parts[0] + "|" + parts[3] + "|" + parts[4] + "|" + parts[5]


def lsv_pvalues_without_coverage_calls(null: pd.DataFrame, level: float = .05) -> pd.DataFrame:
    """Per LSV, MAJIQ's smallest raw edge p-value over the edges the coverage-only null does not call (null p >= level).

    LSVs whose every edge is called are absent. Columns lsv, p_value, edges, kept_edges.
    """
    table = null.dropna(subset=["feature_id"]).assign(lsv=lambda frame: lsv_key(frame.contrast_id, frame.feature_id), p_null=null_pvalue(null))
    table["kept"] = table.p_null >= level
    summary = table.groupby("lsv").agg(edges=("kept", "size"), kept_edges=("kept", "sum"))
    pvalues = table.loc[table.kept].groupby("lsv").raw_pvalue.min().rename("p_value")
    return summary.join(pvalues, how="inner").reset_index()
