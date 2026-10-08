"""Coverage and direction diagnostics, distinct from discovery testing."""

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


def ranked_direction_summary(table, cutoffs=(40, 100, 200)):
    """Summarize discrete cumulative agreement without extrapolating short curves.

    The normalized area at K is the arithmetic mean of the cumulative
    agreement at integer ranks 1 through K, not a ROC area.
    """
    rows = []
    for method, local in table.groupby("method", observed=True):
        local = local.sort_values("rank")
        ranks = local["rank"].to_numpy(dtype=float)
        if not np.array_equal(ranks, np.arange(1, len(local) + 1)):
            raise ValueError("agreement ranks must be unique and consecutive from one")
        if local.pooled_replicated.isna().any():
            raise ValueError("ranked directions must all be evaluable")
        direction = local.pooled_replicated.astype(str).str.lower()
        if not direction.isin(["true", "false"]).all():
            raise ValueError("ranked directions must be boolean agreement indicators")
        agreement = direction.eq("true").to_numpy(dtype=float)
        curve = np.cumsum(agreement) / ranks
        for cutoff in cutoffs:
            if int(cutoff) != cutoff or cutoff <= 0:
                raise ValueError("rank cutoffs must be positive integers")
            available = len(local) >= cutoff
            rows.append({"method": str(method), "cutoff": int(cutoff), "n_available": len(local), "n_agree": int(agreement[:cutoff].sum()) if available else np.nan, "agreement": float(curve[cutoff - 1]) if available else np.nan, "normalized_auc": float(curve[:cutoff].mean()) if available else np.nan})
    return rows


def ranked_category_summary(table, category="event_type", cutoffs=(100, 200)):
    """Describe category composition of fixed global rank prefixes.

    Categories are not reranked or selected for favorable agreement. Missing
    labels remain an explicit category, and external-zero nonagreements stay
    in every denominator. This is descriptive, not a standardized AUC.
    """
    if category not in table:
        raise ValueError("rank category column is missing")
    overall = ranked_direction_summary(table, cutoffs)
    rows = []
    for summary in overall:
        method, cutoff = summary["method"], summary["cutoff"]
        local = table.loc[table.method.astype(str).eq(method) & table["rank"].le(cutoff)].copy()
        local[category] = local[category].fillna("unclassified").astype(str)
        for label, group in local.groupby(category, observed=True):
            n_agree = int(group.pooled_replicated.astype(str).str.lower().eq("true").sum())
            rows.append(dict(method=method, cutoff=cutoff, category=label, n_selected=len(local), complete_prefix=len(local) == cutoff, n_tests=len(group), n_agree=n_agree, fraction_of_prefix=len(group) / len(local), agreement=n_agree / len(group)))
    return rows


def reexpress_event_directions(mapped, effects, effect_column, method):
    """Change direction estimates on an EXACT fixed LR-evaluable event family.

    Mapping, depths, p-values, ranks and external zeros are not reselected.
    Missing, nonfinite or zero replacement directions count as nonagreements.
    Input rows require finite nonzero old short-read effects so replicate
    effects can be recovered from the stored dot products. This is a paired
    direction ablation, not remapping or a PSI estimator from ILR scores.
    """
    keys = ["feature_id", "contrast_id"]
    required = keys + ["short_read_effect", "long_read_effect", "replicate_1_dot_product", "replicate_2_dot_product"]
    if any(column not in mapped for column in required) or effect_column not in effects or any(column not in effects for column in keys):
        raise ValueError("event direction replacement requires aligned effects and mapped endpoints")
    old = mapped.short_read_effect.to_numpy(float)
    if not np.isfinite(old).all() or (old == 0).any() or mapped.duplicated(keys).any() or effects.duplicated(keys).any():
        raise ValueError("fixed mapped event family requires unique identities and finite nonzero original directions")
    lookup = effects[keys + [effect_column]].rename(columns={effect_column: "replacement_direction"})
    result = mapped.merge(lookup, on=keys, how="left", validate="one_to_one", indicator=True, sort=False)
    if not result._merge.eq("both").all() or len(result) != len(mapped):
        raise ValueError("every fixed mapped event needs its requested replacement estimate")
    replacement = result.replacement_direction.to_numpy(float)
    valid = np.isfinite(replacement) & (replacement != 0)
    replacement = np.where(valid, replacement, np.nan)
    result["short_read_effect"] = replacement
    result["direction_available"] = valid
    result["method"] = method
    # Compare signs rather than multiplying possibly large Gaussian units.
    pooled = result.long_read_effect.to_numpy(float)
    result["pooled_replicated"] = valid & np.isfinite(pooled) & (np.sign(replacement) == np.sign(pooled))
    replicate_agreement = []
    for replicate in (1, 2):
        column = f"replicate_{replicate}_dot_product"
        external = result[column].to_numpy(float) / old
        result[f"replicate_{replicate}_long_read_effect"] = external
        result[column] = external * replacement
        replicate_agreement.append(valid & np.isfinite(external) & (np.sign(replacement) == np.sign(external)))
    result["both_replicates_replicated"] = replicate_agreement[0] & replicate_agreement[1]
    return result.drop(columns=["replacement_direction", "_merge"])


def complete_paired_fits(table):
    """Require every eligible pair to fit, for one row per subject/type.

    N_samples equals twice the number of pairs before optimization. This
    invariant must be checked against input metadata by the calling analysis.
    Reporting completeness is separate from inferential eligibility.
    """
    required = ["n_samples", "n_subjects", "converged"]
    if any(column not in table for column in required):
        raise ValueError("paired completeness requires sample and subject counts")
    samples = pd.to_numeric(table.n_samples, errors="coerce")
    subjects = pd.to_numeric(table.n_subjects, errors="coerce")
    complete = table.converged.astype(str).str.lower().eq("true") & samples.ge(8) & samples.mod(2).eq(0) & samples.eq(2 * subjects)
    report_subjects = pd.to_numeric(table.get("report_n_subjects", table.n_subjects), errors="coerce")
    reporting = complete & samples.eq(2 * report_subjects)
    return complete, reporting


def complete_cluster_fit(result, expected_subjects):
    """Distinguish failed fits from structurally uninformative fitted clusters."""
    fitted = result.get("n_fitted_subjects", result["n_subjects"])
    return bool(result["converged"] and expected_subjects >= 4 and fitted == expected_subjects and result["n_subjects"] >= 4)


def coverage_correlation(pvalues, coverage, controls=None):
    """Correlate p (not -log p) with depth, optionally residualizing ranks.

    Controls are N by C, with one row per hypothesis. Negative rho means
    smaller p-values at greater coverage; it alone does not imply null bias.
    """
    pvalues, coverage = np.asarray(pvalues, float), np.asarray(coverage, float)
    valid = np.isfinite(pvalues) & (pvalues >= 0) & (pvalues <= 1) & np.isfinite(coverage) & (coverage > 0)
    if controls is not None:
        controls = np.asarray(controls, float)
        if controls.ndim == 1:
            controls = controls[:, None]
        valid &= np.isfinite(controls).all(axis=1)
    p, depth = pvalues[valid], coverage[valid]
    rho = float(spearmanr(p, depth).statistic) if len(p) > 2 and np.ptp(p) > 0 and np.ptp(depth) > 0 else np.nan
    partial = np.nan
    if controls is not None and len(p) > controls.shape[1] + 2:
        x = np.column_stack([np.ones(len(p)), *[rankdata(column) for column in controls[valid].T]])
        residuals = []
        for values in (p, depth):
            ranks = rankdata(values)
            residuals.append(ranks - x @ np.linalg.lstsq(x, ranks, rcond=None)[0])
        if all(np.linalg.norm(values) > 1e-10 for values in residuals):
            partial = float(np.corrcoef(residuals)[0, 1])
    return {"n": len(p), "rho_p_coverage": rho, "partial_rho": partial}


def aligned_direction(first, second, first_features=None, second_features=None):
    """Compare effects in the same feature order, rejecting incompatible sets.

    Scalar effects use sign concordance; multivariate effects use positive
    inner product, plus cosine and component concordance. Zero components
    are excluded rather than counted as agreement.
    """
    first, second = np.asarray(first, float), np.asarray(second, float)
    if first.size == 0 or second.size == 0:
        raise ValueError("missing effect direction")
    if first_features is not None:
        if len(set(first_features)) != len(first_features) or len(set(second_features)) != len(second_features):
            raise ValueError("duplicate feature identities")
        if set(first_features) != set(second_features):
            raise ValueError("effect feature sets differ")
        second = second[[second_features.index(feature) for feature in first_features]]
    if first.shape != second.shape or first.ndim != 1:
        raise ValueError("effect vectors must have the same dimension")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("nonfinite effect vector")
    norm = np.linalg.norm(first) * np.linalg.norm(second)
    eligible = (first != 0) & (second != 0)
    dot = float(first @ second)
    return {"direction_agrees": bool(dot > 0) if norm > 0 else np.nan,
            "cosine": dot / norm if norm > 0 else np.nan,
            "nonzero_components": int(eligible.sum()),
            "agreeing_components": int(np.sum(first[eligible] * second[eligible] > 0))}
