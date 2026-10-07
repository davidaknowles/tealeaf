"""Separate event-universe, ranking, direction and EC support differences.

Diagnostic only. No LR agreement defines an inferential family, read filter,
hyperparameter or production choice. Fresh complete mapping is supplied by
the existing event LR mapper, not inferred from a method's sign agreement.
"""

import argparse
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.run_suppa2_tealeaf_hybrid import canonical
from extra_scripts.plot_tilgner_method_replication import _rank_table
from tealeaf.sc.replication_audit import ranked_direction_summary


KEYS = ["contrast_id", "feature_id"]


def ranked_summary(table, ranking, direction, scope):
    """Same eligibility and event universe for each crossed ranking/effect."""
    ranked = rank_events(table, ranking)
    ranked["method"] = f"{ranking}, {direction}"
    ranked["pooled_replicated"] = ranked[direction] * ranked.long_read_effect > 0
    return [{**row, "scope": scope, "ranking": ranking, "direction": direction} for row in ranked_direction_summary(ranked)]


def rank_events(table, ranking):
    """Reuse the published curve's calibrated/raw/statistic/feature ordering."""
    ranked = table.copy()
    ranked["p_value"] = ranked[ranking]
    ranked["raw_p_value"] = ranked.get(f"{ranking}_raw", ranked[ranking])
    ranked["statistic"] = ranked.get(f"{ranking}_statistic", np.nan)
    ranked["method"] = ranking
    # The plotting helper also computes a curve. Its agreement is replaced
    # after ranking, so a placeholder must never become a diagnostic outcome.
    ranked["pooled_replicated"] = False
    return _rank_table(ranked, len(ranked))


def cross_ranking_effects(native, hybrid, *, exclude_zero_lr=False):
    """Validate common external truth before attributing ranking/direction."""
    if native.duplicated(KEYS).any() or hybrid.duplicated(KEYS).any():
        raise ValueError("event-contrast mapping must be unique per method")
    native = native.rename(columns={"p_value": "native_p", "raw_p_value": "native_p_raw", "statistic": "native_p_statistic", "short_read_effect": "native_effect"})
    hybrid = hybrid.rename(columns={"p_value": "hybrid_p", "raw_p_value": "hybrid_p_raw", "statistic": "hybrid_p_statistic", "short_read_effect": "hybrid_effect"})
    hybrid_columns = KEYS + ["hybrid_p", "hybrid_p_raw", "hybrid_effect", "long_read_effect"]
    if "hybrid_p_statistic" in hybrid:
        hybrid_columns.append("hybrid_p_statistic")
    shared = native.merge(hybrid[hybrid_columns], on=KEYS, validate="one_to_one", suffixes=("", "_hybrid"))
    if not np.allclose(shared.long_read_effect, shared.long_read_effect_hybrid, rtol=1e-10, atol=1e-12, equal_nan=True):
        raise ValueError("methods do not have the same LR event definition/truth")
    shared = shared.drop(columns="long_read_effect_hybrid")
    # A common nonzero-direction subset is explicit, never treated as the
    # full native universe. Whole input method universes are reported too.
    native = native.loc[np.isfinite(native.native_effect) & native.native_effect.ne(0) & np.isfinite(native.long_read_effect)].copy()
    hybrid = hybrid.loc[np.isfinite(hybrid.hybrid_effect) & hybrid.hybrid_effect.ne(0) & np.isfinite(hybrid.long_read_effect)].copy()
    shared = shared.loc[np.isfinite(shared.native_effect) & shared.native_effect.ne(0) & np.isfinite(shared.hybrid_effect) & shared.hybrid_effect.ne(0) & np.isfinite(shared.long_read_effect)].copy()
    if exclude_zero_lr:
        native, hybrid, shared = [table.loc[table.long_read_effect.ne(0)].copy() for table in (native, hybrid, shared)]
    rows = ranked_summary(native, "native_p", "native_effect", "native own universe")
    rows += ranked_summary(hybrid, "hybrid_p", "hybrid_effect", "hybrid own universe")
    for ranking in ("native_p", "hybrid_p"):
        for direction in ("native_effect", "hybrid_effect"):
            rows += ranked_summary(shared, ranking, direction, "shared event-contrast universe")
    return pd.DataFrame(rows), shared, native, hybrid


def gene_support(genes, gene_transcripts, gene_ecs, designs, features):
    """Only transcript presence/support, never expression or LR screening."""
    output = {}
    for gene, name in enumerate(genes):
        indices = np.asarray(gene_transcripts[gene], dtype=int)
        ecs = np.asarray(gene_ecs[gene], dtype=int)
        supported = np.zeros(len(indices), dtype=bool)
        for mapping in designs:
            supported |= np.asarray(mapping[ecs][:, indices].sum(axis=0)).ravel() > 0
        output[canonical(name)] = {"all_transcripts": {canonical(features[index]) for index in indices}, "supported_transcripts": {canonical(features[index]) for index in indices[supported]}, "n_ecs": len(ecs), "n_transcripts": len(indices), "n_supported_transcripts": int(supported.sum())}
    return output


def transcript_evidence(designs, features, expression, expression_features):
    """Audit global EC support and source expression without using LR outcomes."""
    global_support = np.zeros(len(features), dtype=bool)
    for design in designs:
        global_support |= np.asarray(design.sum(axis=0)).ravel() > 0
    represented = {canonical(value) for value in features}
    supported = {canonical(value) for value, keep in zip(features, global_support) if keep}
    mass = {}
    for feature, value in zip(expression_features, np.asarray(expression.sum(axis=0)).ravel()):
        name = canonical(feature)
        mass[name] = mass.get(name, 0.) + float(value)
    return represented, supported, mass


def event_evidence(included, excluded, local_supported, represented, supported, mass):
    """Distinguish annotation-only absence from discarded observed support."""
    members = included | excluded
    missing = members - local_supported
    total_mass = sum(mass.get(value, 0.) for value in members)
    missing_mass = sum(mass.get(value, 0.) for value in missing)
    result = {"n_missing": len(missing), "n_missing_not_in_features": len(missing - represented), "n_missing_global_zero_support": len((missing & represented) - supported), "n_missing_global_supported": len(missing & supported), "n_missing_positive_expression": sum(mass.get(value, 0.) > 0 for value in missing), "missing_expression_fraction": missing_mass / total_mass if total_mass > 0 else np.nan}
    for label, group in (("included", included), ("excluded", excluded)):
        total = sum(mass.get(value, 0.) for value in group)
        result[f"{label}_missing_expression_fraction"] = sum(mass.get(value, 0.) for value in group - local_supported) / total if total > 0 else np.nan
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--native-tests", type=Path, required=True)
    parser.add_argument("--hybrid-tests", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--data-cache", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expression-prefix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    native_tests = pd.read_csv(args.native_tests, sep="\t")
    hybrid_tests = pd.read_csv(args.hybrid_tests, sep="\t")
    methods = [native_tests.method.iloc[0], hybrid_tests.method.iloc[0]]
    mapped = pd.read_csv(args.mapping, sep="\t")
    valid = mapped.mapping_complete.astype(str).str.lower().eq("true") & mapped.minimum_pooled_depth.ge(20)
    mapped = mapped.loc[valid].copy()
    summary, shared, native, hybrid = cross_ranking_effects(mapped.loc[mapped.method.eq(methods[0])], mapped.loc[mapped.method.eq(methods[1])])
    summary["LR_zero_policy"] = "legacy event endpoint, exact zero counts as nonagreement"
    sensitivity, _, _, _ = cross_ranking_effects(mapped.loc[mapped.method.eq(methods[0])], mapped.loc[mapped.method.eq(methods[1])], exclude_zero_lr=True)
    sensitivity["LR_zero_policy"] = "diagnostic, exact zero LR effect excluded"
    summary = pd.concat([summary, sensitivity], ignore_index=True)
    summary.to_csv(args.output_dir / "ranking_direction_decomposition.tsv", sep="\t", index=False, na_rep="NA")
    shared.to_csv(args.output_dir / "shared_events.tsv.gz", sep="\t", index=False)
    catalog = pd.read_csv(args.catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    with args.data_cache.open("rb") as handle:
        groups, counts, genes, gene_transcripts, gene_ecs, designs = pickle.load(handle)
    features = args.features.read_text().splitlines()
    support = gene_support(genes, gene_transcripts, gene_ecs, designs, features)
    prefix = str(args.expression_prefix)
    expression = sparse.load_npz(prefix + "pseudo_spliced_TPM.npz").tocsr()
    expression_features = Path(prefix + "pseudo_spliced_cols.txt").read_text().splitlines()
    if expression.shape[1] != len(expression_features):
        raise ValueError("source expression feature dimensions disagree")
    represented, global_supported, expression_mass = transcript_evidence(designs, features, expression, expression_features)
    # Freshly mapped source-native top events, not just the older top-500 table.
    ranked_native = rank_events(native, "native_p")
    ranked_native["native_rank"] = np.arange(1, len(ranked_native) + 1)
    tested_hybrid = set(map(tuple, hybrid_tests[KEYS].to_numpy()))
    eligible_hybrid = set(map(tuple, hybrid[KEYS].to_numpy()))
    provenance = []
    for row in ranked_native.itertuples(index=False):
        event = catalog.loc[row.feature_id]
        included = {canonical(value) for value in event.included.split(",") if value}
        excluded = {canonical(value) for value in event.excluded.split(",") if value}
        local = support.get(canonical(event.gene_id), {"all_transcripts": set(), "supported_transcripts": set(), "n_ecs": 0, "n_transcripts": 0, "n_supported_transcripts": 0})
        missing = sorted((included | excluded) - local["supported_transcripts"])
        if (row.contrast_id, row.feature_id) in tested_hybrid:
            reason = "tested"
        elif not local["n_transcripts"]:
            reason = "gene absent from EC preparation"
        elif local["n_ecs"] > 128:
            reason = "gene exceeds hybrid EC limit"
        elif missing:
            reason = "incomplete supported event transcript sets"
        else:
            reason = "coverage/context filter or failed event fit"
        provenance.append({"contrast_id": row.contrast_id, "feature_id": row.feature_id, "event_type": row.event_type, "gene_id": row.gene_id, "native_rank": row.native_rank, "native_p": row.native_p, "native_effect": row.native_effect, "long_read_effect": row.long_read_effect, "native_agrees": row.native_effect * row.long_read_effect > 0, "hybrid_tested": (row.contrast_id, row.feature_id) in tested_hybrid, "hybrid_eligible": (row.contrast_id, row.feature_id) in eligible_hybrid, "reason": reason, "missing_transcripts": ",".join(missing), **{key: local[key] for key in ("n_ecs", "n_transcripts", "n_supported_transcripts")}, **event_evidence(included, excluded, local["supported_transcripts"], represented, global_supported, expression_mass)})
    provenance = pd.DataFrame(provenance)
    provenance.to_csv(args.output_dir / "native_event_support.tsv.gz", sep="\t", index=False)
    pd.concat([provenance.head(top).assign(top=top).groupby(["top", "reason"]).agg(n=("feature_id", "size"), native_agreement=("native_agrees", "mean"), hybrid_eligible=("hybrid_eligible", "sum")).reset_index() for top in (40, 100, 200)]).to_csv(args.output_dir / "top_native_support.tsv", sep="\t", index=False)
    # Do not infer units merely from a legacy filename containing TPM.
    raw_path = Path(prefix + "pseudo_spliced_count.npz")
    raw = sparse.load_npz(raw_path).tocsr() if raw_path.exists() else None
    expression_audit = {"prefix": prefix, "expression_shape": expression.shape, "expression_nnz": expression.nnz, "expression_sum": float(expression.sum()), "EC_pseudobulk_groups": len(groups), "EC_primers": [{"shape": values.shape, "nnz": values.nnz, "sum": float(values.sum())} for values in counts]}
    if raw is not None:
        expression_audit.update(count_shape=raw.shape, count_nnz=raw.nnz, count_sum=float(raw.sum()), identical_arrays=bool(raw.shape == expression.shape and (raw != expression).nnz == 0))
    (args.output_dir / "expression_input_audit.json").write_text(json.dumps(expression_audit, indent=2) + "\n")
    (args.output_dir / "manifest.json").write_text(json.dumps({"mapping": str(args.mapping), "input_fingerprints": {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for name, path in (("mapping", args.mapping), ("native_tests", args.native_tests), ("hybrid_tests", args.hybrid_tests), ("catalog", args.catalog), ("features", args.features))}, "methods": methods, "selection": "complete finite mapped tests at source depth20, no method FDR cutoff; finite exact-zero LR effects count as nonagreement, exclusion is labeled sensitivity only", "shared_eligible": len(shared), "native_eligible": len(native), "hybrid_eligible": len(hybrid), "event_support": "fixed prepared EC support; source expression mass diagnostic sums all pseudobulk rows, never selects a model", "rank_area": "mean cumulative sign agreement at ranks1..K", "production_changes": False}, indent=2) + "\n")
    print(summary.loc[summary.cutoff.eq(100)].to_string(index=False), flush=True)
    print(pd.read_csv(args.output_dir / "top_native_support.tsv", sep="\t").to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
