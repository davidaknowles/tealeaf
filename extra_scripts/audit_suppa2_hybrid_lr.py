"""Separate event-universe, ranking, and effect-direction differences in the hybrid audit."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def summarize(frame, scope, ranking, effect, event_type="all"):
    columns, ascending = [ranking], [True]
    if ranking == "hybrid_p":
        for column, direction in (("raw_p_value", True), ("statistic", False)):
            if column in frame:
                columns.append(column)
                ascending.append(direction)
    columns.extend(["contrast_id", "feature_id"])
    ascending.extend([True, True])
    ordered = frame.sort_values(columns, ascending=ascending, kind="stable")
    rows = []
    for top in (25, 75, 200, len(ordered)):
        local = ordered.head(top)
        rows.append({"scope": scope, "ranking": ranking, "direction": effect, "event_type": event_type, "top": top, "n": len(local), "agreement": float((local[effect] * local.long_read_effect > 0).mean()), "median_abs_lr": local.long_read_effect.abs().median(), "median_subjects": local.n_subjects.median()})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hybrid-mapping", required=True, type=Path)
    parser.add_argument("--hybrid-tests", required=True, type=Path)
    parser.add_argument("--native-tests", required=True, type=Path)
    parser.add_argument("--native-catalog", type=Path)
    parser.add_argument("--hybrid-catalog", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    keys = ["contrast_id", "feature_id"]
    hybrid = pd.read_csv(args.hybrid_mapping, sep="\t")
    tests = pd.read_csv(args.hybrid_tests, sep="\t")
    hybrid = hybrid.merge(tests[keys + ["n_subjects", "mean_difference_norm"]], on=keys, validate="one_to_one")
    hybrid = hybrid[hybrid.minimum_pooled_depth.ge(20) & hybrid.long_read_effect.notna() & hybrid.short_read_effect.notna() & hybrid.short_read_effect.ne(0)].copy()
    hybrid = hybrid.rename(columns={"p_value": "hybrid_p", "short_read_effect": "hybrid_effect"})
    native = pd.read_csv(args.native_tests, sep="\t", usecols=keys + ["p_value", "effect_size", "n_subjects"])
    native = native.rename(columns={"p_value": "native_p", "effect_size": "native_effect", "n_subjects": "native_subjects"})
    matched = hybrid.merge(native, on=keys, validate="one_to_one").dropna(subset=["native_p", "native_effect"])
    rows = summarize(hybrid, "hybrid universe", "hybrid_p", "hybrid_effect")
    for kind, subset in hybrid.groupby("event_type"):
        rows.extend(summarize(subset, "hybrid universe", "hybrid_p", "hybrid_effect", kind))
    for ranking in ("hybrid_p", "native_p"):
        for effect in ("hybrid_effect", "native_effect"):
            rows.extend(summarize(matched, "shared hypotheses", ranking, effect))
            for kind, subset in matched.groupby("event_type"):
                rows.extend(summarize(subset, "shared hypotheses", ranking, effect, kind))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows).drop_duplicates()
    summary.to_csv(args.output_dir / "ranking_direction_audit.tsv", sep="\t", index=False)
    matched.to_csv(args.output_dir / "matched_events.tsv.gz", sep="\t", index=False)
    hybrid.groupby("event_type").agg(tests=("feature_id", "size"), median_subjects=("n_subjects", "median")).to_csv(args.output_dir / "event_universe.tsv", sep="\t")
    if args.native_catalog and args.hybrid_catalog:
        catalogs = [pd.read_csv(path, sep="\t") for path in (args.native_catalog, args.hybrid_catalog)]
        catalog_rows = []
        for label, catalog in zip(("SUPPA2", "hybrid"), catalogs):
            for kind, group in catalog.groupby("event_type"):
                catalog_rows.append({"catalog": label, "event_type": kind, "events": len(group), "shared_event_ids": int(group.feature_id.isin(catalogs[1 if label == "SUPPA2" else 0].feature_id).sum())})
        pd.DataFrame(catalog_rows).to_csv(args.output_dir / "catalog_universe.tsv", sep="\t", index=False)
    print(summary[summary.event_type.eq("all")].to_string(index=False))
    print("Hybrid top-200 event types", hybrid.nsmallest(200, "hybrid_p").event_type.value_counts().to_dict())
    print("Matched tests", len(matched), "short-read sign concordance", (matched.hybrid_effect * matched.native_effect > 0).mean())


if __name__ == "__main__":
    main()
