"""Assess hybrid measurement ablations on fixed long-read effects and matched SUPPA2 hypotheses."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def ordered(table):
    return table.sort_values(["p_value", "raw_p_value", "statistic", "contrast_id", "feature_id"], ascending=[True, True, False, True, True], kind="stable")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--native-tests", required=True, type=Path)
    parser.add_argument("--variant", nargs=2, action="append", required=True, metavar=("LABEL", "TESTS"))
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    keys = ["contrast_id", "feature_id"]
    source = pd.read_csv(args.mapping, sep="\t")
    source = source[source.minimum_pooled_depth.ge(20) & source.long_read_effect.notna()]
    source = source[keys + ["long_read_effect", "minimum_pooled_depth", "minimum_replicate_depth", "event_type"]]
    native = pd.read_csv(args.native_tests, sep="\t", usecols=keys + ["p_value", "effect_size", "n_subjects"])
    native = native.rename(columns={"p_value": "native_p", "effect_size": "native_effect", "n_subjects": "native_subjects"})
    output, details = [], []
    for label, path in args.variant:
        tests = pd.read_csv(path, sep="\t")
        tests = tests[tests.converged & tests.p_value.notna()].copy()
        tests = tests.drop(columns=["event_type"], errors="ignore")
        local = source.merge(tests, on=keys, validate="one_to_one")
        local = local.merge(native, on=keys, how="left", validate="one_to_one")
        local["variant"] = label
        details.append(local)
        scopes = {"all eligible": local, "shared SUPPA2": local[local.native_p.notna()], "SUPPA2 event categories": local[~local.event_type.isin(["AF", "AL"])], "at least 8 subjects": local[local.n_subjects.ge(8)]}
        if "baseline_event_mass" in local:
            scopes["no nuisance isoforms"] = local[local.baseline_event_mass.ge(.999999)]
            scopes["with nuisance isoforms"] = local[local.baseline_event_mass.lt(.999999)]
        for scope, subset in scopes.items():
            for column in ("effect_size", "test_ilr_effect_size", "report_ilr_effect", "report_psi_effect", "native_effect"):
                if column not in subset:
                    continue
                eligible = subset[subset[column].notna() & subset[column].ne(0)].copy()
                ranked = ordered(eligible)
                for top in dict.fromkeys((25, 75, 200, len(ranked))):
                    group = ranked.head(top)
                    paired = group[group.native_effect.notna() & group.native_effect.ne(0)]
                    output.append({"variant": label, "scope": scope, "effect": column, "rank": top, "n": len(group), "positive": int((group[column] * group.long_read_effect > 0).sum()), "agreement": float((group[column] * group.long_read_effect > 0).mean()), "median_subjects": group.n_subjects.median(), "median_abs_lr": group.long_read_effect.abs().median(), "native_sign_agreement": float((paired[column] * paired.native_effect > 0).mean()) if len(paired) else np.nan})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(output)
    summary.to_csv(args.output_dir / "variant_direction_summary.tsv", sep="\t", index=False)
    pd.concat(details, ignore_index=True).to_csv(args.output_dir / "variant_matched_events.tsv.gz", sep="\t", index=False)
    print(summary[summary.scope.eq("all eligible") & summary["rank"].isin([75, 200])].to_string(index=False))


if __name__ == "__main__":
    main()
