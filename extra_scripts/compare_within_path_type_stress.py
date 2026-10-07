"""Compare frozen real-design nuisance changes with no-change controls."""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    tables = []
    keys = ["test_id", "draw", "strategy"]
    reference_keys = None
    for fitter in ("fixed", "free"):
        control = pd.read_csv(args.results / "controls" / fitter / "tests.tsv.gz", sep="\t")
        shape = control[["test_id", "n_transcripts", "n_multiplet_paths", "n_outside_transcripts"]].drop_duplicates()
        if shape.test_id.duplicated().any():
            raise ValueError("inconsistent transcript structure by test")
        for condition, table in (("no within-path change", control), ("within-path type SD1", pd.read_csv(args.results / fitter / "tests.tsv.gz", sep="\t"))):
            actual = set(map(tuple, table[keys].values))
            if table.duplicated(keys).any() or (reference_keys is not None and actual != reference_keys):
                raise ValueError("stress/control families differ or are duplicated")
            reference_keys = actual
            if "n_multiplet_paths" not in table:
                table = table.merge(shape, on="test_id", validate="many_to_one")
            table["within_path_structure"] = table.n_multiplet_paths.gt(0).map({True: "at least one multiplet path", False: "singleton paths, no nuisance tilt possible"})
            table["condition"], table["fitter"] = condition, fitter
            table["converged"] = table.converged.astype(str).str.lower().eq("true")
            tables.append(table)
    matched = pd.concat(tables, ignore_index=True)
    def summarize(keys):
        return matched.groupby(keys).agg(n_requested=("p_value", "size"), n_hypotheses=("test_id", "nunique"), n_converged=("converged", "sum"), native_reject_0_05=("raw_p_value", lambda values: values.le(.05).sum()), calibrated_reject_0_05=("p_value", lambda values: values.le(.05).sum()), native_reject_rate_0_05=("raw_p_value", lambda values: values.le(.05).mean()), calibrated_reject_rate_0_05=("p_value", lambda values: values.le(.05).mean())).reset_index()
    summary = summarize(["fitter", "condition", "strategy"])
    strata = summarize(["fitter", "condition", "strategy", "within_path_structure"])
    matched.to_csv(args.results / "matched_tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary.to_csv(args.results / "matched_summary.tsv", sep="\t", index=False, na_rep="NA")
    strata.to_csv(args.results / "structure_summary.tsv", sep="\t", index=False, na_rep="NA")
    print(summary.to_string(index=False))
    print(strata.to_string(index=False))


if __name__ == "__main__":
    main()
