"""Compare generic/scalar real-data smoke runs, not a full endpoint benchmark."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.compare_count_null_implementations import compare_tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    roots = [args.input_root / name for name in ("generic", "scalar")]
    settings = [json.loads((root / "settings.json").read_text()) for root in roots]
    for setting, fast in zip(settings, (False, True)):
        if setting["arguments"].pop("scalar_fast") != fast:
            raise ValueError("wrong implementation label")
        setting["arguments"].pop("output_dir")
    if settings[0] != settings[1]:
        raise ValueError("smoke run models or inputs differ")
    summaries = [json.loads((root / "summary.json").read_text()) for root in roots]
    if any(value["tests_in_shard"] != 4 or value["completed"] < 1 for value in summaries):
        raise ValueError("smoke must have four declared tests and nonempty fitted results")
    failures = [json.loads((root / "failures.json").read_text()) for root in roots]
    if failures[0] != failures[1]:
        raise ValueError("failure retention differs")
    comparison = []
    for name, keys, columns in (("paired_path.tsv", ["test_id"], ["p_value", "statistic", "effect_size", "test_ilr_effect_size", "mean_difference_norm", "n_subjects", "n_fitted_subjects", "report_n_subjects", "converged"]), ("paired_path_null.tsv.gz", ["test_id", "replicate"], ["p_value"])):
        first, second = [pd.read_csv(root / name, sep="\t") for root in roots]
        comparison.append({"table": name, **compare_tables(first, second, keys, columns)})
    first, second = [pd.read_csv(root / "subject_scores.tsv.gz", sep="\t").set_index(["test_id", "subject"], verify_integrity=True).sort_index() for root in roots]
    if not first.index.equals(second.index):
        raise ValueError("upstream subject families differ")
    arrays_equal = all(np.array_equal(first[column], second[column], equal_nan=True) for column in ("score", "information", "reference_information", "biological_shape", "report_inclusion_a", "report_inclusion_b"))
    passed = all(all(value for key, value in row.items() if key.endswith("_matches")) and all(value == 0 for key, value in row.items() if key.startswith("p_decisions_differ")) for row in comparison)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(comparison).to_csv(args.output_dir / "equivalence.tsv", sep="\t", index=False)
    replay = (roots[0] / "reassessment.json").exists()
    times = [summary["elapsed_seconds"] if not replay or "timing_scope" in summary else None for summary in summaries]
    scope = "archive replay, not count fitting; time unavailable for older untimed replay" if replay else "fitting loop excluding initial loading/screening"
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(passed=passed, upstream_arrays_bit_identical=arrays_equal, requested_tests=4, completed=summaries[0]["completed"], generic_processing_seconds=times[0], scalar_processing_seconds=times[1], timing_scope=scope, limitation="four-test pipeline smoke only, archive replay timing is not new count fitting, not full split or LR performance"), indent=2) + "\n")
    print(f"four-test pipeline smoke, passed={passed}, identical upstream arrays={arrays_equal}", flush=True)
    if not passed:
        raise ValueError("pipeline smoke does not meet numerical tolerances")


if __name__ == "__main__":
    main()
