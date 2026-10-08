#!/usr/bin/env python3
"""Describe pooled sign-null tails without changing calibration or selection."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--null", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    observed = pd.read_csv(args.observed, sep="\t", usecols=["test_id", "calibration_stratum", "raw_p_value", "p_value"], float_precision="round_trip")
    null = pd.read_csv(args.null, sep="\t", usecols=["test_id", "calibration_stratum", "raw_p_value", "replicate"], float_precision="round_trip")
    if observed.test_id.duplicated().any() or null.duplicated(["test_id", "replicate"]).any():
        raise ValueError("unique observed and sign-null identities required")
    lookup = observed.set_index("test_id")
    parent = null.test_id.map(lookup.raw_p_value)
    if not np.isfinite(parent).all() or not null.calibration_stratum.eq(null.test_id.map(lookup.calibration_stratum)).all():
        raise ValueError("null rows must match their declared observed parents and strata")
    for values in (observed.raw_p_value, observed.p_value, null.raw_p_value):
        if not np.isfinite(values).all() or not values.between(0., 1.).all():
            raise ValueError("finite probabilities required for complete-family tail audit")
    null["parent_raw_p_value"] = parent
    records = []
    for stratum, local in null.groupby("calibration_stratum", observed=True):
        tests = observed.loc[observed.calibration_stratum.eq(stratum)]
        for threshold in (.05, .01, .001, .0001, .00001, .000001):
            tail = local.raw_p_value.le(threshold)
            small_parent = local.parent_raw_p_value.le(threshold)
            n_tail = int(tail.sum())
            from_small = int((tail & small_parent).sum())
            records.append(dict(calibration_stratum=stratum, threshold=threshold, observed_tests=len(tests), raw_observed_calls=int(tests.raw_p_value.le(threshold).sum()), calibrated_observed_calls=int(tests.p_value.le(threshold).sum()), training_draws=len(local), training_parents=local.test_id.nunique(), raw_null_tail_draws=n_tail, raw_null_tail_fraction=n_tail / len(local), tail_draws_from_small_observed_raw_parents=from_small, tail_fraction_from_small_observed_raw_parents=from_small / n_tail if n_tail else np.nan, minimum_calibrated_p_value=tests.p_value.min()))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.output_dir / "tail_contributions.tsv", sep="\t", index=False, na_rep="NA")
    receipt = dict(observed_tests=len(observed), training_draws=len(null), input_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (args.observed, args.null)}, scope="supplied pooled sign-null training distribution; cohort completeness is established by the upstream split assessor, not this diagnostic; not independent count-null assessment or split/LR replication", interpretation="small observed raw p-values are descriptive parent labels, not proof of biological alternatives; own draws remain excluded by the unchanged calibration implementation", changes="no changed testing family, fitted model, calibration, ranking or production default")
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(records).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
