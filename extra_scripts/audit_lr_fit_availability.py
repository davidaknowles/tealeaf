#!/usr/bin/env python3
"""Inspect complete candidate fits at an unchanged comparator LR prefix.

This is an availability/direction audit, never the candidate's own ranking.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash


def audit_prefix(prefix, fits):
    keys = ["feature_id", "contrast_id"]
    if prefix.duplicated(keys).any() or fits.duplicated(keys).any():
        raise ValueError("unique event/contrast identities required")
    if set(prefix['rank']) != set(range(1, 201)):
        raise ValueError("unchanged complete top-200 native prefix required")
    columns = [*keys, "converged", "complete_subject_fits", "complete_reporting_fits", "n_subjects", "n_expected_subjects", "n_fitted_subjects", "p_value", "raw_p_value", "effect_size", "test_ilr_effect_size"]
    table = prefix.merge(fits[columns], on=keys, how="left", validate="one_to_one", indicator=True)
    table["requested_in_completed_candidate"] = table._merge.eq("both")
    for column in ("converged", "complete_subject_fits", "complete_reporting_fits"):
        table[column] = table[column].astype(str).str.lower().eq("true")
    table["finite_usage_direction"] = np.isfinite(table.effect_size) & table.effect_size.ne(0.) & table.complete_reporting_fits
    table["finite_score_direction"] = np.isfinite(table.test_ilr_effect_size) & table.test_ilr_effect_size.ne(0.) & table.complete_subject_fits
    return table.drop(columns="_merge")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--fits", type=Path, required=True, help="Complete guarded merged table including failed requests.")
    parser.add_argument("--completion-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("new audit output required")
    receipt = json.loads(args.completion_manifest.read_text())
    cohorts = receipt["cohorts"]
    if len(cohorts) != 3 or {row['fold'] for row in cohorts} != {0, 1, "full"}:
        raise ValueError("completed three-cohort family required")
    full = next(row for row in cohorts if row['fold'] == "full")
    if len(full['shards']) != 64 or Path(full['merged']).resolve() != args.fits.parent.resolve() or any(row['completed'] + row['failures'] != row['tests_in_shard'] for row in full['shards']):
        raise ValueError("all 64 complete compatible full-data shards required")
    prefix = pd.read_csv(args.prefix, sep="\t")
    fits = pd.read_csv(args.fits, sep="\t", low_memory=False)
    if len(fits) != sum(row['tests_in_shard'] for row in full['shards']):
        raise ValueError("all original requested fits, including failures, required")
    table = audit_prefix(prefix, fits)
    summaries = []
    for cutoff in (100, 200):
        local = table.loc[table['rank'].le(cutoff)]
        summaries.append(dict(cutoff=cutoff, requested_in_prefit_audit=int(local.eligibility.eq("requested before fitting").sum()), requested_in_completed_candidate=int(local.requested_in_completed_candidate.sum()), complete_subject_fits=int(local.complete_subject_fits.sum()), complete_reporting_fits=int(local.complete_reporting_fits.sum()), finite_usage_direction=int(local.finite_usage_direction.sum()), finite_score_direction=int(local.finite_score_direction.sum())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "rank_prefix.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(input_hashes={str(path): file_hash(path) for path in (args.prefix, args.fits, args.completion_manifest)}, requested_full_fits=len(fits), scope="completed-fit availability at frozen SUPPA2 LR prefix, not candidate own-ranked LR agreement, power or calibration", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
