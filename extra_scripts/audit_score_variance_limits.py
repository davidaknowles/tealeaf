#!/usr/bin/env python3
"""All-variance influence bounds on the existing frozen subject-audit panels."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_shared_score_variance import prepare_panel
from extra_scripts.reassess_event_score_archive import read_table
from tealeaf.sc.score_variance_pooling import scalar_precision_share_bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", type=Path, required=True)
    parser.add_argument("--diagnostic-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=64)
    args = parser.parse_args()
    if args.output_dir.exists() or args.shard_count < 1:
        raise ValueError("new output and positive shard count required")
    original = read_table(args.diagnostic_dir / "diagnostics.tsv.gz").set_index("test_id", verify_integrity=True)
    if not original.status.eq("ok").all():
        raise ValueError("complete original archive diagnostics required")
    requested, records = set(original.index), {}
    for index in range(args.shard_count):
        archive = read_table(args.cohort_root / f"shard_{index}" / "subject_scores.tsv.gz")
        for key, local in archive.loc[archive.test_id.isin(requested)].groupby("test_id", sort=False):
            if key in records:
                raise ValueError("duplicate archive identity")
            records[key] = local
    if set(records) != requested:
        raise ValueError("every original diagnostic request must remain")
    tests, subjects = [], []
    for key in sorted(records):
        panel, table = prepare_panel(records, [key])
        local = table.iloc[panel.source_positions]
        if panel.n_subjects[0] != original.at[key, "n_informative_subjects"]:
            raise ValueError("original diagnostic information rank changed")
        lower, upper = scalar_precision_share_bounds(local.information, local.biological_shape)
        equivalent = local.information.to_numpy() * local.biological_shape.to_numpy()
        if not np.isfinite(equivalent).all():
            raise ValueError("information-equivalent count exceeds finite range")
        subject_ids = local.subject.astype(str).to_numpy()
        dominant = np.flatnonzero(subject_ids == str(original.at[key, "dominant_subject"]))
        if len(dominant) != 1:
            raise ValueError("original dominant subject is missing or duplicated")
        tests.append(dict(test_id=key, panel=original.at[key, "panel"], n_subjects=len(local), original_maximum_precision_share=original.at[key, "maximum_precision_share"], largest_uniform_lower_bound=float(lower.max()), original_dominant_uniform_lower_bound=float(lower[dominant[0]]), unavoidable_dominant_subject=str(subject_ids[np.argmax(lower)]), n_information_equivalent_counts_ge_1=int((equivalent >= 1).sum()), n_information_equivalent_counts_ge_5=int((equivalent >= 5).sum()), n_information_equivalent_counts_ge_10=int((equivalent >= 10).sum()), median_information_equivalent_count=float(np.median(equivalent))))
        subjects.extend(dict(test_id=key, subject=subject, minimum_share_bound=float(lo), maximum_share_bound=float(hi), information_equivalent_count=float(count)) for subject, lo, hi, count in zip(subject_ids, lower, upper, equivalent))
    table = pd.DataFrame(tests)
    summaries = []
    for name, local in table.groupby("panel"):
        summaries.append(dict(panel=name, requested=len(local), original_share_gt_90=int(local.original_maximum_precision_share.gt(.9).sum()), unavoidable_share_gt_90=int(local.largest_uniform_lower_bound.gt(.9).sum()), original_dominant_unavoidable_share_gt_90=int(local.original_dominant_uniform_lower_bound.gt(.9).sum()), unavoidable_share_gt_half=int(local.largest_uniform_lower_bound.gt(.5).sum()), at_least_four_equivalent_counts_ge_5=int(local.n_information_equivalent_counts_ge_5.ge(4).sum())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "diagnostics.tsv.gz", sep="\t", index=False)
    pd.DataFrame(subjects).to_csv(args.output_dir / "subject_bounds.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    receipt = dict(source=str(args.cohort_root), diagnostic_source=str(args.diagnostic_dir), requested=len(table), proof="w_j(t)/w_i(t) lies between I_j/I_i and B_i/B_j; reciprocal sums of pairwise maxima/minima bound share_i for all t>=0 and infinity", count_proxy="I*B, measurement-equivalent paired path count; equals per-type count for ideal balanced two-path counts, not observed junction/exon counts", scope="same original selected panels and original rank mask; descriptive all-variance influence bound, not a new eligibility filter, direction estimator, statistical test or full-family replication result", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
