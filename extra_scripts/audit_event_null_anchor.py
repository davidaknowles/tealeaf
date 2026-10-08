#!/usr/bin/env python3
"""Separate rare fitted null paths from nuisance-information loss, descriptively."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.audit_shared_score_variance import prepare_panel
from extra_scripts.reassess_event_score_archive import read_table
from tealeaf.sc.path_score_mixed import binary_null_minor_usage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-root", type=Path, required=True)
    parser.add_argument("--diagnostic-dir", type=Path, required=True)
    parser.add_argument("--bound-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=64)
    args = parser.parse_args()
    if args.output_dir.exists() or args.shard_count < 1:
        raise ValueError("new output and positive source shard count required")
    source = read_table(args.diagnostic_dir / "diagnostics.tsv.gz").set_index("test_id", verify_integrity=True)
    bounds = read_table(args.bound_dir / "diagnostics.tsv.gz").set_index("test_id", verify_integrity=True)
    if set(source.index) != set(bounds.index) or not source.status.eq("ok").all():
        raise ValueError("identical complete influence/bound diagnostic family required")
    records, receipts = {}, []
    for index in range(args.shard_count):
        path = args.cohort_root / f"shard_{index}" / "subject_scores.tsv.gz"
        table = read_table(path)
        receipts.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        for key, local in table.loc[table.test_id.isin(source.index)].groupby("test_id", sort=False):
            if key in records:
                raise ValueError("duplicate archive identity")
            records[key] = local
    if set(records) != set(source.index):
        raise ValueError("every frozen diagnostic case must remain")
    tests, subjects = [], []
    for key in sorted(records):
        panel, table = prepare_panel(records, [key])
        local = table.iloc[panel.source_positions].copy()
        if len(local) != source.at[key, "n_informative_subjects"]:
            raise ValueError("original informative subject rank changed")
        minor = binary_null_minor_usage(local.biological_shape)
        retained_fraction = local.information.to_numpy() / local.reference_information.to_numpy()
        equivalent = local.information.to_numpy() * local.biological_shape.to_numpy()
        dominant = np.flatnonzero(local.subject.astype(str).eq(str(source.at[key, "dominant_subject"])))
        if len(dominant) != 1:
            raise ValueError("original dominant subject identity changed")
        tests.append(dict(test_id=key, panel=source.at[key, "panel"], n_subjects=len(local), unavoidable_share_gt_90=bool(bounds.at[key, "largest_uniform_lower_bound"] > .9), n_null_minor_usage_le_1e6=int((minor <= 1e-6).sum()), n_null_minor_usage_ge_01=int((minor >= .01).sum()), median_null_minor_usage=float(np.median(minor)), n_profiled_fraction_ge_10pct=int((retained_fraction >= .1).sum()), median_profiled_information_fraction=float(np.median(retained_fraction)), n_equivalent_counts_ge_5=int((equivalent >= 5).sum()), dominant_null_minor_usage=float(minor[dominant[0]]), dominant_profiled_information_fraction=float(retained_fraction[dominant[0]]), dominant_equivalent_count=float(equivalent[dominant[0]])))
        subjects.extend(dict(test_id=key, subject=str(subject), null_minor_usage=float(usage), profiled_information_fraction=float(fraction), measurement_equivalent_count=float(count), report_inclusion_a=float(report_a), report_inclusion_b=float(report_b)) for subject, usage, fraction, count, report_a, report_b in zip(local.subject, minor, retained_fraction, equivalent, local.report_inclusion_a, local.report_inclusion_b))
    table = pd.DataFrame(tests)
    summaries = []
    for (name, unavoidable), local in table.groupby(["panel", "unavoidable_share_gt_90"]):
        summaries.append(dict(panel=name, unavoidable_share_gt_90=bool(unavoidable), requested=len(local), at_least_four_nonrare_null_subjects=int(local.n_null_minor_usage_ge_01.ge(4).sum()), at_least_four_profiled_fractions_ge_10pct=int(local.n_profiled_fraction_ge_10pct.ge(4).sum()), at_least_four_equivalent_counts_ge_5=int(local.n_equivalent_counts_ge_5.ge(4).sum()), median_null_minor_usage=local.median_null_minor_usage.median(), median_profiled_information_fraction=local.median_profiled_information_fraction.median(), median_dominant_null_minor_usage=local.dominant_null_minor_usage.median(), median_dominant_equivalent_count=local.dominant_equivalent_count.median()))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "diagnostics.tsv.gz", sep="\t", index=False)
    pd.DataFrame(subjects).to_csv(args.output_dir / "subject_geometry.tsv.gz", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(source=receipts, selected_diagnostics_sha256=hashlib.sha256((args.diagnostic_dir / "diagnostics.tsv.gz").read_bytes()).hexdigest(), bound_sha256=hashlib.sha256((args.bound_dir / "diagnostics.tsv.gz").read_bytes()).hexdigest(), requested=len(table), definition="minor null usage is the small root of p*(1-p)=1/B in original binary Helmert ILR units; information fraction is profiled J over unprofiled R", selection="same prior strong-tail and informative-subject-count-matched weak-control panels, no LR selection", scope="diagnostic fitted geometry, not observed read support, fitted heterogeneity, an eligibility rule or new inference; independent A1 reports may differ from constrained null compositions", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
