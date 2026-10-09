#!/usr/bin/env python3
"""Whole selected-panel fitting assessment, not discovery or replication rates."""

import argparse
import json
from pathlib import Path

import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.pilot_library_read_models import MODELS, VARIANTS


def validate_family(table, cases):
    expected = {(row.fold, row.test_id, model, variant) for row in cases.itertuples(index=False) for model in MODELS for variant in VARIANTS}
    columns = ["fold", "test_id", "model", "variant"]
    if table.duplicated(columns).any() or set(table[columns].itertuples(index=False, name=None)) != expected:
        raise ValueError("every original selected case/model/marker identity required")
    if not table.p_value.between(0., 1.).all() or not table.converged.isin([True, False]).all() or not table.loc[~table.converged, "p_value"].eq(1.).all():
        raise ValueError("finite native p-values and failures retained at p=1 required")
    if not table.groupby(["fold", "test_id", "variant"]).counts_sha256.nunique().eq(1).all():
        raise ValueError("both models must receive identical counts within each marker variant")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise ValueError("preserve earlier assessments, new output required")
    source_files = [args.source / name for name in ("manifest.json", "diagnostics.tsv.gz", "subject_support.tsv.gz", "run_support_signatures.tsv.gz")]
    hashes = {str(path): file_hash(path) for path in source_files}
    cases = pd.read_csv(source_files[1], sep="\t").sort_values(["fold", "test_id"]).reset_index(drop=True)
    tables, receipts = [], []
    integration_policy = None
    for index in range(args.shard_count):
        folder = args.pilot_root / f"shard_{index}"
        receipt = json.loads((folder / "manifest.json").read_text())
        adaptive = receipt.get('adaptive_integration', False)
        if not isinstance(adaptive, bool) or integration_policy is not None and integration_policy != adaptive:
            raise ValueError('quadrature integration recipe must agree across all shards')
        integration_policy = adaptive
        requested = len(cases.iloc[index::args.shard_count]) * len(MODELS) * len(VARIANTS)
        if receipt["shard_index"] != index or receipt["shard_count"] != args.shard_count or receipt["source_hashes"] != hashes or receipt["whole_selected_cases"] != len(cases) or receipt["requested_fits"] != requested or receipt["completed_fits"] != requested or receipt["production_changes"] is not False:
            raise ValueError("complete compatible selected-panel shards required")
        table = pd.read_csv(folder / "tests.tsv.gz", sep="\t")
        if len(table) != requested:
            raise ValueError("declared shard output count must match receipt")
        tables.append(table)
        receipts.append(receipt)
    table = pd.concat(tables, ignore_index=True)
    validate_family(table, cases)
    summaries = []
    for keys, local in table.groupby(["fold", "panel", "model", "variant"]):
        summaries.append(dict(zip(["fold", "panel", "model", "variant"], keys), requested=len(local), usable=int(local.converged.sum()), median_runtime_seconds=float(local.runtime_seconds.median()), median_requested_subjects=float(local.requested_subjects.median()), exception_trials=int(local.error.notna().sum())))
    args.output_dir.mkdir(parents=True)
    table.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    manifest = dict(source_hashes=hashes, shards=receipts, requested_cases=len(cases), requested_fits=len(table), adaptive_integration=integration_policy, scope="numerical/local-marker diagnostic on an ascertained strong/weak panel, not full-family FDR, unbiased replication or own-ranked LR performance", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
