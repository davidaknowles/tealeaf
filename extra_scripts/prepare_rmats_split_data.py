#!/usr/bin/env python3
"""Prepare validated rMATS count indices for the production subject halves."""

import argparse
import json
from pathlib import Path

import pandas as pd

from tealeaf.sc.junction_benchmark import index_subject_paired_contrasts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-manifest", type=Path, required=True)
    parser.add_argument("--split-manifest", type=Path, action="append", required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--count-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if len(args.split_manifest) != 2:
        raise ValueError("exactly two subject-fold manifests are required")
    mapping = {}
    for contrast in json.loads(args.full_manifest.read_text()):
        for side in ("a", "b"):
            for sample, index in zip(contrast[f"samples_{side}"], contrast[f"indices_{side}"], strict=True):
                if sample in mapping and mapping[sample] != index:
                    raise ValueError("full-data manifest has inconsistent count indices")
                mapping[sample] = int(index)
    expected_size = max(mapping.values()) + 1
    for event in ("SE", "A3SS", "A5SS", "MXE", "RI"):
        with (args.count_dir / f"JCEC.raw.input.{event}.txt").open() as handle:
            header = next(handle).strip().split("\t")
            row = dict(zip(header, next(handle).strip().split("\t"), strict=True))
        size = sum(len(row[field].split(",")) if row[field] else 0 for field in ("IJC_SAMPLE_1", "IJC_SAMPLE_2"))
        if size != expected_size:
            raise ValueError(f"{event} count-pass size {size} != manifest size {expected_size}")
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    folds = [index_subject_paired_contrasts(json.loads(path.read_text()), samples, mapping) for path in args.split_manifest]
    subjects = [{subject for contrast in fold for subject in contrast["paired_subjects"]} for fold in folds]
    if subjects[0] & subjects[1]:
        raise ValueError("subject halves overlap")
    if {item["contrast_id"] for item in folds[0]} != {item["contrast_id"] for item in folds[1]}:
        raise ValueError("fold contrast catalogues differ")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    offset = 0
    for fold, contrasts in enumerate(folds):
        for index, contrast in enumerate(contrasts):
            contrast["index"] = offset + index
        (args.output_dir / f"fold{fold}_contrasts.json").write_text(json.dumps(contrasts, indent=2) + "\n")
        offset += len(contrasts)
    summary = {"fold_contrasts": list(map(len, folds)), "fold_subjects": list(map(len, subjects)), "array_tasks": offset, "count_pass_samples": expected_size, "sample_indices_recovered": len(mapping)}
    (args.output_dir / "manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
