#!/usr/bin/env python3
"""Validate and merge all paired rMATS split contrasts, never partial arrays."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifests = [json.loads((args.input_dir / f"fold{k}_contrasts.json").read_text()) for k in (0, 1)]
    missing = [contrast["index"] for fold in manifests for contrast in fold if not (args.input_dir / f"results/{contrast['index']}.tsv.gz").exists()]
    if missing:
        raise ValueError(f"missing rMATS tasks, {missing}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for fold, contrasts in enumerate(manifests):
        output = args.output_dir / f"split_data_fold{fold}_tests.tsv.gz"
        with gzip.open(output, "wt") as handle:
            for number, contrast in enumerate(contrasts):
                table = pd.read_csv(args.input_dir / f"results/{contrast['index']}.tsv.gz", sep="\t", low_memory=False)
                for column, expected in (("fold", fold), ("contrast_id", contrast["contrast_id"]), ("level_a", contrast["level_a"]), ("level_b", contrast["level_b"]), ("statistical_model", "PAIRADISE"), ("effect_orientation", "a_minus_b"), ("n_subjects", len(contrast["paired_subjects"]))):
                    if not table[column].eq(expected).all():
                        raise ValueError(f"rMATS task {contrast['index']} has inconsistent {column}")
                if table.feature_id.duplicated().any() or set(table.event_type) != {"SE", "A3SS", "A5SS", "MXE", "RI"}:
                    raise ValueError("rMATS event catalogue is incomplete or duplicated")
                finite = np.isfinite(table.p_value)
                if not finite.any() or not table.loc[finite, "p_value"].between(0, 1).all():
                    raise ValueError("rMATS task has no valid p-values or values outside [0, 1]")
                summaries.append({"fold": fold, "contrast_id": contrast["contrast_id"], "n_subjects": len(contrast["paired_subjects"]), "n_catalogue_events": len(table), "n_finite_tests": int(finite.sum()), "n_native_bh": int(table.q_value.le(.05).sum()), "n_nonzero_effects": int((np.isfinite(table.effect_size) & table.effect_size.ne(0)).sum())})
                table.to_csv(handle, sep="\t", index=False, header=number == 0, na_rep="NA")
    pd.DataFrame(summaries).to_csv(args.output_dir / "split_data_contrast_summary.tsv", sep="\t", index=False)
    summary = {"statistical_model": "PAIRADISE", "count_type": "JCEC", "read_handling": "single_end", "effect_orientation": "a_minus_b", "fold_contrasts": list(map(len, manifests)), "fold_subjects": [len({subject for contrast in fold for subject in contrast["paired_subjects"]}) for fold in manifests], "source": "retained corrected all-sample count pass; split-specific statistical refits"}
    (args.output_dir / "split_data_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(pd.DataFrame(summaries).groupby("fold")[["n_catalogue_events", "n_finite_tests", "n_native_bh"]].sum().to_string(), flush=True)


if __name__ == "__main__":
    main()
