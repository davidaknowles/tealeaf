#!/usr/bin/env python3
"""Compare cached primer PSI, nuisance fits and paired effects to an archive."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--psi-dir", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--repro-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalog = pd.read_csv(args.psi_dir / "event_catalog.tsv.gz", sep="\t").set_index("feature_id")
    previous_catalog = pd.read_csv(args.archive / "event_catalog.tsv.gz", sep="\t").set_index("feature_id")
    shared = catalog.index.intersection(previous_catalog.index)
    for column in catalog.columns:
        if column.startswith(("primer_offset", "primer_dispersion")):
            error = np.abs(catalog.loc[shared, column] - previous_catalog.loc[shared, column])
            print(f"{column} shared={len(shared)} max_error={error.max()} median_error={error.median()}", flush=True)
    print(f"catalogues new={len(catalog)} old={len(previous_catalog)}", flush=True)
    lookup = {sample: i for i, sample in enumerate(json.loads((args.psi_dir / "samples.json").read_text()))}
    summaries, mismatches = [], []
    for fold in (0, 1):
        psi = np.load(args.psi_dir / f"psi_fold{fold}.npz")["psi"]
        old = pd.read_csv(args.archive / f"split_data_fold{fold}_tests.tsv.gz", sep="\t")
        print(f"fold={fold} archived_rows={len(old)} archived_fold_labels={old.fold.unique().tolist()} archived_subject_counts={old.n_subjects.value_counts().sort_index().to_dict()} cache_shape={psi.shape}", flush=True)
        contrasts = json.loads((args.repro_dir / f"fold{fold}/contrasts.json").read_text())
        for contrast in contrasts:
            table = old.loc[old.contrast_id.eq(contrast["contrast_id"])].set_index("feature_id").reindex(catalog.index)
            pairs = [(lookup[a], lookup[b]) for a, b in zip(contrast["samples_a"], contrast["samples_b"]) if a in lookup and b in lookup]
            if not pairs:
                continue
            first, second = psi[:, [a for a, _ in pairs]], psi[:, [b for _, b in pairs]]
            valid = np.isfinite(first) & np.isfinite(second)
            counts = valid.sum(axis=1)
            effects = np.divide(np.where(valid, second - first, 0).sum(axis=1), counts, out=np.full(len(catalog), np.nan), where=counts > 0)
            finite = table.effect_size.notna().to_numpy() & np.isfinite(effects)
            error = np.abs(effects - table.effect_size.to_numpy())
            count_error = counts - table.n_subjects.to_numpy()
            bad = finite & ((error > 1e-10) | (count_error != 0))
            summaries.append({"fold": fold, "contrast_id": contrast["contrast_id"], "compared_events": int(finite.sum()), "mismatched_events": int(bad.sum()), "maximum_effect_error": float(np.max(error[finite])) if finite.any() else np.nan, "maximum_count_error": float(np.max(np.abs(count_error[finite]))) if finite.any() else np.nan})
            for index in np.flatnonzero(bad)[:20]:
                mismatches.append({"fold": fold, "contrast_id": contrast["contrast_id"], "feature_id": catalog.index[index], "new_effect": effects[index], "old_effect": table.effect_size.iloc[index], "new_subjects": counts[index], "old_subjects": table.n_subjects.iloc[index], "new_psi_a": json.dumps(first[index].tolist()), "new_psi_b": json.dumps(second[index].tolist())})
        print(f"fold={fold} mismatched={sum(row['mismatched_events'] for row in summaries if row['fold'] == fold)}", flush=True)
    pd.DataFrame(summaries).to_csv(args.output_dir / "contrast_checks.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(mismatches).to_csv(args.output_dir / "mismatch_examples.tsv.gz", sep="\t", index=False)
    print(pd.DataFrame(mismatches).head(10).drop(columns=["new_psi_a", "new_psi_b"]).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
