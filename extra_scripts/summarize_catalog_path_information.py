"""Collate the complete catalog information audit with explicit denominators."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def summary(group):
    weights = group.n_contrasts.to_numpy(int)
    ratios = group.retained_information_fraction.to_numpy(float)
    finite = np.isfinite(ratios)
    row = dict(events=len(group), genes=group.gene_id.nunique(), event_contrasts=int(weights.sum()), pooled_fit_converged_events=int(group.pooled_baseline_converged.sum()), ratio_undefined_events=int((~finite).sum()), ratio_undefined_contrast_fraction=float(np.average(~finite, weights=weights)))
    if finite.any():
        weighted = np.repeat(ratios[finite], weights[finite])
        for probability in (.1, .25, .5, .75, .9):
            row[f"ratio_q{probability}"] = float(np.quantile(ratios[finite], probability))
            row[f"contrast_weighted_ratio_q{probability}"] = float(np.quantile(weighted, probability))
    for cutoff in (1e-10, .01, .1, .5, .9):
        row[f"ratio_below_{cutoff}_events"] = int((finite & (ratios < cutoff)).sum())
        row[f"ratio_below_{cutoff}_contrast_fraction"] = float(np.average(finite & (ratios < cutoff), weights=weights))
    for cutoff in (1, 5, 10, 25):
        row[f"expected_path_origins_below_{cutoff}_events"] = int(group.minimum_expected_path_origins.lt(cutoff).sum())
        row[f"expected_path_origins_below_{cutoff}_contrast_fraction"] = float(np.average(group.minimum_expected_path_origins.lt(cutoff), weights=weights))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    frames, summaries, metadata = [], [], []
    for source in ("parsimony_binary", "original_binary"):
        directories = [args.input_root / source / f"shard_{index}" for index in range(8)]
        settings = [json.loads((path / "settings.json").read_text()) for path in directories]
        if any(value != settings[0] for value in settings):
            raise ValueError("inconsistent catalog-audit settings")
        table = pd.concat([pd.read_csv(path / "events.tsv.gz", sep="\t") for path in directories], ignore_index=True)
        observed = set(map(tuple, table[["feature_id", "anchor"]].to_numpy()))
        expected = {(feature, anchor) for feature in settings[0]["requested_features"] for anchor in settings[0]["anchors"]}
        if len(table) != len(observed) or observed != expected or not table.source.eq(source).all():
            raise ValueError("incomplete, duplicated or wrong-source audited event family")
        frames.append(table)
        metadata.append(settings[0])
        for anchor, group in table.groupby("anchor", sort=False):
            summaries.append(dict(source=source, anchor=anchor, event_type="all", **summary(group)))
            for event_type, subset in group.groupby("event_type", sort=False):
                summaries.append(dict(source=source, anchor=anchor, event_type=event_type, **summary(subset)))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(args.output_dir / "events.tsv.gz", sep="\t", index=False)
    result = pd.DataFrame(summaries)
    result.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(source_settings=metadata, audit="free type-specific transcript shares versus fixed within-class shares, both retain free event mass", scope="complete frozen supported event definitions, four common-null anchors", denominators="event counts and frozen gene-contrast workload weights, no eligibility/significance/LR subset", caution="label-independent pooled anchor uses all gene pseudobulks, not each actual contrast; expected path origins include ambiguous reads and are model-derived, not observed uniquely informative support", inference="geometry only, no count-null calibration, biological power, data-split replication or own-ranked LR claim", production_changes=False), indent=2) + "\n")
    print(result.loc[result.event_type.eq("all"), ["source", "anchor", "events", "event_contrasts", "ratio_q0.5", "contrast_weighted_ratio_q0.5", "ratio_below_0.1_contrast_fraction", "expected_path_origins_below_5_contrast_fraction"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
