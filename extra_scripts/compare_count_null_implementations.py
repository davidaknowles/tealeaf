"""Check a numerical optimization against identical complete count-null trials."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.summarize_ec_count_null import validate_requested_trials


def normalize_settings(settings):
    """Known legacy defaults only, never discard different model settings."""
    settings = dict(settings)
    for key, default in dict(count_likelihood="multinomial", ec_opportunity_scale=0., kernel_units="prepared", simulation_kernel_units="analysis", scalar_fast=False).items():
        settings.setdefault(key, default)
    return settings


def load(root, name):
    shards = [root / f"shard_{index}" for index in range(16)]
    settings = [json.loads((shard / "settings.json").read_text()) for shard in shards]
    if any(value != settings[0] for value in settings):
        raise ValueError("incomplete or inconsistent implementation shards")
    frames = [pd.read_csv(shard / name, sep="\t") for shard in shards]
    table = pd.concat(frames, ignore_index=True)
    if name == "observed.tsv":
        validate_requested_trials(table, settings[0])
    return table, settings[0]


def compare_tables(original, candidate, keys, columns):
    first, second = original.set_index(keys, verify_integrity=True).sort_index(), candidate.set_index(keys, verify_integrity=True).sort_index()
    if not first.index.equals(second.index):
        raise ValueError("numerical comparison changed the declared trial family")
    row = dict(trials=len(first))
    for column in columns:
        x, y = first[column].to_numpy(float), second[column].to_numpy(float)
        finite = np.isfinite(x) & np.isfinite(y)
        if not np.array_equal(np.isfinite(x), np.isfinite(y)):
            raise ValueError("optimization changed numerical availability")
        row[column + "_matches"] = bool(np.allclose(x, y, rtol=3e-6, atol=1e-8, equal_nan=True))
        difference = np.abs(x[finite] - y[finite])
        row[column + "_max_absolute_difference"] = float(difference.max()) if finite.any() else 0.
        row[column + "_max_scaled_difference"] = float(np.max(difference / np.maximum(np.abs(x[finite]), 1e-8))) if finite.any() else 0.
    for threshold in (.05, .01, .001):
        row[f"p_decisions_differ_at_{threshold}"] = int(np.count_nonzero((first.p_value < threshold) != (second.p_value < threshold)))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    upstream = []
    for scenario in ("strict", "biological", "nuisance"):
        for name, extra, columns in (("observed.tsv", [], ["p_value", "statistic", "mean_difference_norm", "converged", "n_subjects"]), ("null.tsv.gz", ["replicate"], ["p_value", "n_subjects"])):
            first, old = load(args.original / scenario, name)
            second, new = load(args.candidate / scenario, name)
            old, new = normalize_settings(old), normalize_settings(new)
            if old.pop("scalar_fast"):
                raise ValueError("original must use the generic scalar reference")
            new_fast = new.pop("scalar_fast", False)
            if old != new or not new_fast:
                raise ValueError("only the explicitly enabled scalar implementation may differ")
            row = compare_tables(first, second, ["test_id", "draw", "strategy", *extra], columns)
            rows.append(dict(scenario=scenario, table=name, **row))
        first, _ = load(args.original / scenario, "subject_null_diagnostics.tsv.gz")
        second, _ = load(args.candidate / scenario, "subject_null_diagnostics.tsv.gz")
        keys = ["test_id", "draw", "subject"]
        first, second = first.set_index(keys, verify_integrity=True).sort_index(), second.set_index(keys, verify_integrity=True).sort_index()
        if not first.index.equals(second.index):
            raise ValueError("upstream subject fit identities changed")
        for column in ("score", "information", "null_inclusion", "null_event_mass_min", "null_event_mass_max", "pooled_refit_inclusion", "pooled_refit_event_mass"):
            x, y = first[column].to_numpy(float), second[column].to_numpy(float)
            delta = np.abs(x - y)
            upstream.append(dict(scenario=scenario, component=column, subject_fits=len(first), not_bit_identical=int(np.count_nonzero(x != y)), max_absolute_difference=float(delta.max()), max_scaled_difference=float(np.max(delta / np.maximum(np.abs(x), 1e-8)))))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_dir / "equivalence.tsv", sep="\t", index=False)
    pd.DataFrame(upstream).to_csv(args.output_dir / "upstream_fits.tsv", sep="\t", index=False)
    passed = all(all(value for key, value in row.items() if key.endswith("_matches")) and all(value == 0 for key, value in row.items() if key.startswith("p_decisions_differ")) for row in rows)
    (args.output_dir / "manifest.json").write_text(json.dumps(dict(original=str(args.original), candidate=str(args.candidate), passed=passed, tolerance=dict(relative=3e-6, absolute=1e-8), model_settings="identical except scalar implementation; all trial identities and availability preserved", limitation="fresh upstream count/null fitting, not fixed score arrays; consult upstream fit differences and separate same-score controls before attribution; not new power or LR performance"), indent=2) + "\n")
    print(f"complete count-null implementation equivalence, passed={passed}", flush=True)
    if not passed:
        raise ValueError("end-to-end rerun differs; inspect upstream fitting and same-score controls before attributing differences")


if __name__ == "__main__":
    main()
