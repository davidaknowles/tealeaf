"""Collate the declared real-count diagnostics, never a selected-family benchmark."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def collate_source(root):
    settings = json.loads((root / "shard_0/settings.json").read_text())
    frames = []
    for shard in range(settings["shard_count"]):
        directory = root / f"shard_{shard}"
        if json.loads((directory / "settings.json").read_text()) != settings:
            raise ValueError(f"inconsistent settings in {directory}")
        frames.append(pd.read_csv(directory / "tests.tsv", sep="\t"))
    tests = pd.concat(frames, ignore_index=True)
    if len(tests) != settings["n_requested"] or tests.duplicated(["scope", "test_id"]).any():
        raise ValueError("missing or duplicate declared diagnostic hypotheses")
    if not np.isfinite(tests.p_value).all() or not tests.p_value.between(0, 1).all():
        raise ValueError("invalid diagnostic p-values")
    if not tests.loc[~tests.converged, "p_value"].eq(1).all():
        raise ValueError("incomplete fits must remain p1")
    return tests, settings


def summarize(tests):
    rows = []
    for (source, scope), frame in tests.groupby(["source", "scope"], sort=False):
        reporting = frame.report_complete.eq(True)
        inference = frame.converged.astype(bool)
        row = dict(source=source, scope=scope, n_requested=len(frame), n_complete_inference=int(inference.sum()), n_complete_reporting=int(reporting.sum()), n_raw_f_p05=int((frame.p_value < .05).sum()), n_raw_chi_p05=int((inference & (frame.chi_square_p_value < .05)).sum()), n_boundary_pooled=int(((frame.pooled_inclusion < .01) | (frame.pooled_inclusion > .99)).sum()), n_negligible_pooled_event_mass=int((frame.pooled_event_mass < .01).sum()))
        if frame.long_read_effect.notna().any():
            # External exact zeros remain nonagreements. Failures remain in
            # the declared denominator and are also reported separately.
            row.update(n_report_agrees=int((reporting & (frame.report_effect * frame.long_read_effect > 0)).sum()), n_score_agrees=int((inference & (frame.score_effect * frame.long_read_effect > 0)).sum()), n_external_nonzero=int(frame.long_read_effect.ne(0).sum()))
        rows.append(row)
    return pd.DataFrame(rows)


def validate_source_identities(parts, allow_source_specific_random_panels=False):
    """Require matched native panels, explicitly label different random panels."""
    identities = [set(zip(part.scope, part.test_id)) for part in parts]
    if identities[0] != identities[1]:
        if not allow_source_specific_random_panels:
            raise ValueError("source controls do not have the same declared requests")
        native = [{key for key in identity if key[0] == "native top100 diagnostic, not an inference family"} for identity in identities]
        random = [{key for key in identity if key[0] == "fixed random real-data diagnostic"} for identity in identities]
        if native[0] != native[1] or len(native[0]) != 100 or any(len(panel) != 32 for panel in random) or any(len(identity) != 132 for identity in identities):
            raise ValueError("source-specific diagnostics still require 32 frozen random and 100 identical native requests")
    return identities[0] == identities[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-source-specific-random-panels", action="store_true", help="Keep different frozen random diagnostic panels explicit; native-leading identities must still match.")
    args = parser.parse_args()
    sources = ("parsimony_binary", "original_binary")
    parts, settings = [], {}
    for source in sources:
        table, settings[source] = collate_source(args.input_root / source)
        parts.append(table)
    tests = pd.concat(parts, ignore_index=True)
    same_requests = validate_source_identities(parts, args.allow_source_specific_random_panels)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tests.to_csv(args.output_dir / "tests.tsv.gz", sep="\t", index=False, na_rep="NA")
    summary = summarize(tests)
    summary.to_csv(args.output_dir / "summary.tsv", sep="\t", index=False, na_rep="NA")
    failures = tests.loc[~tests.converged].groupby(["source", "scope", "error"], dropna=False).size().reset_index(name="n")
    failures.to_csv(args.output_dir / "failure_summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({"settings": settings, "selection": "132 diagnostic requests/source, 32 fixed random identities and the 100 native-SUPPA2 LR-leading associations; not a full screened hypothesis family", "same_random_panel_across_sources": same_requests, "source_comparison": "native panel identities match; different random panels are not a matched source-only ablation" if not same_requests else "same declared diagnostic identities", "ranking_claim": "none; no A100 or selected-family FDR inference", "production_changes": False}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    print(failures.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
