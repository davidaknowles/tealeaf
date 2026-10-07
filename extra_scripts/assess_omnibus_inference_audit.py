#!/usr/bin/env python3
"""Matched split and all-tested LR assessment of experimental EC omnibus tests."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.summarize_path_reporting_omnibus import load_shards, calibrate_omnibus, split_omnibus_summary, lr_omnibus_summary
from tealeaf.sc.ds_benchmark import benjamini_hochberg


def complete_failed_tests(frame, reference):
    """Keep failed tests at p=1 on a fixed eligible universe, not complete cases."""
    output = []
    for strategy, local in frame.groupby("strategy"):
        retained = local.loc[local.block_id.isin(reference.block_id)].copy()
        missing = reference.loc[~reference.block_id.isin(retained.block_id)].copy()
        missing["strategy"] = strategy
        for field in ("p_value", "raw_p_value", "fdr"):
            missing[field] = 1.
        missing["statistic"] = 0.
        missing["converged"] = False
        missing["adjusted_effects"] = missing.adjusted_effects.map(lambda value: json.dumps(np.full(np.asarray(json.loads(value)).shape, np.nan).tolist()))
        retained["fit_available"] = retained.converged.astype(str).str.lower().eq("true") if "converged" in retained else True
        retained.loc[~retained.fit_available, ["p_value", "raw_p_value"]] = 1.
        retained.loc[~retained.fit_available, "statistic"] = 0.
        retained.loc[~retained.fit_available, "adjusted_effects"] = retained.loc[~retained.fit_available, "adjusted_effects"].map(lambda value: json.dumps(np.full(np.asarray(json.loads(value)).shape, np.nan).tolist()))
        missing["fit_available"] = False
        completed = pd.concat([retained, missing], ignore_index=True)
        completed["fdr"] = benjamini_hochberg(completed.p_value.to_numpy(float))
        output.append(completed)
    if not output:
        raise ValueError("no completed omnibus strategies")
    return pd.concat(output, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--control-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--control-fold-shards", type=int, default=16)
    parser.add_argument("--split-only", action="store_true", help="Assess completed subject halves before full-data fits finish.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    folds, audits, provenance = [], [], []
    full, full_control = None, None
    for fold in ((0, 1) if args.split_only else (0, 1, "full")):
        folder = "omnibus_full" if fold == "full" else f"omnibus_fold{fold}"
        experimental = args.cache / folder
        if fold == "full" and not experimental.exists():
            experimental = args.cache / "full"
        roots = [experimental, args.control_cache / folder]
        settings = []
        for index, root in enumerate(roots):
            files = sorted(root.glob("shard_*/settings.json"))
            expected = args.control_fold_shards if index == 1 and fold != "full" else 32
            if len(files) != expected or {path.parent.name for path in files} != {f"shard_{number}" for number in range(expected)}:
                raise ValueError(f"expected {expected} completed contiguous shards in {root}")
            local = [json.loads(path.read_text()) for path in files]
            if any(item["candidate_settings"] != local[0]["candidate_settings"] for item in local):
                raise ValueError("candidate settings differ within an omnibus run")
            if local[0]["candidate_settings"].get("min_gene_umis") != 25:
                raise ValueError("omnibus audit must use the production minimum coverage")
            settings.append(local[0])
        if settings[0]["candidate_settings"] != settings[1]["candidate_settings"]:
            raise ValueError(f"different omnibus candidate cohorts in {folder}")
        frames = []
        for index, root in enumerate(roots):
            expected = args.control_fold_shards if index == 1 and fold != "full" else 32
            observed, null = load_shards(root, expected=expected), load_shards(root, "null.tsv.gz", expected=expected)
            # Older controls used the documented default 32; full-data controls
            # are rerun at 64 so both the prior center and score comparisons match.
            strengths = settings[len(frames)].get("test_concentration", 32.)
            joint_dm = settings[len(frames)].get("joint_dm_only", False) or settings[len(frames)].get("cox_reid_only", False)
            if not joint_dm and strengths != (64. if fold == "full" else 32.):
                raise ValueError("testing concentrations differ from the prespecified matched control")
            calibrated, _, audit = calibrate_omnibus(observed, null)
            frames.append(calibrated)
            audit["variant"] = args.name if len(frames) == 1 else "uniform control"
            audit["fold"] = fold
            audits.append(audit)
        new, control = frames
        reference = control.loc[control.strategy.str.startswith("null-variance Wald")].drop_duplicates("block_id")
        # Split eligibility is established by the existing control, before looking
        # at experimental p-values. Count nulls separately screen calibration.
        new = complete_failed_tests(new, reference)
        control = complete_failed_tests(control, reference)
        new["strategy"] += f", {args.name}"
        control["strategy"] += ", uniform control"
        combined = pd.concat([new, control], ignore_index=True)
        combined.to_csv(args.output_dir / f"omnibus_{fold}_tests.tsv.gz", sep="\t", index=False, na_rep="NA")
        provenance.append({"fold": fold, "candidate_settings": settings[0]["candidate_settings"], "test_concentration": settings[0].get("test_concentration", 32), "quantification_concentrations": settings[0].get("quantification_concentrations"), "control_testing_concentration": 64 if fold == "full" else 32, "n_control_blocks": len(reference), "n_missing_by_strategy": new.groupby("strategy").fit_available.apply(lambda values: int((~values).sum())).to_dict()})
        if fold == "full":
            full, full_control = new, reference
        else:
            folds.append(combined)
    summary = split_omnibus_summary(folds, args.output_dir)
    pd.concat(audits, ignore_index=True).to_csv(args.output_dir / "held_null_summary.tsv", sep="\t", index=False, na_rep="NA")
    if args.split_only:
        (args.output_dir / "manifest.json").write_text(json.dumps({"variant": args.name, "cohorts": provenance, "selection": "fixed eligible control block universe; missing experimental fits p=1", "split_only": True, "long_read": "not assessed by this invocation", "production_changes": False}, indent=2) + "\n")
        print(summary.to_string(index=False), flush=True)
        return
    # No pooled reporting or old paired discovery list restricts the LR universe.
    # Select the largest source effect among supported types before LR checks.
    if settings[0].get("joint_dm_only", False) or settings[0].get("cox_reid_only", False):
        report = full[["block_id", "levels", "path_signatures", "adjusted_effects", "strategy", "converged"]].copy()
        report["strategy"] = "standardized means, " + report.strategy
    else:
        report = full.drop_duplicates("block_id")[["block_id", "levels", "path_signatures", "adjusted_effects"]].assign(strategy="experimental subject-blocked A1", converged=True)
    # The reusable mapper historically labels its reference "archived
    # production". Here it is an explicitly refitted, matched uniform control.
    lr_omnibus_summary(full, report, full_control, args.run_root / "differential/gencode_vM32_splice_blocks.json.gz", args.matrix_dir, args.gtf, args.output_dir)
    for filename in ("lr_omnibus_summary.tsv", "lr_omnibus_rank.tsv.gz"):
        path = args.output_dir / filename
        table = pd.read_csv(path, sep="\t", low_memory=False)
        table["strategy"] = table.strategy.replace({"archived production omnibus": "refitted uniform control omnibus"})
        if "method" in table:
            table["method"] = table.method.str.replace("archived production omnibus", "refitted uniform control omnibus", regex=False)
        table.to_csv(path, sep="\t", index=False, na_rep="NA")
    reporting_description = "each variant's standardized model means" if settings[0].get("joint_dm_only", False) or settings[0].get("cox_reid_only", False) else "source A1 means"
    (args.output_dir / "manifest.json").write_text(json.dumps({"variant": args.name, "cohorts": provenance, "selection": "fixed eligible control block universe; missing experimental fits have p=1 and unavailable direction", "LR": f"all tested control-universe blocks, remapped without a paired-discovery screen; largest effect contrast using {reporting_description}", "interpretation": "exploratory conditional count-null, split and LR audits, not certification of FDR or an adopted production change"}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
