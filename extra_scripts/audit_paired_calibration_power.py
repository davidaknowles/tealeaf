"""Diagnose native-versus-pooled calibration power without adopting native tails."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_paired_inference_audit import fit_table, split_assessment
from tealeaf.sc.ds_benchmark import benjamini_hochberg
from tealeaf.sc.replication_audit import ranked_direction_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--models", nargs="+", default=["score", "score_free", "baseline", "subject"])
    parser.add_argument("--assessment-root", type=Path)
    parser.add_argument("--native-p-column", choices=("raw_p_value", "chi_square_p_value"), default="raw_p_value", help="Diagnostic analytic reference; this does not promote its tails to production.")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    rows = []
    assessment_root = args.assessment_root or repo / "analyses/paired_inference_production"
    for model in args.models:
        assessment = json.loads((assessment_root / model / "manifest.json").read_text())
        expected = {str(item["fold"]): item["candidate_settings"] for item in assessment["cohort_checks"]}
        for fold in (0, 1, "full"):
            summaries = list((args.cache / f"{model}_{fold}").glob("shard_*/summary.json"))
            if len(summaries) != 32 or any(json.loads(path.read_text()).get("candidate_settings") != expected[str(fold)] for path in summaries):
                raise ValueError("calibration audit cache differs from the complete production assessment cohort")
        folds = []
        for fold in (0, 1):
            table = fit_table(args.cache / f"{model}_{fold}" / "merged" / "paired_path.tsv")
            archived_calibrated_BH = int(table.fdr.le(.05).sum())
            table["fdr"] = benjamini_hochberg(table.p_value.to_numpy(float))
            good = table.converged.astype(str).str.lower().eq("true") & table.n_subjects.ge(4)
            native = table.copy()
            native["p_value"] = np.where(good, native[args.native_p_column], 1.)
            native["fdr"] = benjamini_hochberg(native.p_value.to_numpy(float))
            folds.append(native)
            rows.append({"model": model, "fold": fold, "n_requested": len(table), "n_converged": int(good.sum()), "native_BH": int(native.fdr.le(.05).sum()), "calibrated_BH": int(table.fdr.le(.05).sum()), "archived_calibrated_BH": archived_calibrated_BH, "minimum_native_p": native.p_value.min(), "minimum_calibrated_p": table.p_value.min(), "native_below_minimum_calibrated": int(native.p_value.lt(table.p_value.min()).sum())})
        output = args.output_dir / model
        output.mkdir(parents=True, exist_ok=True)
        split_assessment(folds, repo, output, f"{model}, native diagnostic")
        # Reuse the complete existing LR mapping, not the historical discoveries.
        mapped = pd.read_csv(assessment_root / model / "lr_mapping.tsv.gz", sep="\t")
        full = pd.read_csv(args.cache / f"{model}_full" / "merged" / "paired_path.tsv", sep="\t")
        paired_p = mapped[["test_id", "p_value"]].merge(full[["test_id", "p_value"]], on="test_id", validate="one_to_one", suffixes=("_mapped", "_fitted"))
        if len(paired_p) != len(mapped) or not np.allclose(paired_p.p_value_mapped, paired_p.p_value_fitted, rtol=1e-10, atol=1e-14):
            raise ValueError("LR mapping and fitted statistics are not the same cohort/model")
        good = full.converged.astype(str).str.lower().eq("true") & full.n_subjects.ge(4)
        full["p_value"] = np.where(good, full[args.native_p_column], 1.)
        full["fdr"] = benjamini_hochberg(full.p_value.to_numpy(float))
        mapped = mapped.drop(columns=["p_value", "fdr"]).merge(full[["test_id", "p_value", "fdr"]], on="test_id", validate="one_to_one")
        summaries = []
        for scope, subset in (("all tested", mapped), ("BH discoveries", mapped.loc[mapped.fdr.le(.05)])):
            ranked = subset.loc[subset.eligible].sort_values(["p_value", "statistic", "test_id"], ascending=[True, False, True], kind="stable").copy()
            ranked["rank"] = np.arange(1, len(ranked) + 1)
            ranked["method"] = model
            summaries.extend({**row, "scope": scope, "ranking": "native diagnostic"} for row in ranked_direction_summary(ranked))
        pd.DataFrame(summaries).to_csv(output / "lr_rank_summary.tsv", sep="\t", index=False, na_rep="NA")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(args.output_dir / "calibration_power_summary.tsv", sep="\t", index=False)
    (args.output_dir / "manifest.json").write_text(json.dumps({"scope": "Native-tail diagnostic using unchanged fitted statistics and complete mapped LR families", "native_p_column": args.native_p_column, "cohort_guard": "Every inference shard exactly matches its prior production-cohort assessment; mapped calibrated p-values match the same fitted table", "production_changes": False, "limitation": "Count-null nominal calibration does not validate extreme tails, biological overdispersion or family FDR", "selection": "Frozen comparator-matched gene/pair universes and existing all-tested LR mappings"}, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
