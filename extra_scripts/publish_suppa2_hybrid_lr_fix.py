"""Publish the nuisance-profiled hybrid while retaining its fixed-mass reference."""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def reporting_table(table):
    """Use descriptive PSI effects while preserving the ILR testing estimate."""
    table = table.copy()
    if "test_ilr_effect_size" not in table:
        table["test_ilr_effect_size"] = table.effect_size
    table["effect_size"] = table.report_psi_effect
    table["effect_coordinate"] = "psi"
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--debug-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--native-tests", type=Path, default=Path("analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz"))
    parser.add_argument("--initial-audit-dir", type=Path, default=Path(".cache/suppa2_statistics/hybrid_lr_debug/initial"))
    args = parser.parse_args()
    base = args.debug_dir
    audit = args.output_dir / "lr_model_audit"
    reference = audit / "fixed_mass_reference"
    reference.mkdir(parents=True, exist_ok=True)
    full = args.output_dir / "full_data"
    for name in ("full_data_tests.tsv.gz", "full_data_tilgner.tsv.gz", "full_data_tilgner_summary.tsv", "full_data_tilgner_manifest.json", "full_data_summary.tsv", "null_families.tsv"):
        if not (reference / name).exists():
            shutil.copyfile(full / name, reference / name)
    old = pd.read_csv(base / "fixed64_testonly/merged/paired_path.tsv", sep="\t").set_index("test_id")
    control = pd.read_csv(base / "fixed64_report1/merged/paired_path.tsv", sep="\t").set_index("test_id").loc[old.index]
    tolerances = {"raw_p_value": (1e-8, 1e-10), "statistic": (.01, 1e-4), "mean_difference_norm": (1e-8, 1e-7)}
    checks = []
    for column, (relative, absolute) in tolerances.items():
        np.testing.assert_allclose(old[column], control[column], rtol=relative, atol=absolute)
        checks.append({"field": column, "maximum_absolute_difference": float((old[column] - control[column]).abs().max()), "relative_tolerance": relative, "absolute_tolerance": absolute})
    # A floating-point tie can move a calibrated tail by one null-pool count.
    stratum_size = control.groupby("calibration_stratum")["p_value"].transform("size")
    quantum = 1 / (1 + 32 * (stratum_size - 1))
    difference = (old.p_value - control.p_value).abs()
    assert (difference <= quantum + 1e-12).all(), "reporting control changed calibrated tails beyond one empirical count"
    print("Calibrated control tails differing by one empirical count", int(difference.gt(1e-12).sum()))
    pd.DataFrame(checks).to_csv(audit / "reporting_control_precision.tsv", sep="\t", index=False)
    np.testing.assert_array_equal(old.fdr.le(.05), control.fdr.le(.05))
    ranking = ["p_value", "raw_p_value", "statistic", "test_id"]
    order = [True, True, False, True]
    np.testing.assert_array_equal(old.reset_index().sort_values(ranking, ascending=order).head(200).test_id, control.reset_index().sort_values(ranking, ascending=order).head(200).test_id)
    updated = reporting_table(pd.read_csv(base / "profile64_report1/merged/paired_path.tsv", sep="\t"))
    usage = pd.concat([pd.read_csv(path, sep="\t") for path in (base / "profile64_report1/shards").glob("shard_*/path_usage.tsv.gz")], ignore_index=True)
    report_n = usage.groupby("test_id").subject.nunique()
    updated["report_n_subjects"] = updated.test_id.map(report_n).fillna(0).astype(int)
    print("Reporting fits with different subject count", int(updated.report_n_subjects.ne(updated.n_subjects).sum()))
    updated.to_csv(full / "full_data_tests.tsv.gz", sep="\t", index=False)
    for source, target in (("summary.tsv", "full_data_summary.tsv"), ("null_families.tsv", "null_families.tsv")):
        shutil.copyfile(base / "profile64_report1/merged" / source, full / target)
    split = args.output_dir / "split_data"
    split.mkdir(exist_ok=True)
    for fold in (0, 1):
        reporting_table(pd.read_csv(base / f"profile32_split/fold{fold}_merged/paired_path.tsv", sep="\t")).to_csv(split / f"fold{fold}_tests.tsv.gz", sep="\t", index=False)
    for name in ("split_reproducibility.tsv", "event_diagnostics.tsv", "gene_pvalues.tsv.gz"):
        source = base / "profile32_split/evaluation" / name
        if name == "event_diagnostics.tsv":
            pd.read_csv(source, sep="\t").to_csv(split / name, sep="\t", index=False, na_rep="NA")
        else:
            shutil.copyfile(source, split / name)
    command = [sys.executable, str(Path(__file__).with_name("evaluate_hybrid_lr_variants.py")), "--mapping", str(reference / "full_data_tilgner.tsv.gz"), "--native-tests", str(args.native_tests), "--output-dir", str(audit), "--variant", "original", str(reference / "full_data_tests.tsv.gz")]
    for tag in ("fixed64_testonly", "fixed64_report1", "profile64_report1", "fixed4_report1", "profile4_report1", "fixed64_report1_psi", "profile64_report1_psi"):
        command.extend(["--variant", tag, str(base / tag / "merged/paired_path.tsv")])
    subprocess.run(command, check=True)
    for name in ("ranking_direction_audit.tsv", "event_universe.tsv", "catalog_universe.tsv"):
        shutil.copyfile(args.initial_audit_dir / name, audit / name)


if __name__ == "__main__":
    main()
