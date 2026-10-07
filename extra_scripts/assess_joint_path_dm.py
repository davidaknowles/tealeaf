#!/usr/bin/env python3
"""Assess joint local-path DM tests on frozen comparator-matched families."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from extra_scripts.assess_paired_inference_audit import split_assessment, long_read_assessment
from extra_scripts.evaluate_suppa2_statistics import normalize_pairs
from extra_scripts.summarize_path_reporting_omnibus import load_shards, calibrate_omnibus
from tealeaf.sc.replication_audit import coverage_correlation


def paired_table(table):
    """Use each model's standardized usage effects, not another variant's fits."""
    table = normalize_pairs(table)
    table.loc[~table.converged.astype(str).str.lower().eq("true"), "p_value"] = 1.
    table["converged"] = table.converged.astype(str).str.lower().eq("true")
    table["feature_id"] = table.block_id
    table["coverage"] = table.median_gene_umis
    table["published_q"] = table.fdr
    table["effect_vector"] = table.adjusted_effects.map(lambda value: (np.asarray(json.loads(value))[1] - np.asarray(json.loads(value))[0]).tolist())
    table["effect_features"] = table.path_signatures.map(lambda value: [json.dumps(item, sort_keys=True) for item in json.loads(value)])
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--split-only", action="store_true", help="Assess completed subject halves before the full-data LR fits finish.")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tables, audit, provenance = [], [], []
    backend = None
    for fold in ((0, 1) if args.split_only else (0, 1, "full")):
        root = args.cache / ("pairwise_full" if fold == "full" else f"pairwise_fold{fold}")
        paths = sorted(root.glob("shard_*/settings.json"))
        if len(paths) != 32 or {path.parent.name for path in paths} != {f"shard_{i}" for i in range(32)}:
            raise ValueError("expected 32 completed contiguous shards")
        settings = [json.loads(path.read_text()) for path in paths]
        method = "cox_reid" if settings[0].get("cox_reid_only", False) else "ml"
        if backend is not None and backend != method:
            raise ValueError("different dispersion methods across subject folds")
        backend = method
        if any(item["candidate_settings"] != settings[0]["candidate_settings"] or not (item.get("joint_dm_only") or item.get("cox_reid_only")) or item.get("mode") != "pairwise" or item.get("cox_reid_only", False) != (method == "cox_reid") for item in settings):
            raise ValueError("inconsistent joint DM paired cohort")
        cohort = settings[0]["candidate_settings"]
        if cohort.get("min_gene_umis") != 25 or cohort.get("subject_fold") != (None if fold == "full" else fold):
            raise ValueError("must assess matched production cohorts at minimum 25 gene counts")
        observed = load_shards(root, expected=32)
        null = load_shards(root, "null.tsv.gz", expected=32)
        tested, _, held = calibrate_omnibus(observed, null)
        tables.append(paired_table(tested))
        audit.append(held.assign(fold=fold))
        provenance.append({"fold": fold, "candidate_settings": cohort, "failures": sum(item["n_failures"] for item in settings)})
    coverage = []
    for strategy in tables[0].strategy.unique():
        output = args.output_dir / strategy.replace(" ", "_")
        output.mkdir(parents=True, exist_ok=True)
        selected = [frame.loc[frame.strategy.eq(strategy)].copy() for frame in tables]
        split_assessment(selected[:2], repo, output, strategy)
        if not args.split_only:
            effects = selected[2].set_index("test_id").effect_vector.to_dict()
            long_read_assessment(selected[2], None, args.matrix_dir, args.gtf, args.block_cache, output, strategy, effect_vectors=effects)
        for fold, frame in zip((0, 1) if args.split_only else (0, 1, "full"), selected):
            frame.to_csv(output / f"tests_{fold}.tsv.gz", sep="\t", index=False, na_rep="NA")
            for column in ("p_value", "raw_p_value"):
                coverage.append({"strategy": strategy, "fold": fold, "p_column": column, **coverage_correlation(frame[column], frame.coverage)})
    pd.concat(audit, ignore_index=True).to_csv(args.output_dir / "held_null_summary.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(args.output_dir / "coverage_correlations.tsv", sep="\t", index=False)
    testing = "joint subject-blocked DM LRT; concentration re-estimated under null and alternative" if backend == "ml" else "joint subject-blocked Cox-Reid grid-selected fixed-precision DM, chi-square and heuristic F tails; precision reselected under every permutation design"
    (args.output_dir / "manifest.json").write_text(json.dumps({"cohorts": provenance, "testing": testing + "; fractional covariance-matched path counts", "calibration": "32 within-subject label permutations, leave-own-test-out pooled calibration; another 32 families held for assessment", "reporting": "each model's means standardized over the same subjects, no cross-variant substitution", "split": "frozen Table 1 comparator-matched gene/pair universes with the original Simes/conjunction/BH; failed tests p=1", "split_only": args.split_only, "long_read": "not assessed by this invocation" if args.split_only else "all tested pairs freshly remapped, no significance or historical discovery-list restriction", "production_changes": False}, indent=2) + "\n")


if __name__ == "__main__":
    main()
