#!/usr/bin/env python3
"""Audit event tests on the original matched SUPPA2 split universe.

Observed effect sizes identify the original TPM input rather than changing
quantification while comparing statistics. Moderation is estimated over all
eligible catalogue events within a contrast, before restricting the benchmark.
Subject-synchronized sign flips retain dependence across events and contrasts.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from extra_scripts.audit_suppa2_event_null import synchronized_sign
from extra_scripts.compare_suppa2_primer_aware_tealeaf import load_merged_tealeaf
from extra_scripts.evaluate_suppa2_statistics import normalize_pairs
from extra_scripts.run_suppa2_full_data_comparison import TEST_METHODS, event_test_differences, event_test_pvalues
from extra_scripts.run_suppa2_split_data_comparison import event_psi, load_catalog
from tealeaf.sc import differential
from tealeaf.sc.event_tests import paired_signed_rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--null-replicates", type=int, default=32)
    parser.add_argument("--primer-psi-dir", type=Path, help="Use cached primer-aware PSI instead of transcript TPM.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    repro = args.run_root / "junction_benchmark/reproducibility"
    archive = args.repo_root / "analyses/comparator_suppa_rmats" / ("suppa2_primer_aware" if args.primer_psi_dir else "suppa2")
    catalog = load_catalog((args.primer_psi_dir or archive) / "event_catalog.tsv.gz")
    old = [normalize_pairs(pd.read_csv(archive / f"split_data_fold{k}_tests.tsv.gz", sep="\t")) for k in (0, 1)]
    references = []
    for k in (0, 1):
        table = load_merged_tealeaf(repro / f"fold{k}/tealeaf_paired_path_total_a32_production/paired_path.tsv")
        references.append(normalize_pairs(table))
    shared = set.intersection(*(set(zip(table.gene_id, table.pair_id)) for table in old + references))
    old_fold = old[args.fold].set_index(["contrast_id", "feature_id"])
    contrasts = json.loads((repro / f"fold{args.fold}/contrasts.json").read_text())
    # Verify identical PSI effects and subject counts against the saved run.
    candidates = ("cluster_dx_mouse_min20_umi1e5_em", "cluster_dx_mouse", "weighted_parsimony_min20_umi1e5_em")
    best = None
    checks = []
    if args.primer_psi_dir:
        samples = json.loads((args.primer_psi_dir / "samples.json").read_text())
        lookup = {sample: i for i, sample in enumerate(samples)}
        psi = np.load(args.primer_psi_dir / f"psi_fold{args.fold}.npz")["psi"]
        best = (0., "primer_aware", psi, lookup)
        candidates = ()
    for prefix in candidates:
        matrix_path = args.run_root / f"{prefix}_pseudo_spliced_TPM.npz"
        row_path = args.run_root / f"{prefix}_pseudo_rows.txt"
        col_path = args.run_root / f"{prefix}_pseudo_spliced_cols.txt"
        if not all(path.exists() for path in (matrix_path, row_path, col_path)):
            continue
        rows = row_path.read_text().splitlines()
        lookup = {row.split("__")[-1] + "__" + row.split("__")[0]: i for i, row in enumerate(rows)}
        psi = event_psi(catalog, sparse.load_npz(matrix_path), col_path.read_text().splitlines())
        errors = []
        count_errors = []
        for contrast in contrasts[:5]:
            pairs = [(lookup[a], lookup[b]) for a, b in zip(contrast["samples_a"], contrast["samples_b"]) if a in lookup and b in lookup]
            if not pairs:
                continue
            first, second = psi[:, [a for a, _ in pairs]], psi[:, [b for _, b in pairs]]
            valid = np.isfinite(first) & np.isfinite(second)
            effects = np.divide(np.where(valid, second - first, 0).sum(axis=1), valid.sum(axis=1), out=np.full(len(catalog), np.nan), where=valid.any(axis=1))
            previous = old_fold.loc[contrast["contrast_id"]].reindex(catalog.feature_id)
            finite = np.isfinite(previous.effect_size) & np.isfinite(effects)
            errors.extend(np.abs(previous.effect_size.to_numpy()[finite] - effects[finite]))
            count_errors.extend(np.abs(previous.n_subjects.to_numpy()[finite] - valid.sum(axis=1)[finite]))
        error = float(np.max(errors)) if len(errors) else np.inf
        checks.append({"input": prefix, "maximum_effect_error": error, "maximum_subject_count_error": float(np.max(count_errors)) if len(count_errors) else np.inf})
        if best is None or error < best[0]:
            best = (error, prefix, psi, lookup)
    pd.DataFrame(checks).to_csv(args.output_dir / f"fold{args.fold}_input_checks.tsv", sep="\t", index=False)
    if best is None or best[0] > 1e-10:
        raise ValueError(f"cannot reproduce archived PSI, checks={checks}")
    _, prefix, psi, lookup = best
    observed, nulls, diagnostics = [], [], []
    genes = catalog.gene_id.str.split(".").str[0].to_numpy()
    for contrast in contrasts:
        pair_id = "||".join(sorted((contrast["level_a"], contrast["level_b"])))
        matched = np.array([(gene, pair_id) in shared for gene in genes])
        if not matched.any():
            continue
        pairs = [(lookup[a], lookup[b], a.split("__", 1)[0]) for a, b in zip(contrast["samples_a"], contrast["samples_b"]) if a in lookup and b in lookup]
        if len(pairs) < 8:
            continue
        first, second = psi[:, [a for a, _, _ in pairs]], psi[:, [b for _, b, _ in pairs]]
        valid = np.isfinite(first) & np.isfinite(second)
        enough = valid.sum(axis=1) >= 8
        archived_features = set(old_fold.loc[contrast["contrast_id"]].index)
        archived = catalog.feature_id.isin(archived_features).to_numpy()
        selected = np.flatnonzero(matched & enough & archived)
        if not len(selected):
            continue
        differences = second - first
        for index in selected:
            previous = old_fold.loc[(contrast["contrast_id"], catalog.iloc[index].feature_id)]
            effect = differences[index, valid[index]].mean()
            if abs(effect - previous.effect_size) > 1e-10 or valid[index].sum() != previous.n_subjects:
                raise ValueError("quantification changed relative to the original split comparison")
            diagnostics.append({"fold": args.fold, "contrast_id": contrast["contrast_id"], "gene_id": genes[index], "feature_id": catalog.iloc[index].feature_id, "n_subjects": valid[index].sum(), "nonzero_subjects": np.sum(valid[index] & (differences[index] != 0)), "effect_size": effect})
        for method in TEST_METHODS:
            # Ineligible events must not contribute to the moderation prior.
            test_valid = valid & enough[:, None]
            if method.startswith("wilcoxon"):
                values = event_test_pvalues(first[selected], second[selected], test_valid[selected], method)
            else:
                values = event_test_pvalues(first, second, test_valid, method)[selected]
            meta = catalog.iloc[selected][["gene_id", "feature_id", "event_type"]].copy()
            meta["gene_id"] = genes[selected]
            meta["level_a"], meta["level_b"] = contrast["level_a"], contrast["level_b"]
            meta["contrast_id"], meta["method"], meta["fold"] = contrast["contrast_id"], method, args.fold
            meta["p_value"] = values
            observed.append(meta.copy())
            transformed = event_test_differences(first, second, method)
            for replicate in range(args.null_replicates):
                signs = np.array([synchronized_sign(20261002, replicate, subject) for _, _, subject in pairs])
                flipped = transformed * signs
                if method.startswith("wilcoxon"):
                    values = paired_signed_rank(flipped[selected], test_valid[selected], exact=method == "wilcoxon_exact")
                else:
                    values = differential.vectorized_paired_t_pvalues(flipped, test_valid, moderate=method.endswith("moderated_t"))[selected]
                local = meta.copy()
                local["replicate"], local["p_value"] = replicate, values
                nulls.append(local)
        print(f"fold={args.fold} contrast={contrast['contrast_id']} matched_events={len(selected)}", flush=True)
    pd.concat(observed).to_csv(args.output_dir / f"fold{args.fold}_tests.tsv.gz", sep="\t", index=False)
    pd.concat(nulls).to_csv(args.output_dir / f"fold{args.fold}_null.tsv.gz", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(args.output_dir / f"fold{args.fold}_support.tsv.gz", sep="\t", index=False)
    print(f"matched_pairs={len(shared)} reproduced_input={prefix}", flush=True)


if __name__ == "__main__":
    main()
