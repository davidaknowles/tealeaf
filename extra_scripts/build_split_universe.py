#!/usr/bin/env python3
"""Comparator-defined split universes for the subject-split benchmark.

For each comparator (LeafCutter, MAJIQ Heterogen, scQuint, SUPPA2, primer-aware
SUPPA2, the SUPPA2/Tealeaf hybrid and paired rMATS), the universe is the set of
(gene, cell-type pair) keys with at least one finite feature p-value in both
subject halves, taken from the comparator's own fold outputs. No Tealeaf result
enters. The output holds every comparator feature test on its universe, one
row per (comparison, fold, feature, pair), with columns gene_id, pair_id,
feature_id, p_value, comparison, method and fold; SUPPA2 tests are the exact
paired signed-rank tests on the full event catalogue. Its gene_id and pair_id
columns also serve as the --universe of run_local_path_tests.py.
"""

import argparse
from pathlib import Path

import pandas as pd

from extra_scripts.audit_split_coverage_direction import load_comparators


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--suppa2-split", default="split_data_exact", help="SUPPA2 fold outputs <name>_fold{k}_tests.tsv.gz in analyses/comparator_suppa_rmats/{suppa2,suppa2_primer_aware}; split_data_exact holds exact paired signed-rank tests on the full catalogue")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    repro = args.run_root / "junction_benchmark/reproducibility"
    folds = [load_comparators(repro, repo, fold, suppa2_split=args.suppa2_split) for fold in (0, 1)]
    rows, summary = [], []
    for comparison in sorted(set(folds[0]) & set(folds[1])):
        tables = [fold[comparison].loc[lambda t: t.gene_id.notna() & t.gene_id.ne("nan")] for fold in folds]
        keys = set(zip(tables[0].gene_id, tables[0].pair_id)) & set(zip(tables[1].gene_id, tables[1].pair_id))
        for fold, table in enumerate(tables):
            local = table.loc[[key in keys for key in zip(table.gene_id, table.pair_id)], ["gene_id", "pair_id", "feature_id", "p_value"]]
            rows.append(local.assign(comparison=comparison, method=comparison, fold=fold))
        summary.append({"comparison": comparison, "pairs": len(keys), "genes": len({gene for gene, _ in keys}), "cell_type_pairs": len({pair for _, pair in keys})})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(args.output, sep="\t", index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(args.output.with_name("summary.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
