#!/bin/bash
# Paper figures with every method's coverage-null-called tests removed (targeted rule; SUPPA2 not simulated).
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
source ~/venv/kinnex_py/bin/activate
cd ~/projects/tealeaf; export PYTHONPATH=$PWD
R=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark
F=analyses/coverage_null_filter/brain
T=analyses/tilgner_junction_replication
COMMON="--tealeaf-events analyses/native_align/events_v4_groupfix/test_events.tsv.gz --suppa2 analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz --suppa2-primer-aware analyses/comparator_suppa_rmats/suppa2_primer_aware/full_data_tests.tsv.gz --rmats analyses/comparator_suppa_rmats/rmats_paired_full_tests.tsv.gz --majiq-tests $R/majiq_build_min3_cov3/full_tests/normalized --majiq-coverage-null $F/majiq_coverage_null_targeted.parquet --coverage-flagged $F/flagged_full.tsv.gz"
if [ "${1:-}" != "lr" ]; then  # pass "lr" to redo only the long-read curves
python extra_scripts/plot_discoveries_by_event.py $COMMON --junction-table $T/all_tested_comprehensive/junction_replication.tsv.gz --output docs/figures/discoveries_by_event.pdf
python extra_scripts/plot_discoveries_by_event.py $COMMON --junction-table $R/comparison_corrected/all_tests.tsv.gz --leafcutter-counts $R/leafcutter/clustering/benchmark_perind_numers.counts.gz --junction-annotation $R/pseudobulk_junctions_junctions.tsv.gz --contrasts all --output docs/figures/discoveries_by_event_all.pdf
fi
# long-read curves: filtered ranks of coverage_filter_eval.py, SUPPA2 rows unchanged (not simulated)
python - <<'PY'
import pandas as pd
F = "analyses/coverage_null_filter/brain"
ranks = pd.read_csv(f"{F}/lr_ranks.tsv.gz", sep="\t", low_memory=False)
ranks = ranks.loc[ranks.table.eq("filtered") & ranks["rank"].le(200)]
base = pd.read_csv("analyses/tilgner_junction_replication/sign_agreement_by_rank_comprehensive.tsv", sep="\t", low_memory=False)
comparators = pd.concat([base.loc[base.method.str.startswith("SUPPA2") | base.method.str.contains("hybrid")], ranks.loc[ranks.method.ne("Tealeaf")].reindex(columns=base.columns)], ignore_index=True)
comparators.to_csv(f"{F}/sign_agreement_by_rank_filtered.tsv", sep="\t", index=False, na_rep="NA")
ranks.loc[ranks.method.eq("Tealeaf")].assign(scope="all tested", ranking="continuous raw").to_csv(f"{F}/tealeaf_lr_rank_filtered.tsv.gz", sep="\t", index=False)
print(comparators.groupby("method").size())
PY
python extra_scripts/plot_lr_agreement_curves.py --tealeaf $F/tealeaf_lr_rank_filtered.tsv.gz --comparators $F/sign_agreement_by_rank_filtered.tsv --output docs/figures/lr_agreement.pdf
