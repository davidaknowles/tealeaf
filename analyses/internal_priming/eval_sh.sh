#!/bin/bash
# Coverage-filtered evaluation and discoveries with the shrunk-offset Tealeaf fits (brain v4sh, PBMC _sh) and LeafCutter min_coverage 5.
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD LEAF5=1 SUFFIX=_sh BRAIN=v4sh LR_MAPPING=.cache/internal_priming/compare/v4sh/lr_mapping.tsv.gz PBMC_LR_MAPPING=analyses/pbmc/GSE320507_lr_sh/lr_mapping.tsv.gz
PY=~/venv/cpa_blackwell/bin/python; F=analyses/coverage_null_filter_sh
$PY .cache/adhoc/coverage_flags_full.py brain $F/brain
$PY .cache/adhoc/coverage_flags_full.py pbmc $F/pbmc
$PY .cache/adhoc/coverage_filter_eval.py brain $F/brain
$PY .cache/adhoc/coverage_filter_eval.py pbmc $F/pbmc --skip-lr
$PY .cache/adhoc/pbmc_lr_filter.py $F/pbmc
source ~/venv/kinnex_py/bin/activate
R=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark; B=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark
A=analyses/comparator_suppa_rmats; C=$B/comparator_suppa_rmats
BR="--tealeaf-events analyses/native_align/events_v4sh/test_events.tsv.gz --suppa2 $A/suppa2/full_data_exact_tests.tsv.gz --suppa2-primer-aware $A/suppa2_primer_aware/full_data_exact_tests.tsv.gz --rmats $A/rmats_paired_full_tests.tsv.gz --majiq-tests $R/majiq_build_min3_cov3/full_tests/normalized --majiq-coverage-null $F/brain/majiq_coverage_null_targeted.parquet --coverage-flagged $F/brain/flagged_full.tsv.gz"
python extra_scripts/plot_discoveries_by_event.py $BR --junction-table analyses/tilgner_junction_replication/all_tested_comprehensive_leaf5/junction_replication_merged.tsv.gz --output docs/figures/discoveries_by_event.pdf
python extra_scripts/plot_discoveries_by_event.py $BR --junction-table $R/comparison_corrected_leaf5/all_tests.tsv.gz --leafcutter-counts $R/leafcutter/clustering/benchmark_perind_numers.counts.gz --junction-annotation $R/pseudobulk_junctions_junctions.tsv.gz --contrasts all --output docs/figures/discoveries_by_event_all.pdf
python extra_scripts/plot_discoveries_by_event.py --tealeaf-events analyses/pbmc/GSE320507_events_sh/test_events.tsv.gz --suppa2 $C/suppa2/full_data_tests.tsv.gz --suppa2-primer-aware $C/suppa2_primer_aware/full_data_tests.tsv.gz --rmats $C/rmats/full_data_tests.tsv.gz --majiq-tests $B/majiq/tests/normalized --junction-table $B/comparison_leaf5/all_tests.tsv.gz --leafcutter-counts $B/leafcutter/clustering/benchmark_perind_numers.counts.gz --junction-annotation $B/pseudobulk_junctions_junctions.tsv.gz --contrasts common --majiq-coverage-null $F/pbmc/majiq_coverage_null_targeted.parquet --coverage-flagged $F/pbmc/flagged_full.tsv.gz --output docs/figures/discoveries_by_event_GSE320507.pdf
echo done
