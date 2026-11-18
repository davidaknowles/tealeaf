#!/bin/bash
# Evaluation under the consistent local-anchor rule (assemble_localrule.py: brain v4lr, blood _lr), as eval_sh.sh and
# pbmc_post.sh: event classes, coverage-filtered split and long reads,
# discoveries, dataset replication and statistics. Figures go to analyses/local_anchors/figures (the paper's are untouched).
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1
PY=~/venv/cpa_blackwell/bin/python; O=analyses/local_anchors; F=analyses/coverage_null_filter_lr; mkdir -p $O/figures_lr
M=.cache/pbmc/kinnex_parse/molecules/molecules.parquet; L=.cache/pbmc/kinnex_parse/celltypist_labels.tsv.gz
# event classes
$PY extra_scripts/categorize_block_tests.py --blocks .cache/native_v4/collection_localrule/all/blocks.json.gz --cache .cache/binary_grid --internal-model internal_native_v4lr --terminal-model tj_native_v4lr --assessment .cache/local_anchor/lr_brain_v4lr --output-dir analyses/native_align/events_v4lr
$PY extra_scripts/categorize_block_tests.py --blocks .cache/pbmc/GSE320507/collection_localrule/all/blocks.json.gz --cache .cache/binary_grid --internal-model pbmc_GSE320507_internal_lr --terminal-model pbmc_GSE320507_terminal_lr --assessment analyses/pbmc/GSE320507_lr_lr --output-dir analyses/pbmc/GSE320507_events_lr
# coverage-filtered evaluation
export LEAF5=1 SUFFIX=_lr BRAIN=v4lr BRAIN_BLOCKS=.cache/native_v4/collection_localrule/all/blocks.json.gz LR_MAPPING=.cache/local_anchor/lr_brain_v4lr/lr_mapping.tsv.gz PBMC_LR_MAPPING=analyses/pbmc/GSE320507_lr_lr/lr_mapping.tsv.gz
$PY .cache/adhoc/coverage_flags_full.py brain $F/brain
$PY .cache/adhoc/coverage_flags_full.py pbmc $F/pbmc
$PY .cache/adhoc/coverage_filter_eval.py brain $F/brain
$PY .cache/adhoc/coverage_filter_eval.py pbmc $F/pbmc --skip-lr
$PY .cache/adhoc/pbmc_lr_filter.py $F/pbmc
LR_DEPTH=1000 $PY .cache/adhoc/pbmc_lr_filter.py $F/pbmc
# dataset replication and statistics
for g in GSE320507 GSE327907 GSE334840; do $PY .cache/adhoc/pbmc/pbmc_split.py $g analyses/pbmc/${g}_lr; done
$PY .cache/adhoc/pbmc/pbmc_split.py GSE320507 analyses/pbmc/GSE320507_disease_lr --model pbmc_GSE320507_disease
$PY .cache/adhoc/pbmc/paper_stats.py > analyses/pbmc/paper_stats_lr.txt
# discoveries by event class
source ~/venv/kinnex_py/bin/activate
R=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/junction_benchmark; B=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/GSE320507/junction_benchmark
A=analyses/comparator_suppa_rmats; C=$B/comparator_suppa_rmats
BR="--tealeaf-events analyses/native_align/events_v4lr/test_events.tsv.gz --suppa2 $A/suppa2/full_data_exact_tests.tsv.gz --suppa2-primer-aware $A/suppa2_primer_aware/full_data_exact_tests.tsv.gz --rmats $A/rmats_paired_full_tests.tsv.gz --majiq-tests $R/majiq_build_min3_cov3/full_tests/normalized --majiq-coverage-null $F/brain/majiq_coverage_null_targeted.parquet --coverage-flagged $F/brain/flagged_full.tsv.gz"
python extra_scripts/plot_discoveries_by_event.py $BR --junction-table analyses/tilgner_junction_replication/all_tested_comprehensive_leaf5/junction_replication_merged.tsv.gz --output $O/figures_lr/discoveries_by_event.pdf
python extra_scripts/plot_discoveries_by_event.py $BR --junction-table $R/comparison_corrected_leaf5/all_tests.tsv.gz --leafcutter-counts $R/leafcutter/clustering/benchmark_perind_numers.counts.gz --junction-annotation $R/pseudobulk_junctions_junctions.tsv.gz --contrasts all --output $O/figures_lr/discoveries_by_event_all.pdf
python extra_scripts/plot_discoveries_by_event.py --tealeaf-events analyses/pbmc/GSE320507_events_lr/test_events.tsv.gz --suppa2 $C/suppa2/full_data_tests.tsv.gz --suppa2-primer-aware $C/suppa2_primer_aware/full_data_tests.tsv.gz --rmats $C/rmats/full_data_tests.tsv.gz --majiq-tests $B/majiq/tests/normalized --junction-table $B/comparison_leaf5/all_tests.tsv.gz --leafcutter-counts $B/leafcutter/clustering/benchmark_perind_numers.counts.gz --junction-annotation $B/pseudobulk_junctions_junctions.tsv.gz --contrasts common --majiq-coverage-null $F/pbmc/majiq_coverage_null_targeted.parquet --coverage-flagged $F/pbmc/flagged_full.tsv.gz --output $O/figures_lr/discoveries_by_event_GSE320507.pdf
echo done
