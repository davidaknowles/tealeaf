#!/bin/bash
# Current model (v4) vs primer-offset shrinkage (v4sh) vs internal priming + shrinkage (v4ip): split replication,
# long-read direction and discoveries, all scored by the same code without the coverage-only filter.
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1
PY=~/venv/cpa_blackwell/bin/python
C=.cache/native_v4/collection_parity_p2_v4; R=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential
IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor
for V in v4sh v4ip; do
  for c in 0 1 full; do
    $PY extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort $c --output-dir .cache/binary_grid/internal_native_${V}_$c --shard-count 160 --merge
    $PY extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort $c --output-dir .cache/binary_grid/tj_native_${V}_$c --shard-count 40 --merge
  done
done
for V in v4 v4sh v4ip; do
  $PY -W ignore extra_scripts/assess_block_structure_filter.py --cache .cache/binary_grid --model internal_native_$V --terminal-model tj_native_$V --local-reads $C/all --filter internal \
    --output-dir .cache/internal_priming/compare/$V --matrix-dir .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --gtf $IDX/gencode.vM32.basic.annotation.gtf \
    --long-read-gtf .cache/reference/gencode.vM32.annotation.gtf.gz --chain-mapping --block-cache $R/gencode_vM32_splice_blocks.json.gz --split-universe analyses/split_universe_leaf5/comparator_tests.tsv.gz > /dev/null
done
$PY .cache/internal_priming/compare_summary.py
