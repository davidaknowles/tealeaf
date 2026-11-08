#!/bin/bash
# Shrunk-offset (_sh) PBMC fits: merge, long-read direction, event classes, split replication, paper statistics, disease examples.
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1 SUFFIX=_sh BRAIN=v4sh
PY=~/venv/cpa_blackwell/bin/python
merge() { for c in 0 1 full; do $PY extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort $c --output-dir .cache/binary_grid/$1_$c --shard-count $2 --merge; done; }
for g in GSE320507 GSE334840; do merge pbmc_${g}_internal_sh 32; merge pbmc_${g}_terminal_sh 32; done
merge pbmc_GSE327907_internal_sh 64; merge pbmc_GSE327907_terminal_sh 64
merge pbmc_GSE320507_disease_internal_sh 128; merge pbmc_GSE320507_disease_terminal_sh 128
M=.cache/pbmc/kinnex_parse/molecules/molecules.parquet; L=.cache/pbmc/kinnex_parse/celltypist_labels.tsv.gz
$PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE320507_internal_sh,pbmc_GSE320507_terminal_sh --blocks .cache/pbmc/GSE320507/collection/all/blocks.json.gz --molecules $M --labels $L --output-dir analyses/pbmc/GSE320507_lr_sh
$PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE327907_internal_sh,pbmc_GSE327907_terminal_sh --blocks .cache/pbmc/GSE327907/collection/all/blocks.json.gz --molecules $M --labels $L --label-column Immune_All_High_majority --merge Macrophages=Monocytes --suffix "|Baseline" --output-dir analyses/pbmc/GSE327907_lr_sh
$PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE334840_internal_sh,pbmc_GSE334840_terminal_sh --blocks .cache/pbmc/GSE334840/collection/all/blocks.json.gz --molecules $M --labels $L --label-column Immune_All_High_majority --merge Macrophages=Monocytes --suffix "|DMSO" --output-dir analyses/pbmc/GSE334840_lr_sh
$PY extra_scripts/categorize_block_tests.py --blocks .cache/pbmc/GSE320507/collection/all/blocks.json.gz --cache .cache/binary_grid --internal-model pbmc_GSE320507_internal_sh --terminal-model pbmc_GSE320507_terminal_sh --assessment analyses/pbmc/GSE320507_lr_sh --output-dir analyses/pbmc/GSE320507_events_sh
for g in GSE320507 GSE327907 GSE334840; do $PY .cache/adhoc/pbmc/pbmc_split.py $g analyses/pbmc/${g}_sh; done
$PY .cache/adhoc/pbmc/pbmc_split.py GSE320507 analyses/pbmc/GSE320507_disease_sh --model pbmc_GSE320507_disease
$PY .cache/adhoc/pbmc/paper_stats.py > analyses/pbmc/paper_stats_sh.txt
D=.cache/pbmc/GSE320507
$PY extra_scripts/export_block_examples.py --cache .cache/binary_grid --internal-model pbmc_GSE320507_disease_internal_sh --terminal-model pbmc_GSE320507_disease_terminal_sh --internal-reads $D/disease/all --terminal-reads $D/disease/junction --gtf ~/tealeaf_runs/pbmc/ref/gencode.v49.basic.annotation.gtf --groups LTCS HC --max-fdr .05 --min-gene-umis 50 --unpaired --any-structure --max-window 20000 --output-dir analyses/pbmc/GSE320507_disease/examples_sh
echo done
