#!/bin/bash
# Long-read direction of the local-anchor block tests: brain (scISOr-Seq, as compare.sh) and GSE320507 (Kinnex molecules).
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1
PY=~/venv/cpa_blackwell/bin/python
IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor
$PY -W ignore extra_scripts/assess_block_structure_filter.py --cache .cache/binary_grid --model internal_native_v4sh_L --terminal-model tj_native_v4sh_L --local-reads .cache/native_v4/collection_local/all --filter internal \
  --output-dir .cache/local_anchor/lr_brain_L --matrix-dir .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --gtf $IDX/gencode.vM32.basic.annotation.gtf \
  --long-read-gtf .cache/reference/gencode.vM32.annotation.gtf.gz --chain-mapping --block-cache .cache/local_anchor/brain_local_blocks.json.gz --split-universe analyses/split_universe_leaf5/comparator_tests.tsv.gz > .cache/local_anchor/lr_brain_L.log 2>&1 || echo "brain LR failed, see log"
M=.cache/pbmc/kinnex_parse/molecules/molecules.parquet; L=.cache/pbmc/kinnex_parse/celltypist_labels.tsv.gz
$PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE320507_internal_sh_L,pbmc_GSE320507_terminal_sh_L --blocks .cache/pbmc/GSE320507/blocks_local.json.gz --molecules $M --labels $L --output-dir analyses/pbmc/GSE320507_lr_sh_L
echo done
