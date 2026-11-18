#!/bin/bash
# Long-read direction of the new local-rule block tests (_N models). Usage: lr_new.sh <brain|GSE320507|GSE327907|GSE334840>
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD OMP_NUM_THREADS=1
PY=~/venv/cpa_blackwell/bin/python
M=.cache/pbmc/kinnex_parse/molecules/molecules.parquet; L=.cache/pbmc/kinnex_parse/celltypist_labels.tsv.gz
case $1 in
  brain) IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor
    $PY -W ignore extra_scripts/assess_block_structure_filter.py --cache .cache/binary_grid --model internal_native_v4sh_N --terminal-model tj_native_v4sh_N --local-reads .cache/native_v4/collection_lrnew/all --filter internal \
      --output-dir .cache/local_anchor/lr_brain_N --matrix-dir .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --gtf $IDX/gencode.vM32.basic.annotation.gtf \
      --long-read-gtf .cache/reference/gencode.vM32.annotation.gtf.gz --chain-mapping --block-cache .cache/local_anchor/brain_localrule_blocks.json.gz --split-universe analyses/split_universe_leaf5/comparator_tests.tsv.gz ;;
  GSE320507) $PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE320507_internal_sh_N,pbmc_GSE320507_terminal_sh_N --blocks .cache/pbmc/GSE320507/blocks_localrule_new.json.gz --molecules $M --labels $L --output-dir analyses/pbmc/GSE320507_lr_sh_N ;;
  GSE327907) $PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE327907_internal_sh_N,pbmc_GSE327907_terminal_sh_N --blocks .cache/pbmc/GSE327907/blocks_localrule_new.json.gz --molecules $M --labels $L --label-column Immune_All_High_majority --merge Macrophages=Monocytes --suffix "|Baseline" --output-dir analyses/pbmc/GSE327907_lr_sh_N ;;
  GSE334840) $PY extra_scripts/assess_longread_molecules.py --cache .cache/binary_grid --models pbmc_GSE334840_internal_sh_N,pbmc_GSE334840_terminal_sh_N --blocks .cache/pbmc/GSE334840/blocks_localrule_new.json.gz --molecules $M --labels $L --label-column Immune_All_High_majority --merge Macrophages=Monocytes --suffix "|DMSO" --output-dir analyses/pbmc/GSE334840_lr_sh_N ;;
esac
