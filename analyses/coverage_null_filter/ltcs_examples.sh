#!/bin/bash
# LTCS vs HC example tables and sashimi plots (GSE320507 disease contrasts)
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD
P=~/venv/cpa_blackwell/bin/python
D=.cache/pbmc/GSE320507; O=analyses/pbmc/GSE320507_disease/examples; GTF=~/tealeaf_runs/pbmc/ref/gencode.v49.basic.annotation.gtf
$P extra_scripts/export_block_examples.py --cache .cache/binary_grid --internal-model pbmc_GSE320507_disease_internal --terminal-model pbmc_GSE320507_disease_terminal --internal-reads $D/disease/all --terminal-reads $D/disease/junction --gtf $GTF --groups LTCS HC --max-fdr .05 --min-gene-umis 50 --unpaired --any-structure --max-window 20000 --output-dir $O
$P extra_scripts/plot_block_sashimi.py --examples $O/examples.tsv --blocks $D/disease/all/blocks.json.gz --native-dir $D --recipe $D/recipe.json --gtf $GTF --subject-groups $D/subject_folds.tsv --group-column condition --ggsashimi tealeaf/ggsashimi/tealeaf_ggsashimi.py --output-dir $O/sashimi --plot-setup "module purge; module load R/4.4-custom"
