#!/bin/bash
# Example panels under the hierarchical rule: brain usage panels (examples_v4hl), sashimi plots of the new examples
# (Dock9 L1 brain; INPP4B L9 donor groups; PTPRC L1 classical monocytes vs naive B cells in GSE320507).
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD
PY=~/venv/cpa_blackwell/bin/python; KP=~/venv/kinnex_py/bin/python
IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor; GTF=~/tealeaf_runs/pbmc/ref/gencode.v49.basic.annotation.gtf
D=.cache/pbmc/GSE320507; SETUP_R='module purge; module load R/4.4-custom'
LABELS=(--level-label "ASC=Astrocytes" --level-label "EX_cortical=Cortical excitatory" --level-label "ODC=Oligodendrocytes" --level-label "EX_hippocampus_granule=DG granule" --level-label "INH_medium_spiny=Medium spiny")
# brain usage panels, with and without the path structure
$KP extra_scripts/plot_block_examples.py --tables analyses/native_align/examples_v4hl --output-dir analyses/native_align/examples_v4hl/usage_paper --width 3.2 --height 3.1 --base-size 7 "${LABELS[@]}"
$KP extra_scripts/plot_block_examples.py --tables analyses/native_align/examples_v4hl --output-dir analyses/native_align/examples_v4hl/usage_paper_nopaths --width 3.2 --height 2.1 --base-size 7 --no-structure "${LABELS[@]}"
# Dock9 sashimi (block and gene views come from the same call)
awk -F'\t' 'NR == 1 || $2 ~ /^Dock9/' analyses/native_align/examples_v4hl/examples.tsv > analyses/native_align/examples_v4hl/examples_dock9.tsv
$PY extra_scripts/plot_block_sashimi.py --examples analyses/native_align/examples_v4hl/examples_dock9.tsv --blocks .cache/native_v4/collection_hier/all/blocks.json.gz --native-dir .cache/native_v3 --recipe .cache/native_cells/recipe.json --gtf $IDX/gencode.vM32.basic.annotation.gtf --long-read-gtf .cache/reference/gencode.vM32.annotation.gtf.gz --long-read-matrix .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --ggsashimi tealeaf/ggsashimi/tealeaf_ggsashimi.py --output-dir analyses/native_align/examples_v4hl/sashimi_paper --width 5.6 --track-height 0.36 --base-size 7 --y-breaks 2 --level-label "EX_cortical=Cortical excitatory" --level-label "ODC=Oligodendrocytes" --plot-setup "$SETUP_R"
# INPP4B donor-group sashimi
awk -F'\t' 'NR == 1 || $2 ~ /^INPP4B/' analyses/pbmc/GSE320507_disease/examples_hl/examples.tsv > analyses/pbmc/GSE320507_disease/examples_hl/examples_inpp4b.tsv
$PY extra_scripts/plot_block_sashimi.py --examples analyses/pbmc/GSE320507_disease/examples_hl/examples_inpp4b.tsv --blocks $D/disease_hier/all/blocks.json.gz --native-dir $D --recipe $D/recipe.json --gtf $GTF --subject-groups $D/subject_folds.tsv --group-column condition --ggsashimi tealeaf/ggsashimi/tealeaf_ggsashimi.py --output-dir analyses/pbmc/GSE320507_disease/examples_hl/sashimi_paper --width 5.6 --track-height 0.36 --base-size 7 --y-breaks 2 --plot-setup "$SETUP_R"
# PTPRC exon 3-7 block, classical monocytes vs naive B cells
$PY extra_scripts/export_block_examples.py --assessment analyses/pbmc/GSE320507_lr_hl --cache .cache/binary_grid --internal-model pbmc_GSE320507_internal_hl --terminal-model pbmc_GSE320507_terminal_hl --internal-reads $D/collection_hier/all --terminal-reads $D/collection_hier/junction --gtf $GTF --output-dir analyses/pbmc/GSE320507_ptprc_hl --test-ids "ENSG00000081237.23:L1|cell_type|Classical_monocytes|Naive_B_cells"
$KP extra_scripts/plot_block_examples.py --tables analyses/pbmc/GSE320507_ptprc_hl --output-dir analyses/pbmc/GSE320507_ptprc_hl/usage_paper --width 3.2 --height 3.1 --base-size 7 --level-label "Classical_monocytes=Classical monocytes" --level-label "Naive_B_cells=Naive B cells"
$PY extra_scripts/plot_block_sashimi.py --examples analyses/pbmc/GSE320507_ptprc_hl/examples.tsv --blocks $D/collection_hier/junction/blocks.json.gz --native-dir $D --recipe $D/recipe.json --gtf $GTF --ggsashimi tealeaf/ggsashimi/tealeaf_ggsashimi.py --output-dir analyses/pbmc/GSE320507_ptprc_hl/sashimi_paper --width 5.6 --track-height 0.36 --base-size 7 --y-breaks 2 --level-label "Classical_monocytes=Classical monocytes" --level-label "Naive_B_cells=Naive B cells" --plot-setup "$SETUP_R"
echo done
