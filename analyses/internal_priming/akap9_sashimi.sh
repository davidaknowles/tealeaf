#!/bin/bash
set -euo pipefail
unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done
module load Python/3.12.3-GCCcore-13.3.0
cd ~/projects/tealeaf; export PYTHONPATH=$PWD
IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor
C=.cache/native_v4/collection_parity_p2_v4
~/venv/cpa_blackwell/bin/python extra_scripts/plot_block_sashimi.py --examples .cache/internal_priming/examples_akap9.tsv --blocks $C/all/blocks.json.gz --native-dir .cache/native_v3 --recipe .cache/native_cells/recipe.json --gtf $IDX/gencode.vM32.basic.annotation.gtf --long-read-gtf .cache/reference/gencode.vM32.annotation.gtf.gz --long-read-matrix .cache/tilgner_scisorseq/LR_annotTranscripts_byCT --ggsashimi tealeaf/ggsashimi/tealeaf_ggsashimi.py --output-dir analyses/native_align/examples_v4/sashimi_paper --width 5.6 --track-height 0.36 --base-size 7 --y-breaks 2 --level-label "INH_medium_spiny=Medium spiny" --level-label "ODC=Oligodendrocytes" --plot-setup "module purge; module load R/4.4-custom"
