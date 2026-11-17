#!/bin/bash
# Block sets of the consistent local-anchor rule for each dataset (annotated blocks, novel paths and new blocks from
# unannotated junctions), from the existing alignments and junction support: <run>/blocks_localrule.json.gz.
set -euo pipefail
cd ~/projects/tealeaf
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD'
IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor; R=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/ref
N=.cache/native_v4
sbatch --parsable -p cpu -c 4 --mem=64G -t 6:00:00 -J LR_blocks_brain -o $N/logs_local/localrule_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/native_collect.py --stage junctions --output-dir $N --recipe .cache/native_cells/recipe.json --known-junctions $N/junctions.tsv.gz --all-multiexon-genes --novel-terminal-length 151 --candidate-cache /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential/ec_block_full_pairwise_candidates.pkl --block-cache .cache/local_anchor/brain_localrule_blocks.json.gz --gtf $IDX/gencode.vM32.basic.annotation.gtf --blocks $N/blocks_localrule.json.gz --no-junction-table'"
for spec in GSE320507:136 GSE327907:150 GSE334840:66; do g=${spec%:*}; rl=${spec#*:}; P=.cache/pbmc/$g
  sbatch --parsable -p cpu -c 4 --mem=96G -t 6:00:00 -J LR_blocks_$g -o $P/logs_local/localrule_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/native_collect.py --stage junctions --output-dir $P --recipe $P/recipe.json --known-junctions $P/junctions.tsv.gz --all-multiexon-genes --novel-terminal-length $rl --block-cache .cache/local_anchor/pbmc_localrule_blocks.json.gz --gtf $R/gencode.v49.basic.annotation.gtf.gz --blocks $P/blocks_localrule.json.gz --no-junction-table'"
done
