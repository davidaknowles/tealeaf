#!/bin/bash
# Local-anchor blocks for one dataset, added to an existing run: blocks-only junction stage on the local-anchor block
# cache (existing alignments and junction support), collect, collate into collection_local, then fits of the new blocks
# with the settings of the current shrunk-offset models (new models with suffix _L) and merges.
# Usage: local_chain.sh <brain|GSE320507|GSE327907|GSE334840>
set -euo pipefail
cd ~/projects/tealeaf
DS=$1
SETUP='unset LOADEDMODULES; for v in $(compgen -A variable _ModuleTable); do unset "$v"; done; module load Python/3.12.3-GCCcore-13.3.0; source ~/venv/cpa_blackwell/bin/activate; cd ~/projects/tealeaf; export PYTHONPATH=$PWD NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK'
RUST="--rust-library rust/native_align/target/release/libnative_align.so"
SB="sbatch --parsable -p cpu --exclude=ne1dc4-039"
case $DS in
  brain)
    N=.cache/native_v4; RECIPE=.cache/native_cells/recipe.json; NLIB=8; RL=151; COLLECT_EXTRA=--drop-unfit-reads; MEM=64G
    IDX=/gpfs/commons/home/daknowles/knowles_lab/index/kallisto/mus_musculus/with_precursor
    BLOCKS="--all-multiexon-genes --novel-terminal-length 151 --candidate-cache /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential/ec_block_full_pairwise_candidates.pkl --block-cache .cache/local_anchor/brain_local_blocks.json.gz --gtf $IDX/gencode.vM32.basic.annotation.gtf"
    INTERNAL=internal_native_v4sh_L; TERMINAL=tj_native_v4sh_L; FIT_ENV="UNIVERSE=analyses/split_universe/comparator_tests.tsv.gz"; LEVEL="" ;;
  GSE*)
    N=.cache/pbmc/$DS; RECIPE=$N/recipe.json; NLIB=$(python3 -c "import json; print(len(json.load(open('$N/recipe.json'))['bams']))"); COLLECT_EXTRA=""; MEM=96G
    case $DS in GSE320507) RL=136;; GSE327907) RL=150;; GSE334840) RL=66;; esac
    R=/gpfs/commons/home/daknowles/tealeaf_runs/pbmc/ref
    BLOCKS="--all-multiexon-genes --novel-terminal-length $RL --block-cache .cache/local_anchor/pbmc_local_blocks.json.gz --gtf $R/gencode.v49.basic.annotation.gtf.gz"
    INTERNAL=pbmc_${DS}_internal_sh_L; TERMINAL=pbmc_${DS}_terminal_sh_L; FIT_ENV="UNIVERSE=none,LR_TYPES=none,FOLDS=$N/subject_folds.tsv"
    LEVEL=""; [ -f $N/level_pairs.tsv ] && [ $DS != GSE320507 ] && LEVEL="--level-pairs $N/level_pairs.tsv" ;;
esac
COMMON="--output-dir $N --recipe $RECIPE --known-junctions $N/junctions.tsv.gz"
B=$N/blocks_local.json.gz; C=$N/collection_local; L=$N/logs_local; mkdir -p $L
J=$($SB -c 4 --mem=$MEM -t 6:00:00 -J ${DS}_Ljunctions -o $L/junctions_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/native_collect.py --stage junctions $COMMON $BLOCKS --blocks $B --blocks-only'")
K=$($SB --dependency=afterok:$J -c 4 --mem=$MEM -t 12:00:00 -a 0-$((NLIB - 1)) -J ${DS}_Lcollect -o $L/collect_%A_%a.log --wrap "bash -lc '$SETUP; python extra_scripts/native_collect.py --stage collect $COMMON --blocks $B --collection-dir $C $COLLECT_EXTRA --library-index \$SLURM_ARRAY_TASK_ID'")
T=$($SB --dependency=afterok:$K -c 16 --mem=128G -t 12:00:00 -J ${DS}_Lcollate -o $L/collate_%j.log --wrap "bash -lc '$SETUP; python extra_scripts/native_collect.py --stage collate $COMMON $RUST --blocks $B --collection-dir $C --read-length $RL'")
echo "$DS junctions $J collect $K collate $T"
S=48
for v in "$INTERNAL|$C/all|--structure internal" "$TERMINAL|$C/junction|--structure terminal --group-same-junctions"; do
  IFS='|' read -r model reads flags <<< "$v"
  F=$(sbatch --parsable --exclude=ne1dc4-039 -a 0-$((3 * S - 1)) --dependency=afterok:$T -J $model --export=ALL,SHARDS=$S,$FIT_ENV,MODEL=$model,LOCAL_READS=$reads,EXTRA_ARGS="--primer-offset --primer-offset-prior-sd 1 $flags --min-share .01 $LEVEL" extra_scripts/run_local_path_tests.sbatch)
  M=$($SB --dependency=afterok:$F -c 1 --mem=16G -t 2:00:00 -J ${model}_m -o $L/merge_${model}_%j.log --wrap "bash -lc '$SETUP; for c in 0 1 full; do python extra_scripts/run_local_path_tests.py --local-reads x --subject-folds x --cohort \$c --output-dir .cache/binary_grid/${model}_\$c --shard-count $S --merge; done'")
  echo "$model fits $F merge $M"
done
