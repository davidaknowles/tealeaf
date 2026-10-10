#!/bin/bash
#SBATCH --job-name=binary_ec_glmm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-95
#SBATCH --output=.cache/binary_grid/logs/fit_%A_%a.log
# Array task t fits shard t%32 of cohort fold0/fold1/full; MODEL names the output prefix.
set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/cpa_blackwell/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
differential=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/differential
cohorts=(0 1 full)
cohort=${cohorts[$((SLURM_ARRAY_TASK_ID / 32))]}
shard=$((SLURM_ARRAY_TASK_ID % 32))
if [ "$cohort" = full ]; then candidates=$differential/ec_block_full_pairwise_candidates.pkl; else candidates=$differential/ec_block_reproducibility_pairwise_fold${cohort}_candidates.pkl; fi
python -W ignore extra_scripts/run_binary_ec_glmm.py --data-cache $differential/ec_glmm_prepared.pkl --candidate-cache "$candidates" --output-dir ".cache/binary_grid/${MODEL:?}_${cohort}/shard_${shard}" --shard-index "$shard" --shard-count 32 ${EXTRA_ARGS:-}
