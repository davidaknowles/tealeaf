#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-7%8
#SBATCH --job-name=lr_reporting
#SBATCH --output=.cache/lr_ranking_20261006/effects_%A_%a.out
#SBATCH --error=.cache/lr_ranking_20261006/effects_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/audit_tealeaf_lr_effects.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --replication "${REPO_ROOT}/.cache/lr_ranking_20261006/pairwise_replication_precise.tsv" \
  --output "${REPO_ROOT}/.cache/lr_ranking_20261006/effects/shard_${SLURM_ARRAY_TASK_ID}.tsv" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 8
