#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --job-name=paired_path

set -euo pipefail

source "${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

SHARD_COUNT=${SHARD_COUNT:-96}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}

python "${REPO_ROOT}/extra_scripts/run_paired_path_test.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE:?}" \
  --output-dir "${OUTPUT_ROOT:?}/shard_${SHARD_INDEX}" \
  --null-replicates "${NULL_REPLICATES:-32}" \
  --path-pseudocount "${PATH_PSEUDOCOUNT:-0}" \
  --seed "${SEED:-20260830}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"
