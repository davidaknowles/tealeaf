#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=merge_hybrid_ec
#SBATCH --output=/scratch/daknowles/merge_hybrid_ec_%j.out
#SBATCH --error=/scratch/daknowles/merge_hybrid_ec_%j.err

set -euo pipefail

: "${CONFIG:?set CONFIG}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"
: "${BATCHABILITY:?set BATCHABILITY}"
: "${OUTPUT_TAG:?set OUTPUT_TAG}"

source "${CONFIG}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

command=(
  python "${REPO_ROOT}/analyses/merge_hybrid_batched_ec_glmm.py"
  --null-shards "${RUN_ROOT}/differential/${OUTPUT_TAG}_null_*_shard*"
  --alternative-shards "${RUN_ROOT}/differential/${OUTPUT_TAG}_alternative_*_shard*"
  --candidate-cache "${CANDIDATE_CACHE}"
  --batchability "${BATCHABILITY}"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_merged"
)
if [[ -n "${SCALAR_RESULTS:-}" ]]; then
  command+=(--scalar-results "${SCALAR_RESULTS}")
fi
"${command[@]}"
