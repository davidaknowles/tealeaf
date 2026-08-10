#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --job-name=hybrid_ec_glmm
#SBATCH --output=/scratch/daknowles/hybrid_ec_glmm_%A_%a.out
#SBATCH --error=/scratch/daknowles/hybrid_ec_glmm_%A_%a.err

set -euo pipefail

: "${CONFIG:?set CONFIG}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"
: "${BATCHABILITY:?set BATCHABILITY}"
: "${MODEL:?set MODEL}"
: "${OUTPUT_TAG:?set OUTPUT_TAG}"
: "${SHARD_COUNT:?set SHARD_COUNT}"

ROUTE=${ROUTE:-all}
VENV_PATH=${VENV_PATH:-"${HOME}/venv/cpa_blackwell"}

source "${CONFIG}"
module load "${PYTHON_MODULE}"
source "${VENV_PATH}/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${SLURM_CPUS_PER_TASK}"

command=(
  python "${REPO_ROOT}/analyses/run_hybrid_batched_ec_glmm.py"
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl"
  --candidate-cache "${CANDIDATE_CACHE}"
  --batchability "${BATCHABILITY}"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}_${MODEL}_${ROUTE}_shard${SLURM_ARRAY_TASK_ID}"
  --model "${MODEL}"
  --route "${ROUTE}"
  --batch-size "${BATCH_SIZE:-8}"
  --minimum-batch-size "${MINIMUM_BATCH_SIZE:-4}"
  --maximum-padding-ratio "${MAXIMUM_PADDING_RATIO:-2.0}"
  --maximum-coefficients "${MAXIMUM_COEFFICIENTS:-128}"
  --max-iter "${MAX_ITER:-300}"
  --mode-steps "${MODE_STEPS:-12}"
  --continuation-mode-steps "${CONTINUATION_MODE_STEPS:-30}"
  --retries "${RETRIES:-2}"
  --laplace-cache-size "${LAPLACE_CACHE_SIZE:-16}"
  --shard-index "${SLURM_ARRAY_TASK_ID}"
  --shard-count "${SHARD_COUNT}"
)
if [[ -n "${NULL_FIT_SHARDS:-}" ]]; then
  command+=(--null-fit-shards "${NULL_FIT_SHARDS}")
fi
if [[ -n "${LABEL_PERMUTATION_SEED:-}" ]]; then
  command+=(--label-permutation-seed "${LABEL_PERMUTATION_SEED}")
fi
"${command[@]}"
