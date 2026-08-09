#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=ec_glmm_concurrency

set -euo pipefail

: "${TARGET_SCRIPT:?set TARGET_SCRIPT}"
: "${OUTPUT_TAG:?set OUTPUT_TAG}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"

run_shard() {
  local index=$1
  SLURM_CPUS_PER_TASK=1 \
  SHARD_COUNT=2 \
  SHARD_INDEX=${index} \
  OUTPUT_TAG=${OUTPUT_TAG} \
  CANDIDATE_CACHE=${CANDIDATE_CACHE} \
  CALIBRATION=lrt_bic \
  METHODS=laplace_multinomial \
  NULL_REPLICATES=0 \
  LAPLACE_OPTIMIZER=lbfgs \
  NULL_WARM_START=default \
  bash "${TARGET_SCRIPT}"
}

run_shard 0 &
first=$!
run_shard 1 &
second=$!
status=0
wait "${first}" || status=$?
wait "${second}" || status=$?
exit "${status}"
