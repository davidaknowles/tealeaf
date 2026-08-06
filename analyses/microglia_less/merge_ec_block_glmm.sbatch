#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"
: "${OUTPUT_TAG:?set OUTPUT_TAG}"
: "${MERGE_OUTPUT:?set MERGE_OUTPUT}"
: "${CANDIDATE_CACHE:?set CANDIDATE_CACHE}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

SHARD_COUNT=${SHARD_COUNT:-96}
shards=()
for ((index = 0; index < SHARD_COUNT; index++)); do
  shards+=("${RUN_ROOT}/differential/${OUTPUT_TAG}_shard${index}")
done

python "${REPO_ROOT}/extra_scripts/merge_ec_block_asymptotic.py" \
  --shards "${shards[@]}" \
  --candidate-cache "${CANDIDATE_CACHE}" \
  --output-dir "${MERGE_OUTPUT}"
