#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=microglia_ec_block_packed
#SBATCH --output=/scratch/daknowles/tealeaf_ec_block_packed_%j.out
#SBATCH --error=/scratch/daknowles/tealeaf_ec_block_packed_%j.err

set -euo pipefail

: "${SHARD_COUNT:?set SHARD_COUNT}"
: "${OUTPUT_TAG:?set OUTPUT_TAG}"
: "${CONFIG:?set CONFIG}"

source "${CONFIG}"
worker_count=${PACKED_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}
shard_start=${SHARD_START:-0}
shard_stop=${SHARD_STOP:-${SHARD_COUNT}}
run_root=${RUN_ROOT}
runner=${REPO_ROOT}/analyses/microglia_less/run_ec_block_glmm.sbatch
if (( shard_start < 0 || shard_stop > SHARD_COUNT || shard_start >= shard_stop )); then
  echo "invalid shard interval [${shard_start}, ${shard_stop}) for ${SHARD_COUNT} shards" >&2
  exit 2
fi

run_shard() {
  local shard=$1
  local summary=${run_root}/differential/${OUTPUT_TAG}_shard${shard}/summary.json
  if [[ -s "${summary}" ]]; then
    return
  fi
  SHARD_INDEX=${shard} SLURM_CPUS_PER_TASK=1 bash "${runner}"
}

for ((shard = shard_start; shard < shard_stop; shard++)); do
  while (( $(jobs -pr | wc -l) >= worker_count )); do
    wait -n
  done
  run_shard "${shard}" &
done
wait
