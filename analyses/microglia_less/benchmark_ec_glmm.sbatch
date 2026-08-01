#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-23
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --job-name=ec_glmm_sim

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

SHARD_COUNT=${SHARD_COUNT:-24}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID}}
REPLICATES=${REPLICATES:-96}
OUTPUT_TAG=${OUTPUT_TAG:-ec_glmm_simulation}
FAMILIES=${FAMILIES:-multinomial;dirichlet_multinomial}
ISOFORMS=${ISOFORMS:-2}
INCLUDE_FULL_COVARIANCE=${INCLUDE_FULL_COVARIANCE:-0}
MAX_ITER=${MAX_ITER:-150}
VI_SAMPLES=${VI_SAMPLES:-128}
IFS=';' read -r -a family_args <<< "${FAMILIES}"
full_covariance_args=()
if [[ "${INCLUDE_FULL_COVARIANCE}" == "1" ]]; then
  full_covariance_args+=(--include-full-covariance)
fi

python "${REPO_ROOT}/extra_scripts/benchmark_ec_glmm.py" \
  --output "${RUN_ROOT}/differential/${OUTPUT_TAG}/shard${SHARD_INDEX}.tsv" \
  --replicates "${REPLICATES}" \
  --families "${family_args[@]}" \
  --isoforms "${ISOFORMS}" \
  --max-iter "${MAX_ITER}" \
  --vi-samples "${VI_SAMPLES}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  "${full_covariance_args[@]}"
