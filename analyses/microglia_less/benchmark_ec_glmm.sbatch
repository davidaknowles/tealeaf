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

python "${REPO_ROOT}/extra_scripts/benchmark_ec_glmm.py" \
  --output "${RUN_ROOT}/differential/ec_glmm_simulation/shard${SHARD_INDEX}.tsv" \
  --replicates "${REPLICATES}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}"
