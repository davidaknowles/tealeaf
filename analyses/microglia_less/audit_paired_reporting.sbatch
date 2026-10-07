#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-31%16
#SBATCH --job-name=paired_reporting
#SBATCH --output=.cache/paired_reporting_20261006/%x_%A_%a.out
#SBATCH --error=.cache/paired_reporting_20261006/%x_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${FREE_ISOFORMS:-0}" == 1 ]]; then
  ARGS+=(--free-isoforms)
fi
python "${REPO_ROOT}/extra_scripts/audit_paired_reporting.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" --reference "${REFERENCE}" \
  --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count "${SHARD_COUNT:-32}" \
  --concentrations "${CONCENTRATIONS:-1;4;16;32;64}" "${ARGS[@]}"
