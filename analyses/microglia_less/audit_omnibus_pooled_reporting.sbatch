#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-15%12
#SBATCH --job-name=omnibus_reporting
#SBATCH --output=.cache/path_reporting_20261006/%x_%A_%a.out
#SBATCH --error=.cache/path_reporting_20261006/%x_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${PROFILE_MASS:-0}" == 1 ]]; then
  ARGS+=(--profile-mass)
fi
python "${REPO_ROOT}/extra_scripts/audit_omnibus_pooled_reporting.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_omnibus_min4_candidates.pkl" \
  --output "${AUDIT_CACHE:-${REPO_ROOT}/.cache/path_reporting_20261006}/omnibus_pooled/shard_${SLURM_ARRAY_TASK_ID}/observed.tsv" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 16 "${ARGS[@]}"
