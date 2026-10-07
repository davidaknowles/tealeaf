#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-31%24
#SBATCH --job-name=path_split_audit
#SBATCH --output=.cache/path_reporting_20261006/%x_%A_%a.out
#SBATCH --error=.cache/path_reporting_20261006/%x_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${INCLUDE_SCORE:-0}" == 1 ]]; then
  ARGS+=(--include-score)
fi
if [[ "${SCORE_ONLY:-0}" == 1 ]]; then
  ARGS+=(--score-only)
fi
if [[ -n "${REFERENCE:-}" ]]; then
  ARGS+=(--reference "${REFERENCE}")
fi
if [[ "${PROFILE_MASS:-0}" == 1 ]]; then
  ARGS+=(--profile-mass)
fi
python "${REPO_ROOT}/extra_scripts/audit_path_reporting_omnibus.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" --mode "${MODE}" \
  --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count "${SHARD_COUNT:-32}" \
  --null-replicates "${NULL_REPLICATES:-64}" --test-concentration "${TEST_CONCENTRATION:-32}" \
  --prior-center "${PRIOR_CENTER:-uniform}" "${ARGS[@]}"
