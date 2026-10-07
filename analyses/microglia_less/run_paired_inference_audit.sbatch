#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-31%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=paired_inference
#SBATCH --output=.cache/paired_reporting_20261006/inference_%x_%A_%a.out
#SBATCH --error=.cache/paired_reporting_20261006/inference_%x_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${SCORE_FREE_NULL:-0}" == 1 ]]; then
  ARGS+=(--score-free-null)
fi
python "${REPO_ROOT}/extra_scripts/run_paired_path_test.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count "${SHARD_COUNT:-32}" \
  --max-iter "${MAX_ITER:-100}" \
  --null-replicates 64 --path-pseudocount "${CONCENTRATION:-32}" \
  --path-prior-center "${PRIOR_CENTER:-uniform}" --path-pseudocount-scaling total \
  --paired-inference "${BACKEND:-ec-score}" --score-null-concentration "${NULL_CONCENTRATION:-1}" \
  --export-path-usage "${ARGS[@]}"
