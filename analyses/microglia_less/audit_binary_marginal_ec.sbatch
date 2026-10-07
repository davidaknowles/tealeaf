#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-3%4
#SBATCH --job-name=binary_marginal_ec
#SBATCH --output=.cache/paired_reporting_20261006/marginal_ec_%A_%a.out
#SBATCH --error=.cache/paired_reporting_20261006/marginal_ec_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
panel_args=()
if [[ "${ONLY_NULL:-0}" == 1 ]]; then panel_args+=(--only-null); fi
if [[ "${DIAGNOSTICS:-0}" == 1 ]]; then panel_args+=(--diagnostics); fi
python -m "${EC_PILOT_MODULE:-extra_scripts.audit_binary_marginal_ec}" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" --reference "${REFERENCE}" \
  --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --blocks "${BLOCKS:-4}" --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count "${SHARD_COUNT:-4}" \
  --subject-nodes "${SUBJECT_NODES:-48}" --path-nodes "${PATH_NODES:-32}" \
  --max-iter "${MAX_ITER:-100}" --draws "${DRAWS:-1}" "${panel_args[@]}"
