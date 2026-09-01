#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-95
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --job-name=path_smoothing

set -euo pipefail

source "${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ARGS=(
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl"
  --candidate-cache "${CANDIDATE_CACHE:?}"
  --output "${OUTPUT_ROOT:?}/evidence_${SLURM_ARRAY_TASK_ID}.tsv"
  --alpha-grid "${ALPHA_GRID:-0.25;0.5;1;2;4;8;16;32;64;128}"
  --path-pseudocount-scaling "${PATH_PSEUDOCOUNT_SCALING:-per_path}"
  --path-prior-center "${PATH_PRIOR_CENTER:-uniform}"
  --shard-index "${SLURM_ARRAY_TASK_ID}"
  --shard-count "${SHARD_COUNT:-96}"
)
if [[ -n "${MAX_CANDIDATES:-}" ]]; then
  ARGS+=(--max-candidates "${MAX_CANDIDATES}")
fi

python "${REPO_ROOT}/extra_scripts/estimate_path_smoothing.py" "${ARGS[@]}"
