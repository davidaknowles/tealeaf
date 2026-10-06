#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-15%12
#SBATCH --job-name=wild_omnibus
#SBATCH --output=.cache/path_reporting_20261006/%x_%A_%a.out
#SBATCH --error=.cache/path_reporting_20261006/%x_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/audit_wild_omnibus.py" \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" --candidate-cache "${CANDIDATE_CACHE}" \
  --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 16 --test-concentration "${TEST_CONCENTRATION:-32}"
