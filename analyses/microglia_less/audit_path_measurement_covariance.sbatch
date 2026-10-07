#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-15%8
#SBATCH --job-name=path_measurement_covariance
#SBATCH --output=.cache/paired_reporting_20261006/measurement_%A_%a.out
#SBATCH --error=.cache/paired_reporting_20261006/measurement_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m extra_scripts.audit_path_measurement_covariance \
  --data-cache "${RUN_ROOT}/differential/ec_glmm_prepared.pkl" \
  --candidate-cache "${CANDIDATE_CACHE}" --reference "${REFERENCE}" \
  --output-dir "${OUTPUT_ROOT}/shard_${SLURM_ARRAY_TASK_ID}" \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 16 \
  --blocks "${BLOCKS:-32}" --bootstrap "${BOOTSTRAP:-32}" \
  --biological-concentration "${BIOLOGICAL_CONCENTRATION:-20}" \
  --concentration "${QUANT_CONCENTRATION:-.25}" --subject-scale "${SUBJECT_SCALE:-.5}"
