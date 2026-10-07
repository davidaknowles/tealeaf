#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --job-name=merge_paired_inference
#SBATCH --output=.cache/paired_reporting_20261006/merge_%j.out
#SBATCH --error=.cache/paired_reporting_20261006/merge_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for FOLDER in "${OUTPUT_ROOT}"/shard_*; do
  [[ -f "${FOLDER}/paired_path.tsv" ]] || exit 1
done
ARGS=()
if [[ "${RETAIN_FAILED_FAMILY:-0}" == 1 ]]; then ARGS+=(--retain-failed-family); fi
python "${REPO_ROOT}/extra_scripts/merge_paired_path_test.py" \
  --shards "${OUTPUT_ROOT}"/shard_* --output-dir "${OUTPUT_ROOT}/merged" \
  --moderate-variances --max-null-replicates 32 "${ARGS[@]}"
