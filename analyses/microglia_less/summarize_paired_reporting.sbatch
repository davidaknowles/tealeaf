#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --job-name=paired_reporting_summary
#SBATCH --output=.cache/paired_reporting_20261006/summary_%j.out
#SBATCH --error=.cache/paired_reporting_20261006/summary_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ -n "${FREE_CACHE:-}" ]]; then
  ARGS+=(--free-cache "${FREE_CACHE}")
fi
python "${REPO_ROOT}/extra_scripts/summarize_paired_reporting.py" \
  --cache "${AUDIT_CACHE:-${REPO_ROOT}/.cache/paired_reporting_20261006}" \
  --previous-cache "${REPO_ROOT}/.cache/path_reporting_20261006" \
  --run-root "${RUN_ROOT}" --output-dir "${AUDIT_OUTPUT:-${REPO_ROOT}/analyses/paired_reporting_assessment}" "${ARGS[@]}"
