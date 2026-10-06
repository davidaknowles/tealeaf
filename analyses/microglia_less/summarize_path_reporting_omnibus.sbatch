#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --job-name=path_audit_summary
#SBATCH --output=.cache/path_reporting_20261006/summary_%j.out
#SBATCH --error=.cache/path_reporting_20261006/summary_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${REPORTING_ONLY:-0}" == 1 ]]; then
  ARGS+=(--reporting-only)
fi
if [[ "${INCLUDE_WILD:-0}" == 1 ]]; then
  ARGS+=(--include-wild)
fi
python "${REPO_ROOT}/extra_scripts/summarize_path_reporting_omnibus.py" \
  --run-root "${RUN_ROOT}" --cache "${AUDIT_CACHE:-${REPO_ROOT}/.cache/path_reporting_20261006}" \
  --output-dir "${AUDIT_OUTPUT:-${REPO_ROOT}/analyses/path_reporting_omnibus}" \
  --matrix-dir "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" \
  --gtf "${GENCODE_GTF}" "${ARGS[@]}"
