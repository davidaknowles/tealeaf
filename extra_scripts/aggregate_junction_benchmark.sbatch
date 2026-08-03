#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${RUN_LIST:?tab-separated RUN and STAR output directory}"
: "${CELL_METADATA:?set CELL_METADATA}"
: "${OUTPUT_PREFIX:?set OUTPUT_PREFIX}"
: "${LEAFCUTTER_DIR:?set LEAFCUTTER_DIR}"
: "${GTF:?set GTF}"

module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
run_args=()
while IFS=$'\t' read -r run path; do
  run_args+=(--run "${run}=${path}")
done < "${RUN_LIST}"
python "${REPO_ROOT}/extra_scripts/prepare_junction_benchmark.py" aggregate \
  "${run_args[@]}" --cell-metadata "${CELL_METADATA}" \
  --output-prefix "${OUTPUT_PREFIX}" --leafcutter-dir "${LEAFCUTTER_DIR}" \
  --gtf "${GTF}"
