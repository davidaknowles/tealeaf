#!/bin/bash -l
set -euo pipefail
source "$(dirname "$0")/config.env"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

root="${RUN_ROOT}/junction_benchmark"
mkdir -p "${root}"
python "${REPO_ROOT}/extra_scripts/prepare_junction_benchmark.py" whitelists \
  --barcodes "${PARSE_PERMIT_LIST}" \
  --output-dir "${root}/whitelists"
python "${REPO_ROOT}/extra_scripts/prepare_junction_cell_metadata.py" \
  --combined-groups "${RUN_ROOT}/barcodes_to_cluster_dx_mouse.csv" \
  --output "${root}/cell_metadata.tsv" \
  --combined-output "${root}/barcode_groups.csv"
