#!/bin/bash
#SBATCH --job-name=gse233208_preprocessing_gate
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/validate_parse_preprocessing.py" \
  --manifest "${DATA_ROOT}/manifest/ena_fastqs.tsv" \
  --processed-root "${DATA_ROOT}/processed" \
  --output "${RUN_ROOT}/preprocessing_validation.json"
