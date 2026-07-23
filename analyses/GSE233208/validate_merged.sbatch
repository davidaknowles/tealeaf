#!/bin/bash
#SBATCH --job-name=gse233208_validate
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/validate_alevin_quantification.py" \
  --alevin-dir "${DATA_ROOT}/processed/merged_alevin" \
  --manifest "${DATA_ROOT}/manifest/ena_fastqs.tsv" \
  --reference-labels "${DATA_ROOT}/metadata/reference_annotation.csv" \
  --primer-pairs "${DATA_ROOT}/processed/merged_alevin/primer_pairs.tsv" \
  --min-cell-umis 500 \
  --min-reference-overlap 1000 \
  --output "${RUN_ROOT}/merged_validation.json"
