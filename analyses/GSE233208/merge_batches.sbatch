#!/bin/bash
#SBATCH --job-name=gse233208_merge
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

set -euo pipefail
source "$(dirname "$0")/config.env"
module load python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${RUN_ROOT}/logs"
inputs=()
for batch in "${BATCHES[@]}"; do
  inputs+=(--input "${batch}=${DATA_ROOT}/processed/${batch}/alevin_quant")
done
python "${REPO_ROOT}/extra_scripts/merge_alevin_batches.py" \
  "${inputs[@]}" \
  --output "${DATA_ROOT}/processed/merged_alevin"
python "${REPO_ROOT}/extra_scripts/write_parse_primer_pairs.py" \
  --cell-rows "${DATA_ROOT}/processed/merged_alevin/quants_mat_rows.txt" \
  --rt-barcodes "${PARSE_RT_BARCODES}" \
  --output "${DATA_ROOT}/processed/merged_alevin/primer_pairs.tsv"
