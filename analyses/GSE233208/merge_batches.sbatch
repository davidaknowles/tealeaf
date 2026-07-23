#!/bin/bash
#SBATCH --job-name=gse233208_merge
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load "${PYTHON_MODULE}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${RUN_ROOT}/logs"
inputs=()
while IFS=$'\t' read -r batch run_accession; do
  inputs+=(
    --input
    "${run_accession}=${DATA_ROOT}/processed/${batch}/${run_accession}/alevin_quant"
  )
done < <(awk -F '\t' 'NR > 1 && $7 == 1 {print $1 "\t" $2}' \
  "${DATA_ROOT}/manifest/ena_fastqs.tsv")
"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/merge_alevin_batches.py" \
  "${inputs[@]}" \
  --output "${DATA_ROOT}/processed/merged_alevin"
"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/write_parse_primer_pairs.py" \
  --cell-rows "${DATA_ROOT}/processed/merged_alevin/quants_mat_rows.txt" \
  --rt-barcodes "${PARSE_RT_BARCODES}" \
  --output "${DATA_ROOT}/processed/merged_alevin/primer_pairs.tsv"
