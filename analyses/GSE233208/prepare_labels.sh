#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"

module load python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

metadata="${DATA_ROOT}/metadata/GSE233208_Human_snRNA-Seq_ADDS_metadata.tsv.gz"
"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/prepare_reference_labels.py" \
  --metadata "${metadata}" \
  --label-column annotation \
  --group-column CaseNum \
  --batch-map batch1=Batch1 \
  --batch-map batch2=Batch4 \
  --batch-map batch3=Batch2 \
  --batch-map batch4=Batch3 \
  --batch-map batch5=Batch5 \
  --ena-manifest "${DATA_ROOT}/manifest/ena_fastqs.tsv" \
  --parse-rt-barcodes "${PARSE_RT_BARCODES}" \
  --labels-output "${DATA_ROOT}/metadata/reference_annotation.csv" \
  --groups-output "${DATA_ROOT}/metadata/reference_donor.csv"
