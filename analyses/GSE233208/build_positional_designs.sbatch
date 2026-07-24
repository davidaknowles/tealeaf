#!/bin/bash
#SBATCH --job-name=gse233208_build_pos
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load "${PYTHON_MODULE}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

MANIFEST=${MANIFEST:-"${DATA_ROOT}/manifest/ena_fastqs.tsv"}
ALEVIN_DIR=${ALEVIN_DIR:-"${DATA_ROOT}/processed/merged_alevin"}
OUTPUT_ROOT=${OUTPUT_ROOT:-"${DATA_ROOT}/processed/primer_positional_bias"}
quant_args=()
validation_args=()
while IFS=$'\t' read -r batch run _sample _title _experiment _experiment_title \
    read _url _bytes _md5 _filename; do
  [[ "${batch}" == batch || "${read}" != 1 ]] && continue
  root="${OUTPUT_ROOT}/by_run/${run}"
  for primer in polydt ranhex; do
    test -s "${root}/salmon_${primer}/.tealeaf_complete"
  done
  "${PYTHON_BIN}" \
    "${REPO_ROOT}/extra_scripts/validate_primer_positional_quantification.py" \
    --run-root "${root}" \
    --source-meta "${DATA_ROOT}/processed/${batch}/${run}/salmon_rad/aux_info/meta_info.json" \
    --expected-targets 243927 \
    --output "${root}/validation.json"
  validation_args+=(--validation "${root}/validation.json")
  quant_args+=(
    --polydt-quant "${root}/salmon_polydt"
    --ranhex-quant "${root}/salmon_ranhex"
  )
done < "${MANIFEST}"

"${PYTHON_BIN}" \
  "${REPO_ROOT}/extra_scripts/aggregate_primer_positional_validations.py" \
  "${validation_args[@]}" \
  --output "${OUTPUT_ROOT}/validation_summary.json"

"${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/build_salmon_primer_designs.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  "${quant_args[@]}" \
  --primer-pairs "${ALEVIN_DIR}/primer_pairs.tsv" \
  --output-report "${OUTPUT_ROOT}/design_report.json"
