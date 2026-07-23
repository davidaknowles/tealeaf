#!/bin/bash
#SBATCH --job-name=parse_alevin_fry
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4-00:00:00

set -euo pipefail

: "${SALMON_DIR:?set SALMON_DIR}"
: "${INDEX_DIR:?set INDEX_DIR}"
: "${PERMIT_LIST:?set PERMIT_LIST}"
: "${OUT_DIR:?set OUT_DIR}"
: "${ALEVIN_FRY_BIN:?set ALEVIN_FRY_BIN}"
MAPPING=${MAPPING:-t2t_dedup_v2.tsv}
PERMIT_DIR=${PERMIT_DIR:-"${SALMON_DIR}/permit_known"}
MIN_READS=${MIN_READS:-100}
THREADS=${SLURM_CPUS_PER_TASK:-8}
COMPLETE_MARKER="${OUT_DIR}/.tealeaf_complete"
VALIDATION_REPORT="${OUT_DIR}/validation.json"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
PYTHON_BIN=${PYTHON_BIN:-$(command -v python3 || command -v python || true)}
if [[ -z "${PYTHON_BIN}" ]]; then
  echo "Could not resolve a Python interpreter for validation" >&2
  exit 1
fi

quantification_files_exist() {
  [[
    -s "${OUT_DIR}/alevin/quants_mat_rows.txt"
    && -s "${OUT_DIR}/alevin/quants_mat_cols.txt"
    && -s "${OUT_DIR}/alevin/geqc_counts.mtx"
    && -s "${OUT_DIR}/alevin/gene_eqclass.txt.gz"
    && -s "${OUT_DIR}/alevin/gene_eqclass_probs.tsv.gz"
  ]]
}

if [[ ! -s "${PERMIT_DIR}/permit_freq.bin" ]]; then
  "${ALEVIN_FRY_BIN}" generate-permit-list \
    -d both -i "${SALMON_DIR}" -o "${PERMIT_DIR}" \
    --unfiltered-pl "${PERMIT_LIST}" --min-reads "${MIN_READS}"
fi
if [[ ! -s "${PERMIT_DIR}/map.collated.rad" ]]; then
  "${ALEVIN_FRY_BIN}" collate \
    -r "${SALMON_DIR}" -t "${THREADS}" -i "${PERMIT_DIR}"
fi
if quantification_files_exist; then
  if [[ ! -s "${VALIDATION_REPORT}" ]]; then
    PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" "${PYTHON_BIN}" \
      "${SCRIPT_DIR}/validate_alevin_quantification.py" \
      --alevin-dir "${OUT_DIR}" \
      --min-cell-umis 1 \
      --salmon-meta "${SALMON_DIR}/aux_info/meta_info.json" \
      --output "${VALIDATION_REPORT}"
  fi
  if [[ ! -s "${COMPLETE_MARKER}" ]]; then
    printf 'alevin_fry_complete\n' > "${COMPLETE_MARKER}"
  fi
  echo "Reusing validated alevin-fry output: ${OUT_DIR}"
  exit 0
fi

if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing to overwrite existing output: ${OUT_DIR}" >&2
  exit 1
fi
"${ALEVIN_FRY_BIN}" quant \
  -m "${INDEX_DIR}/${MAPPING}" \
  -i "${PERMIT_DIR}" \
  -o "${OUT_DIR}" \
  -t "${THREADS}" \
  -r parsimony-em \
  --use-mtx \
  --dump-eqclasses
if ! quantification_files_exist; then
  echo "alevin-fry completed without required matrix/probability outputs" >&2
  exit 1
fi
PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" "${PYTHON_BIN}" \
  "${SCRIPT_DIR}/validate_alevin_quantification.py" \
  --alevin-dir "${OUT_DIR}" \
  --min-cell-umis 1 \
  --salmon-meta "${SALMON_DIR}/aux_info/meta_info.json" \
  --output "${VALIDATION_REPORT}"
printf 'alevin_fry_complete\n' > "${COMPLETE_MARKER}"
