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

quantification_files_exist() {
  [[
    -s "${OUT_DIR}/alevin/quants_mat_rows.txt"
    && -s "${OUT_DIR}/alevin/quants_mat_cols.txt"
    && -s "${OUT_DIR}/alevin/geqc_counts.mtx"
    && -s "${OUT_DIR}/alevin/gene_eqclass.txt.gz"
    && -s "${OUT_DIR}/alevin/gene_eqclass_probs.tsv.gz"
  ]]
}

if [[ -s "${COMPLETE_MARKER}" ]] && quantification_files_exist; then
  echo "Reusing completed alevin-fry output: ${OUT_DIR}"
  exit 0
fi

if [[ ! -s "${PERMIT_DIR}/permit_freq.bin" ]]; then
  "${ALEVIN_FRY_BIN}" generate-permit-list \
    -d both -i "${SALMON_DIR}" -o "${PERMIT_DIR}" \
    --unfiltered-pl "${PERMIT_LIST}" --min-reads "${MIN_READS}"
fi
if [[ ! -s "${PERMIT_DIR}/map.collated.rad" ]]; then
  "${ALEVIN_FRY_BIN}" collate \
    -r "${SALMON_DIR}" -t "${THREADS}" -i "${PERMIT_DIR}"
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
printf 'alevin_fry_complete\n' > "${COMPLETE_MARKER}"
