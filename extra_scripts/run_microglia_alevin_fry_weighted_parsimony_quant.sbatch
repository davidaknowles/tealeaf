#!/bin/bash
#SBATCH --job-name=microglia_weighted_fry
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G
#SBATCH --time=4-00:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/alevin_fry_weighted_parsimony-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/alevin_fry_weighted_parsimony-%j.err

set -euo pipefail

DATA_DIR=${DATA_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice}
SALMON_DIR=${SALMON_DIR:-salmon_spliceu_weighted_rad}
INDEX_DIR=${INDEX_DIR:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu}
MAPPING=${MAPPING:-t2t_dedup_v2.tsv}
PERMIT_LIST=${PERMIT_LIST:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta/Parse_expanded_barcodes.txt}
PERMIT_DIR=${PERMIT_DIR:-${SALMON_DIR}/out_permit_known_weighted}
OUT_DIR=${OUT_DIR:-${PERMIT_DIR}/quant_t2t_dedup_parsimony_em}
ALEVIN_FRY_BIN=${ALEVIN_FRY_BIN:-/gpfs/commons/home/daknowles/projects/alevin-fry/target/release/alevin-fry}
THREADS=${SLURM_CPUS_PER_TASK:-20}
MIN_READS=${MIN_READS:-100}
RESOLUTION=${RESOLUTION:-parsimony-em}

cd "${DATA_DIR}"

if [ ! -s "${PERMIT_DIR}/permit_freq.bin" ]; then
  "${ALEVIN_FRY_BIN}" generate-permit-list \
    -d both \
    -i "${SALMON_DIR}" \
    --output-dir "${PERMIT_DIR}" \
    --unfiltered-pl "${PERMIT_LIST}" \
    --min-reads "${MIN_READS}"
else
  echo "Using existing permit list in ${PERMIT_DIR}"
fi

if [ ! -s "${PERMIT_DIR}/map.collated.rad" ]; then
  "${ALEVIN_FRY_BIN}" collate \
    -r "${SALMON_DIR}" \
    -t "${THREADS}" \
    -i "${PERMIT_DIR}"
else
  echo "Using existing collated RAD in ${PERMIT_DIR}"
fi

if [ -e "${OUT_DIR}" ]; then
  echo "Refusing to overwrite existing quant output: ${OUT_DIR}" >&2
  exit 1
fi

"${ALEVIN_FRY_BIN}" quant \
  -m "${INDEX_DIR}/${MAPPING}" \
  -i "${PERMIT_DIR}" \
  -o "${OUT_DIR}" \
  -t "${THREADS}" \
  -r "${RESOLUTION}" \
  --use-mtx \
  --dump-eqclasses
