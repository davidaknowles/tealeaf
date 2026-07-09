#!/bin/bash
#SBATCH --job-name=microglia_fry_t2t
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/alevin_fry_t2t-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/alevin_fry_t2t-%j.err

set -euo pipefail

DATA_DIR=${DATA_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice}
SALMON_DIR=${SALMON_DIR:-salmon_spliceu}
INDEX_DIR=${INDEX_DIR:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu}
MAPPING=${MAPPING:-t2t_dedup_v2.tsv}
PERMIT_LIST=${PERMIT_LIST:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/Parse_meta/Parse_expanded_barcodes.txt}
PERMIT_DIR=${PERMIT_DIR:-${SALMON_DIR}/out_permit_known}
OUT_DIR=${OUT_DIR:-${PERMIT_DIR}/quant_t2t_dedup}
THREADS=${SLURM_CPUS_PER_TASK:-20}
MIN_READS=${MIN_READS:-100}

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate alevinfry

cd "${DATA_DIR}"

alevin-fry generate-permit-list \
  -d both \
  -i "${SALMON_DIR}" \
  --output-dir "${PERMIT_DIR}" \
  --unfiltered-pl "${PERMIT_LIST}" \
  --min-reads "${MIN_READS}"

alevin-fry collate \
  -r "${SALMON_DIR}" \
  -t "${THREADS}" \
  -i "${PERMIT_DIR}"

alevin-fry quant \
  -m "${INDEX_DIR}/${MAPPING}" \
  -i "${PERMIT_DIR}" \
  -o "${OUT_DIR}" \
  -t "${THREADS}" \
  -r cr-like-em \
  --use-mtx \
  --dump-eqclasses
