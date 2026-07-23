#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
batch=${BATCHES[${SLURM_ARRAY_TASK_ID:?array index is required}]}
stage=${STAGE:?set STAGE to salmon or fry}
mkdir -p "${RUN_ROOT}/logs" "${DATA_ROOT}/processed/${batch}"

if [[ "${stage}" == salmon ]]; then
  export FASTQ_DIR="${DATA_ROOT}/fastq/${batch}"
  export INDEX="${SALMON_REF}/index_dir"
  export OUT_DIR="${DATA_ROOT}/processed/${batch}/salmon_rad"
  export SALMON_BIN
  bash "${REPO_ROOT}/extra_scripts/run_parse_salmon_rad.sbatch"
elif [[ "${stage}" == fry ]]; then
  export SALMON_DIR="${DATA_ROOT}/processed/${batch}/salmon_rad"
  export INDEX_DIR="${SALMON_REF}"
  export PERMIT_LIST="${PARSE_PERMIT_LIST}"
  export OUT_DIR="${DATA_ROOT}/processed/${batch}/alevin_quant"
  export ALEVIN_FRY_BIN
  bash "${REPO_ROOT}/extra_scripts/run_parse_alevin_fry.sbatch"
else
  echo "Unknown stage: ${stage}" >&2
  exit 1
fi
