#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
task=${SLURM_ARRAY_TASK_ID:?array index is required}
manifest="${DATA_ROOT}/manifest/ena_fastqs.tsv"
row=$(awk -F '\t' -v target="$((task + 1))" '
  NR > 1 && $7 == 1 {seen += 1; if (seen == target) {print $1 "\t" $2; exit}}
' "${manifest}")
if [[ -z "${row}" ]]; then
  echo "No manifest run for array index ${task}" >&2
  exit 1
fi
IFS=$'\t' read -r batch run_accession <<< "${row}"
stage=${STAGE:?set STAGE to salmon or fry}
run_dir="${DATA_ROOT}/processed/${batch}/${run_accession}"
mkdir -p "${RUN_ROOT}/logs" "${run_dir}"

if [[ "${stage}" == salmon ]]; then
  export FASTQ_R1="${DATA_ROOT}/fastq/${batch}/${run_accession}_1.fastq.gz"
  export FASTQ_R2="${DATA_ROOT}/fastq/${batch}/${run_accession}_2.fastq.gz"
  export INDEX="${SALMON_REF}/index_dir"
  export OUT_DIR="${run_dir}/salmon_rad"
  export SALMON_BIN
  bash "${REPO_ROOT}/extra_scripts/run_parse_salmon_rad.sbatch"
elif [[ "${stage}" == fry ]]; then
  export SALMON_DIR="${run_dir}/salmon_rad"
  export INDEX_DIR="${SALMON_REF}"
  export PERMIT_LIST="${PARSE_PERMIT_LIST}"
  export OUT_DIR="${run_dir}/alevin_quant"
  export ALEVIN_FRY_BIN
  bash "${REPO_ROOT}/extra_scripts/run_parse_alevin_fry.sbatch"
else
  echo "Unknown stage: ${stage}" >&2
  exit 1
fi
