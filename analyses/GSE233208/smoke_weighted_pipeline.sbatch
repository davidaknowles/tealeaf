#!/bin/bash
#SBATCH --job-name=gse233208_real_smoke
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load SeqKit/2.9.0
module load python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

RUN_ACCESSION=${RUN_ACCESSION:-SRR24710695}
BATCH=${BATCH:-batch1}
READ_PAIRS=${READ_PAIRS:-1000000}
SMOKE_ROOT=${SMOKE_ROOT:-"${RUN_ROOT}/real_data_smoke"}
mkdir -p "${SMOKE_ROOT}/fastq" "${RUN_ROOT}/logs"

for read in 1 2; do
  source_fastq="${DATA_ROOT}/fastq/${BATCH}/${RUN_ACCESSION}_${read}.fastq.gz"
  smoke_fastq="${SMOKE_ROOT}/fastq/${RUN_ACCESSION}_${read}.fastq.gz"
  if [[ ! -s "${smoke_fastq}" ]]; then
    seqkit head -n "${READ_PAIRS}" "${source_fastq}" -o "${smoke_fastq}"
  fi
done

export FASTQ_R1="${SMOKE_ROOT}/fastq/${RUN_ACCESSION}_1.fastq.gz"
export FASTQ_R2="${SMOKE_ROOT}/fastq/${RUN_ACCESSION}_2.fastq.gz"
export INDEX="${SALMON_REF}/index_dir"
export OUT_DIR="${SMOKE_ROOT}/salmon_rad"
export SALMON_BIN
bash "${REPO_ROOT}/extra_scripts/run_parse_salmon_rad.sbatch"

export SALMON_DIR="${SMOKE_ROOT}/salmon_rad"
export INDEX_DIR="${SALMON_REF}"
export PERMIT_LIST="${PARSE_PERMIT_LIST}"
export OUT_DIR="${SMOKE_ROOT}/alevin_quant"
export ALEVIN_FRY_BIN
export MIN_READS=5
bash "${REPO_ROOT}/extra_scripts/run_parse_alevin_fry.sbatch"

python "${REPO_ROOT}/extra_scripts/validate_alevin_quantification.py" \
  --alevin-dir "${SMOKE_ROOT}/alevin_quant" \
  --min-cell-umis 1 \
  --max-total-molecules "${READ_PAIRS}" \
  --output "${SMOKE_ROOT}/validation.json"
