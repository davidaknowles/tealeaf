#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail
: "${GTF:?set GTF}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

module purge
module load Python/3.12.3-GCCcore-13.3.0 HTSlib/1.21 zstd/1.5.6-GCCcore-13.3.0
source "${MAJIQ_VENV:-${HOME}/venv/majiq_cluster}/bin/activate"
mkdir -p "${OUTPUT_DIR}"
gff3="${OUTPUT_DIR}/annotation.gff3"
"${GFFREAD:-${HOME}/knowles_lab/software/bin/gffread}" \
  -F --keep-genes "${GTF}" -o "${gff3}"
majiq-v3 gff3 "${gff3}" "${OUTPUT_DIR}/annotated.splicegraph" \
  --license "${MAJIQ_LICENSE_FILE}" --overwrite
