#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00

set -euo pipefail
: "${BAM_MANIFEST:?set BAM_MANIFEST}"
: "${ANNOTATED_SPLICEGRAPH:?set ANNOTATED_SPLICEGRAPH}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

module purge
module load Python/3.12.3-GCCcore-13.3.0 HTSlib/1.21 zstd/1.5.6-GCCcore-13.3.0
source "${MAJIQ_VENV:-${HOME}/venv/majiq_cluster}/bin/activate"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${BAM_MANIFEST}")
IFS=$'\t' read -r sample bam <<< "${line}"
mkdir -p "${OUTPUT_DIR}"
majiq-v3 sj "${bam}" "${ANNOTATED_SPLICEGRAPH}" "${OUTPUT_DIR}/${sample}.sj" \
  --prefix "${sample}" --strandness NONE --nthreads "${SLURM_CPUS_PER_TASK:-4}" \
  --license "${MAJIQ_LICENSE_FILE}" --overwrite
