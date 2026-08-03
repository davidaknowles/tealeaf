#!/bin/bash
set -euo pipefail

: "${INPUT_BAM:?set INPUT_BAM to a pseudobulk-RG-tagged BAM}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"

module load "${SAMTOOLS_MODULE:-SAMtools/1.23.1}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "/scratch/${USER}"
samtools split -@ "${SLURM_CPUS_PER_TASK:-4}" \
  -f "${OUTPUT_DIR}/%!.unsorted.bam" "${INPUT_BAM}"
for bam in "${OUTPUT_DIR}"/*.unsorted.bam; do
  output="${bam%.unsorted.bam}.bam"
  samtools sort -@ "${SLURM_CPUS_PER_TASK:-4}" -o "${output}" "${bam}"
  samtools index -@ "${SLURM_CPUS_PER_TASK:-4}" "${output}"
  mv "${bam}" "/scratch/${USER}/$(basename "${bam}").$$"
done
