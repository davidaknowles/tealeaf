#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail
: "${BAM_LIST:?file containing one RG-tagged spliced BAM per line}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"

module load "${SAMTOOLS_MODULE:-SAMtools/1.23.1}"
mkdir -p "${OUTPUT_DIR}" "/scratch/${USER}"
combined="${OUTPUT_DIR}/all_pseudobulks.unsorted.bam"
mapfile -t bams < "${BAM_LIST}"
samtools merge -c -@ "${SLURM_CPUS_PER_TASK:-8}" -u -o "${combined}" "${bams[@]}"
samtools split -@ "${SLURM_CPUS_PER_TASK:-8}" \
  -f "${OUTPUT_DIR}/%!.unsorted.bam" "${combined}"
export OUTPUT_DIR
find "${OUTPUT_DIR}" -maxdepth 1 -name '*.unsorted.bam' ! -name 'all_pseudobulks.unsorted.bam' -print0 | \
  xargs -0 -P "${SLURM_CPUS_PER_TASK:-8}" -n 1 bash -c '
    bam=$1
    output=${bam%.unsorted.bam}.bam
    samtools sort -@ 1 -o "${output}" "${bam}"
    samtools index "${output}"
    mv "${bam}" "/scratch/${USER}/$(basename "${bam}").$$"
  ' _
mv "${combined}" "/scratch/${USER}/$(basename "${combined}").$$"
