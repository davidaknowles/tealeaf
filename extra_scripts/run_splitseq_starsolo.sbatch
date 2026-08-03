#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00

set -euo pipefail

: "${FASTQ_R1:?set FASTQ_R1}"
: "${FASTQ_R2:?set FASTQ_R2}"
: "${STAR_INDEX:?set STAR_INDEX}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${WHITELIST_1:?set WHITELIST_1}"
: "${WHITELIST_2:?set WHITELIST_2}"
: "${WHITELIST_3:?set WHITELIST_3}"

module load "${STAR_MODULE:-STAR/2.7.11b-GCC-13.2.0}"

mkdir -p "${OUTPUT_DIR}"
if [[ -s "${OUTPUT_DIR}/.tealeaf_complete" ]]; then
  exit 0
fi

read_limit_args=()
if [[ -n "${READ_LIMIT:-}" ]]; then
  read_limit_args+=(--readMapNumber "${READ_LIMIT}")
fi

STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK:-12}" \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${FASTQ_R1}" "${FASTQ_R2}" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Complex \
  --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
  --soloUMIposition 0_0_0_9 \
  --soloCBwhitelist "${WHITELIST_1}" "${WHITELIST_2}" "${WHITELIST_3}" \
  --soloCBmatchWLtype EditDist_2 \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene SJ \
  --soloCellFilter None \
  --soloUMIdedup 1MM_All \
  --soloMultiMappers Unique \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM CR CY UR UY CB UB GX GN sS sQ \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM "${BAM_SORT_RAM:-48000000000}" \
  --outFileNamePrefix "${OUTPUT_DIR}/" \
  "${read_limit_args[@]}"

test -s "${OUTPUT_DIR}/Solo.out/SJ/raw/matrix.mtx"
test -s "${OUTPUT_DIR}/Solo.out/SJ/raw/barcodes.tsv"
test -s "${OUTPUT_DIR}/Aligned.sortedByCoord.out.bam"
printf 'complete\n' > "${OUTPUT_DIR}/.tealeaf_complete"
