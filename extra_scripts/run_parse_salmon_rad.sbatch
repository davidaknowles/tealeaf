#!/bin/bash
#SBATCH --job-name=parse_salmon_rad
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00

set -euo pipefail

: "${INDEX:?set INDEX}"
: "${OUT_DIR:?set OUT_DIR}"
: "${SALMON_BIN:?set SALMON_BIN}"
THREADS=${SLURM_CPUS_PER_TASK:-20}
COMPLETE_MARKER="${OUT_DIR}/.tealeaf_complete"

if [[ -s "${COMPLETE_MARKER}" && -s "${OUT_DIR}/map.rad" ]]; then
  echo "Reusing completed Salmon RAD output: ${OUT_DIR}"
  exit 0
fi
if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing incomplete or unmarked Salmon output: ${OUT_DIR}" >&2
  exit 1
fi

if [[ -n "${FASTQ_R1:-}" || -n "${FASTQ_R2:-}" ]]; then
  : "${FASTQ_R1:?set both FASTQ_R1 and FASTQ_R2}"
  : "${FASTQ_R2:?set both FASTQ_R1 and FASTQ_R2}"
  READ1=("${FASTQ_R1}")
  READ2=("${FASTQ_R2}")
else
  : "${FASTQ_DIR:?set FASTQ_DIR or explicit FASTQ_R1/FASTQ_R2}"
  mapfile -t READ1 < <(find "${FASTQ_DIR}" -maxdepth 1 -type f \( -name '*_1.fastq.gz' -o -name '*_1.fq.gz' \) | sort)
  mapfile -t READ2 < <(find "${FASTQ_DIR}" -maxdepth 1 -type f \( -name '*_2.fastq.gz' -o -name '*_2.fq.gz' \) | sort)
fi
if (( ${#READ1[@]} == 0 || ${#READ1[@]} != ${#READ2[@]} )); then
  echo "Expected one or more matching paired FASTQs" >&2
  exit 1
fi

"${SALMON_BIN}" alevin \
  -i "${INDEX}" \
  -l A \
  -1 "${READ1[@]}" \
  -2 "${READ2[@]}" \
  -p "${THREADS}" \
  --splitseqV2 \
  --rad \
  -o "${OUT_DIR}"
if [[ ! -s "${OUT_DIR}/map.rad" ]]; then
  echo "Salmon completed without a nonempty map.rad" >&2
  exit 1
fi
printf 'salmon_rad_complete\n' > "${COMPLETE_MARKER}"
