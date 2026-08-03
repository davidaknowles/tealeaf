#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --time=08:00:00

set -euo pipefail
: "${GENOME_FASTA:?set GENOME_FASTA}"
: "${GTF:?set GTF}"
: "${STAR_INDEX:?set STAR_INDEX}"
: "${SJDB_OVERHANG:?set SJDB_OVERHANG}"

module load "${STAR_MODULE:-STAR/2.7.11b-GCC-13.2.0}"
mkdir -p "${STAR_INDEX}"
if [[ -s "${STAR_INDEX}/Genome" && -s "${STAR_INDEX}/SA" ]]; then
  exit 0
fi
STAR --runMode genomeGenerate \
  --runThreadN "${SLURM_CPUS_PER_TASK:-20}" \
  --genomeDir "${STAR_INDEX}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  --limitGenomeGenerateRAM "${INDEX_RAM:-170000000000}"
