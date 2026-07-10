#!/bin/bash
#SBATCH --job-name=microglia_weighted_rad
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_weighted_rad-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_weighted_rad-%j.err

set -euo pipefail

FASTQ_DIR=${FASTQ_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA}
INDEX=${INDEX:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/index_dir}
OUT_DIR=${OUT_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad}
SALMON_BIN=${SALMON_BIN:-/gpfs/commons/home/daknowles/projects/salmon-v1.10.3/build/src/salmon}
THREADS=${SLURM_CPUS_PER_TASK:-20}

cd "${FASTQ_DIR}"

"${SALMON_BIN}" alevin \
  -i "${INDEX}" \
  -l A \
  -1 *_1.fastq.gz \
  -2 *_2.fastq.gz \
  -p "${THREADS}" \
  --splitseqV2 \
  --rad \
  -o "${OUT_DIR}"
