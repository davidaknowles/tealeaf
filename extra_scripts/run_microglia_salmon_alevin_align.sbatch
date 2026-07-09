#!/bin/bash
#SBATCH --job-name=microglia_alevin_align
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_alevin_align-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_alevin_align-%j.err

set -euo pipefail

FASTQ_DIR=${FASTQ_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA}
INDEX=${INDEX:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/index_dir}
OUT_DIR=${OUT_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu}
THREADS=${SLURM_CPUS_PER_TASK:-20}

cd "${FASTQ_DIR}"

salmon alevin \
  -i "${INDEX}" \
  -l A \
  -1 *_1.fastq.gz \
  -2 *_2.fastq.gz \
  -p "${THREADS}" \
  --splitseqV2 \
  --sketch \
  -o "${OUT_DIR}"
