#!/bin/bash
#SBATCH --job-name=microglia_eqweights
#SBATCH --cpus-per-task=20
#SBATCH --mem=164G
#SBATCH --time=4-00:00:00
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_eqweights-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_eqweights-%j.err

set -euo pipefail

FASTQ_DIR=${FASTQ_DIR:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/SRA}
INDEX=${INDEX:-/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/index_dir}
OUT_DIR=${OUT_DIR:-/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/salmon_eqweights}
THREADS=${SLURM_CPUS_PER_TASK:-20}
RANGE_BINS=${RANGE_BINS:-4}

mkdir -p "${OUT_DIR}"
cd "${FASTQ_DIR}"

# This produces aux_info/eq_classes.txt with conditional transcript weights and
# quant.sf with Salmon effective lengths. It is a bulk Salmon pass over the
# same FASTQs used for the single-cell Alevin RAD generation.
salmon quant \
  -i "${INDEX}" \
  -l A \
  -1 *_1.fastq.gz \
  -2 *_2.fastq.gz \
  -p "${THREADS}" \
  --validateMappings \
  --rangeFactorizationBins "${RANGE_BINS}" \
  --dumpEq \
  --dumpEqWeights \
  -o "${OUT_DIR}"
