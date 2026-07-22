#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 04:00:00
#SBATCH -J tealeaf_polish_bench
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

DESIGN=${DESIGN:?Set DESIGN to binary or weighted}
RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
ALEVIN=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin

python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/benchmark_factorized_polish.py \
  --alevin-dir "${ALEVIN}" \
  --salmon-ref /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa \
  --primer-pairs "${RUN}/primer_pairs.tsv" \
  --design "${DESIGN}" \
  --initial-prefix "${RUN}/paired_primer_${DESIGN}_factorized_theta_" \
  --output "${RUN}/paired_primer_${DESIGN}_polish_benchmark.json" \
  --iterations 32 \
  --inner-steps 1 \
  --inner-steps 2 \
  --inner-steps 4 \
  --inner-steps 8 \
  --inner-steps 16 \
  --inner-steps 32 \
  --inner-steps 64 \
  --inner-steps 128
