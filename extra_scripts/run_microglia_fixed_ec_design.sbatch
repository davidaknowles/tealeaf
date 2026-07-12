#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -J tealeaf_ec_design
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/fixed_ec_design-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/fixed_ec_design-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

ALEVIN=${ALEVIN:-/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin}

python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/build_fixed_ec_design.py \
  --alevin-dir "${ALEVIN}"
