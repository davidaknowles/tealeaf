#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -J tealeaf_primer_design
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
ALEVIN=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin
python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/build_primer_ec_designs.py \
  --alevin-dir "${ALEVIN}" \
  --pairs "${RUN}/primer_pairs.tsv"
