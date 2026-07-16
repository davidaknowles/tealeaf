#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -J tealeaf_sc_score
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run

python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/score_single_cell_glm_labels.py \
  --labels "${RUN}/barcodes_to_cluster_name.csv" \
  --groups "${RUN}/barcodes_to_cluster_dx_mouse.csv" \
  --output-dir "${RUN}/single_cell_label_scores" \
  --fit "binary_admm_factorized=${RUN}/single_cell_binary_admm_factorized_theta_" \
  --fit "binary_frank_wolfe=${RUN}/single_cell_binary_frank_wolfe_theta_" \
  --fit "binary_factorized=${RUN}/single_cell_binary_factorized_theta_" \
  --fit "weighted_admm_factorized=${RUN}/single_cell_weighted_admm_factorized_theta_" \
  --fit "weighted_frank_wolfe=${RUN}/single_cell_weighted_frank_wolfe_theta_" \
  --fit "weighted_factorized=${RUN}/single_cell_weighted_factorized_theta_"
