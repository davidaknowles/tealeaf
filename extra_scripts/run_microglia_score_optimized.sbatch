#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -J tealeaf_score_optimized
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
  --transcript-to-gene /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu_t2g.tsv \
  --standard-h5ad /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_per_sublib/t2g_permit_known_combined.h5ad \
  --output-dir "${RUN}/single_cell_log_gene_pca_scores_minibatch" \
  --device cuda \
  --fit "binary_factorized_cv_minibatch=${RUN}/single_cell_binary_factorized_cv_minibatch_theta_" \
  --fit "weighted_factorized_cv_minibatch=${RUN}/single_cell_weighted_factorized_cv_minibatch_theta_" \
  --fit "binary_admm_factorized_cv_minibatch=${RUN}/single_cell_binary_admm_factorized_cv_minibatch_theta_" \
  --fit "weighted_admm_factorized_cv_minibatch=${RUN}/single_cell_weighted_admm_factorized_cv_minibatch_theta_"
