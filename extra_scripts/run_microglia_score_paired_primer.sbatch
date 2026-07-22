#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -J tealeaf_score_paired
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/score_single_cell_glm_labels.py \
  --labels "${RUN}/paired_primer_labels.csv" \
  --groups "${RUN}/paired_primer_groups.csv" \
  --transcript-to-gene /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu_t2g.tsv \
  --standard-h5ad /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_per_sublib/t2g_permit_known_combined.h5ad \
  --output-dir "${RUN}/paired_primer_log_gene_pca_scores" \
  --device cuda \
  --fit "paired_binary_accelerated=${RUN}/paired_primer_binary_factorized_accelerated_theta_" \
  --fit "paired_weighted_accelerated=${RUN}/paired_primer_weighted_factorized_accelerated_theta_"
