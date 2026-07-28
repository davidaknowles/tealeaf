#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=microglia_count_score

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

OUTPUT_TAG=${OUTPUT_TAG:-weighted_rank32}
COUNT_TAG=${COUNT_TAG:-weighted_nb10_rank32}

python "${REPO_ROOT}/extra_scripts/score_single_cell_glm_labels.py" \
  --labels "${RUN_ROOT}/paired_primer_labels.csv" \
  --groups "${RUN_ROOT}/paired_primer_groups.csv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --standard-h5ad "${STANDARD_H5AD}" \
  --output-dir "${RUN_ROOT}/count_hierarchical_factor_scores" \
  --representation factor \
  --device cuda \
  --fit "count_pca=${RUN_ROOT}/count_glm/${COUNT_TAG}_pca_" \
  --fit "count_gene=${RUN_ROOT}/count_glm/${COUNT_TAG}_gene_" \
  --fit "count_joint=${RUN_ROOT}/count_glm/${COUNT_TAG}_joint_" \
  --fit "hierarchical_pca=${RUN_ROOT}/hierarchical/${OUTPUT_TAG}_gene_pca_" \
  --fit "hierarchical_gene=${RUN_ROOT}/hierarchical/${OUTPUT_TAG}_gene_" \
  --fit "hierarchical_joint=${RUN_ROOT}/hierarchical/${OUTPUT_TAG}_joint_"
