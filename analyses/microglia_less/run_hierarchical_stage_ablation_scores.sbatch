#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=microglia_hier_stage_score

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

python "${REPO_ROOT}/extra_scripts/score_single_cell_glm_labels.py" \
  --labels "${RUN_ROOT}/paired_primer_labels.csv" \
  --groups "${RUN_ROOT}/paired_primer_groups.csv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --standard-h5ad "${STANDARD_H5AD}" \
  --output-dir "${RUN_ROOT}/hierarchical_stage_ablation_scores" \
  --representation factor \
  --device cuda \
  --fit "three_stage=${RUN_ROOT}/hierarchical/weighted_rank32_joint_" \
  --fit "gene_joint=${RUN_ROOT}/hierarchical/ablation_gene10_joint20_joint_" \
  --fit "isoform_joint=${RUN_ROOT}/hierarchical/ablation_isoform10_joint20_joint_" \
  --fit "joint_lr1e3=${RUN_ROOT}/hierarchical/ablation_joint30_lr1e3_joint_" \
  --fit "joint_lr3e3=${RUN_ROOT}/hierarchical/ablation_joint30_lr3e3_joint_" \
  --fit "joint_lr1e2=${RUN_ROOT}/hierarchical/ablation_joint30_lr1e2_joint_"
