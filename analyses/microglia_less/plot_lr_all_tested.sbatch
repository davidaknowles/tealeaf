#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=lr_rank_plot
#SBATCH --output=.cache/lr_ranking_20261006/plot_%j.out
#SBATCH --error=.cache/lr_ranking_20261006/plot_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/plot_tilgner_method_replication.py" \
  --summary "${REPO_ROOT}/analyses/tilgner_junction_replication/method_comparison_summary.tsv" \
  --output "${REPO_ROOT}/.cache/lr_ranking_20261006/discovery_comparison.pdf" \
  --rank-input "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.tsv" \
  --rank-replacement "${REPO_ROOT}/analyses/tilgner_junction_replication/all_tested/junction_replication.tsv.gz" \
  --rank-replacement "${REPO_ROOT}/analyses/tilgner_junction_replication/all_tested/majiq/junction_replication.tsv.gz" \
  --rank-replacement "${REPO_ROOT}/analyses/tilgner_long_read_replication/ranking_audit/production_rank.tsv" \
  --rank-table "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.tsv" \
  --rank-summary "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_rank_summary.tsv" \
  --rank-output "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.pdf"
