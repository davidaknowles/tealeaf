#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=eval_hybrid_tilgner
#SBATCH --output=.cache/suppa2_statistics/hybrid-tilgner-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-tilgner-%j.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

BASE="${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/full_data"
python "${REPO_ROOT}/extra_scripts/assess_event_tilgner_replication.py" \
  --tests "${BASE}/full_data_tests.tsv.gz" \
  --events "${BASE}/event_catalog.tsv.gz" \
  --tilgner-matrix "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" \
  --gtf "${GENCODE_GTF}" \
  --output "${BASE}/full_data_tilgner.tsv.gz" \
  --summary "${BASE}/full_data_tilgner_summary.tsv" \
  --top-per-contrast 0 \
  --minimum-depth 20

python "${REPO_ROOT}/extra_scripts/plot_tilgner_method_replication.py" \
  --summary "${REPO_ROOT}/analyses/tilgner_junction_replication/method_comparison_summary.tsv" \
  --output "${REPO_ROOT}/.cache/suppa2_statistics/hybrid-method-comparison.pdf" \
  --rank-input "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.tsv" \
  --hybrid-replication "${BASE}/full_data_tilgner.tsv.gz" \
  --rank-table "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.tsv" \
  --rank-output "${REPO_ROOT}/analyses/tilgner_junction_replication/sign_agreement_by_rank.pdf" \
  --max-rank 200
