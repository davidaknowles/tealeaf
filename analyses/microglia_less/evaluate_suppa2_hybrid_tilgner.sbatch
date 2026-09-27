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
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

HYBRID_TAG=${HYBRID_TAG:-tealeaf_hybrid_full}
OUTPUT_TAG=${OUTPUT_TAG:-tilgner_rank}
python "${REPO_ROOT}/extra_scripts/compare_suppa2_hybrid_tilgner_sign.py" \
  --native-tests "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz" \
  --hybrid-tests "${REPO_ROOT}/.cache/suppa2_statistics/${HYBRID_TAG}/merged/paired_path.tsv" \
  --tilgner-audit "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2/full_data_tilgner.tsv.gz" \
  --output-dir "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/${OUTPUT_TAG}" \
  --top-rank 200 \
  --minimum-depth 20
