#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=merge_suppa_tealeaf_full
#SBATCH --output=.cache/suppa2_statistics/hybrid-full-merge-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-full-merge-%j.err

set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

OUTPUT_TAG=${OUTPUT_TAG:-tealeaf_hybrid_full}
BASE="${REPO_ROOT}/.cache/suppa2_statistics/${OUTPUT_TAG}"
python "${REPO_ROOT}/extra_scripts/merge_paired_path_test.py" \
  --shards "${BASE}/shards"/shard_* \
  --output-dir "${BASE}/merged" \
  --moderate-variances \
  --calibration empirical
