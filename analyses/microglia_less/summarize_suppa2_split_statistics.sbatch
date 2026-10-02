#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=suppa_stats_summary
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/split_summary_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/split_summary_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/summarize_suppa2_split_statistics.py" --run-root "${RUN_ROOT}" --repo-root "${REPO_ROOT}" --audit-dir "${REPO_ROOT}/.cache/suppa2_statistics/split_audit_20261002" --output-dir "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2/statistics_audit"
