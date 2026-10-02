#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --job-name=suppa_split_stats
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/split_audit_%A_%a.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/split_audit_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m pytest "${REPO_ROOT}/tests/test_suppa2_primer_aware.py" -q
python "${REPO_ROOT}/extra_scripts/audit_suppa2_split_statistics.py" --run-root "${RUN_ROOT}" --repo-root "${REPO_ROOT}" --fold "${SLURM_ARRAY_TASK_ID}" --output-dir "${REPO_ROOT}/.cache/suppa2_statistics/split_audit_20261002"
