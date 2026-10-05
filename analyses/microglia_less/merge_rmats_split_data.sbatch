#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=rmats_split_merge
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/merge_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/merge_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/merge_rmats_split_data.py" \
  --input-dir "${REPO_ROOT}/.cache/rmats_split_20261005" \
  --output-dir "${REPO_ROOT}/analyses/comparator_suppa_rmats/rmats"
