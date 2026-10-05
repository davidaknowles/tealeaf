#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=split_coverage
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/split_coverage_audit/audit_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/split_coverage_audit/audit_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m pytest "${REPO_ROOT}/tests/test_replication_audit.py" -q
ARGS=()
if [[ "${WITH_DIRECTIONS:-0}" == 1 ]]; then
  ARGS+=(--directions-dir "${REPO_ROOT}/.cache/split_coverage_audit/archived_directions")
fi
python "${REPO_ROOT}/extra_scripts/audit_split_coverage_direction.py" \
  --run-root "${RUN_ROOT}" --repo-root "${REPO_ROOT}" \
  --output-dir "${REPO_ROOT}/analyses/split_coverage_direction" --junction-directions "${ARGS[@]}"
