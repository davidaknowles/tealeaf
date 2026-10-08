#!/bin/bash
#SBATCH --job-name=conditional_read_null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=.cache/conditional_read_null_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
pytest -q tests/test_conditional_read_odds.py
python extra_scripts/audit_conditional_read_null.py --scenario "${SLURM_ARRAY_TASK_ID:?array scenario required}" --draws 64 --nodes 21 --output-dir "analyses/conditional_local_read_null/scenario$SLURM_ARRAY_TASK_ID"
