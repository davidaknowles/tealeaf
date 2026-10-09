#!/bin/bash
#SBATCH --job-name=unconditional_read_null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=.cache/unconditional_read_null_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
task_integration_args=()
if [[ "${READ_ADAPTIVE_INTEGRATION:-0}" == 1 ]]; then
    task_integration_args+=(--adaptive-integration)
fi
python extra_scripts/audit_unconditional_read_null.py --scenario "${SLURM_ARRAY_TASK_ID:?array scenario required}" --draws "${READ_NULL_DRAWS:-64}" --nodes 11 --output-dir "${READ_NULL_ROOT:-analyses/unconditional_local_read_null}/scenario$SLURM_ARRAY_TASK_ID" "${task_integration_args[@]}"
