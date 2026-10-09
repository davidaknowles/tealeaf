#!/bin/bash
#SBATCH --job-name=library_read_model_pilot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=.cache/library_read_model_pilot_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
task_integration_args=()
if [[ "${READ_ADAPTIVE_INTEGRATION:-0}" == 1 ]]; then
    task_integration_args+=(--adaptive-integration)
fi
python extra_scripts/pilot_library_read_models.py --source "${READ_SUPPORT_PUBLIC:-analyses/event_sequence_local_read_support_library}" --shard-count 16 --shard-index "${SLURM_ARRAY_TASK_ID:?pilot shard required}" --output-root "${READ_MODEL_PILOT_ROOT:-.cache/event_sequence_local_read_model_pilot}" "${task_integration_args[@]}"
