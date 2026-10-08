#!/bin/bash
#SBATCH --job-name=event_local_read_support
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=.cache/event_local_read_support_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/cpa_blackwell/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python extra_scripts/audit_event_local_read_support.py collect --output-dir .cache/event_sequence_local_read_support --shard-index "${SLURM_ARRAY_TASK_ID:?array task required}" --threads "$SLURM_CPUS_PER_TASK"
