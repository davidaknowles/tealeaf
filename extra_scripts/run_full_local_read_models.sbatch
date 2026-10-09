#!/bin/bash
#SBATCH --job-name=full_local_read_models
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --output=.cache/full_local_read_models_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/cpa_blackwell/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
task_index="${SLURM_ARRAY_TASK_ID:?cohort/shard array index required}"
task_cohorts=(fold0 fold1 full)
task_cohort="${task_cohorts[$((task_index / 64))]}"
task_shard="$((task_index % 64))"
python extra_scripts/full_local_read_models.py fit --output-root "${READ_MODEL_FULL_ROOT:?frozen full-model family required}" --cohort "$task_cohort" --shard-index "$task_shard" --marker-variant "${READ_MARKER_VARIANT:-all}"
