#!/bin/bash
#SBATCH --job-name=event_geometry_calibration
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=.cache/event_geometry_calibration_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
task_cohort="${COHORT:-fold${SLURM_ARRAY_TASK_ID:-0}}"
task_source=".cache/event_sequence_controls_L150_E001"
task_output=".cache/event_sequence_geometry_calibration"
if [[ "$task_cohort" == full ]]; then
    task_merged="$task_source/assess_mixed_score_reference/original_binary/full/merged"
elif [[ "$task_cohort" == fold0 || "$task_cohort" == fold1 ]]; then
    task_fold="${task_cohort#fold}"
    task_merged="$task_source/assess_mixed_score_reference_split_only/original_binary/$task_fold/merged"
else
    echo "Expected full, fold0 or fold1" >&2
    exit 2
fi
python extra_scripts/assess_event_geometry_calibration.py cohort --cohort-root "$task_source/mixed_score_reference/original_binary/$task_cohort" --merged-dir "$task_merged" --output-dir "$task_output/$task_cohort" --shard-count 64
