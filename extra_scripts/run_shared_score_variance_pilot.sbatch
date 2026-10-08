#!/bin/bash
#SBATCH --job-name=shared_score_variance_pilot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=.cache/shared_score_variance_pilot_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
task_fold="${SLURM_ARRAY_TASK_ID:?array task required}"
[[ "$task_fold" == 0 || "$task_fold" == 1 ]] || exit 2
task_source=.cache/event_sequence_controls_L150_E001
python extra_scripts/audit_shared_score_variance.py --cohort-root "$task_source/mixed_score_reference/original_binary/fold$task_fold" --merged-dir "$task_source/assess_mixed_score_reference_split_only/original_binary/$task_fold/merged" --diagnostic-dir "analyses/event_sequence_subject_influence/fold$task_fold" --output-dir "analyses/event_sequence_shared_variance_pilot/fold$task_fold" --shard-count 64
