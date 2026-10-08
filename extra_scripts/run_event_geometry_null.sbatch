#!/bin/bash
#SBATCH --job-name=event_geometry_count_null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=.cache/event_geometry_count_null_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
case "${SLURM_ARRAY_TASK_ID:?array task required}" in
    0) task_panel=strict; task_root=.cache/event_sequence_controls_L150_E001/count_null_extended/original_binary/strict ;;
    1) task_panel=biological; task_root=.cache/event_sequence_controls_L150_E001/count_null_extended/original_binary/biological ;;
    2) task_panel=nuisance; task_root=.cache/event_sequence_controls_L150_E001/count_null_extended/original_binary/nuisance ;;
    3) task_panel=mass_only; task_root=.cache/event_sequence_event_mass_null/original_binary/mass_only ;;
    4) task_panel=mass_and_within; task_root=.cache/event_sequence_event_mass_null/original_binary/mass_and_within ;;
    *) echo "Expected array task 0--4" >&2; exit 2 ;;
esac
python extra_scripts/assess_event_geometry_calibration.py count-null --training-dir .cache/event_sequence_geometry_calibration/fold0 --null-root "$task_root" --output-dir "analyses/event_sequence_geometry_count_null/$task_panel" --shard-count 16
