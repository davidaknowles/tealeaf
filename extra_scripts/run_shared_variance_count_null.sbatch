#!/bin/bash
#SBATCH --job-name=shared_variance_count_null
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=.cache/shared_variance_count_null_%A_%a.log

set -euo pipefail
module load Python/3.12.3-GCCcore-13.3.0
source /gpfs/commons/home/daknowles/venv/jax/bin/activate
cd /gpfs/commons/home/daknowles/projects/tealeaf
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
case "${SLURM_ARRAY_TASK_ID:?array task required}" in
    0) task_panel=strict ;;
    1) task_panel=biological ;;
    2) task_panel=nuisance ;;
    3) task_panel=mass_only ;;
    4) task_panel=mass_and_within ;;
    *) exit 2 ;;
esac
python extra_scripts/audit_shared_variance_count_null.py --source-assessment "analyses/event_sequence_geometry_count_null/$task_panel" --output-dir "analyses/event_sequence_shared_variance_null/$task_panel"
