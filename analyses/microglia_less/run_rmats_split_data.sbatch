#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-379%16
#SBATCH --job-name=rmats_split_paired
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/fit_%A_%a.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/fit_%A_%a.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load rmats/4.1.2
BASE="${REPO_ROOT}/.cache/rmats_split_20261005"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ ! -e /scratch/rmats_paired_env ]]; then
  ln -s /gpfs/commons/home/daknowles/rmats_paired_env /scratch/rmats_paired_env || test -e /scratch/rmats_paired_env
fi
export R_LIBS_USER=/gpfs/commons/home/daknowles/R/x86_64-pc-linux-gnu-library/4.3
/gpfs/commons/home/daknowles/rmats_paired_env/bin/python "${REPO_ROOT}/extra_scripts/run_rmats_comparison.py" \
  --contrasts "${BASE}/fold0_contrasts.json" --contrasts "${BASE}/fold1_contrasts.json" \
  --contrast-index "${SLURM_ARRAY_TASK_ID}" --output "${BASE}/results/${SLURM_ARRAY_TASK_ID}.tsv.gz" \
  --work-dir "${BASE}/work/${SLURM_ARRAY_TASK_ID}" \
  --old-output /gpfs/commons/home/daknowles/rmats_paired_data/rmats_base_single/output \
  --jobs 1 --threads-per-job "${SLURM_CPUS_PER_TASK}" --paired-stats --paired-jcec-only \
  --rscript /gpfs/commons/home/daknowles/rscript_bin/Rscript --r-libs-user "${R_LIBS_USER}"
