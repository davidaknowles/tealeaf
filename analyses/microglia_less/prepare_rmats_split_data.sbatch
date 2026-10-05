#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=rmats_split_prepare
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/prepare_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/rmats_split_20261005/logs/prepare_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m pytest "${REPO_ROOT}/tests/test_rmats_split_data.py" -q
if [[ ! -e /scratch/rmats_paired_env ]]; then
  ln -s /gpfs/commons/home/daknowles/rmats_paired_env /scratch/rmats_paired_env || test -e /scratch/rmats_paired_env
fi
export R_LIBS_USER=/gpfs/commons/home/daknowles/R/x86_64-pc-linux-gnu-library/4.3
/gpfs/commons/home/daknowles/rmats_paired_env/bin/Rscript --vanilla -e 'library(PAIRADISE); cat("PAIRADISE", as.character(packageVersion("PAIRADISE")), "\n")'
python "${REPO_ROOT}/extra_scripts/prepare_rmats_split_data.py" \
  --full-manifest /gpfs/commons/home/daknowles/rmats_paired_data/rmats_celltype_contrasts.json \
  --split-manifest "${RUN_ROOT}/junction_benchmark/reproducibility/fold0/contrasts.json" \
  --split-manifest "${RUN_ROOT}/junction_benchmark/reproducibility/fold1/contrasts.json" \
  --samples "${RUN_ROOT}/junction_benchmark/pseudobulk_junctions_samples.tsv" \
  --count-dir /gpfs/commons/home/daknowles/rmats_paired_data/rmats_base_single/output \
  --output-dir "${REPO_ROOT}/.cache/rmats_split_20261005"
