#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --job-name=assess_joint_path_dm
#SBATCH --output=.cache/paired_reporting_20261006/joint_paired_assess_%j.out
#SBATCH --error=.cache/paired_reporting_20261006/joint_paired_assess_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m extra_scripts.assess_joint_path_dm \
  --cache "${AUDIT_CACHE}" --output-dir "${AUDIT_OUTPUT}" \
  --matrix-dir "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" --gtf "${GENCODE_GTF}" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz"
