#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --job-name=omnibus_inference_assessment
#SBATCH --output=.cache/paired_reporting_20261006/omnibus_assess_%j.out
#SBATCH --error=.cache/paired_reporting_20261006/omnibus_assess_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=()
if [[ "${SPLIT_ONLY:-0}" == 1 ]]; then
  ARGS+=(--split-only)
fi
python -m extra_scripts.assess_omnibus_inference_audit \
  --cache "${AUDIT_CACHE}" --control-cache "${CONTROL_CACHE}" --name "${VARIANT}" \
  --output-dir "${AUDIT_OUTPUT}" --run-root "${RUN_ROOT}" \
  --matrix-dir "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" --gtf "${GENCODE_GTF}" "${ARGS[@]}"
