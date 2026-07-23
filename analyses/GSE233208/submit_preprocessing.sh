#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
module load "${PYTHON_MODULE}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${DATA_ROOT}/manifest" "${RUN_ROOT}/logs"

"${PYTHON_BIN}" "${REPO_ROOT}/analyses/GSE233208/prepare_manifest.py" \
  --output "${DATA_ROOT}/manifest/ena_fastqs.tsv"

download=${DOWNLOAD_JOB:-}
if [[ -z "${download}" ]]; then
  download=$(sbatch --parsable \
    -p io --job-name=gse233208_download --cpus-per-task=8 --mem=16G \
    --time=3-00:00:00 \
    -o "${RUN_ROOT}/logs/download-%j.out" \
    -e "${RUN_ROOT}/logs/download-%j.err" \
    --wrap "bash -lc 'export PYTHONPATH=${REPO_ROOT} && cd ${REPO_ROOT} && ${PYTHON_BIN} extra_scripts/fetch_ena_fastqs.py --manifest ${DATA_ROOT}/manifest/ena_fastqs.tsv --destination ${DATA_ROOT}/fastq --workers 8'")
fi
salmon=$(sbatch --parsable --array=0-39 --dependency="afterok:${download}" \
  --job-name=gse233208_salmon \
  -o "${RUN_ROOT}/logs/salmon-%A_%a.out" \
  -e "${RUN_ROOT}/logs/salmon-%A_%a.err" \
  --export=ALL,STAGE=salmon \
  "${REPO_ROOT}/analyses/GSE233208/run_batch_stage.sbatch")
fry=$(sbatch --parsable --array=0-39 --dependency="aftercorr:${salmon}" \
  --job-name=gse233208_fry --mem=40G --time=02:00:00 \
  -o "${RUN_ROOT}/logs/fry-%A_%a.out" \
  -e "${RUN_ROOT}/logs/fry-%A_%a.err" \
  --export=ALL,STAGE=fry \
  "${REPO_ROOT}/analyses/GSE233208/run_batch_stage.sbatch")
gate=$(sbatch --parsable --dependency="afterany:${fry}" \
  -o "${RUN_ROOT}/logs/preprocessing-gate-%j.out" \
  -e "${RUN_ROOT}/logs/preprocessing-gate-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/validate_preprocessing.sbatch")
merge=$(sbatch --parsable --dependency="afterok:${gate}" \
  -o "${RUN_ROOT}/logs/merge-%j.out" \
  -e "${RUN_ROOT}/logs/merge-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/merge_batches.sbatch")

printf 'download=%s\nsalmon=%s\nfry=%s\ngate=%s\nmerge=%s\n' \
  "${download}" "${salmon}" "${fry}" "${gate}" "${merge}"
