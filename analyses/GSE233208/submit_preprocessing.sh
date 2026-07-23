#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
module load python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${DATA_ROOT}/manifest" "${RUN_ROOT}/logs"

python "${REPO_ROOT}/analyses/GSE233208/prepare_manifest.py" \
  --output "${DATA_ROOT}/manifest/ena_fastqs.tsv"

download=$(sbatch --parsable \
  -p io --job-name=gse233208_download --cpus-per-task=8 --mem=16G \
  --time=3-00:00:00 \
  -o "${RUN_ROOT}/logs/download-%j.out" \
  -e "${RUN_ROOT}/logs/download-%j.err" \
  --wrap "bash -lc 'module load python/3.10.8-GCCcore-12.2.0 && export PYTHONPATH=${REPO_ROOT} && cd ${REPO_ROOT} && python extra_scripts/fetch_ena_fastqs.py --manifest ${DATA_ROOT}/manifest/ena_fastqs.tsv --destination ${DATA_ROOT}/fastq --workers 8'")
salmon=$(sbatch --parsable --array=0-4 --dependency="afterok:${download}" \
  --job-name=gse233208_salmon \
  -o "${RUN_ROOT}/logs/salmon-%A_%a.out" \
  -e "${RUN_ROOT}/logs/salmon-%A_%a.err" \
  --export=ALL,STAGE=salmon \
  "${REPO_ROOT}/analyses/GSE233208/run_batch_stage.sbatch")
fry=$(sbatch --parsable --array=0-4 --dependency="afterok:${salmon}" \
  --job-name=gse233208_fry --mem=96G \
  -o "${RUN_ROOT}/logs/fry-%A_%a.out" \
  -e "${RUN_ROOT}/logs/fry-%A_%a.err" \
  --export=ALL,STAGE=fry \
  "${REPO_ROOT}/analyses/GSE233208/run_batch_stage.sbatch")
merge=$(sbatch --parsable --dependency="afterok:${fry}" \
  -o "${RUN_ROOT}/logs/merge-%j.out" \
  -e "${RUN_ROOT}/logs/merge-%j.err" \
  "${REPO_ROOT}/analyses/GSE233208/merge_batches.sbatch")

printf 'download=%s\nsalmon=%s\nfry=%s\nmerge=%s\n' \
  "${download}" "${salmon}" "${fry}" "${merge}"
