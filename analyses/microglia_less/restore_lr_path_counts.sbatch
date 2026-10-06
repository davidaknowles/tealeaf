#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=lr_count_precision
#SBATCH --output=.cache/lr_ranking_20261006/precision_%j.out
#SBATCH --error=.cache/lr_ranking_20261006/precision_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python "${REPO_ROOT}/extra_scripts/restore_tilgner_path_counts.py" \
  --input "${REPO_ROOT}/analyses/tilgner_long_read_replication/pairwise_replication.tsv" \
  --candidate-cache "${RUN_ROOT}/differential/ec_block_full_pairwise_candidates.pkl" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
  --tilgner-matrix "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" \
  --gtf "${GENCODE_GTF}" \
  --output "${REPO_ROOT}/.cache/lr_ranking_20261006/pairwise_replication_precise.tsv"
