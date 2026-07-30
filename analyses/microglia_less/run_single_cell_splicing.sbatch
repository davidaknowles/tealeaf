#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=microglia_sc_splice

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

FIT_TAG=${FIT_TAG:-ablation_joint30_lr1e2_joint}
SOURCE_TAG=${SOURCE_TAG:-weighted_hierarchical_subisoform_min25}
OUTPUT_TAG=${OUTPUT_TAG:-single_cell_hierarchical_matched}
PERMUTATIONS=${PERMUTATIONS:-20}
MAX_BLOCKS=${MAX_BLOCKS:-}
SHARD_COUNT=${SHARD_COUNT:-1}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}

if (( SHARD_COUNT > 1 )); then
  OUTPUT_TAG="${OUTPUT_TAG}_shard${SHARD_INDEX}"
fi

max_blocks_args=()
if [[ -n "${MAX_BLOCKS}" ]]; then
  max_blocks_args+=(--max-blocks "${MAX_BLOCKS}")
fi

python "${REPO_ROOT}/extra_scripts/run_single_cell_splicing.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --fit-prefix "${RUN_ROOT}/hierarchical/${FIT_TAG}_" \
  --gtf "${GENCODE_GTF}" \
  --barcode-groups "${RUN_ROOT}/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
  --target-tests "${RUN_ROOT}/differential/${SOURCE_TAG}/differential_splicing.tsv" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}" \
  --min-half-umis 500 \
  --min-cells 20 \
  --min-pseudobulk-umis 100000 \
  --min-gene-umis 25 \
  --max-paths 8 \
  --max-logratio-variance 0.125 \
  --min-path-proportion 0.15 \
  --permutations "${PERMUTATIONS}" \
  --shard-index "${SHARD_INDEX}" \
  --shard-count "${SHARD_COUNT}" \
  "${max_blocks_args[@]}"
