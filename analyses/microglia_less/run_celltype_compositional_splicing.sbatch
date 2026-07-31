#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=microglia_ct_comp

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

SOURCE_TAG=${SOURCE_TAG:-weighted_hierarchical_subisoform_min25}
OUTPUT_TAG=${OUTPUT_TAG:-celltype_compositional_identifiable}
MIN_PATH_PROPORTION=${MIN_PATH_PROPORTION:-0.05}
MAX_PATHS=${MAX_PATHS:-0}
MAX_LOGRATIO_VARIANCE=${MAX_LOGRATIO_VARIANCE:-inf}
PERMUTATIONS=${PERMUTATIONS:-500}
CELLTYPE_TEST=${CELLTYPE_TEST:-omnibus}
LIKELIHOOD=${LIKELIHOOD:-dirichlet-multinomial}
SHARD_COUNT=${SHARD_COUNT:-1}
SHARD_INDEX=${SHARD_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}

if (( SHARD_COUNT > 1 )); then
  OUTPUT_TAG="${OUTPUT_TAG}_shard${SHARD_INDEX}"
fi

ARGS=(
  --estimates "${RUN_ROOT}/differential/${SOURCE_TAG}/subisoform_estimates.jsonl.gz"
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz"
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}"
  --min-celltype-mice 3
  --min-path-proportion "${MIN_PATH_PROPORTION}"
  --max-logratio-variance "${MAX_LOGRATIO_VARIANCE}"
  --permutations "${PERMUTATIONS}"
  --celltype-test "${CELLTYPE_TEST}"
  --likelihood "${LIKELIHOOD}"
  --shard-index "${SHARD_INDEX}"
  --shard-count "${SHARD_COUNT}"
)
if (( MAX_PATHS > 0 )); then
  ARGS+=(--max-paths "${MAX_PATHS}")
fi

python "${REPO_ROOT}/extra_scripts/run_celltype_compositional_splicing.py" \
  "${ARGS[@]}"
