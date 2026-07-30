#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=microglia_comp_splice

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

SOURCE_TAG=${SOURCE_TAG:-weighted_hierarchical_subisoform_min25}
OUTPUT_TAG=${OUTPUT_TAG:-compositional_subisoform_min25}
PERMUTATIONS=${PERMUTATIONS:-20}

python "${REPO_ROOT}/extra_scripts/run_compositional_splicing.py" \
  --estimates "${RUN_ROOT}/differential/${SOURCE_TAG}/subisoform_estimates.jsonl.gz" \
  --gls-tests "${RUN_ROOT}/differential/${SOURCE_TAG}/differential_splicing.tsv" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}" \
  --min-condition-replicates 3 \
  --permutations "${PERMUTATIONS}"
