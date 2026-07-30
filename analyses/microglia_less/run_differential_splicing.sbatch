#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=microglia_diff_splice

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

FIT_TAG=${FIT_TAG:-ablation_joint30_lr1e2_joint}
OUTPUT_TAG=${OUTPUT_TAG:-weighted_hierarchical_subisoform}
MIN_GENE_UMIS=${MIN_GENE_UMIS:-25}
MAX_PATHS=${MAX_PATHS:-8}
MAX_LOGRATIO_VARIANCE=${MAX_LOGRATIO_VARIANCE:-0.125}
MIN_PATH_PROPORTION=${MIN_PATH_PROPORTION:-0.15}
MIN_CELLTYPE_MICE=${MIN_CELLTYPE_MICE:-3}
PERMUTATIONS=${PERMUTATIONS:-20}
BOOTSTRAP_BLOCKS=${BOOTSTRAP_BLOCKS:-20}
BOOTSTRAP_REPS=${BOOTSTRAP_REPS:-100}
MAX_BLOCKS=${MAX_BLOCKS:-}

max_blocks_args=()
if [[ -n "${MAX_BLOCKS}" ]]; then
  max_blocks_args+=(--max-blocks "${MAX_BLOCKS}")
fi

python "${REPO_ROOT}/extra_scripts/run_differential_splicing.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --fit-prefix "${RUN_ROOT}/hierarchical/${FIT_TAG}_" \
  --gtf "${GENCODE_GTF}" \
  --barcode-groups "${RUN_ROOT}/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv" \
  --block-cache "${RUN_ROOT}/differential/gencode_vM32_splice_blocks.json.gz" \
  --output-dir "${RUN_ROOT}/differential/${OUTPUT_TAG}" \
  --min-half-umis 500 \
  --min-cells 20 \
  --min-pseudobulk-umis 100000 \
  --min-gene-umis "${MIN_GENE_UMIS}" \
  --max-paths "${MAX_PATHS}" \
  --max-logratio-variance "${MAX_LOGRATIO_VARIANCE}" \
  --min-path-proportion "${MIN_PATH_PROPORTION}" \
  --min-celltype-mice "${MIN_CELLTYPE_MICE}" \
  --permutations "${PERMUTATIONS}" \
  --bootstrap-blocks "${BOOTSTRAP_BLOCKS}" \
  --bootstrap-reps "${BOOTSTRAP_REPS}" \
  "${max_blocks_args[@]}"
