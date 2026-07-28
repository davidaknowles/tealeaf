#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=microglia_hierarchical

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

OUTPUT_TAG=${OUTPUT_TAG:-weighted_rank32}
RANK=${RANK:-32}
BATCH_CELLS=${BATCH_CELLS:-256}
GENE_EPOCHS=${GENE_EPOCHS:-10}
ISOFORM_EPOCHS=${ISOFORM_EPOCHS:-10}
JOINT_EPOCHS=${JOINT_EPOCHS:-10}
MAX_CELLS=${MAX_CELLS:-}

max_cells_args=()
if [[ -n "${MAX_CELLS}" ]]; then
  max_cells_args+=(--max-cells "${MAX_CELLS}")
fi

mkdir -p "${RUN_ROOT}/hierarchical"
python "${REPO_ROOT}/extra_scripts/fit_hierarchical.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --output-prefix "${RUN_ROOT}/hierarchical/${OUTPUT_TAG}_" \
  --ec-design weighted \
  --primer-sampling-model oligodt_tpm \
  --min-eq 5 \
  --min-half-umis 500 \
  --min-cell-umis 500 \
  --rank "${RANK}" \
  --batch-cells "${BATCH_CELLS}" \
  --gene-epochs "${GENE_EPOCHS}" \
  --isoform-epochs "${ISOFORM_EPOCHS}" \
  --joint-epochs "${JOINT_EPOCHS}" \
  --write-pca-initialization \
  --device cuda \
  "${max_cells_args[@]}"
