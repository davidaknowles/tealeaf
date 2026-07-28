#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

OUTPUT_TAG=${OUTPUT_TAG:-hierarchical_rank32}
RANK=${RANK:-32}
BATCH_CELLS=${BATCH_CELLS:-256}
GENE_EPOCHS=${GENE_EPOCHS:-10}
ISOFORM_EPOCHS=${ISOFORM_EPOCHS:-10}
JOINT_EPOCHS=${JOINT_EPOCHS:-10}
JOINT_GENE_WEIGHT=${JOINT_GENE_WEIGHT:-1}
MAX_CELLS=${MAX_CELLS:-}

max_cells_args=()
if [[ -n "${MAX_CELLS}" ]]; then
  max_cells_args+=(--max-cells "${MAX_CELLS}")
fi

python "${REPO_ROOT}/extra_scripts/fit_hierarchical.py" \
  --alevin-dir "${DATA_ROOT}/processed/merged_alevin" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${DATA_ROOT}/processed/merged_alevin/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --output-prefix "${RUN_ROOT}/hierarchical/${OUTPUT_TAG}_" \
  --ec-design positional \
  --primer-sampling-model oligodt_tpm \
  --min-eq 5 \
  --min-half-umis 500 \
  --min-cell-umis 500 \
  --rank "${RANK}" \
  --batch-cells "${BATCH_CELLS}" \
  --gene-epochs "${GENE_EPOCHS}" \
  --isoform-epochs "${ISOFORM_EPOCHS}" \
  --joint-epochs "${JOINT_EPOCHS}" \
  --joint-gene-weight "${JOINT_GENE_WEIGHT}" \
  --write-pca-initialization \
  --device cuda \
  "${max_cells_args[@]}"
