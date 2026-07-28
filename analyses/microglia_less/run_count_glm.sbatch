#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=microglia_count_glm

set -euo pipefail

CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

OUTPUT_TAG=${OUTPUT_TAG:-weighted_nb10_rank32}
RANK=${RANK:-32}
HVG=${HVG:-2000}
MIN_FEATURE_CELLS=${MIN_FEATURE_CELLS:-50}
CONCENTRATION=${CONCENTRATION:-10}
GENE_FAMILY=${GENE_FAMILY:-negative_binomial}
BATCH_CELLS=${BATCH_CELLS:-256}
GENE_EPOCHS=${GENE_EPOCHS:-10}
EC_EPOCHS=${EC_EPOCHS:-10}
JOINT_EPOCHS=${JOINT_EPOCHS:-10}
MAX_CELLS=${MAX_CELLS:-}

max_cells_args=()
if [[ -n "${MAX_CELLS}" ]]; then
  max_cells_args+=(--max-cells "${MAX_CELLS}")
fi

mkdir -p "${RUN_ROOT}/count_glm"
python "${REPO_ROOT}/extra_scripts/fit_count_glm.py" \
  --alevin-dir "${ALEVIN_DIR}" \
  --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --transcript-to-gene "${SALMON_REF}/spliceu_t2g.tsv" \
  --output-prefix "${RUN_ROOT}/count_glm/${OUTPUT_TAG}_" \
  --ec-design weighted \
  --primer-sampling-model oligodt_tpm \
  --rank "${RANK}" \
  --hvg "${HVG}" \
  --min-feature-cells "${MIN_FEATURE_CELLS}" \
  --gene-family "${GENE_FAMILY}" \
  --concentration "${CONCENTRATION}" \
  --batch-cells "${BATCH_CELLS}" \
  --gene-epochs "${GENE_EPOCHS}" \
  --ec-epochs "${EC_EPOCHS}" \
  --joint-epochs "${JOINT_EPOCHS}" \
  --device cuda \
  "${max_cells_args[@]}"
