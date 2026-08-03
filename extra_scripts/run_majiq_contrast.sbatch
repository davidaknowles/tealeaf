#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${CONTRASTS:?set CONTRASTS}"
: "${PSICOV:?set PSICOV}"
: "${SPLICEGRAPH:?set SPLICEGRAPH}"
: "${ANNOTATED_SPLICEGRAPH:?set ANNOTATED_SPLICEGRAPH}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

module purge
module load Python/3.12.3-GCCcore-13.3.0 HTSlib/1.21 zstd/1.5.6-GCCcore-13.3.0
source "${MAJIQ_VENV:-${HOME}/venv/majiq_cluster}/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
mkdir -p "${OUTPUT_DIR}/raw" "${OUTPUT_DIR}/normalized"
mapfile -t fields < <(python - "${CONTRASTS}" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json, sys
x = json.load(open(sys.argv[1]))[int(sys.argv[2])]
print(x["contrast_id"])
print(x["level_a"])
print(x["level_b"])
print(" ".join(x["samples_a"]))
print(" ".join(x["samples_b"]))
PY
)
contrast=${fields[0]}
read -ra group_a <<< "${fields[3]}"
read -ra group_b <<< "${fields[4]}"
raw="${OUTPUT_DIR}/raw/${contrast}.tsv"
majiq-v3 heterogen -psi1 "${PSICOV}" -psi2 "${PSICOV}" \
  --select-grp1-prefixes "${group_a[@]}" \
  --select-grp2-prefixes "${group_b[@]}" \
  -n "${fields[1]}" "${fields[2]}" --min-experiments 0.5 \
  --stats mannwhitneyu --psisamples 100 --pvalue-quantiles 0.95 \
  --splicegraph "${SPLICEGRAPH}" --annotated "${ANNOTATED_SPLICEGRAPH}" \
  --output-tsv "${raw}" --nthreads "${SLURM_CPUS_PER_TASK:-4}" \
  --license "${MAJIQ_LICENSE_FILE}"
python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" normalize-majiq \
  --contrasts "${CONTRASTS}" --contrast-index "${SLURM_ARRAY_TASK_ID}" \
  --input "${raw}" --output "${OUTPUT_DIR}/normalized/${contrast}.tsv"
