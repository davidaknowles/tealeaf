#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:b6k:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=4-00:00:00

set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${ALEVIN_DIR:?set ALEVIN_DIR}"
: "${SALMON_REF:?set SALMON_REF}"
: "${RUN_DIR:?set RUN_DIR}"
: "${DESIGN:?set DESIGN to binary or weighted}"
: "${METHOD:?set METHOD}"
module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

ANALYSIS_TAG=${ANALYSIS_TAG:-standard}
REPORT="${RUN_DIR}/cv_${ANALYSIS_TAG}_${DESIGN}_${METHOD}_theta.json"
VALUE=$(python -c '
import json, sys
report = json.load(open(sys.argv[1]))
if report.get("best_on_boundary"):
    raise SystemExit("CV optimum remains on an open boundary")
if not report.get("convergence_requirement_met", True):
    raise SystemExit("no converged CV candidate")
if not report.get("nondegeneracy_requirement_met", True):
    raise SystemExit("no nondegenerate CV candidate")
print(report["selected_rank"] if report["method"] == "factorized"
      else report["selected_full_value"])
' "${REPORT}")

MIN_CELL_UMIS=${MIN_CELL_UMIS:-500}
MIN_HALF_UMIS=${MIN_HALF_UMIS:-500}
BATCH_CELLS=${BATCH_CELLS:-16384}
EXACT_INNER_STEPS=${EXACT_INNER_STEPS:-32}
method_args=()
glm_rank=64
if [[ "${METHOD}" == factorized ]]; then
  glm_rank=${VALUE}
  method_args+=(--nucnorm_max_iter 8192)
  method_args+=(--glm_no_minibatch --glm_exact_inner_steps "${EXACT_INNER_STEPS}")
  full_tol=1e-5
elif [[ "${METHOD}" == admm_factorized ]]; then
  method_args+=(--nucnorm_lambda "${VALUE}" --nucnorm_max_iter 8192)
  full_tol=1e-6
elif [[ "${METHOD}" == frank_wolfe_penalized ]]; then
  method_args+=(--nucnorm_tau "${VALUE}" --nucnorm_max_iter 4096)
  method_args+=(--fw_max_atoms 4096 --fw_power_iter 10)
  full_tol=1e-4
else
  echo "Unsupported scalable method: ${METHOD}" >&2
  exit 1
fi

primer_args=()
if [[ -n "${PRIMER_PAIRS:-}" ]]; then
  primer_args+=(--primer_pairs "${PRIMER_PAIRS}" --min_half_umis "${MIN_HALF_UMIS}")
fi

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir "${ALEVIN_DIR}" \
  --salmon_ref "${SALMON_REF}" \
  --ref_dir "${RUN_DIR}" \
  --cell_mode single_cell \
  --quant_method "${METHOD}" \
  --regularization_target theta \
  --ec_design "${DESIGN}" \
  --glm_device cuda \
  --glm_rank "${glm_rank}" \
  --glm_batch_cells "${BATCH_CELLS}" \
  --glm_data_backend auto \
  --glm_no_write_chunks \
  --min_cell_umis "${MIN_CELL_UMIS}" \
  --nucnorm_min_iter 50 \
  --nucnorm_patience 10 \
  --nucnorm_tol "${full_tol}" \
  --admm_rho 1 \
  --admm_rho_update_interval 10 \
  --admm_rho_balance 10 \
  --admm_rho_scale 2 \
  --min_eq 5 \
  "${primer_args[@]}" \
  "${method_args[@]}" \
  -o "${RUN_DIR}/${ANALYSIS_TAG}_${DESIGN}_${METHOD}_theta_"
