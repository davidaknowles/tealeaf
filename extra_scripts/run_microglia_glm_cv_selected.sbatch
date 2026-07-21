#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_glm_cv_fit
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

DESIGN=${DESIGN:?Set DESIGN to binary or weighted}
METHOD=${METHOD:?Set METHOD to factorized, admm_factorized, or frank_wolfe_penalized}
REPORT_TAG=${REPORT_TAG:-_minumi500_all_minibatch}
BATCH_CELLS=${BATCH_CELLS:-16384}
DATA_BACKEND=${DATA_BACKEND:-auto}
MIN_CELL_UMIS=${MIN_CELL_UMIS:-500}
RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
REPORT=${RUN}/cv_${DESIGN}_${METHOD}_theta${REPORT_TAG}.json
VALUE=$(python -c '
import json, sys
report = json.load(open(sys.argv[1]))
if report.get("best_on_boundary"):
    raise SystemExit("CV optimum remains on an open boundary after automatic grid expansion")
if not report.get("convergence_requirement_met", True):
    raise SystemExit("No CV candidate satisfied the convergence requirement")
if not report.get("nondegeneracy_requirement_met", True):
    raise SystemExit("No CV candidate satisfied the nondegeneracy requirement")
if report["method"] == "factorized":
    print(report["selected_rank"])
else:
    print(report["selected_full_value"])
' "${REPORT}")

method_args=()
glm_rank=64
if [[ "${METHOD}" == admm_factorized ]]; then
  method_args+=(--nucnorm_lambda "${VALUE}" --nucnorm_max_iter 4096)
  full_tol=1e-6
elif [[ "${METHOD}" == factorized ]]; then
  glm_rank=${VALUE}
  method_args+=(--nucnorm_max_iter 2048)
  full_tol=1e-5
else
  method_args+=(--nucnorm_tau "${VALUE}" --nucnorm_max_iter 2048)
  method_args+=(--fw_max_atoms 2048 --fw_power_iter 10 --fw_nonnegative_penalty 1.0)
  full_tol=1e-4
fi

python -m tealeaf.sc.tealeaf_sc \
  --alevin_dir /gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin \
  --salmon_ref /gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa \
  --ref_dir /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/ref \
  --cell_mode single_cell \
  --quant_method "${METHOD}" \
  --regularization_target theta \
  --ec_design "${DESIGN}" \
  --glm_device cuda \
  --glm_rank "${glm_rank}" \
  --glm_batch_cells "${BATCH_CELLS}" \
  --glm_data_backend "${DATA_BACKEND}" \
  --glm_polish_max_iter 32 \
  --glm_no_write_chunks \
  --min_cell_umis "${MIN_CELL_UMIS}" \
  --nucnorm_min_iter 50 \
  --nucnorm_patience 10 \
  --nucnorm_tol "${full_tol}" \
  --admm_rho 1.0 \
  --admm_rho_update_interval 10 \
  --admm_rho_balance 10 \
  --admm_rho_scale 2 \
  --min_eq 5 \
  "${method_args[@]}" \
  -o "${RUN}/single_cell_${DESIGN}_${METHOD}_cv_minibatch_theta_"
