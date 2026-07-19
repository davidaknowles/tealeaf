#!/bin/bash -l
#SBATCH -p gpu
#SBATCH --gres gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH -t 2-00:00:00
#SBATCH -J tealeaf_glm_cv
#SBATCH -o /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.out
#SBATCH -e /gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run/%x-%j.err

set -euo pipefail

module load Python/3.12.3-GCCcore-13.3.0
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH=/gpfs/commons/home/daknowles/projects/tealeaf:${PYTHONPATH:-}

DESIGN=${DESIGN:?Set DESIGN to binary or weighted}
METHOD=${METHOD:?Set METHOD to admm_factorized or frank_wolfe_penalized}
RUN=/gpfs/commons/home/daknowles/tealeaf_runs/microglia_less/run
ALEVIN=/gpfs/commons/groups/knowles_lab/data/sc/splitpool/microglia_less_mice/salmon_spliceu_weighted_rad/out_permit_known_weighted/quant_t2t_dedup_parsimony_em/alevin
REF=/gpfs/commons/home/daknowles/knowles_lab/index/salmon/mus_spliceu/spliceu.fa

if [[ "${METHOD}" == admm_factorized ]]; then
  multipliers=(0 1e-9 1e-8 1e-7 1e-6 1e-5)
  max_iter=512
  expansion_factor=10
  selection_rule=minimum
  cv_tol=1e-5
else
  multipliers=(0 0.25 1 4 16 64 256)
  max_iter=1024
  expansion_factor=4
  selection_rule=one_standard_error
  cv_tol=1e-4
fi
multiplier_args=()
for value in "${multipliers[@]}"; do
  multiplier_args+=(--multiplier "${value}")
done

python /gpfs/commons/home/daknowles/projects/tealeaf/extra_scripts/tune_single_cell_glm.py \
  --alevin-dir "${ALEVIN}" \
  --salmon-ref "${REF}" \
  --design "${DESIGN}" \
  --method "${METHOD}" \
  --output "${RUN}/cv_${DESIGN}_${METHOD}_theta.json" \
  --cells 5000 \
  --folds 3 \
  --max-iter "${max_iter}" \
  --rank 64 \
  --power-iter 5 \
  --scale-power-iter 10 \
  --batch-cells 2048 \
  --device cuda \
  --max-grid-expansions 4 \
  --grid-expansion-factor "${expansion_factor}" \
  --selection-rule "${selection_rule}" \
  --require-converged \
  --tol "${cv_tol}" \
  "${multiplier_args[@]}"
