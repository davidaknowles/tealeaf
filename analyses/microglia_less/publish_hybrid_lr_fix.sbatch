#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=hybrid_publish_fix
#SBATCH --output=.cache/suppa2_statistics/hybrid-publish-fix-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-publish-fix-%j.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python extra_scripts/publish_suppa2_hybrid_lr_fix.py \
  --debug-dir "${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}" \
  --output-dir analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid
bash analyses/microglia_less/evaluate_suppa2_hybrid_tilgner.sbatch
