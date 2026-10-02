#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --job-name=audit_hybrid_lr
#SBATCH --output=.cache/suppa2_statistics/hybrid-lr-audit-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-lr-audit-%j.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python extra_scripts/audit_suppa2_hybrid_lr.py \
  --hybrid-mapping analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/lr_model_audit/fixed_mass_reference/full_data_tilgner.tsv.gz \
  --hybrid-tests analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/lr_model_audit/fixed_mass_reference/full_data_tests.tsv.gz \
  --native-tests analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz \
  --native-catalog analyses/comparator_suppa_rmats/suppa2/event_catalog.tsv.gz \
  --hybrid-catalog analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/full_data/event_catalog.tsv.gz \
  --output-dir .cache/suppa2_statistics/hybrid_lr_debug/initial
