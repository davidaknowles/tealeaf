#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=hybrid_psi_audit
#SBATCH --output=.cache/suppa2_statistics/hybrid-psi-audit-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-psi-audit-%j.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
BASE=${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}
for TAG in fixed64_report1 profile64_report1; do
  python extra_scripts/retest_suppa2_hybrid_psi.py --shards "${BASE}/${TAG}/shards" --output-dir "${BASE}/${TAG}_psi/raw"
  python extra_scripts/merge_paired_path_test.py --shards "${BASE}/${TAG}_psi/raw" --output-dir "${BASE}/${TAG}_psi/merged" --moderate-variances --calibration empirical
done
python extra_scripts/evaluate_hybrid_lr_variants.py \
  --mapping analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/lr_model_audit/fixed_mass_reference/full_data_tilgner.tsv.gz \
  --native-tests analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz \
  --variant fixed_psi "${BASE}/fixed64_report1_psi/merged/paired_path.tsv" \
  --variant profile_psi "${BASE}/profile64_report1_psi/merged/paired_path.tsv" \
  --output-dir "${BASE}/evaluation_psi"
