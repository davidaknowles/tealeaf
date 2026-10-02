#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=merge_hybrid_ablation
#SBATCH --output=.cache/suppa2_statistics/hybrid-lr-merge-%j.out
#SBATCH --error=.cache/suppa2_statistics/hybrid-lr-merge-%j.err
set -euo pipefail
source analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
BASE=${DEBUG_ROOT:-.cache/suppa2_statistics/hybrid_lr_analytic}
TEST_ZETA=${TEST_ZETA:-64}
if [[ "${REPORT_ZETA:-1}" == none ]]; then
  python extra_scripts/merge_paired_path_test.py --shards "${BASE}/fixed${TEST_ZETA}_testonly/shards"/shard_* --output-dir "${BASE}/fixed${TEST_ZETA}_testonly/merged" --moderate-variances --calibration empirical
  exit 0
fi
for TAG in fixed${TEST_ZETA}_report1 profile${TEST_ZETA}_report1; do
  python extra_scripts/merge_paired_path_test.py --shards "${BASE}/${TAG}/shards"/shard_* --output-dir "${BASE}/${TAG}/merged" --moderate-variances --calibration empirical
done
python extra_scripts/evaluate_hybrid_lr_variants.py \
  --mapping analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/lr_model_audit/fixed_mass_reference/full_data_tilgner.tsv.gz \
  --native-tests analyses/comparator_suppa_rmats/suppa2/full_data_tests.tsv.gz \
  --variant original analyses/comparator_suppa_rmats/suppa2_tealeaf_hybrid/lr_model_audit/fixed_mass_reference/full_data_tests.tsv.gz \
  --variant fixed${TEST_ZETA}_report1 "${BASE}/fixed${TEST_ZETA}_report1/merged/paired_path.tsv" \
  --variant profile${TEST_ZETA}_report1 "${BASE}/profile${TEST_ZETA}_report1/merged/paired_path.tsv" \
  --output-dir "${BASE}/evaluation${TEST_ZETA}"
