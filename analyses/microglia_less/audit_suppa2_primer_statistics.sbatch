#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=suppa_primer_stats
#SBATCH --output=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/primer_audit_%j.out
#SBATCH --error=/gpfs/commons/home/daknowles/projects/tealeaf/.cache/suppa2_statistics/primer_audit_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BASE="${REPO_ROOT}/.cache/suppa2_statistics/primer_split_audit_20261002"
REPRO="${RUN_ROOT}/junction_benchmark/reproducibility"
python "${REPO_ROOT}/extra_scripts/run_suppa2_primer_aware.py" \
  --alevin-dir "${ALEVIN_DIR}" --salmon-ref "${SALMON_REF}/spliceu.fa" \
  --primer-pairs "${RUN_ROOT}/primer_pairs.tsv" \
  --barcode-groups "${RUN_ROOT}/barcodes_to_cluster_dx_mouse.filtered_min20_umi1e5.csv" \
  --event-catalog "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2/event_catalog.tsv.gz" \
  --contrasts "${REPRO}/fold0/contrasts.json" --contrasts "${REPRO}/fold1/contrasts.json" \
  --output "${BASE}/unused_tests.tsv.gz" --event-output "${BASE}/event_catalog.tsv.gz" \
  --work-dir "${BASE}/work" --psi-cache-dir "${BASE}/psi" --prepare-only
for fold in 0 1; do
  python "${REPO_ROOT}/extra_scripts/audit_suppa2_split_statistics.py" --run-root "${RUN_ROOT}" --repo-root "${REPO_ROOT}" --fold "${fold}" --output-dir "${BASE}" --primer-psi-dir "${BASE}/psi"
done
python "${REPO_ROOT}/extra_scripts/summarize_suppa2_split_statistics.py" --run-root "${RUN_ROOT}" --repo-root "${REPO_ROOT}" --audit-dir "${BASE}" --output-dir "${REPO_ROOT}/analyses/comparator_suppa_rmats/suppa2_primer_aware/statistics_audit" --primer-aware
