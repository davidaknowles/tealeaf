#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --job-name=lr_all_tested
#SBATCH --output=.cache/lr_ranking_20261006/all_tested_%j.out
#SBATCH --error=.cache/lr_ranking_20261006/all_tested_%j.err
set -euo pipefail
source /gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env
module load "${PYTHON_MODULE}"
source "${HOME}/venv/jax/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
IFS=';' read -r -a methods <<< "${AUDIT_METHODS:-LeafCutter;scQuint}"
extra=()
if [[ " ${methods[*]} " == *" MAJIQ Heterogen "* ]]; then
  extra+=(--majiq-tests "${RUN_ROOT}/junction_benchmark/majiq_build_min3_cov3/tilgner_tests/normalized" --majiq-build-min-experiments 3 --majiq-psi-minimum-reads 3 --majiq-psi-minimum-bins 3 --majiq-heterogen-minimum-experiment-fraction .5)
fi
python "${REPO_ROOT}/extra_scripts/assess_tilgner_junction_replication.py" \
  --comparator-tests "${RUN_ROOT}/junction_benchmark/comparison_corrected/all_tests.tsv.gz" \
  --junction-bundle "${RUN_ROOT}/junction_benchmark/pseudobulk_junctions" \
  --leafcutter-counts "${RUN_ROOT}/junction_benchmark/leafcutter/clustering/benchmark_perind_numers.counts.gz" \
  --majiq-raw "${RUN_ROOT}/junction_benchmark/majiq_build_min3_cov3/tilgner_tests/raw" \
  --tilgner-matrix "${REPO_ROOT}/.cache/tilgner_scisorseq/LR_annotTranscripts_byCT" \
  --gtf "${GENCODE_GTF}" \
  --output-dir "${AUDIT_DIR:-${REPO_ROOT}/analyses/tilgner_junction_replication/all_tested}" \
  --all-tested --methods "${methods[@]}" --skip-paired-clr "${extra[@]}"
