#!/bin/bash
set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${BUNDLE:?set BUNDLE}"
: "${BENCHMARK_ROOT:?set BENCHMARK_ROOT}"
: "${PYTHON_MODULE:=Python/3.12.3-GCCcore-13.3.0}"
: "${LEAFCUTTER_ROOT:=${HOME}/projects/leafcutter}"
: "${SCQUINT_ROOT:=${HOME}/projects/scquint}"

mkdir -p "${BENCHMARK_ROOT}/logs"
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
contrasts="${BENCHMARK_ROOT}/contrasts.json"
python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" plan \
  --bundle "${BUNDLE}" --min-subjects 4 --output "${contrasts}"
n=$(python -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' "${contrasts}")
if [[ "${n}" -eq 0 ]]; then
  echo "No replicate-supported contrasts" >&2
  exit 1
fi

leaf_job=$(sbatch --parsable \
  --job-name=leafcutter-cluster \
  --output="${BENCHMARK_ROOT}/logs/leafcutter-cluster.%j.out" \
  --export="ALL,JUNCTION_LIST=${BENCHMARK_ROOT}/leafcutter/junction_files.txt,OUTPUT_DIR=${BENCHMARK_ROOT}/leafcutter/clustering,LEAFCUTTER_ROOT=${LEAFCUTTER_ROOT},PYTHON_MODULE=${PYTHON_MODULE}" \
  "${REPO_ROOT}/extra_scripts/run_leafcutter_cluster.sbatch")
leaf_array=$(sbatch --parsable --dependency="afterok:${leaf_job}" --array="0-$((n - 1))%32" \
  --job-name=leafcutter-ds \
  --output="${BENCHMARK_ROOT}/logs/leafcutter-ds.%A_%a.out" \
  --export="ALL,REPO_ROOT=${REPO_ROOT},BUNDLE=${BUNDLE},CONTRASTS=${contrasts},COUNTS=${BENCHMARK_ROOT}/leafcutter/clustering/benchmark_perind_numers.counts.gz,OUTPUT_DIR=${BENCHMARK_ROOT}/leafcutter/tests,LEAFCUTTER_ROOT=${LEAFCUTTER_ROOT},PYTHON_MODULE=${PYTHON_MODULE}" \
  "${REPO_ROOT}/extra_scripts/run_leafcutter_contrast.sbatch")
scquint_array=$(sbatch --parsable --array="0-$((n - 1))%32" \
  --job-name=scquint-ds \
  --output="${BENCHMARK_ROOT}/logs/scquint-ds.%A_%a.out" \
  --export="ALL,REPO_ROOT=${REPO_ROOT},BUNDLE=${BUNDLE},CONTRASTS=${contrasts},OUTPUT_DIR=${BENCHMARK_ROOT}/scquint/tests,SCQUINT_ROOT=${SCQUINT_ROOT},PYTHON_MODULE=${PYTHON_MODULE}" \
  "${REPO_ROOT}/extra_scripts/run_scquint_contrast.sbatch")
printf 'leaf_cluster=%s\nleaf_tests=%s\nscquint_tests=%s\n' \
  "${leaf_job}" "${leaf_array}" "${scquint_array}" | tee "${BENCHMARK_ROOT}/comparator_jobs.txt"
