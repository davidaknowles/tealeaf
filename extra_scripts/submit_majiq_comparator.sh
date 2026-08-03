#!/bin/bash
set -euo pipefail

: "${REPO_ROOT:?set REPO_ROOT}"
: "${BENCHMARK_ROOT:?set BENCHMARK_ROOT}"
: "${GTF:?set GTF}"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

mkdir -p "${BENCHMARK_ROOT}/majiq/logs"
if [[ -z "${BUNDLE:-}" ]]; then
  manifest_json=$(find "${BENCHMARK_ROOT}" -maxdepth 2 -name '*_manifest.json' -print -quit)
  if [[ -z "${manifest_json}" ]]; then
    echo "No junction bundle manifest found under ${BENCHMARK_ROOT}" >&2
    exit 1
  fi
  BUNDLE=${manifest_json%_manifest.json}
fi
manifest="${BENCHMARK_ROOT}/majiq/bams.tsv"
search_root=${BAM_DIR:-${BENCHMARK_ROOT}}
find "${search_root}" -type f -name '*.bam' ! -name '*.unsorted.bam' -printf '%f\t%p\n' | \
  sed 's/\.bam\t/\t/' | sort > "${manifest}"
n_samples=$(wc -l < "${manifest}")
if [[ "${n_samples}" -lt 8 ]]; then
  echo "MAJIQ requires at least eight pseudobulk BAMs; found ${n_samples}" >&2
  exit 1
fi
contrasts="${BENCHMARK_ROOT}/contrasts.json"
if [[ ! -s "${contrasts}" ]]; then
  module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
  source "${HOME}/venv/cpa_blackwell/bin/activate"
  export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
  python "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" plan \
    --bundle "${BUNDLE}" --min-subjects 4 --output "${contrasts}"
fi
n_contrasts=$(python -c 'import json,sys; print(len(json.load(open(sys.argv[1]))))' "${contrasts}")
common="ALL,MAJIQ_LICENSE_FILE=${MAJIQ_LICENSE_FILE}"
annotation_job=$(sbatch --parsable --job-name=majiq-annotation \
  --output="${BENCHMARK_ROOT}/majiq/logs/annotation.%j.out" \
  --export="${common},GTF=${GTF},OUTPUT_DIR=${BENCHMARK_ROOT}/majiq" \
  "${REPO_ROOT}/extra_scripts/run_majiq_annotation.sbatch")
sj_job=$(sbatch --parsable --dependency="afterok:${annotation_job}" \
  --array="0-$((n_samples - 1))%32" --job-name=majiq-sj \
  --output="${BENCHMARK_ROOT}/majiq/logs/sj.%A_%a.out" \
  --export="${common},BAM_MANIFEST=${manifest},ANNOTATED_SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/annotated.splicegraph,OUTPUT_DIR=${BENCHMARK_ROOT}/majiq/sj" \
  "${REPO_ROOT}/extra_scripts/run_majiq_sj.sbatch")
build_job=$(sbatch --parsable --dependency="afterok:${sj_job}" --job-name=majiq-build \
  --output="${BENCHMARK_ROOT}/majiq/logs/build.%j.out" \
  --export="${common},BAM_MANIFEST=${manifest},ANNOTATED_SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/annotated.splicegraph,OUTPUT_DIR=${BENCHMARK_ROOT}/majiq" \
  "${REPO_ROOT}/extra_scripts/run_majiq_build.sbatch")
test_job=$(sbatch --parsable --dependency="afterok:${build_job}" \
  --array="0-$((n_contrasts - 1))%16" --job-name=majiq-ds \
  --output="${BENCHMARK_ROOT}/majiq/logs/test.%A_%a.out" \
  --export="${common},REPO_ROOT=${REPO_ROOT},CONTRASTS=${contrasts},PSICOV=${BENCHMARK_ROOT}/majiq/all.psicov,SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/combined.splicegraph,ANNOTATED_SPLICEGRAPH=${BENCHMARK_ROOT}/majiq/annotated.splicegraph,OUTPUT_DIR=${BENCHMARK_ROOT}/majiq/tests" \
  "${REPO_ROOT}/extra_scripts/run_majiq_contrast.sbatch")
comparator_jobs="${BENCHMARK_ROOT}/comparator_jobs.txt"
for attempt in $(seq 1 60); do
  [[ -s "${comparator_jobs}" ]] && break
  sleep 5
done
summary_dependencies="${test_job}"
if [[ -s "${comparator_jobs}" ]]; then
  leaf_tests=$(awk -F= '$1 == "leaf_tests" {print $2}' "${comparator_jobs}")
  scquint_tests=$(awk -F= '$1 == "scquint_tests" {print $2}' "${comparator_jobs}")
  summary_dependencies+="${leaf_tests:+:${leaf_tests}}${scquint_tests:+:${scquint_tests}}"
fi
summary_job=$(sbatch --parsable --dependency="afterok:${summary_dependencies}" \
  --job-name=junction-summary \
  --output="${BENCHMARK_ROOT}/majiq/logs/summary.%j.out" \
  --export="ALL,REPO_ROOT=${REPO_ROOT},BENCHMARK_ROOT=${BENCHMARK_ROOT}" \
  "${REPO_ROOT}/extra_scripts/run_junction_summary.sbatch")
printf 'annotation=%s\nsj=%s\nbuild=%s\ntests=%s\nsummary=%s\n' \
  "${annotation_job}" "${sj_job}" "${build_job}" "${test_job}" "${summary_job}" | \
  tee "${BENCHMARK_ROOT}/majiq/jobs.txt"
