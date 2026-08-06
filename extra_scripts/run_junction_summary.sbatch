#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${BENCHMARK_ROOT:?set BENCHMARK_ROOT}"

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
output="${BENCHMARK_ROOT}/comparison"
mkdir -p "${output}"
"${HOME}/venv/cpa_blackwell/bin/python" \
  "${REPO_ROOT}/extra_scripts/run_junction_comparator.py" summarize \
  --input-dir "${BENCHMARK_ROOT}/leafcutter/tests/normalized" \
  --input-dir "${BENCHMARK_ROOT}/scquint/tests" \
  --input-dir "${BENCHMARK_ROOT}/majiq/tests/normalized" \
  --output "${output}/all_tests.tsv.gz" \
  --omnibus-output "${output}/cell_type_simes.tsv.gz" \
  --summary-output "${output}/summary.tsv"
