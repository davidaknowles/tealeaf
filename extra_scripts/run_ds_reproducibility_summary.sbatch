#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail
: "${REPO_ROOT:?set REPO_ROOT}"
: "${BENCHMARK_ROOT:?set BENCHMARK_ROOT}"

unset LOADEDMODULES
for variable in $(compgen -A variable _ModuleTable); do
  unset "${variable}"
done
module load "${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python_bin="${PYTHON_BIN:-${HOME}/venv/cpa_blackwell/bin/python}"
pair_args=()
if [[ "${MATCH_PAIRS:-0}" == "1" ]]; then
  pair_args+=(--match-pairs)
fi

"${python_bin}" \
  "${REPO_ROOT}/extra_scripts/score_ds_reproducibility.py" \
  --fold-dir "${BENCHMARK_ROOT}/reproducibility/fold0" \
  --fold-dir "${BENCHMARK_ROOT}/reproducibility/fold1" \
  --junction-table "${BENCHMARK_ROOT}/pseudobulk_junctions_junctions.tsv.gz" \
  --leafcutter-counts \
    "${BENCHMARK_ROOT}/leafcutter/clustering/benchmark_perind_numers.counts.gz" \
  --leafcutter-map \
    "${BENCHMARK_ROOT}/reproducibility/leafcutter_cluster_gene.tsv.gz" \
  --tealeaf-subdir "${TEALEAF_SUBDIR:-tealeaf}" \
  --tealeaf-label "${TEALEAF_LABEL:-Tealeaf EC GLMM}" \
  --tealeaf-fit-method \
    "${TEALEAF_FIT_METHOD:-laplace_multinomial}" \
  --min-median-gene-umis "${MIN_MEDIAN_GENE_UMIS:-0}" \
  --min-paired-subjects "${MIN_PAIRED_SUBJECTS:-0}" \
  "${pair_args[@]}" \
  --output-dir \
    "${BENCHMARK_ROOT}/reproducibility/${COMPARISON_SUBDIR:-comparison}"
