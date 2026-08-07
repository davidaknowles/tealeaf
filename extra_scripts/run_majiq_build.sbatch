#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00

set -euo pipefail
: "${BAM_MANIFEST:?set BAM_MANIFEST}"
: "${ANNOTATED_SPLICEGRAPH:?set ANNOTATED_SPLICEGRAPH}"
: "${OUTPUT_DIR:?set OUTPUT_DIR}"
: "${MAJIQ_LICENSE_FILE:?set MAJIQ_LICENSE_FILE}"

module purge
module load Python/3.12.3-GCCcore-13.3.0 HTSlib/1.21 zstd/1.5.6-GCCcore-13.3.0
source "${MAJIQ_VENV:-${HOME}/venv/majiq_cluster}/bin/activate"
mkdir -p "${OUTPUT_DIR}"
sj_dir=${MAJIQ_SJ_DIR:-${OUTPUT_DIR}/sj}
min_experiments=${MAJIQ_BUILD_MIN_EXPERIMENTS:-0.5}
quantify_minreads=${MAJIQ_QUANTIFY_MINREADS:-10}
quantify_minbins=${MAJIQ_QUANTIFY_MINBINS:-3}
sjs=()
prefixes=()
while IFS=$'\t' read -r sample bam; do
  sjs+=("${sj_dir}/${sample}.sj")
  prefixes+=("${sample}")
done < "${BAM_MANIFEST}"
majiq-v3 build "${ANNOTATED_SPLICEGRAPH}" "${OUTPUT_DIR}/combined.splicegraph" \
  --sjs "${sjs[@]}" --no-introns \
  --min-experiments "${min_experiments}" \
  --nthreads "${SLURM_CPUS_PER_TASK:-8}" \
  --license "${MAJIQ_LICENSE_FILE}" --overwrite
majiq-v3 psi-coverage "${OUTPUT_DIR}/combined.splicegraph" "${OUTPUT_DIR}/all.psicov" \
  "${sjs[@]}" --prefixes "${prefixes[@]}" \
  --minreads "${quantify_minreads}" --minbins "${quantify_minbins}" \
  --nthreads "${SLURM_CPUS_PER_TASK:-8}" \
  --license "${MAJIQ_LICENSE_FILE}" --overwrite
