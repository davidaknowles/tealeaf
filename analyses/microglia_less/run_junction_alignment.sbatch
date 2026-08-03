#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-15
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=20:00:00
#SBATCH --job-name=microglia_sj

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/microglia_less/config.env}
source "${CONFIG}"

mapfile -t reads < <(find "${DATA_ROOT}/SRA" -maxdepth 1 -name '*_1.fastq.gz' | sort)
read1=${reads[${SLURM_ARRAY_TASK_ID}]}
run=$(basename "${read1}" _1.fastq.gz)
benchmark_root="${RUN_ROOT}/junction_benchmark"
output_dir="${benchmark_root}/starsolo/${run}"
whitelists="${benchmark_root}/whitelists"

export FASTQ_R1="${read1}"
export FASTQ_R2="${read1%_1.fastq.gz}_2.fastq.gz"
export STAR_INDEX=${STAR_INDEX:-${HOME}/knowles_lab/index/STAR/tealeaf_GRCm39_gencodevM32}
export OUTPUT_DIR="${output_dir}"
export WHITELIST_1="${whitelists}/barcode_segment_1.txt"
export WHITELIST_2="${whitelists}/barcode_segment_2.txt"
export WHITELIST_3="${whitelists}/barcode_segment_3.txt"
export BAM_SORT_RAM=${BAM_SORT_RAM:-80000000000}
bash "${REPO_ROOT}/extra_scripts/run_splitseq_starsolo.sbatch"

module purge
module load "${PYTHON_MODULE}"
source "${HOME}/venv/cpa_blackwell/bin/activate"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
python "${REPO_ROOT}/extra_scripts/tag_spliced_pseudobulk_bam.py" \
  --input "${output_dir}/Aligned.sortedByCoord.out.bam" \
  --cell-metadata "${benchmark_root}/cell_metadata.tsv" \
  --run "${run}" \
  --output "${output_dir}/spliced_pseudobulk.bam"
