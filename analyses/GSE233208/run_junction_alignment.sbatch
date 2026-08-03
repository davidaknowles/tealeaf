#!/bin/bash -l
#SBATCH -p cpu
#SBATCH --array=0-39
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --job-name=gse233208_sj

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"

mapfile -t records < <(awk -F '\t' 'NR > 1 && $7 == 1 {print $1 "\t" $2}' "${DATA_ROOT}/manifest/ena_fastqs.tsv")
record=${records[${SLURM_ARRAY_TASK_ID}]}
batch=${record%%$'\t'*}
run=${record#*$'\t'}
benchmark_root="${RUN_ROOT}/junction_benchmark"
output_dir="${benchmark_root}/starsolo/${run}"
whitelists="${benchmark_root}/whitelists"

export FASTQ_R1="${DATA_ROOT}/fastq/${batch}/${run}_1.fastq.gz"
export FASTQ_R2="${DATA_ROOT}/fastq/${batch}/${run}_2.fastq.gz"
export STAR_INDEX=${STAR_INDEX:-${HOME}/knowles_lab/index/STAR/tealeaf_hg38_gencode30_74}
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
