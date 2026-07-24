#!/bin/bash
#SBATCH --job-name=gse233208_posbias
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=16:00:00

set -euo pipefail
CONFIG=${CONFIG:-/gpfs/commons/home/daknowles/projects/tealeaf/analyses/GSE233208/config.env}
source "${CONFIG}"
module load "${PYTHON_MODULE}"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
PIGZ_BIN=${PIGZ_BIN:-/nfs/sw/easybuild/software/pigz/2.7-GCCcore-12.2.0/bin/pigz}

MANIFEST=${MANIFEST:-"${DATA_ROOT}/manifest/ena_fastqs.tsv"}
ALEVIN_DIR=${ALEVIN_DIR:-"${DATA_ROOT}/processed/merged_alevin"}
MAX_READS=${MAX_READS:-}
SALMON_THREADS=${SALMON_THREADS:-6}
PIGZ_THREADS=${PIGZ_THREADS:-2}
BALANCED_PREFIX_READS=${BALANCED_PREFIX_READS:-0}
LIBTYPE=${LIBTYPE:-U}

read1=()
read2=()
while IFS=$'\t' read -r batch run _sample _title _experiment _experiment_title \
    read _url _bytes _md5 _filename; do
  [[ "${batch}" == batch ]] && continue
  path="${DATA_ROOT}/fastq/${batch}/${run}_${read}.fastq.gz"
  if [[ "${read}" == 1 ]]; then
    read1+=("${path}")
  elif [[ "${read}" == 2 ]]; then
    read2+=("${path}")
  fi
done < "${MANIFEST}"
if (( ${#read1[@]} == 0 || ${#read1[@]} != ${#read2[@]} )); then
  echo "Manifest does not contain matching Parse read pairs" >&2
  exit 1
fi
run_index=${RUN_INDEX:-${SLURM_ARRAY_TASK_ID:-}}
if [[ -n "${run_index}" ]]; then
  if (( run_index < 0 || run_index >= ${#read1[@]} )); then
    echo "RUN_INDEX is outside the manifest" >&2
    exit 1
  fi
  run_accession=$(basename "${read1[run_index]}" _1.fastq.gz)
  read1=("${read1[run_index]}")
  read2=("${read2[run_index]}")
  OUTPUT_ROOT=${OUTPUT_ROOT:-"${DATA_ROOT}/processed/primer_positional_bias/by_run/${run_accession}"}
  BUILD_DESIGNS=${BUILD_DESIGNS:-0}
else
  OUTPUT_ROOT=${OUTPUT_ROOT:-"${DATA_ROOT}/processed/primer_positional_bias"}
  BUILD_DESIGNS=${BUILD_DESIGNS:-1}
fi
STATS=${STATS:-"${OUTPUT_ROOT}/demultiplex_stats.json"}
REPORT=${REPORT:-"${OUTPUT_ROOT}/design_report.json"}
mkdir -p "${OUTPUT_ROOT}" "${RUN_ROOT}/logs"

polydt_quant="${OUTPUT_ROOT}/salmon_polydt"
ranhex_quant="${OUTPUT_ROOT}/salmon_ranhex"
if [[ -s "${polydt_quant}/.tealeaf_complete" &&
      -s "${ranhex_quant}/.tealeaf_complete" ]]; then
  for quant in "${polydt_quant}" "${ranhex_quant}"; do
    observed_type=$("${PYTHON_BIN}" -c '
import json, sys
types = json.load(open(sys.argv[1])).get("library_types", [])
print(",".join(types))
' "${quant}/aux_info/meta_info.json")
    if [[ "${observed_type}" != "${LIBTYPE}" ]]; then
      echo "Refusing completed Salmon output with library type ${observed_type}; expected ${LIBTYPE}: ${quant}" >&2
      exit 1
    fi
  done
  echo "Reusing completed primer-specific Salmon quantifications"
else
  if [[ -e "${polydt_quant}" || -e "${ranhex_quant}" ]]; then
    echo "Refusing partial primer-specific Salmon output" >&2
    exit 1
  fi
  mkdir -p "/scratch/${USER}"
  fifo_root=$(mktemp -d "/scratch/${USER}/gse233208-posbias.XXXXXX")
  polydt_fifo="${fifo_root}/polydt.fastq"
  ranhex_fifo="${fifo_root}/ranhex.fastq"
  mkfifo "${polydt_fifo}" "${ranhex_fifo}"

  salmon_common=(
    quant -i "${SALMON_REF}/index_dir" -l "${LIBTYPE}"
    --validateMappings --posBias --dumpEq --dumpEqWeights
    --fldMean 250 --fldSD 25 --rangeFactorizationBins 4
    -p "${SALMON_THREADS}"
  )
  "${SALMON_BIN}" "${salmon_common[@]}" -r "${polydt_fifo}" \
    -o "${polydt_quant}" > "${OUTPUT_ROOT}/salmon_polydt.log" 2>&1 &
  polydt_pid=$!
  "${SALMON_BIN}" "${salmon_common[@]}" -r "${ranhex_fifo}" \
    -o "${ranhex_quant}" > "${OUTPUT_ROOT}/salmon_ranhex.log" 2>&1 &
  ranhex_pid=$!
  trap 'kill "${polydt_pid}" "${ranhex_pid}" 2>/dev/null || true' EXIT

  demux_args=()
  for index in "${!read1[@]}"; do
    demux_args+=(--read1 "${read1[index]}" --read2 "${read2[index]}")
  done
  max_args=()
  [[ -n "${MAX_READS}" ]] && max_args+=(--max-reads "${MAX_READS}")
  prefix_args=()
  if [[ -z "${MAX_READS}" && "${BALANCED_PREFIX_READS}" -gt 0 ]]; then
    prefix_args+=(--balanced-prefix-reads "${BALANCED_PREFIX_READS}")
  fi
  "${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/demultiplex_parse_primers.py" \
    "${demux_args[@]}" \
    --rt-barcodes "${PARSE_RT_BARCODES}" \
    --polydt-output "${polydt_fifo}" \
    --ranhex-output "${ranhex_fifo}" \
    --stats "${STATS}" \
    --pigz "${PIGZ_BIN}" \
    --pigz-threads "${PIGZ_THREADS}" \
    "${max_args[@]}" \
    "${prefix_args[@]}"
  wait "${polydt_pid}"
  wait "${ranhex_pid}"
  trap - EXIT
  for quant in "${polydt_quant}" "${ranhex_quant}"; do
    test -s "${quant}/quant.sf"
    test -s "${quant}/aux_info/eq_classes.txt.gz"
    test -s "${quant}/aux_info/obs5_pos.gz"
    test -s "${quant}/aux_info/exp5_pos.gz"
    printf 'salmon_positional_bias_complete\n' > "${quant}/.tealeaf_complete"
  done
fi

if [[ "${BUILD_DESIGNS}" == 1 ]]; then
  "${PYTHON_BIN}" "${REPO_ROOT}/extra_scripts/build_salmon_primer_designs.py" \
    --alevin-dir "${ALEVIN_DIR}" \
    --salmon-ref "${SALMON_REF}/spliceu.fa" \
    --polydt-quant "${polydt_quant}" \
    --ranhex-quant "${ranhex_quant}" \
    --primer-pairs "${ALEVIN_DIR}/primer_pairs.tsv" \
    --output-report "${REPORT}"
fi
