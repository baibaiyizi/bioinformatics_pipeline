#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Spatial transcriptomics pipeline
#
# Usage:
#   bash run_spatial.sh
#   bash run_spatial.sh --steps 3,5-8
#   bash run_spatial.sh --list
#
# The runner mirrors the local wgbs_sperm style:
# - one official entry point
# - selected-step execution
# - active-step result cleanup only
# - auditable logs under result/logs
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SPATIAL_CONFIG_FILE:-${ROOT_DIR}/config/spatial_config.yaml}"
MANIFEST_FILE="${SPATIAL_MANIFEST_FILE:-${ROOT_DIR}/01rawdata/spatial_info.tsv}"
PYTHON_BIN="${SPATIAL_PYTHON_BIN:-python3}"
LOG_DIR="${SPATIAL_LOG_DIR:-${ROOT_DIR}/result/logs}"
AUTO_VALIDATE="${SPATIAL_AUTO_VALIDATE:-1}"

STEPS=(
  "01_validate_inputs"
  "02_build_objects"
  "03_qc_filter"
  "04_normalize_features"
  "05_integrate_cluster"
  "06_spatial_neighbors_stats"
  "07_spatial_domains"
  "08_cell_type_annotation"
  "09_deconvolution"
  "10_segment_objects"
  "11_programs_differential"
  "12_neighborhood_communication"
  "13_multi_sample_summary"
  "14_report_index"
)

usage() {
  cat <<EOF
Usage:
  bash run_spatial.sh                 # run all steps
  bash run_spatial.sh --list          # list steps
  bash run_spatial.sh --steps 3,5-8   # run selected steps

Environment overrides:
  SPATIAL_CONFIG_FILE      default: config/spatial_config.yaml
  SPATIAL_MANIFEST_FILE    default: 01rawdata/spatial_info.tsv
  SPATIAL_PYTHON_BIN       default: python3
  SPATIAL_AUTO_VALIDATE    default: 1
EOF
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

list_steps() {
  local i
  for i in "${!STEPS[@]}"; do
    printf '%2d  %s\n' "$((i + 1))" "${STEPS[$i]}"
  done
}

split_step_tokens() {
  local token="$1"
  token="${token//,/ }"
  printf '%s\n' ${token}
}

append_step_token() {
  local token="$1"
  local step
  [[ -n "${token}" ]] || return 0
  if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    local start="${BASH_REMATCH[1]}"
    local end="${BASH_REMATCH[2]}"
    if (( start < 1 || end > ${#STEPS[@]} || start > end )); then
      die "Invalid step range: ${token}. Valid range: 1-${#STEPS[@]}"
    fi
    for step in $(seq "${start}" "${end}"); do
      REQUESTED_STEPS+=("${step}")
    done
  elif [[ "${token}" =~ ^[0-9]+$ ]]; then
    if (( token < 1 || token > ${#STEPS[@]} )); then
      die "Invalid step: ${token}. Valid range: 1-${#STEPS[@]}"
    fi
    REQUESTED_STEPS+=("${token}")
  else
    die "Invalid step token: ${token}"
  fi
}

dedupe_steps() {
  local -n input_ref="$1"
  local -a out=()
  local seen=" "
  local item
  for item in "${input_ref[@]}"; do
    if [[ "${seen}" != *" ${item} "* ]]; then
      out+=("${item}")
      seen="${seen}${item} "
    fi
  done
  printf '%s\n' "${out[@]}"
}

guard_result_dir() {
  local step_name="$1"
  local result_dir="${ROOT_DIR}/result/${step_name}"
  case "${result_dir}" in
    "${ROOT_DIR}/result/"*) ;;
    *) die "Refusing to clean unsafe result path: ${result_dir}" ;;
  esac
  rm -rf "${result_dir}"
  mkdir -p "${result_dir}/tables" "${result_dir}/plots"
}

run_step() {
  local step_num="$1"
  local step_name="${STEPS[$((step_num - 1))]}"
  echo "[$(date '+%F %T')] START ${step_num} ${step_name}"
  guard_result_dir "${step_name}"
  "${PYTHON_BIN}" "${ROOT_DIR}/scripts/spatial_pipeline.py" run-step \
    --step "${step_name}" \
    --project-root "${ROOT_DIR}" \
    --config "${CONFIG_FILE}" \
    --manifest "${MANIFEST_FILE}"
  echo "[$(date '+%F %T')] END   ${step_num} ${step_name}"
}

REQUESTED_STEPS=()
if [[ $# -eq 0 ]]; then
  for i in "${!STEPS[@]}"; do
    REQUESTED_STEPS+=("$((i + 1))")
  done
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --list)
        list_steps
        exit 0
        ;;
      --steps)
        [[ $# -ge 2 ]] || die "--steps requires a value"
        while IFS= read -r token; do
          append_step_token "${token}"
        done < <(split_step_tokens "$2")
        shift 2
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      *)
        die "Unknown argument: $1"
        ;;
    esac
  done
fi

[[ -f "${CONFIG_FILE}" ]] || die "Config not found: ${CONFIG_FILE}"
[[ -f "${MANIFEST_FILE}" ]] || die "Manifest not found: ${MANIFEST_FILE}"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/run_spatial_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${LOG_FILE}") 2>&1

if [[ "${AUTO_VALIDATE}" == "1" ]]; then
  has_validate=0
  for step in "${REQUESTED_STEPS[@]}"; do
    [[ "${step}" == "1" ]] && has_validate=1
  done
  if [[ "${has_validate}" == "0" ]]; then
    REQUESTED_STEPS=("1" "${REQUESTED_STEPS[@]}")
  fi
fi

mapfile -t REQUESTED_STEPS < <(dedupe_steps REQUESTED_STEPS)

echo "Project root: ${ROOT_DIR}"
echo "Config:       ${CONFIG_FILE}"
echo "Manifest:     ${MANIFEST_FILE}"
echo "Log:          ${LOG_FILE}"
echo "Steps:        ${REQUESTED_STEPS[*]}"

for step_num in "${REQUESTED_STEPS[@]}"; do
  run_step "${step_num}"
done

echo "[$(date '+%F %T')] Pipeline finished"
