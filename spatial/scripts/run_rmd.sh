#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RMD_DIR="${ROOT_DIR}/rmd"
LOG_DIR="${ROOT_DIR}/result/logs"
mkdir -p "${LOG_DIR}"

MODULES=(
  "02_qc_overview"
  "03_spatial_domains"
  "04_deconvolution_review"
  "05_neighborhood_communication"
  "06_segment_programs"
  "07_multi_sample_summary"
)

usage() {
  cat <<EOF
Usage:
  bash scripts/run_rmd.sh
  bash scripts/run_rmd.sh --list
  bash scripts/run_rmd.sh --modules 2,4-6
EOF
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

list_modules() {
  local i
  for i in "${!MODULES[@]}"; do
    printf '%2d  %s\n' "$((i + 2))" "${MODULES[$i]}"
  done
}

append_module_token() {
  local token="$1"
  local module
  [[ -n "${token}" ]] || return 0
  token="${token//,/ }"
  for module in ${token}; do
    if [[ "${module}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      local start="${BASH_REMATCH[1]}"
      local end="${BASH_REMATCH[2]}"
      for i in $(seq "${start}" "${end}"); do
        REQUESTED+=("${i}")
      done
    elif [[ "${module}" =~ ^[0-9]+$ ]]; then
      REQUESTED+=("${module}")
    else
      die "Invalid module token: ${module}"
    fi
  done
}

render_module() {
  local module_num="$1"
  local idx=$((module_num - 2))
  (( idx >= 0 && idx < ${#MODULES[@]} )) || die "Invalid Rmd module: ${module_num}"
  local module="${MODULES[$idx]}"
  local rmd="${RMD_DIR}/${module}.Rmd"
  [[ -f "${rmd}" ]] || die "Rmd not found: ${rmd}"
  echo "Rendering ${module}"
  Rscript -e "rmarkdown::render('${rmd}', output_dir='${RMD_DIR}', quiet=FALSE)"
}

REQUESTED=()
if [[ $# -eq 0 ]]; then
  for i in "${!MODULES[@]}"; do
    REQUESTED+=("$((i + 2))")
  done
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --list)
        list_modules
        exit 0
        ;;
      --modules)
        [[ $# -ge 2 ]] || die "--modules requires a value"
        append_module_token "$2"
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

LOG_FILE="${LOG_DIR}/run_rmd_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "Project root: ${ROOT_DIR}"
echo "Log: ${LOG_FILE}"
for module in "${REQUESTED[@]}"; do
  render_module "${module}"
done
