#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
CONFIG_FILE="${METABO_CONFIG_FILE:-${ROOT_DIR}/config/project.yml}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
DRY_RUN=0
FORCE=0
THREADS=""
RESULT_DIR=""
REQUESTED=()
METHOD_ARGS=()

declare -A STEP_NAMES=(
  [1]="inventory"
  [2]="xcms"
)
DEFAULT_STEPS=(1 2)

log() { echo "[metabo-first] $*"; }
die() { echo "[metabo-first][error] $*" >&2; exit 1; }

usage() {
  cat <<EOF
Usage:
  bash ${SCRIPT_NAME} [--all]
  bash ${SCRIPT_NAME} --steps 1,2
  bash ${SCRIPT_NAME} --list | --dry-run | --force

Steps:
  1 inventory   inventory mzML/peak-table, metadata and real acquisition order
  2 xcms        modern XCMS detection, alignment, grouping and gap filling

Options:
  --config FILE       YAML configuration (default: ${CONFIG_FILE})
  --result-dir DIR    override result directory
  --steps LIST        comma-separated step ids or names
  --threads N         XCMS worker count override
  --dry-run           show actions without writing results
  --force             recompute selected stages with the same verified inputs
EOF
}

parse_steps() {
  local token part name
  local -a out=()
  for token in "$@"; do
    IFS=',' read -r -a parts <<< "${token}"
    for part in "${parts[@]}"; do
      [[ -z "${part}" ]] && continue
      if [[ "${part}" =~ ^[0-9]+$ ]]; then
        part="$((10#${part}))"
        (( part >= 1 && part <= 2 )) || die "step out of range: ${part}"
        out+=("${part}")
      else
        name="${part,,}"
        case "${name}" in
          inventory) out+=(1) ;;
          xcms) out+=(2) ;;
          *) die "unknown step: ${part}" ;;
        esac
      fi
    done
  done
  ((${#out[@]})) || die "no steps selected"
  printf '%s\n' "${out[@]}" | awk '!seen[$0]++'
}

run_module() {
  local module="$1"
  local -a cmd=("${RSCRIPT_BIN}" "${ROOT_DIR}/scripts/pipeline.R" --config "${CONFIG_FILE}" --module "${module}")
  [[ -n "${RESULT_DIR}" ]] && cmd+=(--result-dir "${RESULT_DIR}")
  [[ -n "${THREADS}" ]] && cmd+=(--threads "${THREADS}")
  (( FORCE )) && cmd+=(--force)
  (( DRY_RUN )) && cmd+=(--dry-run)
  cmd+=("${METHOD_ARGS[@]}")
  log "${cmd[*]}"
  (( DRY_RUN )) || "${cmd[@]}"
}

print_steps() { printf '  1: inventory\n  2: xcms\n'; }

while (($#)); do
  case "$1" in
    --all) REQUESTED=("${DEFAULT_STEPS[@]}"); shift ;;
    --list) print_steps; exit 0 ;;
    --config) [[ $# -ge 2 ]] || die "--config requires a file"; CONFIG_FILE="$2"; shift 2 ;;
    --result-dir) [[ $# -ge 2 ]] || die "--result-dir requires a directory"; RESULT_DIR="$2"; shift 2 ;;
    --steps|-s) shift; [[ $# -ge 1 ]] || die "--steps requires a value"; parsed_steps="$(parse_steps "$1")" || exit $?; mapfile -t REQUESTED <<< "${parsed_steps}"; shift ;;
    --threads) [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || die "--threads requires a positive integer"; THREADS="$2"; shift 2 ;;
    --qc-role|--qc-filter|--drift-correction|--normalization|--run-order-adjustment|--significance|--company-comparison|--qc-rsd-cutoff|--detection-rate|--fc-cutoff|--p-cutoff|--fdr-cutoff)
      [[ $# -ge 2 ]] || die "$1 requires a value"; METHOD_ARGS+=("$1" "$2"); shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    [0-9]*|inventory|xcms) parsed_steps="$(parse_steps "$@")" || exit $?; mapfile -t REQUESTED <<< "${parsed_steps}"; break ;;
    *) die "unknown option: $1" ;;
  esac
done
(( ${#REQUESTED[@]} )) || REQUESTED=("${DEFAULT_STEPS[@]}")

mapfile -t REQUESTED < <(printf '%s\n' "${REQUESTED[@]}" | sort -n -u)
declare -A need=()
for step in "${REQUESTED[@]}"; do
  case "$step" in
    1) need[1]=1 ;;
    2) need[1]=1; need[2]=1 ;;
  esac
done
for step in 1 2; do
  [[ "${need[$step]:-0}" == 1 ]] && run_module "${STEP_NAMES[$step]}"
done
log "upstream finished"
