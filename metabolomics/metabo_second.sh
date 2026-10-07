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
  [1]="audit"
  [2]="preprocess"
  [3]="plots"
  [4]="annotate"
  [5]="stats"
  [6]="pathway"
  [7]="report"
)
DEFAULT_STEPS=(1 2 3 4 5 6 7)

log() { echo "[metabo-second] $*"; }
die() { echo "[metabo-second][error] $*" >&2; exit 1; }

usage() {
  cat <<EOF
Usage:
  bash ${SCRIPT_NAME} [--all]
  bash ${SCRIPT_NAME} --steps 2,4,5
  bash ${SCRIPT_NAME} --list --dry-run --force

Steps:
  1 audit       metadata, order, effective config, software and input hashes
  2 preprocess  pre-gap-fill detection filter, QC role/filter and normalization
  3 plots       QC, PCA, correlation and run-order figures
  4 annotate    MetDNA exact/nearest m/z-RT evidence matching
  5 stats       limma case-control differential model and candidate tiers
  6 pathway     pathway quantitative supplement and optional KEGG ORA
  7 report      Chinese report and output map

Method switches (CLI overrides YAML):
  --qc-role unverified|pooled
  --qc-filter off|rsd
  --qc-rsd-cutoff FLOAT
  --drift-correction none|linear
  --normalization median|pqn|none
  --run-order-adjustment off|on
  --significance dual|fdr|nominal
  --company-comparison on|off
  --detection-rate FLOAT --fc-cutoff FLOAT --p-cutoff FLOAT --fdr-cutoff FLOAT

General options:
  --config FILE --result-dir DIR --threads N --dry-run --force --list
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
        (( part >= 1 && part <= 7 )) || die "step out of range: ${part}"
        out+=("${part}")
      else
        name="${part,,}"
        case "${name}" in
          audit|inventory) out+=(1) ;;
          preprocess|qc) out+=(2) ;;
          plots|plot) out+=(3) ;;
          annotate|annotation) out+=(4) ;;
          stats|stat|differential) out+=(5) ;;
          pathway) out+=(6) ;;
          report) out+=(7) ;;
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
  local -a cmd=("${RSCRIPT_BIN}" "${ROOT_DIR}/pipeline.R" --config "${CONFIG_FILE}" --module "${module}")
  [[ -n "${RESULT_DIR}" ]] && cmd+=(--result-dir "${RESULT_DIR}")
  [[ -n "${THREADS}" ]] && cmd+=(--threads "${THREADS}")
  (( FORCE )) && cmd+=(--force)
  (( DRY_RUN )) && cmd+=(--dry-run)
  cmd+=("${METHOD_ARGS[@]}")
  log "${cmd[*]}"
  (( DRY_RUN )) || "${cmd[@]}"
}

print_steps() {
  for i in $(seq 1 7); do printf '  %d: %s\n' "$i" "${STEP_NAMES[$i]}"; done
}

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
    [0-9]*|audit|inventory|preprocess|qc|plots|plot|annotate|annotation|stats|stat|differential|pathway|report) parsed_steps="$(parse_steps "$1")" || exit $?; mapfile -t REQUESTED <<< "${parsed_steps}"; shift ;;
    *) die "unknown option: $1" ;;
  esac
done
(( ${#REQUESTED[@]} )) || REQUESTED=("${DEFAULT_STEPS[@]}")
mapfile -t REQUESTED < <(printf '%s\n' "${REQUESTED[@]}" | sort -n -u)

# Resolve prerequisites: statistics needs preprocessing and annotation;
# pathway/report consume the complete downstream state. The xcms module is
# the shared upstream prerequisite for every downstream step.
declare -A need=()
for step in "${REQUESTED[@]}"; do
  case "$step" in
    1) need[1]=1 ;;
    2) need[1]=1; need[2]=1 ;;
    3) need[1]=1; need[2]=1; need[3]=1 ;;
    4) need[1]=1; need[2]=1; need[4]=1 ;;
    5) need[1]=1; need[2]=1; need[4]=1; need[5]=1 ;;
    6) need[1]=1; need[2]=1; need[4]=1; need[5]=1; need[6]=1 ;;
    7) need[1]=1; need[2]=1; need[3]=1; need[4]=1; need[5]=1; need[6]=1; need[7]=1 ;;
  esac
done
run_module audit
if [[ "${need[2]:-0}" == 1 || "${need[3]:-0}" == 1 || "${need[4]:-0}" == 1 || "${need[5]:-0}" == 1 || "${need[6]:-0}" == 1 || "${need[7]:-0}" == 1 ]]; then
  run_module xcms
fi
for step in 2 3 4 5 6 7; do
  [[ "${need[$step]:-0}" == 1 ]] && run_module "${STEP_NAMES[$step]}"
done
log "downstream finished"
