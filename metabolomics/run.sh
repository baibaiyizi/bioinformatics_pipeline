#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${METABO_CONFIG_FILE:-${ROOT_DIR}/config/project.yml}"

usage() {
  cat <<EOF
Usage:
  bash run.sh all [options]
  bash run.sh first [options]
  bash run.sh second [options]
  bash run.sh check [options]

Use metabo_first.sh --help or metabo_second.sh --help for step and method switches.
EOF
}

[[ $# -gt 0 ]] || set -- all
MODE="$1"
shift || true

case "${MODE}" in
  all)
    "${RSCRIPT_BIN:-Rscript}" "${ROOT_DIR}/pipeline.R" --config "${CONFIG_FILE}" --module all "$@"
    ;;
  first)
    METABO_CONFIG_FILE="${CONFIG_FILE}" bash "${ROOT_DIR}/metabo_first.sh" "$@"
    ;;
  second)
    METABO_CONFIG_FILE="${CONFIG_FILE}" bash "${ROOT_DIR}/metabo_second.sh" "$@"
    ;;
  check)
    "${RSCRIPT_BIN:-Rscript}" "${ROOT_DIR}/pipeline.R" --config "${CONFIG_FILE}" --check "$@"
    ;;
  -h|--help) usage ;;
  *) echo "Unknown mode: ${MODE}" >&2; usage >&2; exit 1 ;;
esac
