#!/usr/bin/env bash
export RNASEQ_ROOT_DIR="${ROOT_DIR}"
export RNASEQ_CONFIG="${RNASEQ_CONFIG:-${ROOT_DIR}/config/project.env}"
[[ -f "${RNASEQ_CONFIG}" ]] || { echo "MISSING config: ${RNASEQ_CONFIG}" >&2; exit 1; }
set -a
source "${RNASEQ_CONFIG}"
set +a
for key in RNASEQ_SAMPLES_TSV RNASEQ_COUNTS_FILE REFERENCE_DIR HISAT2_INDEX GTF_FILE GENOME_FASTA STRING_INPUT_DIR; do
  value="${!key:-}"
  if [[ -n "${value}" && "${value}" != /* ]]; then printf -v "${key}" '%s/%s' "${ROOT_DIR}" "${value}"; fi
  export "${key}"
done
export RNASEQ_READS_INFO="${RNASEQ_SAMPLES_TSV}"
for argument in "$@"; do
  if [[ "${argument}" == --preflight ]]; then
    "${PYTHON_BIN:-python3}" "${ROOT_DIR}/scripts/input_contract.py" preflight
    exit $?
  fi
done
