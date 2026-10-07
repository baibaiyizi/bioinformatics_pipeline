#!/usr/bin/env bash
# Caller provides ROOT_DIR and optionally SMALLRNA_CONFIG. No writes.
SMALLRNA_CONFIG="${SMALLRNA_CONFIG:-${ROOT_DIR}/config/project.env}"
[[ -f "$SMALLRNA_CONFIG" ]] || { echo "MISSING config: $SMALLRNA_CONFIG" >&2; exit 1; }
set -a
source "$SMALLRNA_CONFIG"
set +a
export SMALLRNA_ROOT_DIR="$ROOT_DIR"
for variable in SAMPLES_TSV SPORTS_RESULT_DIR ANALYSIS_ROOT SPORTS_PL REFERENCE_ROOT DB_ROOT GENOME_FASTA GENOME_BASE MIRNA_BASE RRNA_BASE TRNA_BASE TRNA_MITO_BASE PIRNA_BASE ENSEMBL_BASE RFAM_BASE; do
  value="${!variable:-}"
  [[ -n "$value" ]] || { echo "Missing config variable: $variable" >&2; exit 1; }
  [[ "$value" == /* ]] || printf -v "$variable" '%s/%s' "$ROOT_DIR" "$value"
  export "$variable"
done
if [[ -n "${ENDOSIRNA_REFERENCE_TSV:-}" && "$ENDOSIRNA_REFERENCE_TSV" != /* ]]; then
  export ENDOSIRNA_REFERENCE_TSV="$ROOT_DIR/$ENDOSIRNA_REFERENCE_TSV"
fi
export SMALLRNA_CONFIG
