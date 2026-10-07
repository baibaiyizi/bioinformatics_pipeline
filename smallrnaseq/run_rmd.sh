#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES='1'
MODE=run
RUN_ID="run_$(date +%Y%m%d_%H%M%S)_$$"
while (($#)); do
  case "$1" in
    --help|-h) echo 'Usage: run_rmd.sh [--config project.env] [--modules "1,2,3"|--all] [--run-id ID] [--check|--dry-run|--list]'; echo 'Default: module 1 expression matrices + descriptive summary report. Optional modules 2-8; module 9 requires RUN_ENDOSIRNA=1.'; exit 0 ;;
    --list|--list-steps) printf '%s\n' '1 setup+expression+summary report' '2 per-library QC' '3 cohort PCA (all qualified samples)' '4 global sncRNA' '5 tsRNA' '6 rsRNA' '7 miRNA' '8 piRNA' '9 optional endo-siRNA'; exit 0 ;;
    --config) export SMALLRNA_CONFIG="$2"; shift 2 ;;
    --modules) MODULES="$2"; shift 2 ;;
    --all) MODULES='1,2,3,4,5,6,7,8,9'; shift ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --check) MODE=check; shift ;;
    --dry-run) MODE=dry-run; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid run-id' >&2; exit 2; }
source "$ROOT_DIR/scripts/load_config.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
export SMALLRNA_RESULT_DIR="$ANALYSIS_ROOT/$RUN_ID"
IFS=', ' read -r -a requested <<< "$MODULES"
ordered=(1)
for token in "${requested[@]}"; do
  [[ "$token" =~ ^0?[1-9]$ ]] || { echo "Invalid module: $token" >&2; exit 2; }
  number=$((10#$token))
  [[ " ${ordered[*]} " == *" $number "* ]] || ordered+=("$number")
done
python3 "$ROOT_DIR/scripts/preflight.py" --mode "$( [[ "$MODE" == dry-run ]] && echo config || echo downstream )"
if [[ "$MODE" == dry-run ]]; then printf 'modules=%s\noutput=%s\n' "${ordered[*]}" "$SMALLRNA_RESULT_DIR"; exit 0; fi
"$R_SCRIPT" -e 'for(p in c("ggplot2","dplyr","tidyr","stringr","purrr","tibble","knitr")) if(!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)'
if [[ "$MODE" == check ]]; then exit 0; fi
[[ ! -e "$SMALLRNA_RESULT_DIR" ]] || { echo "Existing run-id rejected: $SMALLRNA_RESULT_DIR" >&2; exit 1; }
mkdir -p "$SMALLRNA_RESULT_DIR/logs" "$SMALLRNA_RESULT_DIR/audit"
cp "$SMALLRNA_CONFIG" "$SMALLRNA_RESULT_DIR/audit/project.env.frozen"
cp "$SAMPLES_TSV" "$SMALLRNA_RESULT_DIR/audit/samples.tsv.frozen"
sha256sum "$SMALLRNA_CONFIG" "$SAMPLES_TSV" "$ROOT_DIR"/rmd/* > "$SMALLRNA_RESULT_DIR/audit/input_code_sha256.txt"
python3 "$ROOT_DIR/scripts/preflight.py" --mode downstream --sha256 > "$SMALLRNA_RESULT_DIR/audit/input_sha256.tsv"
printf 'module\texecution\n' > "$SMALLRNA_RESULT_DIR/audit/module_status.tsv"
ACTIVE_MODULE=prepare
trap 'code=$?; if [[ $code -ne 0 ]]; then printf "%s\tfailed\n" "$ACTIVE_MODULE" >> "$SMALLRNA_RESULT_DIR/audit/module_status.tsv"; fi' EXIT
for number in "${ordered[@]}"; do
  if [[ "$number" == 9 && "$RUN_ENDOSIRNA" != 1 ]]; then printf '9\tdisabled\n' >> "$SMALLRNA_RESULT_DIR/audit/module_status.tsv"; continue; fi
  ACTIVE_MODULE="$number"
  export SMALLRNA_MODULE="$number"
  "$R_SCRIPT" "$ROOT_DIR/scripts/run_module.R" > "$SMALLRNA_RESULT_DIR/logs/module_${number}.log" 2>&1
  printf '%s\tcompleted\n' "$number" >> "$SMALLRNA_RESULT_DIR/audit/module_status.tsv"
done
ACTIVE_MODULE=report
"$R_SCRIPT" "$ROOT_DIR/scripts/render_summary.R" > "$SMALLRNA_RESULT_DIR/logs/report.log" 2>&1
printf 'report\tcompleted\n' >> "$SMALLRNA_RESULT_DIR/audit/module_status.tsv"
echo "$SMALLRNA_RESULT_DIR/report.html"
