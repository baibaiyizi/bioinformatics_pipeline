#!/usr/bin/env bash
set -euo pipefail
METHOD_COMPARISON_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export METHOD_COMPARISON_ROOT
stage=all
check=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --all) stage=all; shift ;;
    --stage) stage="$2"; shift 2 ;;
    --check) check=1; shift ;;
    --help|-h) echo '用法：bash run.sh [--check] [--all | --stage describe|differential|report]'; exit 0 ;;
    *) echo "未知参数：$1" >&2; exit 2 ;;
  esac
done
case "$stage" in all|describe|differential|report) ;; *) echo '无效 stage' >&2;exit 2 ;; esac
source "$METHOD_COMPARISON_ROOT/config/project.sh"
export PANDORA_PROJECT_ROOT="$(cd "$METHOD_COMPARISON_ROOT/$SMALLRNA_SOURCE" && pwd)"
unset PROJECT_SPECIES RUN_ENDOSIRNA DB_ROOT REFERENCE_ROOT SPORTS_PL ENDOSIRNA_REFS
set -a
source "$PANDORA_PROJECT_ROOT/config/project.sh"
set +a
export PANDORA_SHARED_ROOT="$METHOD_COMPARISON_ROOT/shared"
export PANDORA_SPECIES=Homo PANDORA_MODE=without_endo
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
python3 "$METHOD_COMPARISON_ROOT/scripts/check_inputs.py"
[[ "$check" == 0 ]] || exit 0
for current in describe differential report; do
  [[ "$stage" == all || "$stage" == "$current" ]] || continue
  if [[ "$current" == report ]]; then
    python3 "$METHOD_COMPARISON_ROOT/scripts/build_report.py"
  else
    export PANDORA_STAGE="$current" PANDORA_RESULT_DIR="$METHOD_COMPARISON_ROOT/result/$current"
    export PANDORA_SAMPLES_FILE="$PANDORA_PROJECT_ROOT/config/samples.tsv"
    export ROOT_DIR="$PANDORA_PROJECT_ROOT" RESULT_DIR="$PANDORA_RESULT_DIR" SHARED_ROOT="$PANDORA_SHARED_ROOT" SPECIES=Homo STAGE="$current"
    export PROJECT_CONTRAST_BASELINE=1 PROJECT_CONTRAST_CASE=2
    Rscript "$PANDORA_SHARED_ROOT/R/run_analysis.R" --stage "$current"
  fi
done
