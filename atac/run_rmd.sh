#!/usr/bin/env bash
set -euo pipefail

# 独立 ATAC 下游模块入口；配置集中于 config/project.env。默认保留结果。
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"

source "$ROOT_DIR/scripts/load_config.sh"
case "$ATAC_INPUT_MODE" in
 counts) exec bash "$ROOT_DIR/scripts/run_counts.sh" "$@" ;;
 bam) ;;
 *) echo "ATAC_INPUT_MODE 必须为 bam/counts" >&2; exit 2 ;;
esac

# ==============================================================================
# 以下为步骤调度代码，一般不修改
# ==============================================================================

RMD_DIR="${ROOT_DIR}/rmd"

RESULT_DIR="${ROOT_DIR}/result"
LOG_DIR="${RESULT_DIR}/logs"
REPORT_SCRIPT="${RMD_DIR}/generate_result_report.R"

TOTAL_STEPS=16
ACTIVE_STEPS=(1 2 3 4 5 6 7 8 10 12 13 14 15 16)
DRY_RUN=false
RUN_ALL=false

declare -A STEP_FILES=(
  [1]="01_main_setup.R"
  [2]="02_project_inventory.Rmd"
  [3]="03_atac_qc.Rmd"
  [4]="04_peak_catalog_annotation.Rmd"
  [5]="05_differential_accessibility.Rmd"
  [6]="06_genome_tracks.Rmd"
  [7]="07_peak_gene_enrichment.Rmd"
  [8]="08_motif_activity.Rmd"
  [10]="10_footprinting_tf_evidence.Rmd"
  [12]="12_tf_gene_network.Rmd"
  [13]="13_conservation_gc.Rmd"
  [14]="14_machine_learning.Rmd"
  [15]="15_keygenes_locus_reports.Rmd"
  [16]="16_export_summary_manifest.Rmd"
)

declare -A STEP_TITLES=(
  [1]="shared setup, paths and helpers"
  [2]="project inventory and run context"
  [3]="ATAC QC and integrated QC grading"
  [4]="peak catalog, annotation and coverage profiles"
  [5]="differential accessibility with configurable significance threshold"
  [6]="genome tracks"
  [7]="peak-associated gene enrichment"
  [8]="motif enrichment and motif activity"
  [10]="footprinting and TF evidence"
  [12]="TF-gene network"
  [13]="sequence conservation and GC context"
  [14]="machine-learning feature ranking"
  [15]="key genes and locus reports"
  [16]="export, summary and result manifest"
)

declare -A STEP_RESULT_DIRS=(
  [2]="02_project_inventory"
  [3]="03_atac_qc"
  [4]="04_peak_catalog_annotation"
  [5]="05_differential_accessibility"
  [6]="06_genome_tracks"
  [7]="07_peak_gene_enrichment"
  [8]="08_motif_activity"
  [10]="10_footprinting_tf_evidence"
  [12]="12_tf_gene_network"
  [13]="13_sequence_conservation_gc"
  [14]="14_machine_learning"
  [15]="15_keygenes_locus_reports"
  [16]="16_export_summary_manifest"
)

step_abs_path() {
  local step_id="$1"
  printf '%s/%s\n' "${RMD_DIR}" "${STEP_FILES[${step_id}]}"
}

step_result_dir_name() {
  local step_id="$1"
  printf '%s\n' "${STEP_RESULT_DIRS[${step_id}]:-}"
}

print_step_table() {
  local i
  for i in "${ACTIVE_STEPS[@]}"; do
    local result_hint
    result_hint="$(step_result_dir_name "${i}")"
    if [[ -n "${result_hint}" ]]; then
      result_hint=" -> result/${result_hint}"
    fi
    printf '  %02d  %-14s  %s (%s)%s\n' "${i}" "$(step_alias_hint "${i}")" "${STEP_TITLES[${i}]}" "$(step_abs_path "${i}")" "${result_hint}"
  done
}

usage() {
  cat <<EOF
Usage:
  ${SCRIPT_NAME}
  ${SCRIPT_NAME} --all
  ${SCRIPT_NAME} --list
  ${SCRIPT_NAME} --dry-run --steps qc,diff,motif,export
  ${SCRIPT_NAME} --steps 5,8
  ${SCRIPT_NAME} --steps 5 8
  ${SCRIPT_NAME} 5 8

Behavior:
  1) No arguments or --all runs the current ATAC steps (9 and 11 have been removed).
  2) --steps accepts ids, aliases, comma-separated values, and numeric ranges.
  3) Positional step ids/aliases are accepted for quick use.
  4) Requested steps auto-run their required upstream steps in scientific order.
  5) 保留现有输出；配置签名不匹配则拒绝复用。缺省关闭跨组学、footprint、保守性；机器学习停用。
  6) 完整参数见 config/project.env。常用参数：ATAC_RUN_PROFILE, ATAC_REUSE_STEP_TABLES,
     ATAC_STRICT_MODE, ATAC_SAVE_PNG_PREVIEWS, ATAC_PNG_DPI,
     ATAC_SIGNIFICANCE_MODE, ATAC_PRIMARY_REQUIRE_FC, ATAC_PRIMARY_P_CUTOFF,
     ATAC_DIFF_FDR_CUTOFF, ATAC_DIFF_LOG2FC_CUTOFF, ATAC_DIFFBIND_RPKM_FILTER,
     ATAC_REFERENCE_GROUP, ATAC_COMPARISON_GROUP, ATAC_RUN_DEEPTOOLS_PROFILES,
     ATAC_RUN_MEME_TOMTOM, ATAC_TOMTOM_DB.

Steps:
$(print_step_table)
EOF
}

step_alias_hint() {
  case "$1" in
    1) echo "setup" ;;
    2) echo "inventory" ;;
    3) echo "qc" ;;
    4) echo "peaks" ;;
    5) echo "diff" ;;
    6) echo "tracks" ;;
    7) echo "enrichment" ;;
    8) echo "motif" ;;
    10) echo "footprint" ;;
    12) echo "network" ;;
    13) echo "conservation" ;;
    14) echo "ml" ;;
    15) echo "keygenes" ;;
    16) echo "export" ;;
  esac
}

ensure_rscript() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "[ERROR] Rscript not found in PATH" >&2
    exit 1
  fi
}

validate_public_config() {
  local mode="${ATAC_SIGNIFICANCE_MODE:-pvalue}"
  mode="$(printf '%s' "${mode}" | tr '[:upper:]' '[:lower:]')"
  case "${mode}" in
    pvalue|padj|fdr) ;;
    *) echo "[ERROR] ATAC_SIGNIFICANCE_MODE must be pvalue, padj, or fdr; got: ${ATAC_SIGNIFICANCE_MODE}" >&2; exit 2 ;;
  esac

  local profile="${ATAC_RUN_PROFILE:-full}"
  case "$(printf '%s' "${profile}" | tr '[:upper:]' '[:lower:]')" in
    full|quick|report) ;;
    *) echo "[ERROR] ATAC_RUN_PROFILE must be full, quick, or report; got: ${profile}" >&2; exit 2 ;;
  esac

  local boolean_name boolean_value
  for boolean_name in ATAC_PRIMARY_REQUIRE_FC ATAC_REUSE_STEP_TABLES ATAC_STRICT_MODE \
    ATAC_SHOW_CODE ATAC_CACHE_RESULTS ATAC_EMIT_LIGHT_PLOTS ATAC_SAVE_PNG_PREVIEWS \
    ATAC_DEEP_PLOT_AUDIT ATAC_LIGHT_SETUP ATAC_RUN_DEEPTOOLS_PROFILES \
    ATAC_RUN_MEME_TOMTOM ATAC_RUN_PHASTCONS ATAC_ENABLE_LIFTOVER_RELATED_GENE_ANNOTATION; do
    boolean_value="${!boolean_name:-FALSE}"
    case "$(printf '%s' "${boolean_value}" | tr '[:upper:]' '[:lower:]')" in
      0|1|true|false|yes|no|y|n|on|off) ;;
      *) echo "[ERROR] ${boolean_name} must be a boolean value; got: ${boolean_value}" >&2; exit 2 ;;
    esac
  done

  local reference="${ATAC_REFERENCE_GROUP:-}"
  local comparison="${ATAC_COMPARISON_GROUP:-}"
  if [[ -n "${reference}" && -z "${comparison}" ]] || [[ -z "${reference}" && -n "${comparison}" ]]; then
    echo "[ERROR] ATAC_REFERENCE_GROUP and ATAC_COMPARISON_GROUP must be set together." >&2
    exit 2
  fi
  if [[ -n "${reference}" && "${reference}" == "${comparison}" ]]; then
    echo "[ERROR] ATAC_REFERENCE_GROUP and ATAC_COMPARISON_GROUP must be different." >&2
    exit 2
  fi

  local name value
  for name in ATAC_PRIMARY_P_CUTOFF ATAC_DIFF_FDR_CUTOFF ATAC_NOMINAL_P_CUTOFF; do
    value="${!name:-0.05}"
    awk -v v="${value}" 'BEGIN { if (v !~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || v <= 0 || v > 1) exit 1 }' || {
      echo "[ERROR] ${name} must be numeric in (0,1]; got: ${value}" >&2
      exit 2
    }
  done
  for name in ATAC_DIFF_LOG2FC_CUTOFF ATAC_DIFFBIND_RPKM_FILTER; do
    value="${!name:-}"
    [[ -z "${value}" ]] && continue
    awk -v v="${value}" 'BEGIN { if (v !~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || v < 0) exit 1 }' || {
      echo "[ERROR] ${name} must be a non-negative number; got: ${value}" >&2
      exit 2
    }
  done

  case "$(printf '%s' "${ATAC_DIFFBIND_FILTER_FUN}" | tr '[:upper:]' '[:lower:]')" in
    max|mean) ;;
    *) echo "[ERROR] ATAC_DIFFBIND_FILTER_FUN must be max or mean; got: ${ATAC_DIFFBIND_FILTER_FUN}" >&2; exit 2 ;;
  esac

  for name in ATAC_PNG_DPI ATAC_PDF_PREVIEW_MAX_PAGES ATAC_FOOTPRINT_MIN_PRIMARY_PEAKS ATAC_DIFFBIND_MIN_OVERLAP; do
    value="${!name:-}"
    [[ "${value}" =~ ^[0-9]+$ ]] && (( value >= 1 )) || {
      echo "[ERROR] ${name} must be a positive integer; got: ${value}" >&2
      exit 2
    }
  done
  for name in ATAC_DIFFBIND_SUMMITS ATAC_DIFFBIND_MIN_COUNT; do
    value="${!name:-}"
    [[ "${value}" =~ ^[0-9]+$ ]] || {
      echo "[ERROR] ${name} must be a non-negative integer; got: ${value}" >&2
      exit 2
    }
  done
}

print_user_config_summary() {
  local prefix="${1:-[INFO]}"
  echo "${prefix} Config | profile=${ATAC_RUN_PROFILE}, reuse_tables=${ATAC_REUSE_STEP_TABLES}, strict=${ATAC_STRICT_MODE}"
  echo "${prefix} Config | primary=${ATAC_SIGNIFICANCE_MODE}, require_FC=${ATAC_PRIMARY_REQUIRE_FC}, p=${ATAC_PRIMARY_P_CUTOFF}, FDR=${ATAC_DIFF_FDR_CUTOFF}, |log2FC|=${ATAC_DIFF_LOG2FC_CUTOFF}"
  echo "${prefix} Config | DiffBind: minOverlap=${ATAC_DIFFBIND_MIN_OVERLAP}, summits=${ATAC_DIFFBIND_SUMMITS}, RPKM_filter=${ATAC_DIFFBIND_RPKM_FILTER}, filterFun=${ATAC_DIFFBIND_FILTER_FUN}, minCount=${ATAC_DIFFBIND_MIN_COUNT}"
  echo "${prefix} Config | nominal_p=${ATAC_NOMINAL_P_CUTOFF}"
  echo "${prefix} Config | contrast=${ATAC_COMPARISON_GROUP:-AUTO} vs ${ATAC_REFERENCE_GROUP:-AUTO}"
  echo "${prefix} Config | deepTools=${ATAC_RUN_DEEPTOOLS_PROFILES:-AUTO_BY_PROFILE}, MEME/Tomtom=${ATAC_RUN_MEME_TOMTOM:-AUTO_BY_PROFILE}, phastCons=${ATAC_RUN_PHASTCONS}"
}

validate_step_id() {
  local step_id="$1"
  [[ "${step_id}" =~ ^[0-9]+$ ]] || return 1
  [[ -n "${STEP_FILES[${step_id}]:-}" ]]
}

step_id_from_token() {
  local token="$1"
  token="$(printf '%s' "${token}" | tr '[:upper:]' '[:lower:]')"
  token="${token//_/-}"

  if [[ "${token}" =~ ^0*[0-9]+$ ]]; then
    token="$((10#${token}))"
    validate_step_id "${token}" || return 1
    printf '%s\n' "${token}"
    return 0
  fi

  case "${token}" in
    setup|main) echo 1 ;;
    inventory|project|context) echo 2 ;;
    qc|atac-qc|integrated-qc) echo 3 ;;
    peak|peaks|annotation|catalog|peak-annotation) echo 4 ;;
    diff|differential|diffbind|da|accessibility) echo 5 ;;
    track|tracks|gviz|genome-tracks) echo 6 ;;
    enrichment|peak-gene|geneset|go|kegg) echo 7 ;;
    motif|motifs|chromvar|activity) echo 8 ;;
    footprint|footprinting|tobias|tf-evidence) echo 10 ;;
    network|tf-network|tf-gene) echo 12 ;;
    conservation|gc|sequence) echo 13 ;;
    ml|machine-learning|machine|feature) echo 14 ;;
    keygene|keygenes|locus|loci) echo 15 ;;
    export|manifest|summary|delivery) echo 16 ;;
    *) return 1 ;;
  esac
}

parse_step_tokens() {
  local tokens=("$@")
  local token part start end i

  for token in "${tokens[@]}"; do
    token="${token//[[:space:]]/}"
    [[ -z "${token}" ]] && continue
    IFS=',' read -r -a raw_parts <<< "${token}"
    for part in "${raw_parts[@]}"; do
      part="${part//[[:space:]]/}"
      [[ -z "${part}" ]] && continue
      if [[ "${part}" =~ ^([0-9]+)(-|\\.\\.)([0-9]+)$ ]]; then
        start="$(step_id_from_token "${BASH_REMATCH[1]}")" || {
          echo "[ERROR] Invalid step range start: ${BASH_REMATCH[1]}" >&2
          exit 1
        }
        end="$(step_id_from_token "${BASH_REMATCH[3]}")" || {
          echo "[ERROR] Invalid step range end: ${BASH_REMATCH[3]}" >&2
          exit 1
        }
        if (( start <= end )); then
          for i in $(seq "${start}" "${end}"); do validate_step_id "${i}" || { echo "[ERROR] Removed or invalid step: ${i}" >&2; return 1; }; printf '%s\n' "${i}"; done
        else
          for i in $(seq "${start}" -1 "${end}"); do validate_step_id "${i}" || { echo "[ERROR] Removed or invalid step: ${i}" >&2; return 1; }; printf '%s\n' "${i}"; done
        fi
      else
        step_id_from_token "${part}" || {
          echo "[ERROR] Invalid step token: ${part}" >&2
          exit 1
        }
      fi
    done
  done
}

dedupe_preserve_order() {
  local input=("$@")
  local out=()
  local seen=" "
  local x
  for x in "${input[@]}"; do
    if [[ "${seen}" == *" ${x} "* ]]; then
      continue
    fi
    out+=("${x}")
    seen+="${x} "
  done
  printf '%s\n' "${out[@]}"
}

RUN_ORDER=()
declare -A RUN_NEEDED=()

mark_with_dependencies() {
  local step_id="$1"
  RUN_NEEDED["${step_id}"]=1
  case "${step_id}" in
    1)
      ;;
    2|3|4)
      mark_with_dependencies 1
      ;;
    5)
      mark_with_dependencies 1
      mark_with_dependencies 4
      ;;
    6|7|8)
      mark_with_dependencies 1
      mark_with_dependencies 4
      mark_with_dependencies 5
      ;;
    10)
      mark_with_dependencies 1
      mark_with_dependencies 4
      mark_with_dependencies 5
      mark_with_dependencies 8
      ;;
    12)
      mark_with_dependencies 1
      mark_with_dependencies 4
      mark_with_dependencies 5
      mark_with_dependencies 7
      mark_with_dependencies 8
      mark_with_dependencies 10
      ;;
    13|14)
      mark_with_dependencies 1
      mark_with_dependencies 4
      mark_with_dependencies 5
      ;;
    15)
      mark_with_dependencies 1
      mark_with_dependencies 4
      mark_with_dependencies 5
      mark_with_dependencies 7
      mark_with_dependencies 8
      mark_with_dependencies 10
      mark_with_dependencies 12
      mark_with_dependencies 13
      mark_with_dependencies 14
      ;;
    16)
      local i
      for i in "${ACTIVE_STEPS[@]}"; do [[ "$i" == 16 ]] || mark_with_dependencies "${i}"; done
      ;;
  esac
}

resolve_run_order() {
  local requested=("$@")
  local step_id i
  RUN_ORDER=()
  RUN_NEEDED=()
  for step_id in "${requested[@]}"; do
    mark_with_dependencies "${step_id}"
  done
  for i in "${ACTIVE_STEPS[@]}"; do
    if [[ -n "${RUN_NEEDED[${i}]:-}" ]]; then
      RUN_ORDER+=("${i}")
    fi
  done
}

ensure_step_file_exists() {
  local step_id="$1"
  local input_file
  input_file="$(step_abs_path "${step_id}")"
  if [[ ! -f "${input_file}" ]]; then
    echo "[ERROR] Missing step file: ${input_file}" >&2
    exit 1
  fi
}

run_step() {
  local step_id="$1"
  local input_file setup_file

  CURRENT_STEP="$step_id"
  printf "%s\t%s\tSTARTED\n" "$(date -Is)" "$step_id" >> "$RESULT_DIR/execution.tsv"
  input_file="$(step_abs_path "${step_id}")"
  setup_file="$(step_abs_path 1)"
  ensure_step_file_exists "${step_id}"

  echo "=============================================================="
  echo "[INFO] Step ${step_id}: ${STEP_TITLES[${step_id}]}"
  echo "[INFO] File: ${input_file}"
  echo "[INFO] Start: $(date '+%Y-%m-%d %H:%M:%S')"

  Rscript --vanilla "$ROOT_DIR/scripts/run_module.R" "$step_id" "$ROOT_DIR" "$RMD_DIR" "${STEP_FILES[$step_id]}"

  printf "%s\t%s\tCOMPLETED\n" "$(date -Is)" "$step_id" >> "$RESULT_DIR/execution.tsv"
  echo "[INFO] Done step ${step_id} at $(date '+%Y-%m-%d %H:%M:%S')"
}

STEP_ARGS=()

if [[ $# -eq 0 ]]; then
  RUN_ALL=true
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --report-only)
        python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC
        python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC --verify downstream
          RESULT_REPORT_ROOT="$ROOT_DIR" Rscript --vanilla "$ROOT_DIR/rmd/generate_result_report.R"
        exit 0
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      --list)
        print_step_table
        exit 0
        ;;
      --all)
        RUN_ALL=true
        shift
        ;;
      --dry-run|--check)
        DRY_RUN=true
        shift
        ;;
      --steps)
        shift
        if [[ $# -eq 0 ]]; then
          echo "[ERROR] --steps requires step ids or aliases" >&2
          exit 1
        fi
        while [[ $# -gt 0 && "$1" != --* ]]; do
          STEP_ARGS+=("$1")
          shift
        done
        ;;
      --*)
        echo "[ERROR] Unknown option: $1; use --help to see supported --steps syntax." >&2
        exit 2
        ;;
      *)
        STEP_ARGS+=("$1")
        shift
        ;;
    esac
  done
fi

if [[ "${RUN_ALL}" == "true" && ${#STEP_ARGS[@]} -eq 0 ]]; then
  REQUESTED_STEPS=("${ACTIVE_STEPS[@]}")
elif [[ ${#STEP_ARGS[@]} -eq 0 ]]; then
  REQUESTED_STEPS=("${ACTIVE_STEPS[@]}")
else
  REQUESTED_STEPS=($(parse_step_tokens "${STEP_ARGS[@]}"))
  REQUESTED_STEPS=($(dedupe_preserve_order "${REQUESTED_STEPS[@]}"))
fi

# 可选足迹和保守性按显式开关执行；ATAC 候选与网络不依赖跨组学开关。
filtered=()
for step in "${REQUESTED_STEPS[@]}"; do
  [[ "$step" != 14 ]] || { echo "[SKIPPED] 14 machine_learning 已停用" >&2; continue; }
  if [[ "$step" == 10 ]]; then [[ "$ATAC_ENABLE_FOOTPRINT" == 1 ]] || { echo "[SKIPPED] footprint 未启用: $step" >&2; continue; }; fi
  if [[ "$step" == 13 ]]; then [[ "$ATAC_ENABLE_CONSERVATION" == 1 ]] || { echo "[SKIPPED] conservation 未启用: $step" >&2; continue; }; fi
  filtered+=("$step")
done
REQUESTED_STEPS=("${filtered[@]}")
resolve_run_order "${REQUESTED_STEPS[@]}"
filtered=()
for step in "${RUN_ORDER[@]}"; do
  [[ "$step" != 14 ]] || { echo "[SKIPPED] 14 machine_learning 已停用" >&2; continue; }
  if [[ "$step" == 10 ]]; then [[ "$ATAC_ENABLE_FOOTPRINT" == 1 ]] || { echo "[SKIPPED] footprint 未启用: $step" >&2; continue; }; fi
  if [[ "$step" == 13 ]]; then [[ "$ATAC_ENABLE_CONSERVATION" == 1 ]] || { echo "[SKIPPED] conservation 未启用: $step" >&2; continue; }; fi
  filtered+=("$step")
done
RUN_ORDER=("${filtered[@]}")
python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC
validate_public_config

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "[DRY-RUN] Project root: ${ROOT_DIR}"
  echo "[DRY-RUN] Rmd template dir: ${RMD_DIR}"
  echo "[DRY-RUN] Requested steps: ${REQUESTED_STEPS[*]}"
  echo "[DRY-RUN] Final run order: ${RUN_ORDER[*]}"
  print_user_config_summary "[DRY-RUN]"
  for step_id in "${RUN_ORDER[@]}"; do
    ensure_step_file_exists "${step_id}"
    result_hint="$(step_result_dir_name "${step_id}")"
    if [[ -n "${result_hint}" ]]; then
      result_hint=" result/${result_hint}"
    else
      result_hint=" setup-only"
    fi
    printf '[DRY-RUN] %02d %-14s %s%s\n' "${step_id}" "$(step_alias_hint "${step_id}")" "$(step_abs_path "${step_id}")" "${result_hint}"
  done
  exit 0
fi

ensure_rscript
if compgen -G "$ROOT_DIR/02mapping/bam/*.final.bam" >/dev/null; then
    python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC --verify upstream
  fi
  python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix ATAC --record downstream
mkdir -p "${LOG_DIR}" "${ROOT_DIR}/.cache/rmd"
printf 'module\tenabled\n10_footprint\t%s\n13_conservation\t%s\n14_machine_learning\tDISABLED\n' "$ATAC_ENABLE_FOOTPRINT" "$ATAC_ENABLE_CONSERVATION" > "$RESULT_DIR/optional_module_gates.tsv"


[[ -f "$RESULT_DIR/execution.tsv" ]] || printf "time\tmodule\tstatus\n" > "$RESULT_DIR/execution.tsv"
for step in 10 13 14; do
  disabled=0
  case "$step" in
    10) [[ "$ATAC_ENABLE_FOOTPRINT" == 1 ]] || disabled=1 ;;
    13) [[ "$ATAC_ENABLE_CONSERVATION" == 1 ]] || disabled=1 ;;
    14) disabled=1 ;;
  esac
  [[ "$disabled" == 0 ]] || printf '%s\t%s\tSKIPPED_DISABLED\n' "$(date -Is)" "$step" >> "$RESULT_DIR/execution.tsv"
done
CURRENT_STEP=report
RUN_TAG="run_rmd_$(date +%Y%m%d_%H%M%S)"
RUN_LOG_FILE="${LOG_DIR}/${RUN_TAG}.log"
if [[ -e "${RUN_LOG_FILE}" ]]; then
  RUN_LOG_FILE="${LOG_DIR}/${RUN_TAG}_$$.log"
fi

exec > >(tee -a "${RUN_LOG_FILE}") 2>&1

print_failure_footer() {
  local exit_code="$1"
  echo "[ERROR] run_rmd.sh failed with exit code ${exit_code}"
  echo "[ERROR] Log file: ${RUN_LOG_FILE}"
}

on_exit() {
  local exit_code="$?"
  if [[ "${exit_code}" -ne 0 ]]; then
    printf "%s\t%s\tFAILED\n" "$(date -Is)" "$CURRENT_STEP" >> "$RESULT_DIR/execution.tsv"
    print_failure_footer "${exit_code}"
  fi
}
trap on_exit EXIT

echo "[INFO] Log file: ${RUN_LOG_FILE}"
echo "[INFO] Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] Project root: ${ROOT_DIR}"
echo "[INFO] Rmd template dir: ${RMD_DIR}"
echo "[INFO] Requested steps (dedup): ${REQUESTED_STEPS[*]}"
echo "[INFO] Final run order: ${RUN_ORDER[*]}"
print_user_config_summary "[INFO]"

for step_id in "${RUN_ORDER[@]}"; do
  run_step "${step_id}"
done

if [[ ! -f "${REPORT_SCRIPT}" ]]; then
  echo "[ERROR] Missing result report generator: ${REPORT_SCRIPT}" >&2
  exit 1
fi
CURRENT_STEP=report
printf '%s\treport\tSTARTED\n' "$(date -Is)" >> "$RESULT_DIR/execution.tsv"
echo "[INFO] Refresh result report: ${RESULT_DIR}/index.html"
RESULT_REPORT_ROOT="${ROOT_DIR}" Rscript --vanilla "${REPORT_SCRIPT}"
printf '%s\treport\tCOMPLETED\n' "$(date -Is)" >> "$RESULT_DIR/execution.tsv"

echo "[INFO] End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] Completed steps: ${RUN_ORDER[*]}"
