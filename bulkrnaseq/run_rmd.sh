#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
if [[ -n "${RNASEQ_RMD_DIR:-}" ]]; then
  RMD_DIR="$(cd "${RNASEQ_RMD_DIR}" && pwd)"
else
  RMD_DIR="${ROOT_DIR}/rmd"
  ##### 移至rnaseq下用 RMD_DIR="$(cd "${ROOT_DIR}/.." && pwd)/rmd"
fi
RESULT_DIR="${ROOT_DIR}/result"
LOG_DIR="${RESULT_DIR}/logs"
CACHE_DIR="${ROOT_DIR}/.cache/rmd"
CACHE_06_TABLES="${CACHE_DIR}/06_main_enrichment_tables.rds"
CACHE_09_GSVA_STATE="${CACHE_DIR}/09_gsva_state.rds"
CACHE_15_STATE="${CACHE_DIR}/15_splicing_state.rds"
CACHE_16_STATE="${CACHE_DIR}/16_immune_state.rds"

mkdir -p "${LOG_DIR}" "${CACHE_DIR}"

TOTAL_MODULES=18
# 文件编号已经按科学顺序重排，默认/--all 顺序即 01-18。
DEFAULT_MODULES=($(seq 1 "${TOTAL_MODULES}"))

declare -A MODULE_FILES=(
  [1]="01_main_setup.R"
  [2]="02_selection.Rmd"
  [3]="03_qc.Rmd"
  [4]="04_pca.Rmd"
  [5]="05_de.Rmd"
  [6]="06_enrichment.Rmd"
  [7]="07_go.Rmd"
  [8]="08_gsea.Rmd"
  [9]="09_gsva.Rmd"
  [10]="10_tf_regulator.Rmd"
  [11]="11_ppi.Rmd"
  [12]="12_wgcna.Rmd"
  [13]="13_apa.Rmd"
  [14]="14_der.Rmd"
  [15]="15_splicing.Rmd"
  [16]="16_immune.Rmd"
  [17]="17_pathway_activity.Rmd"
  [18]="18_circle.Rmd"
)

declare -A MODULE_TITLES=(
  [1]="main setup"
  [2]="selection"
  [3]="qc"
  [4]="pca"
  [5]="de"
  [6]="enrichment"
  [7]="go"
  [8]="gsea"
  [9]="gsva"
  [10]="tf regulator"
  [11]="ppi"
  [12]="wgcna"
  [13]="apa"
  [14]="der"
  [15]="splicing"
  [16]="immune"
  [17]="pathway activity"
  [18]="circle"
)

declare -A MODULE_RESULT_DIRS=(
  [1]="01_setup"
  [2]="02_selection"
  [3]="03_qc"
  [4]="04_pca"
  [5]="05_de"
  [6]="06_enrichment"
  [7]="07_go"
  [8]="08_gsea"
  [9]="09_gsva"
  [10]="10_tf_regulator"
  [11]="11_ppi"
  [12]="12_wgcna"
  [13]="13_apa"
  [14]="14_der"
  [15]="15_splicing"
  [16]="16_immune"
  [17]="17_pathway_activity"
  [18]="18_circle"
)

print_module_table() {
  local i
  for i in $(seq 1 "${TOTAL_MODULES}"); do
    printf '  %02d: %s (%s)\n' "${i}" "${MODULE_TITLES[${i}]}" "$(module_abs_path "${i}")"
  done
}

print_default_order() {
  printf '%s\n' "${DEFAULT_MODULES[*]}"
}

usage() {
  cat <<EOF
Usage:
  ${SCRIPT_NAME}
  ${SCRIPT_NAME} --all
  ${SCRIPT_NAME} --list
  ${SCRIPT_NAME} --modules 6,8,17
  ${SCRIPT_NAME} --modules 6 8 17
  ${SCRIPT_NAME} 6 8 17
  ${SCRIPT_NAME} 6,8,17

Behavior:
  1) 默认无参数/--all 执行全模块，顺序为: $(print_default_order)
  2) 模块 02-18 自动补 01
  3) 模块 09/10/13/14/15/16/17 自动加载 post context；11/12/18 仅加载 01
  4) 运行 07/08/09/10/13-18 时，仅当 ${CACHE_06_TABLES} 缺失才补跑 06
  5) 运行 10 时，若 ${CACHE_09_GSVA_STATE} 缺失，会先补跑 09
  6) 运行 17 时，若 ${CACHE_15_STATE} 或 ${CACHE_16_STATE} 缺失，会按 13/14/15/16 补齐整合缓存
  7) 可选环境变量 RNASEQ_SPECIES=mouse|human（默认 mouse）
  8) 每个实际执行模块开始时清理 result/<module> 和该模块缓存；可用 RNASEQ_CLEAN_MODULE_RESULT=false 或 RNASEQ_CLEAN_MODULE_CACHE=false 关闭

Modules:
$(print_module_table)
EOF
}

ensure_rscript() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "[ERROR] Rscript not found in PATH" >&2
    exit 1
  fi
}

validate_module_id() {
  local module_id="$1"
  [[ "${module_id}" =~ ^[0-9]+$ ]] || return 1
  (( module_id >= 1 && module_id <= TOTAL_MODULES ))
}

# 去重并保持顺序
# shellcheck disable=SC2120
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

parse_module_tokens() {
  local tokens=("$@")
  local parsed=()
  local token part raw

  for token in "${tokens[@]}"; do
    token="${token//[[:space:]]/}"
    [[ -z "${token}" ]] && continue
    IFS=',' read -r -a raw <<< "${token}"
    for part in "${raw[@]}"; do
      part="${part//[[:space:]]/}"
      [[ -z "${part}" ]] && continue
      if [[ ! "${part}" =~ ^0*[0-9]+$ ]]; then
        echo "[ERROR] Invalid module token: ${part}" >&2
        exit 1
      fi
      # 处理前导 0
      part="$((10#${part}))"
      if ! validate_module_id "${part}"; then
        echo "[ERROR] Unsupported module id: ${part}" >&2
        exit 1
      fi
      parsed+=("${part}")
    done
  done

  if [[ ${#parsed[@]} -eq 0 ]]; then
    echo "[ERROR] No valid module ids parsed." >&2
    exit 1
  fi

  printf '%s\n' "${parsed[@]}"
}

is_main_module() {
  local m="$1"
  (( m >= 2 && m <= 10 ))
}

is_post_module() {
  local m="$1"
  (( m >= 11 && m <= TOTAL_MODULES ))
}

is_knit_module() {
  local m="$1"
  (( m >= 2 && m <= TOTAL_MODULES ))
}

uses_post_context() {
  local m="$1"
  [[ "${m}" == "9" || "${m}" == "10" || "${m}" == "13" || "${m}" == "14" || "${m}" == "15" || "${m}" == "16" || "${m}" == "17" ]]
}

requires_enrichment_cache() {
  local m="$1"
  [[ "${m}" == "7" || "${m}" == "8" || "${m}" == "18" ]] || uses_post_context "${m}"
}

add_once() {
  local module_id="$1"
  local current
  for current in "${RUN_ORDER[@]}"; do
    if [[ "${current}" == "${module_id}" ]]; then
      return 0
    fi
  done
  RUN_ORDER+=("${module_id}")
}

add_with_setup_dep() {
  local module_id="$1"
  if [[ "${module_id}" != "1" ]]; then
    add_once 1
  fi
  add_once "${module_id}"
}

add_module_with_cache_deps() {
  local module_id="$1"
  if requires_enrichment_cache "${module_id}" && [[ ! -f "${CACHE_06_TABLES}" ]]; then
    add_with_setup_dep 6
  fi
  add_with_setup_dep "${module_id}"
}

resolve_run_order() {
  local requested=("$@")
  local m

  RUN_ORDER=()

  for m in "${requested[@]}"; do
    if requires_enrichment_cache "${m}" && [[ ! -f "${CACHE_06_TABLES}" ]]; then
      add_with_setup_dep 6
    fi

    if [[ "${m}" == "10" && ! -f "${CACHE_09_GSVA_STATE}" ]]; then
      add_module_with_cache_deps 9
    fi

    if [[ "${m}" == "17" ]]; then
      if [[ ! -f "${CACHE_15_STATE}" ]]; then
        add_module_with_cache_deps 13
        add_module_with_cache_deps 14
        add_module_with_cache_deps 15
      fi
      if [[ ! -f "${CACHE_16_STATE}" ]]; then
        add_module_with_cache_deps 16
      fi
    fi

    add_with_setup_dep "${m}"
  done
}

module_abs_path() {
  local module_id="$1"
  printf '%s/%s\n' "${RMD_DIR}" "${MODULE_FILES[${module_id}]}"
}

module_result_step() {
  local module_id="$1"
  printf '%s\n' "${MODULE_RESULT_DIRS[${module_id}]}"
}

ensure_module_file_exists() {
  local module_id="$1"
  local f
  f="$(module_abs_path "${module_id}")"
  if [[ ! -f "${f}" ]]; then
    echo "[ERROR] Missing module file: ${f}" >&2
    exit 1
  fi
}

run_module() {
  local module_id="$1"
  local input_file setup_file active_result_step

  input_file="$(module_abs_path "${module_id}")"
  active_result_step="$(module_result_step "${module_id}")"
  ensure_module_file_exists "${module_id}"

  echo "=============================================================="
  echo "[INFO] Module ${module_id}: ${MODULE_TITLES[${module_id}]}"
  echo "[INFO] File: ${input_file}"
  echo "[INFO] Result cleanup target: result/${active_result_step}"
  echo "[INFO] Start: $(date '+%Y-%m-%d %H:%M:%S')"

  if [[ "${module_id}" == "1" ]]; then
    RNASEQ_ACTIVE_RESULT_STEP="${active_result_step}" Rscript -e "root_dir <- '${ROOT_DIR}'; rmd_dir <- '${RMD_DIR}'; Sys.setenv(RNASEQ_ROOT_DIR = root_dir, RNASEQ_RMD_DIR = rmd_dir); setup_file <- '${input_file}'; setwd(root_dir); source(setup_file, local = globalenv())"
  else
    setup_file="$(module_abs_path 1)"

    RNASEQ_ACTIVE_RESULT_STEP="${active_result_step}" Rscript -e "module_id <- as.integer('${module_id}'); root_dir <- '${ROOT_DIR}'; rmd_dir <- '${RMD_DIR}'; Sys.setenv(RNASEQ_ROOT_DIR = root_dir, RNASEQ_RMD_DIR = rmd_dir); input_file <- '${input_file}'; setup_file <- '${setup_file}'; setwd(root_dir); if (!file.exists(setup_file)) stop('Missing setup file: ', setup_file); source(setup_file, local = globalenv()); if (module_id %in% c(9L, 10L, 13L, 14L, 15L, 16L, 17L)) load_post_context(globalenv()); knitr::opts_knit\$set(root.dir = root_dir); temp_output <- tempfile(sprintf('module_%02d_', module_id), fileext = '.md'); on.exit(unlink(temp_output, force = TRUE), add = TRUE); invisible(knitr::knit(input_file, output = temp_output, envir = globalenv(), quiet = TRUE))"
  fi

  echo "[INFO] Done module ${module_id} at $(date '+%Y-%m-%d %H:%M:%S')"
}

MODULE_ARGS=()
RUN_ALL=false

if [[ $# -eq 0 ]]; then
  RUN_ALL=true
else
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --help|-h)
        usage
        exit 0
        ;;
      --list)
        print_module_table
        exit 0
        ;;
      --all)
        RUN_ALL=true
        shift
        ;;
      --modules)
        shift
        if [[ $# -eq 0 ]]; then
          echo "[ERROR] --modules requires module ids" >&2
          exit 1
        fi
        while [[ $# -gt 0 && "$1" != --* ]]; do
          MODULE_ARGS+=("$1")
          shift
        done
        ;;
      *)
        MODULE_ARGS+=("$1")
        shift
        ;;
    esac
  done
fi

ensure_rscript

if [[ "${RUN_ALL}" == "true" && ${#MODULE_ARGS[@]} -eq 0 ]]; then
  REQUESTED_MODULES=("${DEFAULT_MODULES[@]}")
elif [[ ${#MODULE_ARGS[@]} -eq 0 ]]; then
  REQUESTED_MODULES=("${DEFAULT_MODULES[@]}")
else
  REQUESTED_MODULES=($(parse_module_tokens "${MODULE_ARGS[@]}"))
  REQUESTED_MODULES=($(dedupe_preserve_order "${REQUESTED_MODULES[@]}"))
fi

resolve_run_order "${REQUESTED_MODULES[@]}"

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
    print_failure_footer "${exit_code}"
  fi
}
trap on_exit EXIT

echo "[INFO] Log file: ${RUN_LOG_FILE}"
echo "[INFO] Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] Project root: ${ROOT_DIR}"
echo "[INFO] Rmd template dir: ${RMD_DIR}"
echo "[INFO] Requested modules (dedup): ${REQUESTED_MODULES[*]}"
echo "[INFO] Final run order: ${RUN_ORDER[*]}"

for module_id in "${RUN_ORDER[@]}"; do
  run_module "${module_id}"
done

echo "[INFO] End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "[INFO] Completed modules: ${RUN_ORDER[*]}"
