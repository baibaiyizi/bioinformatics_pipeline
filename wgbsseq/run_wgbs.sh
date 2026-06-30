#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# WGBS downstream Rmd pipeline for sperm samples
#
# 用法:
#   bash run_wgbs.sh
#     运行下游全部模块 1-11。
#
#   bash run_wgbs.sh --modules 3,5-8
#     只运行指定模块；默认自动前置 1_setup。这是唯一的选模块入口。
#
#   bash run_wgbs.sh --list
#     列出所有可运行模块。
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RMD_DIR="${ROOT_DIR}/rmd"
READS_INFO="${ROOT_DIR}/01rawdata/reads_info.txt"

# 可选参数集中配置。运行前可用同名环境变量覆盖，例如: WGBS_AUTO_SETUP=0 bash run_wgbs.sh --modules 3
WGBS_LOG_DIR="${WGBS_LOG_DIR:-${ROOT_DIR}/log}" # 下游运行日志目录。
WGBS_PANDOC_FALLBACK_DIR="${WGBS_PANDOC_FALLBACK_DIR:-${HOME}/miniconda3/bin}" # pandoc 不在 PATH 时使用该目录。
WGBS_AUTO_SETUP="${WGBS_AUTO_SETUP:-1}" # 选择部分模块时，是否自动前置 1_setup。
WGBS_RMARKDOWN_QUIET="${WGBS_RMARKDOWN_QUIET:-0}" # 是否静默 rmarkdown::render 输出；默认保留详细日志。

if ! command -v pandoc >/dev/null 2>&1 && [[ -x "${WGBS_PANDOC_FALLBACK_DIR}/pandoc" ]]; then
  export PATH="${WGBS_PANDOC_FALLBACK_DIR}:${PATH}"
  export RSTUDIO_PANDOC="${RSTUDIO_PANDOC:-${WGBS_PANDOC_FALLBACK_DIR}}"
fi

STEPS=(
  "1_setup"
  "2_qc"
  "3_dmr"
  "4_dmroverview"
  "5_enrichment"
  "6_genomecontext"
  "7_regulatory"
  "8_multiomics"
  "9_degfocused"
  "10_browser"
  "11_sperm"
)

trim_carriage_return() {
  local value="$1"
  value="${value%$'\r'}"
  printf '%s' "${value}"
}

reads_info_die() {
  echo "[ERROR] $*" >&2
  exit 1
}

reads_info_contains() {
  local needle="$1"
  shift || true
  local item
  for item in "$@"; do
    if [[ "${item}" == "${needle}" ]]; then
      return 0
    fi
  done
  return 1
}

reads_info_resolve_path() {
  local root_dir="$1"
  local path_value="$2"
  if [[ "${path_value}" = /* ]]; then
    printf '%s' "${path_value}"
  else
    printf '%s/%s' "${root_dir}" "${path_value}"
  fi
}

load_reads_info() {
  local reads_info_file="$1"
  local root_dir="$2"
  local require_fastq_exists="${3:-0}"

  [[ -f "${reads_info_file}" ]] || reads_info_die "未找到 reads_info 文件: ${reads_info_file}"

  local header_line
  IFS= read -r header_line < "${reads_info_file}" || true
  [[ -n "${header_line}" ]] || reads_info_die "reads_info.txt 为空: ${reads_info_file}"

  local -a header_cols=()
  IFS=$'\t' read -r -a header_cols <<< "${header_line}"

  declare -A col_idx=()
  local col i
  for i in "${!header_cols[@]}"; do
    col="$(trim_carriage_return "${header_cols[$i]}")"
    header_cols[$i]="${col}"
    [[ -n "${col}" ]] && col_idx["${col}"]="${i}"
  done

  [[ -n "${col_idx[sample]:-}" ]] || reads_info_die "reads_info.txt 缺少 sample 列"
  [[ -n "${col_idx[group]:-}" ]] || reads_info_die "reads_info.txt 缺少 group 列"

  local raw_fq1_col="raw_fq1"
  local raw_fq2_col="raw_fq2"
  if [[ -z "${col_idx[raw_fq1]:-}" && -n "${col_idx[fq1]:-}" ]]; then
    raw_fq1_col="fq1"
  fi
  if [[ -z "${col_idx[raw_fq2]:-}" && -n "${col_idx[fq2]:-}" ]]; then
    raw_fq2_col="fq2"
  fi

  [[ -n "${col_idx[${raw_fq1_col}]:-}" ]] || reads_info_die "reads_info.txt 缺少 raw_fq1/fq1 列"
  [[ -n "${col_idx[${raw_fq2_col}]:-}" ]] || reads_info_die "reads_info.txt 缺少 raw_fq2/fq2 列"

  declare -ag READS_INFO_SAMPLES=()
  declare -ag READS_INFO_GROUPS=()
  declare -ag READS_INFO_GROUP_STAGES=()
  declare -ag READS_INFO_STAGES=()
  declare -ag READS_INFO_RAW_FQ1=()
  declare -ag READS_INFO_RAW_FQ2=()
  declare -ag READS_INFO_RESOLVED_FQ1=()
  declare -ag READS_INFO_RESOLVED_FQ2=()
  declare -ag READS_INFO_GROUP_LEVELS=()
  declare -Ag READS_INFO_SAMPLE_TO_GROUP=()
  declare -Ag READS_INFO_GROUP_COUNTS=()

  local line_no=1
  local line
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line_no=$((line_no + 1))
    [[ -n "${line//[$'\t\r ']/}" ]] || continue

    local -a fields=()
    IFS=$'\t' read -r -a fields <<< "${line}"

    local sample="${fields[${col_idx[sample]}]:-}"
    local group="${fields[${col_idx[group]}]:-}"
    local raw_fq1="${fields[${col_idx[${raw_fq1_col}]}]:-}"
    local raw_fq2="${fields[${col_idx[${raw_fq2_col}]}]:-}"
    local group_stage=""
    local stage=""

    sample="$(trim_carriage_return "${sample}")"
    group="$(trim_carriage_return "${group}")"
    raw_fq1="$(trim_carriage_return "${raw_fq1}")"
    raw_fq2="$(trim_carriage_return "${raw_fq2}")"

    if [[ -n "${col_idx[group_stage]:-}" ]]; then
      group_stage="$(trim_carriage_return "${fields[${col_idx[group_stage]}]:-}")"
    fi
    if [[ -n "${col_idx[stage]:-}" ]]; then
      stage="$(trim_carriage_return "${fields[${col_idx[stage]}]:-}")"
    fi

    [[ -n "${sample}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 sample 为空"
    [[ -n "${group}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 group 为空: ${sample}"
    [[ -n "${raw_fq1}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 raw_fq1 为空: ${sample}"
    [[ -n "${raw_fq2}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 raw_fq2 为空: ${sample}"

    if [[ -n "${READS_INFO_SAMPLE_TO_GROUP[${sample}]:-}" ]]; then
      reads_info_die "reads_info.txt 存在重复 sample: ${sample}"
    fi

    if [[ -z "${stage}" ]]; then
      if [[ -n "${group_stage}" && "${group_stage}" =~ [Ss]perm$ ]]; then
        stage="sperm"
      elif [[ -n "${group_stage}" ]]; then
        stage="${group_stage}"
      else
        stage="sperm"
      fi
    fi
    [[ -n "${group_stage}" ]] || group_stage="${group}${stage}"

    local resolved_fq1 resolved_fq2
    resolved_fq1="$(reads_info_resolve_path "${root_dir}" "${raw_fq1}")"
    resolved_fq2="$(reads_info_resolve_path "${root_dir}" "${raw_fq2}")"

    if [[ "${require_fastq_exists}" == "1" ]]; then
      [[ -f "${resolved_fq1}" ]] || reads_info_die "未找到 FASTQ 文件: ${resolved_fq1}"
      [[ -f "${resolved_fq2}" ]] || reads_info_die "未找到 FASTQ 文件: ${resolved_fq2}"
    fi

    READS_INFO_SAMPLES+=("${sample}")
    READS_INFO_GROUPS+=("${group}")
    READS_INFO_GROUP_STAGES+=("${group_stage}")
    READS_INFO_STAGES+=("${stage}")
    READS_INFO_RAW_FQ1+=("${raw_fq1}")
    READS_INFO_RAW_FQ2+=("${raw_fq2}")
    READS_INFO_RESOLVED_FQ1+=("${resolved_fq1}")
    READS_INFO_RESOLVED_FQ2+=("${resolved_fq2}")
    READS_INFO_SAMPLE_TO_GROUP["${sample}"]="${group}"

    if ! reads_info_contains "${group}" "${READS_INFO_GROUP_LEVELS[@]}"; then
      READS_INFO_GROUP_LEVELS+=("${group}")
      READS_INFO_GROUP_COUNTS["${group}"]=0
    fi
    READS_INFO_GROUP_COUNTS["${group}"]=$((READS_INFO_GROUP_COUNTS["${group}"] + 1))
  done < <(tail -n +2 "${reads_info_file}")

  (( ${#READS_INFO_SAMPLES[@]} > 0 )) || reads_info_die "reads_info.txt 中没有可用样本"
  (( ${#READS_INFO_GROUP_LEVELS[@]} == 2 )) || reads_info_die "reads_info.txt 必须恰好包含两个分组，当前得到: ${READS_INFO_GROUP_LEVELS[*]}"

  READS_INFO_GROUP1_NAME="${READS_INFO_GROUP_LEVELS[0]}"
  READS_INFO_GROUP2_NAME="${READS_INFO_GROUP_LEVELS[1]}"
  READS_INFO_GROUP1_SIZE="${READS_INFO_GROUP_COUNTS[${READS_INFO_GROUP1_NAME}]}"
  READS_INFO_GROUP2_SIZE="${READS_INFO_GROUP_COUNTS[${READS_INFO_GROUP2_NAME}]}"

  (( READS_INFO_GROUP1_SIZE > 0 )) || reads_info_die "第一个分组没有样本: ${READS_INFO_GROUP1_NAME}"
  (( READS_INFO_GROUP2_SIZE > 0 )) || reads_info_die "第二个分组没有样本: ${READS_INFO_GROUP2_NAME}"
}

usage() {
  cat <<EOF
用法:
  bash run_wgbs.sh                 # 运行下游全部模块 1-11
  bash run_wgbs.sh --list          # 列出模块
  bash run_wgbs.sh --modules 3,5-8 # 只运行指定模块，默认自动前置 1_setup

说明:
  run_wgbs.sh 只负责下游 Rmd 调度，不运行 fastp/Bismark/MethylDackel。
  只保留 --modules 这一种选模块方式；范围和逗号可以组合，例如 --modules 3,5-8。
  上游测序处理请运行 sperm.sh。
EOF
}

list_modules() {
  local i
  for i in "${!STEPS[@]}"; do
    printf '%2d  %s\n' "$((i + 1))" "${STEPS[$i]}"
  done
}

split_module_tokens() {
  local token="$1"
  token="${token//,/ }"
  printf '%s\n' ${token}
}

append_module_token() {
  local token="$1"
  local module

  [[ -n "${token}" ]] || return 0
  if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    local start="${BASH_REMATCH[1]}"
    local end="${BASH_REMATCH[2]}"
    if (( start < 1 || end > ${#STEPS[@]} || start > end )); then
      reads_info_die "无效模块范围: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    for module in $(seq "${start}" "${end}"); do
      REQUESTED_MODULES+=("${module}")
    done
  elif [[ "${token}" =~ ^[0-9]+$ ]]; then
    if (( token < 1 || token > ${#STEPS[@]} )); then
      reads_info_die "无效模块: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    REQUESTED_MODULES+=("${token}")
  else
    reads_info_die "无法识别模块参数: ${token}"
  fi
}

boolean_or_die() {
  local name="$1"
  local value="${!name}"

  case "${value}" in
    0|1) ;;
    *) reads_info_die "${name} 只能是 0/1，当前为: ${value}" ;;
  esac
}

validate_config() {
  boolean_or_die WGBS_AUTO_SETUP
  boolean_or_die WGBS_RMARKDOWN_QUIET
}

select_all_modules() {
  local i

  REQUESTED_MODULES=()
  for i in $(seq 1 "${#STEPS[@]}"); do
    REQUESTED_MODULES+=("${i}")
  done
}

parse_args() {
  REQUESTED_MODULES=()

  if (( $# == 0 )); then
    select_all_modules
  else
    while (( $# > 0 )); do
      case "$1" in
        --help|-h)
          usage
          exit 0
          ;;
        --list)
          list_modules
          exit 0
          ;;
        --modules)
          shift
          (( $# > 0 )) || reads_info_die "--modules 后需要指定模块，例如 --modules 3,5-8"
          [[ ${#REQUESTED_MODULES[@]} -eq 0 ]] || reads_info_die "只能指定一次 --modules"
          local token
          while IFS= read -r token; do
            append_module_token "${token}"
          done < <(split_module_tokens "$1")
          shift
          ;;
        --*)
          reads_info_die "未知参数: $1；选模块只支持 --modules 3,5-8"
          ;;
        *)
          reads_info_die "未知位置参数: $1；选模块只支持 --modules 3,5-8"
          ;;
      esac
    done
  fi

  (( ${#REQUESTED_MODULES[@]} > 0 )) || reads_info_die "未指定要运行的模块"

  local -A seen=()
  local -a sorted=()
  local module
  for module in "${REQUESTED_MODULES[@]}"; do
    seen["${module}"]=1
  done

  for module in $(seq 1 "${#STEPS[@]}"); do
    if [[ -n "${seen[${module}]:-}" ]]; then
      sorted+=("${module}")
    fi
  done

  if [[ "${WGBS_AUTO_SETUP}" == "1" && -z "${seen[1]:-}" ]]; then
    sorted=("1" "${sorted[@]}")
  fi

  REQUESTED_MODULES=("${sorted[@]}")
}

render_modules() {
  local -a modules=("$@")
  local -a step_names=()
  local module

  for module in "${modules[@]}"; do
    step_names+=("${STEPS[$((module - 1))]}")
  done

  mkdir -p "${WGBS_LOG_DIR}"

  local log_file="${WGBS_LOG_DIR}/run_wgbs_$(date '+%Y%m%d_%H%M%S')_$$.log"
  local r_script
  r_script="$(mktemp "${TMPDIR:-/tmp}/run_wgbs.XXXXXX.R")"

  cat > "${r_script}" <<'RSCRIPT'
root_dir <- Sys.getenv("WGBS_ROOT_DIR")
rmd_dir <- Sys.getenv("WGBS_RMD_DIR")
step_env <- Sys.getenv("WGBS_RMD_STEPS")
quiet <- identical(Sys.getenv("WGBS_RMARKDOWN_QUIET"), "1")

if (!nzchar(root_dir) || !dir.exists(root_dir)) {
  stop("WGBS_ROOT_DIR 无效: ", root_dir)
}
if (!nzchar(rmd_dir) || !dir.exists(rmd_dir)) {
  stop("WGBS_RMD_DIR 无效: ", rmd_dir)
}
if (!nzchar(step_env)) {
  stop("WGBS_RMD_STEPS 为空")
}

steps <- strsplit(step_env, "::", fixed = TRUE)[[1]]
steps <- steps[nzchar(steps)]

suppressPackageStartupMessages(library(rmarkdown))
setwd(root_dir)

for (step_name in steps) {
  rmd_file <- file.path(rmd_dir, paste0(step_name, ".Rmd"))
  if (!file.exists(rmd_file)) {
    stop("Rmd 文件不存在: ", rmd_file)
  }

  cat("================================================================\n")
  cat("▶ ", step_name, "\n", sep = "")
  cat("  开始: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("================================================================\n")

  tryCatch(
    {
      rmarkdown::render(
        input = rmd_file,
        knit_root_dir = root_dir,
        envir = .GlobalEnv,
        quiet = quiet
      )
      cat("✅ ", step_name, " 完成: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
    },
    error = function(e) {
      cat("❌ ", step_name, " 失败: ", conditionMessage(e), "\n", sep = "")
      quit(status = 1, save = "no")
    }
  )
}
RSCRIPT

  local steps_joined
  steps_joined="${step_names[0]}"
  local step_name
  for step_name in "${step_names[@]:1}"; do
    steps_joined+="::${step_name}"
  done

  echo "下游模块: ${modules[*]} (${step_names[*]})"
  echo "工作目录: ${ROOT_DIR}"
  echo "reads_info: ${READS_INFO}"
  echo "分组顺序: ${READS_INFO_GROUP1_NAME} -> ${READS_INFO_GROUP2_NAME}"
  echo "样本数量: ${#READS_INFO_SAMPLES[@]}"
  echo "日志: ${log_file}"

  if WGBS_ROOT_DIR="${ROOT_DIR}" \
    WGBS_RMD_DIR="${RMD_DIR}" \
    WGBS_RMD_STEPS="${steps_joined}" \
    WGBS_RMARKDOWN_QUIET="${WGBS_RMARKDOWN_QUIET}" \
    Rscript "${r_script}" > "${log_file}" 2>&1; then
    rm -f "${r_script}"
    echo "下游模块运行完成"
  else
    rm -f "${r_script}"
    echo "下游模块运行失败，日志: ${log_file}" >&2
    tail -n 80 "${log_file}" >&2 || true
    return 1
  fi
}

main() {
  parse_args "$@"
  validate_config

  command -v Rscript >/dev/null 2>&1 || reads_info_die "缺少命令: Rscript"
  [[ -d "${RMD_DIR}" ]] || reads_info_die "Rmd 目录不存在: ${RMD_DIR}"

  load_reads_info "${READS_INFO}" "${ROOT_DIR}" 0
  render_modules "${REQUESTED_MODULES[@]}"
}

main "$@"
