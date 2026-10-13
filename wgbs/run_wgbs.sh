#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# WGBS downstream Rmd pipeline for sperm samples
#
# 用法:
#   bash run_wgbs.sh
#     运行已启用的下游模块。
#
#   bash run_wgbs.sh --steps 3,5-8
#     只运行指定步骤；默认自动前置 01_setup。
#
#   bash run_wgbs.sh --list
#     列出所有可运行模块。
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RMD_DIR="${ROOT_DIR}/rmd"

# 纯绘图必须在 Rmd 调度、依赖补跑及清理之前分支。
if [[ " ${*} " == *" --plots-only "* ]]; then
  plot_args=(--root "${ROOT_DIR}" --assay wgbs)
  while (( $# )); do
    case "$1" in
      --plots-only) shift ;;
      --steps) plot_args+=(--steps "$2"); shift 2 ;;
      --check) plot_args+=(--check); shift ;;
      *) echo "[ERROR] --plots-only 不支持参数: $1" >&2; exit 2 ;;
    esac
  done
  source "${ROOT_DIR}/scripts/load_config.sh"
  exec python3 "${ROOT_DIR}/scripts/replot.py" "${plot_args[@]}"
fi

READS_INFO="${ROOT_DIR}/01rawdata/reads_info.txt"
REPORT_SCRIPT="${ROOT_DIR}/rmd/generate_result_report.R"

source "${ROOT_DIR}/scripts/load_config.sh"
READS_INFO="$WGBS_SAMPLES"
# 内部运行状态，默认 0（正常执行）；传入 --check 后设为 1，只检查配置、样本表和 Rmd，不渲染报告。
WGBS_CHECK_ONLY=0


declare -A STEPS=(
  [1]="01_setup"
  [2]="02_qc"
  [3]="03_dmr"
  [4]="04_dmroverview"
  [5]="05_enrichment"
  [6]="06_genomecontext"
  [7]="07_regulatory"
  [8]="08_gene_methylation"
  [9]="09_browser"
  [10]="10_sperm"
)
ACTIVE_STEPS=(1 2 3 4 5 6 7 8 9 10)

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
  local header_encoded="${header_line//$'\t'/$'\x1f'}"
  IFS=$'\x1f' read -r -a header_cols <<< "${header_encoded}"

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

  [[ -n "${col_idx[${raw_fq1_col}]:-}" ]] || reads_info_die "reads_info.txt 缺少 raw_fq1 列"
  [[ -n "${col_idx[${raw_fq2_col}]:-}" ]] || reads_info_die "reads_info.txt 缺少 raw_fq2 列"

  declare -ag READS_INFO_SAMPLES=()
  declare -ag READS_INFO_GROUPS=()
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
    local line_encoded="${line//$'\t'/$'\x1f'}"
    IFS=$'\x1f' read -r -a fields <<< "${line_encoded}"

    local sample="${fields[${col_idx[sample]}]:-}"
    local group="${fields[${col_idx[group]}]:-}"
    local raw_fq1="${fields[${col_idx[${raw_fq1_col}]}]:-}"
    local raw_fq2="${fields[${col_idx[${raw_fq2_col}]}]:-}"

    sample="$(trim_carriage_return "${sample}")"
    group="$(trim_carriage_return "${group}")"
    raw_fq1="$(trim_carriage_return "${raw_fq1}")"
    raw_fq2="$(trim_carriage_return "${raw_fq2}")"

    [[ -n "${sample}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 sample 为空"
    [[ -n "${group}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 group 为空: ${sample}"
    [[ -n "${raw_fq1}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 raw_fq1 为空: ${sample}"
    [[ -n "${raw_fq2}" ]] || reads_info_die "reads_info.txt 第 ${line_no} 行 raw_fq2 为空: ${sample}"

    if [[ -n "${READS_INFO_SAMPLE_TO_GROUP[${sample}]:-}" ]]; then
      reads_info_die "reads_info.txt 存在重复 sample: ${sample}"
    fi

    local resolved_fq1 resolved_fq2
    resolved_fq1="$(reads_info_resolve_path "${root_dir}" "${raw_fq1}")"
    resolved_fq2="$(reads_info_resolve_path "${root_dir}" "${raw_fq2}")"

    if [[ "${require_fastq_exists}" == "1" ]]; then
      [[ -f "${resolved_fq1}" ]] || reads_info_die "未找到 FASTQ 文件: ${resolved_fq1}"
      [[ -f "${resolved_fq2}" ]] || reads_info_die "未找到 FASTQ 文件: ${resolved_fq2}"
    fi

    READS_INFO_SAMPLES+=("${sample}")
    READS_INFO_GROUPS+=("${group}")
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

  if [[ -n "${WGBS_REFERENCE_GROUP}" ]]; then
    reads_info_contains "${WGBS_REFERENCE_GROUP}" "${READS_INFO_GROUP_LEVELS[@]}" || reads_info_die "reference group 不在 reads_info 中: ${WGBS_REFERENCE_GROUP}"
    reads_info_contains "${WGBS_COMPARISON_GROUP}" "${READS_INFO_GROUP_LEVELS[@]}" || reads_info_die "comparison group 不在 reads_info 中: ${WGBS_COMPARISON_GROUP}"
    [[ "${WGBS_REFERENCE_GROUP}" != "${WGBS_COMPARISON_GROUP}" ]] || reads_info_die "reference/comparison group 不能相同"
    READS_INFO_GROUP1_NAME="${WGBS_REFERENCE_GROUP}"
    READS_INFO_GROUP2_NAME="${WGBS_COMPARISON_GROUP}"
  else
    reads_info_die "必须显式设置 reference 和 comparison"
  fi
  READS_INFO_GROUP1_SIZE="${READS_INFO_GROUP_COUNTS[${READS_INFO_GROUP1_NAME}]}"
  READS_INFO_GROUP2_SIZE="${READS_INFO_GROUP_COUNTS[${READS_INFO_GROUP2_NAME}]}"

  (( READS_INFO_GROUP1_SIZE > 0 )) || reads_info_die "第一个分组没有样本: ${READS_INFO_GROUP1_NAME}"
  (( READS_INFO_GROUP2_SIZE > 0 )) || reads_info_die "第二个分组没有样本: ${READS_INFO_GROUP2_NAME}"
}

usage() {
  cat <<EOF
用法:
  bash run_wgbs.sh                 # 运行已启用的下游模块
  bash run_wgbs.sh --list          # 列出模块
  bash run_wgbs.sh --plots-only    # 从冻结数值重绘，不补跑统计
  bash run_wgbs.sh --steps 3,5-8   # 只运行指定步骤，默认自动前置 01_setup
  bash run_wgbs.sh --check --steps 8,10 # 只检查配置、样本表和 Rmd

说明:
  run_wgbs.sh 只负责下游 Rmd 调度，不运行 fastp/Bismark/MethylDackel。
  --steps 为选步入口。范围和逗号可以组合，例如 --steps 3,5-8。
  WGBS_SIGNIFICANCE_MODE=pvalue|padj 控制 富集筛选口径，默认 pvalue。
  WGBS_REFERENCE_GROUP/WGBS_COMPARISON_GROUP 必须显式固定 group2-group1 的方向。
  上游测序处理请运行 run_upstream.sh。
EOF
}

list_modules() {
  local i
  for i in "${ACTIVE_STEPS[@]}"; do
    printf '%02d  %s\n' "${i}" "${STEPS[$i]}"
  done
}

split_module_tokens() {
  local token="$1"
  token="${token//,/ }"
  local -a parts=()
  read -r -a parts <<< "${token}"
  printf '%s\n' "${parts[@]}"
}

append_module_token() {
  local token="$1"
  local module

  [[ -n "${token}" ]] || return 0
  if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    local start=$((10#${BASH_REMATCH[1]}))
    local end=$((10#${BASH_REMATCH[2]}))
    if (( start < 1 || end > 10 || start > end )); then
      reads_info_die "无效模块范围: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    for module in $(seq "${start}" "${end}"); do
      [[ -n "${STEPS[${module}]:-}" ]] || reads_info_die "模块已移除或不存在: ${module}"
      REQUESTED_MODULES+=("${module}")
    done
  elif [[ "${token}" =~ ^[0-9]+$ ]]; then
    local token_num=$((10#${token}))
    if (( token_num < 1 || token_num > 10 )); then
      reads_info_die "无效模块: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    [[ -n "${STEPS[${token_num}]:-}" ]] || reads_info_die "模块已移除或不存在: ${token_num}"
    REQUESTED_MODULES+=("${token_num}")
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
  boolean_or_die WGBS_RMARKDOWN_QUIET
  case "${WGBS_SIGNIFICANCE_MODE}" in
    pvalue|padj) ;;
    *) reads_info_die "WGBS_SIGNIFICANCE_MODE 只能是 pvalue/padj，当前为: ${WGBS_SIGNIFICANCE_MODE}" ;;
  esac
  if [[ -n "${WGBS_REFERENCE_GROUP}" || -n "${WGBS_COMPARISON_GROUP}" ]]; then
    [[ -n "${WGBS_REFERENCE_GROUP}" && -n "${WGBS_COMPARISON_GROUP}" ]] || reads_info_die "WGBS_REFERENCE_GROUP 与 WGBS_COMPARISON_GROUP 必须同时设置"
  fi
}

select_all_modules() {
  local i

  REQUESTED_MODULES=()
  for i in "${ACTIVE_STEPS[@]}"; do
    if [[ "$i" == 11 ]]; then [[ "$WGBS_ENABLE_SPERM" == 1 ]] || continue; fi
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
        --report-only)
          python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS
          python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS --verify downstream
          RESULT_REPORT_ROOT="$ROOT_DIR" Rscript --vanilla "$ROOT_DIR/rmd/generate_result_report.R"
          exit 0
          ;;
        --help|-h)
          usage
          exit 0
          ;;
        --list)
          list_modules
          exit 0
          ;;
        --check|--dry-run)
          WGBS_CHECK_ONLY=1
          shift
          ;;
        --steps)
          shift
          (( $# > 0 )) || reads_info_die "--steps 后需要指定步骤，例如 --steps 3,5-8"
          [[ ${#REQUESTED_MODULES[@]} -eq 0 ]] || reads_info_die "只能指定一次 --steps"
          local token
          while IFS= read -r token; do
            append_module_token "${token}"
          done < <(split_module_tokens "$1")
          shift
          ;;
        --*)
          reads_info_die "未知参数: $1；选步骤请使用 --steps 3,5-8"
          ;;
        *)
          reads_info_die "未知位置参数: $1；选步骤请使用 --steps 3,5-8"
          ;;
      esac
    done
  fi

  if (( ${#REQUESTED_MODULES[@]} == 0 )) && [[ "${WGBS_CHECK_ONLY}" == "1" ]]; then
    select_all_modules
  fi
  (( ${#REQUESTED_MODULES[@]} > 0 )) || reads_info_die "未指定要运行的模块"

  local -A seen=()
  local -a sorted=()
  local module
  for module in "${REQUESTED_MODULES[@]}"; do
    seen["${module}"]=1
  done

  for module in "${ACTIVE_STEPS[@]}"; do
    if [[ -n "${seen[${module}]:-}" ]]; then
      sorted+=("${module}")
    fi
  done

  if [[ -z "${seen[1]:-}" ]]; then
    sorted=("1" "${sorted[@]}")
  fi

  REQUESTED_MODULES=("${sorted[@]}")
}

render_modules() {
  local -a modules=("$@")
  local -a step_names=()
  local module

  for module in "${modules[@]}"; do
    step_names+=("${STEPS[${module}]}")
  done

  local step_name
  for step_name in "${step_names[@]}"; do
    [[ -f "${RMD_DIR}/${step_name}.Rmd" ]] || reads_info_die "Rmd 文件不存在: ${RMD_DIR}/${step_name}.Rmd"
  done

  mkdir -p "${WGBS_LOG_DIR}" "$ROOT_DIR/result"

  local log_file="${WGBS_LOG_DIR}/run_wgbs_$(date '+%Y%m%d_%H%M%S')_$$.log"
  local r_script
  r_script="$(mktemp "${TMPDIR:-/tmp}/run_wgbs.XXXXXX.R")"
  trap 'rm -f -- "${r_script:-}"' EXIT
  trap 'exit 130' INT
  trap 'exit 143' TERM

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
status_file <- file.path(root_dir, "result", "execution.tsv")
record_status <- function(module, status) {
  row <- data.frame(time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), module = module, status = status)
  write.table(row, status_file, sep = "\t", row.names = FALSE, col.names = !file.exists(status_file), append = file.exists(status_file), quote = FALSE)
}
if(Sys.getenv("WGBS_ENABLE_SPERM") != "1") record_status("10_sperm", "SKIPPED_DISABLED")
for (step_name in steps) {
  rmd_file <- file.path(rmd_dir, paste0(step_name, ".Rmd"))
  if (!file.exists(rmd_file)) {
    stop("Rmd 文件不存在: ", rmd_file)
  }

  if (!identical(step_name, "01_setup")) {
    step_plot_dir <- file.path(root_dir, "result", step_name)
    step_plot_dirs <- list.dirs(step_plot_dir, recursive = TRUE, full.names = TRUE)
    step_plot_dirs <- step_plot_dirs[basename(step_plot_dirs) == "plots"]
    root_result_dir <- normalizePath(file.path(root_dir, "result"), winslash = "/", mustWork = FALSE)
    plot_dir_norm <- normalizePath(step_plot_dir, winslash = "/", mustWork = FALSE)
    if (!startsWith(plot_dir_norm, paste0(root_result_dir, "/"))) {
      stop("Refusing to clean plot directory outside result/: ", step_plot_dir)
    }
    dir.create(step_plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cat("================================================================\n")
  cat("▶ ", step_name, "\n", sep = "")
  cat("  开始: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("================================================================\n")

  record_status(step_name, "STARTED")
  tryCatch(
    {
      rmarkdown::render(
        input = rmd_file,
        knit_root_dir = root_dir,
        envir = .GlobalEnv,
        quiet = quiet
      )
      record_status(step_name, "COMPLETED")
      cat("✅ ", step_name, " 完成: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
    },
    error = function(e) {
      record_status(step_name, "FAILED")
      cat("❌ ", step_name, " 失败: ", conditionMessage(e), "\n", sep = "")
      quit(status = 1, save = "no")
    }
  )
}
RSCRIPT

  local steps_joined
  steps_joined="${step_names[0]}"
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
    WGBS_SIGNIFICANCE_MODE="${WGBS_SIGNIFICANCE_MODE}" \
    WGBS_REFERENCE_GROUP="${WGBS_REFERENCE_GROUP}" \
    WGBS_COMPARISON_GROUP="${WGBS_COMPARISON_GROUP}" \
    Rscript "${r_script}" > "${log_file}" 2>&1; then
    rm -f "${r_script}"
    r_script=""
    trap - EXIT INT TERM
    echo "下游模块运行完成"
  else
    rm -f "${r_script}"
    r_script=""
    trap - EXIT INT TERM
    echo "下游模块运行失败，日志: ${log_file}" >&2
    tail -n 80 "${log_file}" >&2 || true
    return 1
  fi
}

generate_result_report() {
  [[ -f "${REPORT_SCRIPT}" ]] || reads_info_die "结果报告脚本不存在: ${REPORT_SCRIPT}"
  echo "刷新离线结果报告: ${ROOT_DIR}/result/index.html"
  RESULT_REPORT_ROOT="${ROOT_DIR}" Rscript --vanilla "${REPORT_SCRIPT}"
}

normalize_result_pngs() {
  WGBS_PROJECT_ROOT="${ROOT_DIR}" Rscript --vanilla - <<'RSCRIPT'
root_dir <- Sys.getenv("WGBS_PROJECT_ROOT")
source(normalizePath(file.path(root_dir, "rmd", "figure_contract.R"), mustWork = TRUE))
wgbs_normalize_png_widths(file.path(root_dir, "result"))
RSCRIPT
}

check_selection() {
  local module step_name
  local -a step_names=()
  for module in "$@"; do
    step_name="${STEPS[${module}]}"
    [[ -f "${RMD_DIR}/${step_name}.Rmd" ]] || reads_info_die "Rmd 文件不存在: ${RMD_DIR}/${step_name}.Rmd"
    step_names+=("${step_name}")
  done
  echo "检查通过"
  echo "步骤: $* (${step_names[*]})"
  echo "分组: ${READS_INFO_GROUP1_NAME} -> ${READS_INFO_GROUP2_NAME}"
  echo "样本数: ${#READS_INFO_SAMPLES[@]}"
  echo "显著性口径: ${WGBS_SIGNIFICANCE_MODE}"
}

main() {
  parse_args "$@"
  validate_config

  command -v Rscript >/dev/null 2>&1 || reads_info_die "缺少命令: Rscript"
  [[ -d "${RMD_DIR}" ]] || reads_info_die "Rmd 目录不存在: ${RMD_DIR}"

  python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS
  load_reads_info "${READS_INFO}" "${ROOT_DIR}" 0
  for module in "${REQUESTED_MODULES[@]}"; do
    if [[ "$module" == 10 ]]; then [[ "$WGBS_ENABLE_SPERM" == 1 ]] || reads_info_die "精子专题已禁用"; fi
  done
  if [[ "${WGBS_CHECK_ONLY}" == "1" ]]; then
    check_selection "${REQUESTED_MODULES[@]}"
    return 0
  fi
  if [[ -f "$ROOT_DIR/04methykit/CpG.RData" ]]; then
    python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS --verify upstream
  fi
  python3 "$ROOT_DIR/scripts/input_audit.py" --root "$ROOT_DIR" --prefix WGBS --record downstream
  command -v pandoc >/dev/null 2>&1 || reads_info_die "缺少 pandoc：请在 PATH 显式配置工具"
  render_modules "${REQUESTED_MODULES[@]}"
  normalize_result_pngs
  generate_result_report
}

main "$@"
