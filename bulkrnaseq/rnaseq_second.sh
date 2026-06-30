#!/usr/bin/env bash

# Example usage:
#   bash rnaseq_second.sh --list
#   bash rnaseq_second.sh --dry-run --steps 4
#   nohup bash rnaseq_second.sh --all > 2.log 2>&1 &
# Note: use comma/space-separated step ids; ranges like 2-5 are not supported.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
FIRST_SCRIPT="${ROOT_DIR}/rnaseq_first.sh"
MERGE_DIR="${ROOT_DIR}/3.Merge_result"
DE_DIR="${ROOT_DIR}/4.DE_analysis"

SAMPLES_FILE="${DE_DIR}/samples.txt"
CONTRASTS_FILE="${DE_DIR}/contrasts.txt"
MATRIX_FILE="${MERGE_DIR}/genes.counts.matrix"
DE_RESULTS_FILE="${DE_DIR}/DE_results"

THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
FIG_PNG_DPI="${FIG_PNG_DPI:-300}"
RNASEQ_SECOND_SKIP_DE="${RNASEQ_SECOND_SKIP_DE:-false}"

TOTAL_STEPS=5
DEFAULT_STEPS=(1 2 3 4 5)

declare -A STEP_NAMES=(
    [1]="validate"
    [2]="merge_selected"
    [3]="contrast"
    [4]="deseq2"
    [5]="merge_de_results"
)

declare -A STEP_TITLES=(
    [1]="validate selected samples and inputs"
    [2]="rebuild selected-sample matrices"
    [3]="write DE contrast"
    [4]="run DESeq2"
    [5]="merge DESeq2 result files"
)

DE_RUN_DIR=""
DE_SINGLE_RESULT_FILE=""
declare -a REQUESTED_STEPS=()
declare -a RUN_ORDER=()
REQUESTED_STEPS_LABEL=""
RUN_ORDER_LABEL=""

log() {
    echo "[second] $*"
}

die() {
    echo "[second][error] $*" >&2
    exit 1
}

print_default_order() {
    printf '%s\n' "${DEFAULT_STEPS[*]}"
}

format_steps() {
    local step
    local -a out=()
    for step in "$@"; do
        out+=("${step}:${STEP_NAMES[${step}]}")
    done
    printf '%s\n' "${out[*]}"
}

print_step_table() {
    local step
    for step in $(seq 1 "${TOTAL_STEPS}"); do
        printf '  %d: %-16s %s\n' "${step}" "${STEP_NAMES[${step}]}" "${STEP_TITLES[${step}]}"
    done
}

usage() {
    cat <<EOF
Examples:
  bash ${SCRIPT_NAME} --list
  bash ${SCRIPT_NAME} --dry-run --steps 4
  nohup bash ${SCRIPT_NAME} --all > 2.log 2>&1 &

Usage:
  ${SCRIPT_NAME}
  ${SCRIPT_NAME} --all
  ${SCRIPT_NAME} --list
  ${SCRIPT_NAME} --steps 2,4
  ${SCRIPT_NAME} --steps 2 4
  ${SCRIPT_NAME} 2 4

Options:
  --all              Run the full second-stage workflow: $(print_default_order)
  --steps LIST       Run selected steps. LIST can be comma- or space-separated.
  --list             Show available steps.
  --dry-run          Print planned commands without running heavy tools.
  --force            Re-run completed outputs.
  --threads N        Override THREADS for this run.
  --skip-de          Skip DESeq2 and reuse existing ${DE_RESULTS_FILE}.
  --run-de           Run DESeq2 even if RNASEQ_SECOND_SKIP_DE=true in the environment.
  -h, --help         Show this help.

Note:
  Use comma/space-separated step ids; ranges like 2-5 are not supported.

Steps:
$(print_step_table)
EOF
}

validate_step_id() {
    local step="$1"
    [[ "${step}" =~ ^[0-9]+$ ]] || return 1
    (( step >= 1 && step <= TOTAL_STEPS ))
}

step_id_from_token() {
    local token="$1"
    token="${token,,}"
    token="${token//-/_}"

    if [[ "${token}" =~ ^0*[0-9]+$ ]]; then
        token="$((10#${token}))"
        validate_step_id "${token}" || return 1
        printf '%s\n' "${token}"
        return 0
    fi

    case "${token}" in
        validate) printf '1\n' ;;
        merge_selected|merge) printf '2\n' ;;
        contrast|contrasts) printf '3\n' ;;
        deseq2|de) printf '4\n' ;;
        merge_de_results|merge_results) printf '5\n' ;;
        *) return 1 ;;
    esac
}

dedupe_preserve_order() {
    local input=("$@")
    local -a out=()
    local seen=" "
    local x
    for x in "${input[@]}"; do
        [[ "${seen}" == *" ${x} "* ]] && continue
        out+=("${x}")
        seen+="${x} "
    done
    printf '%s\n' "${out[@]}"
}

parse_step_tokens() {
    local tokens=("$@")
    local -a parsed=()
    local token part step_id
    local -a raw

    for token in "${tokens[@]}"; do
        token="${token//[[:space:]]/}"
        [[ -z "${token}" ]] && continue
        IFS=',' read -r -a raw <<< "${token}"
        for part in "${raw[@]}"; do
            part="${part//[[:space:]]/}"
            [[ -z "${part}" ]] && continue
            step_id="$(step_id_from_token "${part}")" || die "不支持的步骤: ${part}。请使用 --list 查看可选步骤。"
            parsed+=("${step_id}")
        done
    done

    [[ ${#parsed[@]} -gt 0 ]] || die "没有解析到有效步骤。"
    dedupe_preserve_order "${parsed[@]}"
}

add_once() {
    local step="$1"
    local current
    for current in "${RUN_ORDER[@]}"; do
        [[ "${current}" == "${step}" ]] && return 0
    done
    RUN_ORDER+=("${step}")
}

step_scheduled() {
    local target="$1"
    local current
    for current in "${RUN_ORDER[@]}"; do
        [[ "${current}" == "${target}" ]] && return 0
    done
    return 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "未找到命令: $1"
}

is_true() {
    case "${1:-}" in
        1|true|TRUE|True|yes|YES|on|ON)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

show_selected_samples() {
    awk 'NF >= 2 && $1 !~ /^#/ {printf("  - %s\t%s\n", $1, $2)}' "${SAMPLES_FILE}"
}

get_group_order() {
    awk 'NF >= 2 && $1 !~ /^#/ {if (!seen[$1]++) print $1}' "${SAMPLES_FILE}"
}

read_selected_sample_names() {
    awk 'NF >= 2 && $1 !~ /^#/ {print $2}' "${SAMPLES_FILE}"
}

matrix_header_matches_samples() {
    local matrix_file="$1"
    [[ -f "${matrix_file}" ]] || return 1

    local matrix_samples=""
    local selected_samples=""
    matrix_samples="$(
        awk -F'\t' 'NR == 1 {
            for (i = 2; i <= NF; i++) {
                print $i
            }
        }' "${matrix_file}"
    )"
    selected_samples="$(read_selected_sample_names)"

    [[ -n "${matrix_samples}" && -n "${selected_samples}" ]] || return 1
    [[ "${matrix_samples}" == "${selected_samples}" ]]
}

selected_merge_outputs_match() {
    [[ -f "${SAMPLES_FILE}" ]] || return 1
    matrix_header_matches_samples "${MERGE_DIR}/genes.counts.matrix" || return 1
    matrix_header_matches_samples "${MERGE_DIR}/genes.DESeq2.normalized_counts.matrix" || return 1
    matrix_header_matches_samples "${MERGE_DIR}/genes.counts.logCPM.matrix" || return 1
    return 0
}

contrasts_complete() {
    [[ -s "${CONTRASTS_FILE}" ]]
}

validate_inputs() {
    [[ -f "${SAMPLES_FILE}" ]] || die "未找到样本文件: ${SAMPLES_FILE}。请先运行 bash rnaseq_first.sh --steps 5"
    [[ -f "${FIRST_SCRIPT}" ]] || die "未找到第一阶段脚本: ${FIRST_SCRIPT}"

    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "samples.txt 中必须恰好包含两个分组，当前得到: ${groups[*]:-<empty>}"

    log "Step 1/5: validate selected samples and inputs"
    log "samples file : ${SAMPLES_FILE}"
    log "matrix file  : ${MATRIX_FILE}"
    log "de results   : ${DE_RESULTS_FILE}"
    log "dry run      : ${DRY_RUN}"
    log "selected samples:"
    show_selected_samples
}

add_step_with_deps() {
    local step="$1"
    case "${step}" in
        1)
            add_once 1
            ;;
        2)
            add_step_with_deps 1
            add_once 2
            ;;
        3)
            if ! selected_merge_outputs_match || step_scheduled 2; then
                add_step_with_deps 2
            else
                add_step_with_deps 1
            fi
            add_once 3
            ;;
        4)
            add_step_with_deps 3
            add_once 4
            ;;
        5)
            if ! is_true "${RNASEQ_SECOND_SKIP_DE}"; then
                add_step_with_deps 4
            else
                add_step_with_deps 1
            fi
            add_once 5
            ;;
        *)
            die "内部错误: 不支持的步骤 ID ${step}"
            ;;
    esac
}

resolve_run_order() {
    local requested=("$@")
    local requested_set=" "
    local step

    RUN_ORDER=()
    for step in "${requested[@]}"; do
        requested_set+="${step} "
    done

    for step in $(seq 1 "${TOTAL_STEPS}"); do
        [[ "${requested_set}" == *" ${step} "* ]] || continue
        add_step_with_deps "${step}"
    done
}

run_step() {
    local step="$1"
    case "${step}" in
        1) validate_inputs ;;
        2) merge_selected_samples ;;
        3) generate_contrasts_file ;;
        4) run_de_analysis ;;
        5) merge_de_results ;;
        *) die "内部错误: 不支持的步骤 ID ${step}" ;;
    esac
}

generate_contrasts_file() {
    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "samples.txt 中必须恰好包含两个分组，当前得到: ${groups[*]:-<empty>}"

    log "Step 3/5: 写入 DE 对比 ${groups[1]} vs ${groups[0]}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ printf '%s %s\\n' ${groups[1]} ${groups[0]} > ${CONTRASTS_FILE}"
        return
    fi

    mkdir -p "${DE_DIR}"
    printf "%s %s\n" "${groups[1]}" "${groups[0]}" > "${CONTRASTS_FILE}"
}

merge_selected_samples() {
    [[ -f "${FIRST_SCRIPT}" ]] || die "未找到第一阶段脚本: ${FIRST_SCRIPT}"

    if [[ "${FORCE}" != "1" ]] && selected_merge_outputs_match; then
        log "当前 3.Merge_result 已与 samples.txt 匹配，直接复用，不再回调 rnaseq_first.sh。"
        return
    fi

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ RNASEQ_SAMPLES_FILE=${SAMPLES_FILE} bash ${FIRST_SCRIPT} --dry-run --steps 6"
        return
    fi

    RNASEQ_SAMPLES_FILE="${SAMPLES_FILE}" \
    FORCE="${FORCE}" \
    DRY_RUN=0 \
    bash "${FIRST_SCRIPT}" --steps 6
}

prepare_de_run_dir() {
    local stamp
    local groups=()
    stamp="$(date +%m%d%H%M%S)"
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "无法为 DESeq2 输出目录解析两个分组"
    DE_RUN_DIR="${DE_DIR}/DESeq2.${stamp}.dir"
    DE_SINGLE_RESULT_FILE="${DE_RUN_DIR}/genes.counts.matrix.${groups[1]}_vs_${groups[0]}.DESeq2.DE_results"
}

run_de_analysis() {
    if is_true "${RNASEQ_SECOND_SKIP_DE}"; then
        log "跳过 DESeq2，沿用现有结果。"
        [[ -f "${DE_RESULTS_FILE}" ]] || log "提醒: 当前未检测到 ${DE_RESULTS_FILE}。"
        return
    fi

    require_cmd Rscript
    prepare_de_run_dir

    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "DESeq2 需要两个分组。"

    log "Step 4/5: 运行内嵌 DESeq2"
    log "  矩阵: ${MATRIX_FILE}"
    log "  样本: ${SAMPLES_FILE}"
    log "  对比: ${groups[1]} vs ${groups[0]}"
    log "  输出: ${DE_RUN_DIR}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline R/DESeq2 -> ${DE_RUN_DIR}"
        return
    fi

    mkdir -p "${DE_RUN_DIR}"

    Rscript --vanilla - "${MATRIX_FILE}" "${SAMPLES_FILE}" "${DE_RUN_DIR}" "${groups[0]}" "${groups[1]}" "${FIG_PNG_DPI}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
matrix_file <- args[[1]]
samples_file <- args[[2]]
out_dir <- args[[3]]
group1 <- args[[4]]
group2 <- args[[5]]
fig_png_dpi <- suppressWarnings(as.numeric(args[[6]]))
if (!is.finite(fig_png_dpi) || fig_png_dpi <= 0) {
  fig_png_dpi <- 300
}

options(scipen = 999)
suppressPackageStartupMessages({
  library(DESeq2)
})

quiet_known_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("S4Vectors:::anyMissing\\(runValue", msg) ||
          grepl("Use 'anyNA\\(\\)' instead\\.", msg) ||
          grepl("makeTxDbFromGRanges\\(\\) has moved from GenomicFeatures to the txdbmaker", msg) ||
          grepl("The \"phase\" metadata column contains non-NA values for features of type", msg, fixed = TRUE) ||
          grepl("genome version information is not available for this TxDb object", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

counts <- read.table(
  matrix_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
counts <- round(data.matrix(counts))

sample_tbl <- read.table(
  samples_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
sample_tbl <- sample_tbl[nzchar(sample_tbl[[1]]) & nzchar(sample_tbl[[2]]), 1:2, drop = FALSE]
colnames(sample_tbl) <- c("group", "sample")

if (!all(sample_tbl$sample %in% colnames(counts))) {
  missing_samples <- sample_tbl$sample[!sample_tbl$sample %in% colnames(counts)]
  stop(sprintf("Samples not found in count matrix: %s", paste(missing_samples, collapse = ", ")))
}

sample_tbl <- sample_tbl[match(sample_tbl$sample, colnames(counts)), , drop = FALSE]
sample_tbl <- sample_tbl[!is.na(sample_tbl$sample), , drop = FALSE]
counts <- counts[, sample_tbl$sample, drop = FALSE]

condition <- factor(sample_tbl$group, levels = c(group1, group2))
col_data <- data.frame(condition = condition, row.names = sample_tbl$sample, check.names = FALSE)

dds <- quiet_known_warnings(DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition))
dds <- quiet_known_warnings(DESeq(dds, quiet = TRUE))
# 主结果保留项目既有口径：不启用 DESeq2 independent filtering，后续仍使用 nominal pvalue < 0.05。
# 辅助结果 res_if 使用 DESeq2 默认 independent filtering，可在正式报告中作为 FDR/过滤敏感性参考。
res <- quiet_known_warnings(results(dds, contrast = c("condition", group2, group1), independentFiltering = FALSE))
res_if <- quiet_known_warnings(results(dds, contrast = c("condition", group2, group1), independentFiltering = TRUE))
norm_counts <- counts(dds, normalized = TRUE)

base_mean_a <- rowMeans(norm_counts[, sample_tbl$group == group2, drop = FALSE])
base_mean_b <- rowMeans(norm_counts[, sample_tbl$group == group1, drop = FALSE])

result_df <- data.frame(
  sampleA = rep(group2, nrow(res)),
  sampleB = rep(group1, nrow(res)),
  baseMeanA = base_mean_a,
  baseMeanB = base_mean_b,
  baseMean = res$baseMean,
  log2FoldChange = res$log2FoldChange,
  lfcSE = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  row.names = rownames(res),
  check.names = FALSE
)

result_df_out <- data.frame(
  gene_id = rownames(result_df),
  result_df,
  check.names = FALSE,
  row.names = NULL
)

prefix <- sprintf("genes.counts.matrix.%s_vs_%s.DESeq2", group2, group1)
result_path <- file.path(out_dir, paste0(prefix, ".DE_results"))
write.table(
  result_df_out,
  file = result_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

result_if_df <- data.frame(
  sampleA = rep(group2, nrow(res_if)),
  sampleB = rep(group1, nrow(res_if)),
  baseMeanA = base_mean_a,
  baseMeanB = base_mean_b,
  baseMean = res_if$baseMean,
  log2FoldChange = res_if$log2FoldChange,
  lfcSE = res_if$lfcSE,
  stat = res_if$stat,
  pvalue = res_if$pvalue,
  padj = res_if$padj,
  row.names = rownames(res_if),
  check.names = FALSE
)
result_if_df_out <- data.frame(
  gene_id = rownames(result_if_df),
  result_if_df,
  check.names = FALSE,
  row.names = NULL
)
write.table(
  result_if_df_out,
  file = file.path(out_dir, paste0(prefix, ".independent_filtering.DE_results")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(norm_counts),
  file = file.path(out_dir, paste0(prefix, ".normalized_counts.tsv")),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  data.frame(sample = names(sizeFactors(dds)), size_factor = sizeFactors(dds), check.names = FALSE),
  file = file.path(out_dir, paste0(prefix, ".size_factors.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ma_pdf <- file.path(out_dir, paste0(prefix, ".MA.pdf"))
ma_png <- sub("\\.pdf$", ".png", ma_pdf)
pdf(ma_pdf, width = 7, height = 6)
plotMA(res, ylim = c(-6, 6), alpha = 0.05, main = sprintf("%s vs %s", group2, group1))
dev.off()
png(ma_png, width = 7, height = 6, units = "in", res = fig_png_dpi, type = "cairo")
plotMA(res, ylim = c(-6, 6), alpha = 0.05, main = sprintf("%s vs %s", group2, group1))
dev.off()

summary_path <- file.path(out_dir, paste0(prefix, ".summary.txt"))
summary_lines <- c(
  sprintf("contrast\t%s_vs_%s", group2, group1),
  sprintf("sampleA\t%s", group2),
  sprintf("sampleB\t%s", group1),
  sprintf("n_genes\t%s", nrow(result_df)),
  sprintf("pvalue_lt_0.05\t%s", sum(!is.na(result_df$pvalue) & result_df$pvalue < 0.05)),
  sprintf("pvalue_lt_0.05_and_abs_log2fc_ge_0.58\t%s", sum(!is.na(result_df$pvalue) & result_df$pvalue < 0.05 & !is.na(result_df$log2FoldChange) & abs(result_df$log2FoldChange) >= 0.58))
)
writeLines(summary_lines, summary_path)

summary_if_path <- file.path(out_dir, paste0(prefix, ".independent_filtering.summary.txt"))
summary_if_lines <- c(
  sprintf("contrast\t%s_vs_%s", group2, group1),
  "note\tDESeq2 default independentFiltering=TRUE; advisory output only",
  sprintf("n_genes\t%s", nrow(result_if_df)),
  sprintf("pvalue_lt_0.05\t%s", sum(!is.na(result_if_df$pvalue) & result_if_df$pvalue < 0.05)),
  sprintf("padj_lt_0.05\t%s", sum(!is.na(result_if_df$padj) & result_if_df$padj < 0.05)),
  sprintf("padj_lt_0.05_and_abs_log2fc_ge_0.58\t%s", sum(!is.na(result_if_df$padj) & result_if_df$padj < 0.05 & !is.na(result_if_df$log2FoldChange) & abs(result_if_df$log2FoldChange) >= 0.58))
)
writeLines(summary_if_lines, summary_if_path)
EOF
}

merge_de_results() {
    if is_true "${RNASEQ_SECOND_SKIP_DE}"; then
        log "Step 5/5: 跳过 DE 合并，沿用既有 ${DE_RESULTS_FILE}"
        return
    fi

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "Step 5/5: 合并 DESeq2 结果"
        log "+ 合并 ${DE_RUN_DIR:-<new_DESeq2_run_dir>}/*.DE_results -> ${DE_RESULTS_FILE}"
        return
    fi

    [[ -n "${DE_RUN_DIR}" && -d "${DE_RUN_DIR}" ]] || die "未找到本次 DESeq2 输出目录"

    local first_result_file=""
    first_result_file="$(find "${DE_RUN_DIR}" -maxdepth 1 -type f -name '*.DESeq2.DE_results' | sort | head -1)"
    [[ -n "${first_result_file}" ]] || die "${DE_RUN_DIR} 中未找到 *.DE_results 文件"

    head -1 "${first_result_file}" > "${DE_RESULTS_FILE}"
    awk 'FNR > 1' "${DE_RUN_DIR}"/*.DESeq2.DE_results >> "${DE_RESULTS_FILE}"

    log "DE: 合并 DE 结果"
    log "  合并结果: ${DE_RESULTS_FILE}"
    log "  总行数  : $(wc -l < "${DE_RESULTS_FILE}")"
}

main() {
    local -a step_args=()
    local run_all=false
    local parsed_steps=""

    if [[ $# -eq 0 ]]; then
        run_all=true
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                usage
                exit 0
                ;;
            --list)
                print_step_table
                exit 0
                ;;
            --all)
                run_all=true
                shift
                ;;
            --steps|--step|-s)
                shift
                [[ $# -gt 0 ]] || die "--steps 需要至少一个步骤编号"
                local step_arg_count_before=${#step_args[@]}
                while [[ $# -gt 0 && "$1" != --* ]]; do
                    step_args+=("$1")
                    shift
                done
                [[ ${#step_args[@]} -gt ${step_arg_count_before} ]] || die "--steps 需要至少一个步骤编号"
                ;;
            --dry-run)
                DRY_RUN=1
                shift
                ;;
            --force)
                FORCE=1
                shift
                ;;
            --threads)
                shift
                [[ $# -gt 0 ]] || die "--threads 需要一个正整数"
                [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "--threads 需要一个正整数，收到: $1"
                THREADS="$1"
                shift
                ;;
            --skip-de)
                RNASEQ_SECOND_SKIP_DE=true
                shift
                ;;
            --run-de)
                RNASEQ_SECOND_SKIP_DE=false
                shift
                ;;
            --)
                shift
                while [[ $# -gt 0 ]]; do
                    step_args+=("$1")
                    shift
                done
                ;;
            --*)
                die "未知参数: $1"
                ;;
            *)
                step_args+=("$1")
                shift
                ;;
        esac
    done

    if [[ "${run_all}" == "true" && ${#step_args[@]} -gt 0 ]]; then
        die "--all 不能和具体步骤混用"
    fi

    if [[ "${run_all}" == "true" || ${#step_args[@]} -eq 0 ]]; then
        REQUESTED_STEPS=("${DEFAULT_STEPS[@]}")
    else
        parsed_steps="$(parse_step_tokens "${step_args[@]}")"
        mapfile -t REQUESTED_STEPS <<< "${parsed_steps}"
    fi

    resolve_run_order "${REQUESTED_STEPS[@]}"
    REQUESTED_STEPS_LABEL="$(format_steps "${REQUESTED_STEPS[@]}")"
    RUN_ORDER_LABEL="$(format_steps "${RUN_ORDER[@]}")"

    log "requested steps : ${REQUESTED_STEPS_LABEL}"
    log "run order       : ${RUN_ORDER_LABEL}"
    log "force           : ${FORCE}"
    log "skip DESeq2     : ${RNASEQ_SECOND_SKIP_DE}"

    local step
    for step in "${RUN_ORDER[@]}"; do
        run_step "${step}"
    done

    log "完成: ${RUN_ORDER_LABEL}"
}

main "$@"
