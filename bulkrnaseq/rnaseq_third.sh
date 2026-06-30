#!/usr/bin/env bash

# Example usage:
#   bash rnaseq_third.sh --list
#   bash rnaseq_third.sh --dry-run --steps 6,9,11
#   nohup bash rnaseq_third.sh --all > 3.log 2>&1 &
# Note: use comma/space-separated step ids; ranges like 6-11 are not supported.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
FIRST_SCRIPT="${ROOT_DIR}/rnaseq_first.sh"
READS_INFO="${ROOT_DIR}/rawdata/reads_info.txt"
FASTQC_DIR="${ROOT_DIR}/fastqc"
MERGE_DIR="${ROOT_DIR}/3.Merge_result"
DE_DIR="${ROOT_DIR}/4.DE_analysis"
POST_TX_DIR="${ROOT_DIR}/5.post_transcriptional_regulation"
SPLICING_DIR="${POST_TX_DIR}/splicing"
RMATS_OUTPUT_DIR="${RMATS_OUTPUT_DIR:-${SPLICING_DIR}/rmats}"
RMATS_TMP_DIR="${RMATS_TMP_DIR:-${POST_TX_DIR}/tmp/rmats}"
EXON_USAGE_DIR="${EXON_USAGE_DIR:-${SPLICING_DIR}/exon_usage}"
ISOFORM_SWITCH_DIR="${ISOFORM_SWITCH_DIR:-${POST_TX_DIR}/isoform_switch}"
APA_DIR="${APA_DIR:-${POST_TX_DIR}/apa}"
APA_REF_DIR="${APA_REF_DIR:-${APA_DIR}/reference_cache}"
DER_DIR="${DER_DIR:-${POST_TX_DIR}/der}"
ISOFORM_REF_DIR="${ISOFORM_REF_DIR:-${ISOFORM_SWITCH_DIR}/reference_cache}"
SALMON_QUANT_DIR="${SALMON_QUANT_DIR:-${ISOFORM_SWITCH_DIR}/salmon_quant}"
ISOFORM_RESULTS_DIR="${ISOFORM_RESULTS_DIR:-${ISOFORM_SWITCH_DIR}/results}"
AS_INTEGRATED_DIR="${AS_INTEGRATED_DIR:-${POST_TX_DIR}/integrated}"

REFERENCE_DIR="${REFERENCE_DIR:-/home/h1028/workspace/reference/GRCm39}"
GTF_FILE="${GTF_FILE:-}"
GENOME_FASTA="${GENOME_FASTA:-}"
SALMON_BIN="${SALMON_BIN:-/home/h1028/miniconda3/bin/salmon}"

SAMPLES_FILE="${DE_DIR}/samples.txt"
CONTRASTS_FILE="${DE_DIR}/contrasts.txt"
MATRIX_FILE="${MERGE_DIR}/genes.counts.matrix"
DE_RESULTS_FILE="${DE_DIR}/DE_results"

THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
RNASEQ_THIRD_SKIP_DE="${RNASEQ_THIRD_SKIP_DE:-${RNASEQ_SECOND_SKIP_DE:-true}}"
RNASEQ_THIRD_SKIP_POST="${RNASEQ_THIRD_SKIP_POST:-${RNASEQ_SECOND_SKIP_AS:-false}}"
RNASEQ_THIRD_METHODS="${RNASEQ_THIRD_METHODS:-${RNASEQ_SECOND_AS_METHODS:-rmats,exon_usage,isoform_switch}}"
RNASEQ_THIRD_EXTRA_METHODS="${RNASEQ_THIRD_EXTRA_METHODS:-${RNASEQ_SECOND_EXTRA_METHODS:-apa,der}}"
FIG_PNG_DPI="${FIG_PNG_DPI:-300}"
THIRD_HEARTBEAT_INTERVAL="${THIRD_HEARTBEAT_INTERVAL:-120}"

RMATS_THREADS="${RMATS_THREADS:-${THREADS}}"
RMATS_READ_LENGTH="${RMATS_READ_LENGTH:-}"
RMATS_LIBTYPE="${RMATS_LIBTYPE:-fr-firststrand}"
RMATS_TASK="${RMATS_TASK:-both}"
RMATS_CSTAT="${RMATS_CSTAT:-0.0001}"
RMATS_VARIABLE_READ_LENGTH="${RMATS_VARIABLE_READ_LENGTH:-true}"
RMATS_ALLOW_CLIPPING="${RMATS_ALLOW_CLIPPING:-false}"
RMATS_NOVELSS="${RMATS_NOVELSS:-false}"
RMATS_INDIVIDUAL_COUNTS="${RMATS_INDIVIDUAL_COUNTS:-true}"
RMATS_EXECUTOR="${RMATS_EXECUTOR:-auto}"
RMATS_IMAGE="${RMATS_IMAGE:-quay.io/biocontainers/rmats:4.3.0--py310ha9d9618_5}"
RMATS_BIND_ROOT="${RMATS_BIND_ROOT:-/home/h1028/workspace}"
TRANSCRIPT_FASTA="${TRANSCRIPT_FASTA:-${ISOFORM_REF_DIR}/transcripts.fa}"
TRANSCRIPT_TX2GENE="${TRANSCRIPT_TX2GENE:-${ISOFORM_REF_DIR}/tx2gene.tsv}"
SALMON_INDEX="${SALMON_INDEX:-${ISOFORM_REF_DIR}/salmon_index}"
SALMON_LIBTYPE="${SALMON_LIBTYPE:-A}"
ISOFORM_SWITCH_ALPHA="${ISOFORM_SWITCH_ALPHA:-0.05}"
ISOFORM_SWITCH_DIF_CUTOFF="${ISOFORM_SWITCH_DIF_CUTOFF:-0.1}"
SASHIMI_TOP_N="${SASHIMI_TOP_N:-8}"
APA_PVALUE="${APA_PVALUE:-0.05}"
APA_DELTA_PAU="${APA_DELTA_PAU:-0.1}"
APA_STRANDTYPE="${APA_STRANDTYPE:-auto}"
APA_SEQTYPE="${APA_SEQTYPE:-ThreeMostPairEnd}"
APA_TEST_METHOD="${APA_TEST_METHOD:-unpaired t-test}"
DER_PVALUE="${DER_PVALUE:-0.05}"
DER_CUTOFF="${DER_CUTOFF:-5}"
DER_GENOME_STYLE="${DER_GENOME_STYLE:-auto}"

COLOR_WT="${COLOR_WT:-#A1C9F4}"
COLOR_PFOS="${COLOR_PFOS:-#FF9F9B}"
COLOR_UP="${COLOR_UP:-#FF9F9B}"
COLOR_DOWN="${COLOR_DOWN:-#A1C9F4}"
COLOR_SIG="${COLOR_SIG:-#FF9F9B}"
COLOR_NS="${COLOR_NS:-#D4D9DE}"

DE_RUN_DIR=""
DE_SINGLE_RESULT_FILE=""
RMATS_B1_FILE=""
RMATS_B2_FILE=""
RMATS_RUN_TMP_DIR=""
RMATS_METADATA_FILE="${RMATS_OUTPUT_DIR}/run_metadata.tsv"
EXON_USAGE_METADATA_FILE="${EXON_USAGE_DIR}/run_metadata.tsv"
ISOFORM_METADATA_FILE="${ISOFORM_RESULTS_DIR}/run_metadata.tsv"
APA_METADATA_FILE="${APA_DIR}/run_metadata.tsv"
DER_METADATA_FILE="${DER_DIR}/run_metadata.tsv"
RMATS_GROUP1=""
RMATS_GROUP2=""
declare -a RMATS_GROUP1_SAMPLES=()
declare -a RMATS_GROUP2_SAMPLES=()
declare -a RMATS_GROUP1_BAMS=()
declare -a RMATS_GROUP2_BAMS=()
declare -a RMATS_CMD=()
declare -a REQUESTED_STEPS=()
declare -a RUN_ORDER=()

TOTAL_STEPS=11
DEFAULT_STEPS=(1 2 3 4 5 6 7 8 9 10 11)
REQUESTED_STEPS_LABEL=""
RUN_ORDER_LABEL=""

declare -A STEP_NAMES=(
    [1]="validate"
    [2]="merge_selected"
    [3]="contrast"
    [4]="deseq2"
    [5]="merge_de_results"
    [6]="rmats"
    [7]="exon_usage"
    [8]="isoform_switch"
    [9]="apa"
    [10]="der"
    [11]="post_tx_integration"
)

declare -A STEP_TITLES=(
    [1]="validate selected samples and inputs"
    [2]="rebuild selected-sample matrices"
    [3]="write DE contrast"
    [4]="run or check DESeq2"
    [5]="merge DESeq2 result files"
    [6]="run rMATS"
    [7]="run exon-usage diffSpliceDGE"
    [8]="run isoform switch analysis"
    [9]="run APA analysis"
    [10]="run DER analysis"
    [11]="integrate post-transcriptional evidence"
)

method_enabled() {
    local method="${1,,}"
    local methods=",${RNASEQ_THIRD_METHODS,,},"
    [[ "${methods}" == *",${method},"* ]]
}

extra_method_enabled() {
    local method="${1,,}"
    local methods=",${RNASEQ_THIRD_EXTRA_METHODS,,},"
    [[ "${methods}" == *",${method},"* ]]
}

log() {
    echo "[third] $*"
}

die() {
    echo "[third][error] $*" >&2
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
        printf '  %02d: %-20s %s\n' "${step}" "${STEP_NAMES[${step}]}" "${STEP_TITLES[${step}]}"
    done
}

usage() {
    cat <<EOF
Examples:
  bash ${SCRIPT_NAME} --list
  bash ${SCRIPT_NAME} --dry-run --steps 6,9,11
  nohup bash ${SCRIPT_NAME} --all > 3.log 2>&1 &

Usage:
  ${SCRIPT_NAME}
  ${SCRIPT_NAME} --all
  ${SCRIPT_NAME} --list
  ${SCRIPT_NAME} --steps 6,9,11
  ${SCRIPT_NAME} --steps 6 9 11
  ${SCRIPT_NAME} 6 9 11

Options:
  --all              Run the full third-stage workflow: $(print_default_order)
  --steps LIST       Run selected steps. LIST can be comma- or space-separated.
  --list             Show available steps.
  --dry-run          Print planned commands without running heavy tools.
  --force            Re-run completed outputs.
  --threads N        Override THREADS/RMATS_THREADS for this run.
  --skip-de          Reuse existing ${DE_RESULTS_FILE}.
  --run-de           Run DESeq2 in steps 4/5.
  --skip-post        Skip post-transcriptional method steps.
  --run-post         Enable post-transcriptional method steps.
  -h, --help         Show this help.

Note:
  Use comma/space-separated step ids; ranges like 6-11 are not supported.

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
        rmats) printf '6\n' ;;
        exon_usage|diffsplice) printf '7\n' ;;
        isoform_switch|isoform) printf '8\n' ;;
        apa) printf '9\n' ;;
        der) printf '10\n' ;;
        post_tx_integration|integration|summary) printf '11\n' ;;
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

format_elapsed() {
    local total_seconds="${1:-0}"
    local hours minutes seconds
    (( total_seconds < 0 )) && total_seconds=0
    hours=$(( total_seconds / 3600 ))
    minutes=$(( (total_seconds % 3600) / 60 ))
    seconds=$(( total_seconds % 60 ))
    if (( hours > 0 )); then
        printf '%02dh:%02dm:%02ds' "${hours}" "${minutes}" "${seconds}"
    else
        printf '%02dm:%02ds' "${minutes}" "${seconds}"
    fi
}

rmats_heartbeat_status() {
    local tmp_dir tmp_files output_files tmp_size output_size
    tmp_dir="${RMATS_RUN_TMP_DIR:-${RMATS_TMP_DIR}}"
    tmp_files="$(find "${tmp_dir}" -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')"
    output_files="$(find "${RMATS_OUTPUT_DIR}" -maxdepth 2 -type f 2>/dev/null | wc -l | tr -d ' ')"
    tmp_size="$(du -sh "${tmp_dir}" 2>/dev/null | awk 'NR==1 {print $1}')"
    output_size="$(du -sh "${RMATS_OUTPUT_DIR}" 2>/dev/null | awk 'NR==1 {print $1}')"
    printf 'tmp=%s(%s files) output=%s(%s files)' \
        "${tmp_size:-0}" "${tmp_files:-0}" "${output_size:-0}" "${output_files:-0}"
}

run_with_heartbeat() {
    local label="$1"
    local status_callback="$2"
    shift 2

    local interval="${THIRD_HEARTBEAT_INTERVAL}"
    local start_ts now elapsed status extra_message
    local worker_pid heartbeat_pid

    start_ts="$(date +%s)"
    # Child long-running jobs like singularity/rMATS can still receive HUP
    # when this wrapper itself is launched under nohup and then backgrounds
    # the worker again. Explicitly nohup the worker to keep it alive.
    nohup "$@" &
    worker_pid=$!

    (
        while kill -0 "${worker_pid}" 2>/dev/null; do
            sleep "${interval}" || true
            kill -0 "${worker_pid}" 2>/dev/null || break
            now="$(date +%s)"
            elapsed="$(( now - start_ts ))"
            extra_message=""
            if [[ -n "${status_callback}" ]] && declare -F "${status_callback}" >/dev/null 2>&1; then
                extra_message="$("${status_callback}" 2>/dev/null || true)"
            fi
            if [[ -n "${extra_message}" ]]; then
                log "  ${label} 仍在运行... elapsed=$(format_elapsed "${elapsed}") ${extra_message}"
            else
                log "  ${label} 仍在运行... elapsed=$(format_elapsed "${elapsed}")"
            fi
        done
    ) &
    heartbeat_pid=$!

    if wait "${worker_pid}"; then
        status=0
    else
        status=$?
    fi

    kill "${heartbeat_pid}" 2>/dev/null || true
    wait "${heartbeat_pid}" 2>/dev/null || true

    now="$(date +%s)"
    elapsed="$(( now - start_ts ))"
    if [[ "${status}" -eq 0 ]]; then
        log "  ${label} 完成，耗时 $(format_elapsed "${elapsed}")"
    else
        log "  ${label} 失败，耗时 $(format_elapsed "${elapsed}")"
    fi
    return "${status}"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "未找到命令: $1"
}

resolve_bin() {
    local current="${1:-}"
    local cmd_name="${2:-}"
    if [[ -n "${current}" && -x "${current}" ]]; then
        printf '%s' "${current}"
        return 0
    fi
    if command -v "${cmd_name}" >/dev/null 2>&1; then
        command -v "${cmd_name}"
        return 0
    fi
    local candidate
    for candidate in \
        "/home/h1028/miniconda3/bin/${cmd_name}" \
        /home/h1028/miniconda3/envs/*/bin/"${cmd_name}" \
        /home/h1028/.conda/envs/*/bin/"${cmd_name}" \
        "/usr/local/bin/${cmd_name}" \
        "/usr/bin/${cmd_name}"
    do
        [[ -x "${candidate}" ]] || continue
        printf '%s' "${candidate}"
        return 0
    done
    return 1
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

normalize_group() {
    printf '%s' "$1" | sed -E 's/([._-]?[Ll]iver)$//'
}

discover_gtf() {
    [[ -d "${REFERENCE_DIR}" ]] || die "未找到参考目录: ${REFERENCE_DIR}"
    if [[ -z "${GTF_FILE}" ]]; then
        GTF_FILE="$(find "${REFERENCE_DIR}" -maxdepth 1 -type f -name 'Mus_musculus.GRCm39*.gtf' | sort | head -1)"
    fi
    [[ -n "${GTF_FILE}" && -f "${GTF_FILE}" ]] || die "未找到 GTF 注释文件，请设置 GTF_FILE 或检查 ${REFERENCE_DIR}"
}

discover_genome_fasta() {
    [[ -d "${REFERENCE_DIR}" ]] || die "未找到参考目录: ${REFERENCE_DIR}"
    if [[ -z "${GENOME_FASTA}" ]]; then
        GENOME_FASTA="$(find "${REFERENCE_DIR}" -maxdepth 1 -type f \( -name 'Mus_musculus.GRCm39*.dna*.fa' -o -name '*.fa' -o -name '*.fasta' \) | sort | head -1)"
    fi
    [[ -n "${GENOME_FASTA}" && -f "${GENOME_FASTA}" ]] || die "未找到 genome FASTA，请设置 GENOME_FASTA 或检查 ${REFERENCE_DIR}"
}

require_r_packages() {
    local label="$1"
    shift
    local pkgs=("$@")
    [[ ${#pkgs[@]} -gt 0 ]] || return 0

    Rscript --vanilla - "${label}" "${pkgs[@]}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
label <- args[[1]]
pkgs <- args[-1]
installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)
if (length(missing) > 0) {
  message(sprintf("%s 缺少 R 包: %s", label, paste(missing, collapse = ", ")))
  quit(status = 1)
}
EOF
}

resolve_clean_fastq() {
    local sample="$1"
    local mate="$2"
    local raw_path=""
    local fallback="${FASTQC_DIR}/${sample}_${mate}.fastq.gz"

    if [[ -f "${READS_INFO}" ]]; then
        raw_path="$(
            awk -F'\t' -v sample="${sample}" -v mate="${mate}" '
                NR == 1 {
                    for (i = 1; i <= NF; i++) {
                        header[$i] = i
                    }
                    next
                }
                header["sample"] > 0 && $header["sample"] == sample {
                    idx = (mate == "1" ? header["raw_fq1"] : header["raw_fq2"])
                    if (idx > 0) {
                        print $idx
                    }
                    exit
                }
            ' "${READS_INFO}"
        )"
    fi

    if [[ -n "${raw_path}" ]]; then
        printf '%s\n' "${FASTQC_DIR}/$(basename "${raw_path}")"
    else
        printf '%s\n' "${fallback}"
    fi
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
    if is_true "${RNASEQ_THIRD_SKIP_DE}" && [[ ! -f "${DE_RESULTS_FILE}" ]] && { step_scheduled 5 || step_scheduled 11; }; then
        die "未找到 DE 结果: ${DE_RESULTS_FILE}。请先运行 bash rnaseq_second.sh，或使用 --run-de。"
    fi

    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "samples.txt 中必须恰好包含两个分组，当前得到: ${groups[*]:-<empty>}"

    log "Step 1/${TOTAL_STEPS}: validate selected samples and inputs"
    log "samples file : ${SAMPLES_FILE}"
    log "matrix file  : ${MATRIX_FILE}"
    log "de results   : ${DE_RESULTS_FILE}"
    log "post-tx dir  : ${POST_TX_DIR}"
    log "rmats dir    : ${RMATS_OUTPUT_DIR}"
    log "methods      : ${RNASEQ_THIRD_METHODS}"
    log "extra methods: ${RNASEQ_THIRD_EXTRA_METHODS}"
    log "dry run      : ${DRY_RUN}"
    log "skip DESeq2  : ${RNASEQ_THIRD_SKIP_DE}"
    log "skip post-tx : ${RNASEQ_THIRD_SKIP_POST}"
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
            if ! is_true "${RNASEQ_THIRD_SKIP_DE}"; then
                add_step_with_deps 3
            else
                add_step_with_deps 1
            fi
            add_once 4
            ;;
        5)
            if ! is_true "${RNASEQ_THIRD_SKIP_DE}"; then
                add_step_with_deps 4
            else
                add_step_with_deps 1
            fi
            add_once 5
            ;;
        6|7|8|9|10)
            add_step_with_deps 1
            add_once "${step}"
            ;;
        11)
            add_step_with_deps 1
            add_once 11
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
        6) run_rmats ;;
        7) run_exon_usage ;;
        8) run_isoform_switch ;;
        9) run_apa ;;
        10) run_der ;;
        11) run_post_tx_summary ;;
        *) die "内部错误: 不支持的步骤 ID ${step}" ;;
    esac
}

generate_contrasts_file() {
    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "samples.txt 中必须恰好包含两个分组，当前得到: ${groups[*]:-<empty>}"

    log "Step 3/${TOTAL_STEPS}: 写入 DE 对比 ${groups[1]} vs ${groups[0]}"

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
    if is_true "${RNASEQ_THIRD_SKIP_DE}"; then
        log "跳过 DESeq2，沿用现有结果。"
        [[ -f "${DE_RESULTS_FILE}" ]] || log "提醒: 当前未检测到 ${DE_RESULTS_FILE}，后续 DE-DAS 重叠会按空集合处理。"
        return
    fi

    require_cmd Rscript
    prepare_de_run_dir

    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "DESeq2 需要两个分组。"

    log "Step 4/${TOTAL_STEPS}: 运行内嵌 DESeq2"
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
      if (grepl("S4Vectors:::anyMissing", msg, fixed = TRUE) ||
          grepl("Use 'anyNA()' instead.", msg, fixed = TRUE) ||
          grepl("makeTxDbFromGRanges", msg, fixed = TRUE) ||
          grepl("txdbmaker", msg, fixed = TRUE) ||
          grepl("The \"phase\" metadata column contains non-NA values for features of type", msg, fixed = TRUE) ||
          grepl("stop_codon", msg, fixed = TRUE) ||
          grepl("genome version information is not available for this TxDb object", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

drop_na_seqnames <- function(gr) {
  if (is.null(gr) || length(gr) == 0) {
    return(gr)
  }
  keep <- !is.na(as.character(GenomicRanges::seqnames(gr))) &
    nzchar(as.character(GenomicRanges::seqnames(gr)))
  keep[is.na(keep)] <- FALSE
  gr[keep]
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
    if is_true "${RNASEQ_THIRD_SKIP_DE}"; then
        log "Step 5/${TOTAL_STEPS}: 跳过 DE 合并，沿用既有 ${DE_RESULTS_FILE}"
        return
    fi

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "Step 5/${TOTAL_STEPS}: 合并 DESeq2 结果"
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

infer_read_length() {
    local candidate=""
    if [[ -f "${READS_INFO}" ]]; then
        candidate="$(
            awk '
                NR == 1 && $1 == "group_stage" && $2 == "sample" {next}
                NR == 1 && $1 == "sample" && $2 == "group" {next}
                NF > 0 {print (NF >= 6 ? $5 : $3); exit}
            ' "${READS_INFO}"
        )"
        candidate="${ROOT_DIR}/rawdata/$(basename "${candidate}")"
    fi
    if [[ -z "${candidate}" || ! -f "${candidate}" ]]; then
        candidate="$(find "${ROOT_DIR}/rawdata" -maxdepth 1 -type f -name '*.fastq.gz' | sort | head -1)"
    fi
    if [[ -n "${candidate}" && -f "${candidate}" ]]; then
        python3 - "${candidate}" <<'EOF'
import gzip
import sys

fq = sys.argv[1]
with gzip.open(fq, "rt") as handle:
    handle.readline()
    seq = handle.readline().strip()
print(len(seq))
EOF
        return
    fi

    local bam_candidate=""
    local one_bam=""
    for one_bam in "${RMATS_GROUP1_BAMS[@]}" "${RMATS_GROUP2_BAMS[@]}"; do
        if [[ -n "${one_bam}" && -f "${one_bam}" ]]; then
            bam_candidate="${one_bam}"
            break
        fi
    done
    if [[ -z "${bam_candidate}" ]]; then
        bam_candidate="$(find "${ROOT_DIR}/1.Mapping" -maxdepth 1 -type f -name '*.bam' | sort | head -1)"
    fi
    [[ -n "${bam_candidate}" && -f "${bam_candidate}" ]] || die "无法自动推断 read length：既没有 rawdata FASTQ，也没有可用 BAM，请设置 RMATS_READ_LENGTH。"

    require_cmd samtools
    local bam_read_length=""
    bam_read_length="$(
        set +o pipefail
        samtools view "${bam_candidate}" | awk '
            length($10) > 0 && $10 != "*" {
                print length($10)
                found = 1
                exit
            }
            END {
                if (!found) exit 1
            }
        '
    )"
    [[ -n "${bam_read_length}" ]] || die "无法从 BAM 推断 read length: ${bam_candidate}"
    printf '%s\n' "${bam_read_length}"
}

collect_rmats_inputs() {
    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "rMATS 需要两个分组。"
    RMATS_GROUP1="${groups[0]}"
    RMATS_GROUP2="${groups[1]}"
    RMATS_GROUP1_SAMPLES=()
    RMATS_GROUP2_SAMPLES=()
    RMATS_GROUP1_BAMS=()
    RMATS_GROUP2_BAMS=()

    while read -r group sample _; do
        [[ -n "${group:-}" && -n "${sample:-}" ]] || continue
        [[ "${group}" =~ ^# ]] && continue
        local bam="${ROOT_DIR}/1.Mapping/${sample}.bam"
        if [[ "${DRY_RUN}" != "1" ]]; then
            [[ -f "${bam}" ]] || die "缺少 BAM 文件: ${bam}"
        fi
        if [[ "${group}" == "${RMATS_GROUP1}" ]]; then
            RMATS_GROUP1_SAMPLES+=("${sample}")
            RMATS_GROUP1_BAMS+=("${bam}")
        elif [[ "${group}" == "${RMATS_GROUP2}" ]]; then
            RMATS_GROUP2_SAMPLES+=("${sample}")
            RMATS_GROUP2_BAMS+=("${bam}")
        else
            die "samples.txt 中存在额外分组: ${group}"
        fi
    done < "${SAMPLES_FILE}"

    [[ ${#RMATS_GROUP1_BAMS[@]} -gt 0 ]] || die "${RMATS_GROUP1} 组没有可用 BAM"
    [[ ${#RMATS_GROUP2_BAMS[@]} -gt 0 ]] || die "${RMATS_GROUP2} 组没有可用 BAM"
}

join_by_comma() {
    local IFS=,
    echo "$*"
}

rmats_outputs_complete() {
    local event mode
    for event in SE A5SS A3SS RI MXE; do
        for mode in JC JCEC; do
            if [[ ! -f "${RMATS_OUTPUT_DIR}/${event}.MATS.${mode}.txt" && ! -f "${RMATS_OUTPUT_DIR}/fromGTF.${event}.MATS.${mode}.txt" ]]; then
                return 1
            fi
        done
    done
    return 0
}

rmats_metadata_matches() {
    [[ -f "${RMATS_METADATA_FILE}" ]] || return 1
    local expected_key
    expected_key="$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
    grep -Fqx "sample_key=${expected_key}" "${RMATS_METADATA_FILE}" || return 1
    grep -Fqx "group1=${RMATS_GROUP1}" "${RMATS_METADATA_FILE}" || return 1
    grep -Fqx "group2=${RMATS_GROUP2}" "${RMATS_METADATA_FILE}" || return 1
}

write_rmats_metadata() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${RMATS_OUTPUT_DIR}"
    {
        echo "group1=${RMATS_GROUP1}"
        echo "group2=${RMATS_GROUP2}"
        echo "group1_samples=$(join_by_comma "${RMATS_GROUP1_SAMPLES[@]}")"
        echo "group2_samples=$(join_by_comma "${RMATS_GROUP2_SAMPLES[@]}")"
        echo "sample_key=$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
        echo "read_length=${RMATS_READ_LENGTH}"
        echo "executor=${RMATS_EXECUTOR}"
        echo "image=${RMATS_IMAGE}"
        echo "task=${RMATS_TASK}"
        echo "tmp_dir=${RMATS_RUN_TMP_DIR:-${RMATS_TMP_DIR}}"
        echo "timestamp=$(date '+%F %T')"
    } > "${RMATS_METADATA_FILE}"
}

exon_usage_outputs_complete() {
    local required=(
        "${EXON_USAGE_DIR}/exon_count_matrix.tsv"
        "${EXON_USAGE_DIR}/exon_annotation.tsv"
        "${EXON_USAGE_DIR}/exon_logCPM.tsv"
        "${EXON_USAGE_DIR}/diffsplice_gene_simes.csv"
        "${EXON_USAGE_DIR}/diffsplice_gene_f.csv"
        "${EXON_USAGE_DIR}/diffsplice_exon.csv"
        "${EXON_USAGE_DIR}/significant_genes.csv"
        "${EXON_USAGE_DIR}/significant_exons.csv"
        "${EXON_USAGE_DIR}/gene_exon_burden.csv"
    )
    local path
    for path in "${required[@]}"; do
        [[ -s "${path}" ]] || return 1
    done
    return 0
}

exon_usage_metadata_matches() {
    [[ -f "${EXON_USAGE_METADATA_FILE}" ]] || return 1
    local expected_key
    expected_key="$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
    grep -Fqx "sample_key=${expected_key}" "${EXON_USAGE_METADATA_FILE}" || return 1
    grep -Fqx "group1=${RMATS_GROUP1}" "${EXON_USAGE_METADATA_FILE}" || return 1
    grep -Fqx "group2=${RMATS_GROUP2}" "${EXON_USAGE_METADATA_FILE}" || return 1
    grep -Fqx "gtf=${GTF_FILE}" "${EXON_USAGE_METADATA_FILE}" || return 1
}

write_exon_usage_metadata() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${EXON_USAGE_DIR}"
    {
        echo "group1=${RMATS_GROUP1}"
        echo "group2=${RMATS_GROUP2}"
        echo "group1_samples=$(join_by_comma "${RMATS_GROUP1_SAMPLES[@]}")"
        echo "group2_samples=$(join_by_comma "${RMATS_GROUP2_SAMPLES[@]}")"
        echo "sample_key=$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
        echo "gtf=${GTF_FILE}"
        echo "threads=${THREADS}"
        echo "method=edgeR::diffSpliceDGE"
        echo "method_note=exon_usage_test_not_gene_level_DESeq2"
        echo "timestamp=$(date '+%F %T')"
    } > "${EXON_USAGE_METADATA_FILE}"
}

isoform_outputs_complete() {
    local required=(
        "${TRANSCRIPT_FASTA}"
        "${TRANSCRIPT_TX2GENE}"
        "${ISOFORM_RESULTS_DIR}/sample_manifest.tsv"
        "${ISOFORM_RESULTS_DIR}/salmon_quant_manifest.tsv"
        "${ISOFORM_RESULTS_DIR}/tximport_counts.tsv"
        "${ISOFORM_RESULTS_DIR}/tximport_abundance.tsv"
        "${ISOFORM_RESULTS_DIR}/isoform_fraction_matrix.tsv"
        "${ISOFORM_RESULTS_DIR}/isoform_feature_table.csv"
        "${ISOFORM_RESULTS_DIR}/significant_switch_isoforms.csv"
        "${ISOFORM_RESULTS_DIR}/significant_switch_genes.csv"
        "${ISOFORM_RESULTS_DIR}/switch_summary.csv"
        "${ISOFORM_RESULTS_DIR}/top_switches.csv"
        "${ISOFORM_RESULTS_DIR}/top_switches_with_consequences.csv"
        "${ISOFORM_RESULTS_DIR}/consequence_summary.csv"
        "${ISOFORM_RESULTS_DIR}/representative_switch_pairs.csv"
    )
    local path
    for path in "${required[@]}"; do
        [[ -s "${path}" ]] || return 1
    done
    [[ -d "${SALMON_INDEX}" ]] || return 1
    return 0
}

isoform_metadata_matches() {
    [[ -f "${ISOFORM_METADATA_FILE}" ]] || return 1
    local expected_key
    expected_key="$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
    grep -Fqx "sample_key=${expected_key}" "${ISOFORM_METADATA_FILE}" || return 1
    grep -Fqx "group1=${RMATS_GROUP1}" "${ISOFORM_METADATA_FILE}" || return 1
    grep -Fqx "group2=${RMATS_GROUP2}" "${ISOFORM_METADATA_FILE}" || return 1
    grep -Fqx "gtf=${GTF_FILE}" "${ISOFORM_METADATA_FILE}" || return 1
    grep -Fqx "genome_fasta=${GENOME_FASTA}" "${ISOFORM_METADATA_FILE}" || return 1
    grep -Fqx "salmon_libtype=${SALMON_LIBTYPE}" "${ISOFORM_METADATA_FILE}" || return 1
}

write_isoform_metadata() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${ISOFORM_RESULTS_DIR}"
    {
        echo "group1=${RMATS_GROUP1}"
        echo "group2=${RMATS_GROUP2}"
        echo "group1_samples=$(join_by_comma "${RMATS_GROUP1_SAMPLES[@]}")"
        echo "group2_samples=$(join_by_comma "${RMATS_GROUP2_SAMPLES[@]}")"
        echo "sample_key=$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
        echo "gtf=${GTF_FILE}"
        echo "genome_fasta=${GENOME_FASTA}"
        echo "transcript_fasta=${TRANSCRIPT_FASTA}"
        echo "transcript_tx2gene=${TRANSCRIPT_TX2GENE}"
        echo "salmon_index=${SALMON_INDEX}"
        echo "salmon_libtype=${SALMON_LIBTYPE}"
        echo "alpha=${ISOFORM_SWITCH_ALPHA}"
        echo "dif_cutoff=${ISOFORM_SWITCH_DIF_CUTOFF}"
        echo "timestamp=$(date '+%F %T')"
    } > "${ISOFORM_METADATA_FILE}"
}

resolve_apa_strandtype() {
    local libtype="${RMATS_LIBTYPE,,}"
    case "${APA_STRANDTYPE}" in
        auto|AUTO|Auto)
            case "${libtype}" in
                fr-firststrand|fr-first|rf|firststrand|first-strand)
                    printf '%s\n' "invert"
                    ;;
                fr-secondstrand|fr-second|fr|secondstrand|second-strand)
                    printf '%s\n' "forward"
                    ;;
                *)
                    printf '%s\n' "NONE"
                    ;;
            esac
            ;;
        *)
            printf '%s\n' "${APA_STRANDTYPE}"
            ;;
    esac
}

apa_outputs_complete() {
    local required=(
        "${APA_DIR}/sample_manifest.tsv"
        "${APA_DIR}/pas_reference.rds"
        "${APA_DIR}/reference_gene_map.tsv"
        "${APA_DIR}/3utr_expression_raw.tsv"
        "${APA_DIR}/ipa_expression_raw.tsv"
        "${APA_DIR}/apa_3utr_events.csv"
        "${APA_DIR}/apa_ipa_events.csv"
        "${APA_DIR}/apa_all_events.csv"
        "${APA_DIR}/significant_apa_events.csv"
        "${APA_DIR}/apa_gene_summary.csv"
        "${APA_DIR}/top_apa_candidates.csv"
    )
    local path
    for path in "${required[@]}"; do
        [[ -s "${path}" ]] || return 1
    done
    return 0
}

apa_metadata_matches() {
    [[ -f "${APA_METADATA_FILE}" ]] || return 1
    local expected_key
    expected_key="$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
    grep -Fqx "sample_key=${expected_key}" "${APA_METADATA_FILE}" || return 1
    grep -Fqx "group1=${RMATS_GROUP1}" "${APA_METADATA_FILE}" || return 1
    grep -Fqx "group2=${RMATS_GROUP2}" "${APA_METADATA_FILE}" || return 1
    grep -Fqx "gtf=${GTF_FILE}" "${APA_METADATA_FILE}" || return 1
    grep -Fqx "genome_fasta=${GENOME_FASTA}" "${APA_METADATA_FILE}" || return 1
}

write_apa_metadata() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${APA_DIR}"
    {
        echo "group1=${RMATS_GROUP1}"
        echo "group2=${RMATS_GROUP2}"
        echo "group1_samples=$(join_by_comma "${RMATS_GROUP1_SAMPLES[@]}")"
        echo "group2_samples=$(join_by_comma "${RMATS_GROUP2_SAMPLES[@]}")"
        echo "sample_key=$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
        echo "gtf=${GTF_FILE}"
        echo "genome_fasta=${GENOME_FASTA}"
        echo "apa_strandtype=$(resolve_apa_strandtype)"
        echo "apa_seqtype=${APA_SEQTYPE}"
        echo "apa_pvalue=${APA_PVALUE}"
        echo "apa_delta_pau=${APA_DELTA_PAU}"
        echo "timestamp=$(date '+%F %T')"
    } > "${APA_METADATA_FILE}"
}

der_outputs_complete() {
    local required=(
        "${DER_DIR}/sample_manifest.tsv"
        "${DER_DIR}/coverage_matrix.tsv"
        "${DER_DIR}/all_ders.csv"
        "${DER_DIR}/significant_ders.csv"
        "${DER_DIR}/annotated_ders.csv"
        "${DER_DIR}/der_gene_summary.csv"
    )
    local path
    for path in "${required[@]}"; do
        [[ -s "${path}" ]] || return 1
    done
    return 0
}

der_metadata_matches() {
    [[ -f "${DER_METADATA_FILE}" ]] || return 1
    local expected_key
    expected_key="$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
    grep -Fqx "sample_key=${expected_key}" "${DER_METADATA_FILE}" || return 1
    grep -Fqx "group1=${RMATS_GROUP1}" "${DER_METADATA_FILE}" || return 1
    grep -Fqx "group2=${RMATS_GROUP2}" "${DER_METADATA_FILE}" || return 1
    grep -Fqx "gtf=${GTF_FILE}" "${DER_METADATA_FILE}" || return 1
    grep -Fqx "genome_fasta=${GENOME_FASTA}" "${DER_METADATA_FILE}" || return 1
}

write_der_metadata() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${DER_DIR}"
    {
        echo "group1=${RMATS_GROUP1}"
        echo "group2=${RMATS_GROUP2}"
        echo "group1_samples=$(join_by_comma "${RMATS_GROUP1_SAMPLES[@]}")"
        echo "group2_samples=$(join_by_comma "${RMATS_GROUP2_SAMPLES[@]}")"
        echo "sample_key=$(printf '%s\n' "${RMATS_GROUP1_SAMPLES[@]}" "${RMATS_GROUP2_SAMPLES[@]}" | sort | paste -sd ',' -)"
        echo "gtf=${GTF_FILE}"
        echo "genome_fasta=${GENOME_FASTA}"
        echo "der_cutoff=${DER_CUTOFF}"
        echo "der_pvalue=${DER_PVALUE}"
        echo "read_length=${RMATS_READ_LENGTH}"
        echo "timestamp=$(date '+%F %T')"
    } > "${DER_METADATA_FILE}"
}

build_rmats_command() {
    local singularity_image_ref
    if [[ "${RMATS_IMAGE}" == docker://* || "${RMATS_IMAGE}" == oras://* || "${RMATS_IMAGE}" == library://* ]]; then
        singularity_image_ref="${RMATS_IMAGE}"
    else
        singularity_image_ref="docker://${RMATS_IMAGE}"
    fi

    case "${RMATS_EXECUTOR}" in
        auto)
            if command -v rmats.py >/dev/null 2>&1; then
                RMATS_CMD=(python "$(command -v rmats.py)")
                return
            fi
            if command -v singularity >/dev/null 2>&1; then
                RMATS_CMD=(
                    singularity exec
                    --bind "${RMATS_BIND_ROOT}:${RMATS_BIND_ROOT}"
                    "${singularity_image_ref}"
                    rmats.py
                )
                return
            fi
            if command -v docker >/dev/null 2>&1; then
                RMATS_CMD=(
                    docker run --rm
                    -u "$(id -u):$(id -g)"
                    -v "${RMATS_BIND_ROOT}:${RMATS_BIND_ROOT}"
                    -w "${ROOT_DIR}"
                    "${RMATS_IMAGE}"
                    rmats.py
                )
                return
            fi
            ;;
        local)
            command -v rmats.py >/dev/null 2>&1 || die "RMATS_EXECUTOR=local，但未找到 rmats.py"
            RMATS_CMD=(python "$(command -v rmats.py)")
            return
            ;;
        singularity)
            command -v singularity >/dev/null 2>&1 || die "RMATS_EXECUTOR=singularity，但未找到 singularity"
            RMATS_CMD=(
                singularity exec
                --bind "${RMATS_BIND_ROOT}:${RMATS_BIND_ROOT}"
                "${singularity_image_ref}"
                rmats.py
            )
            return
            ;;
        docker)
            command -v docker >/dev/null 2>&1 || die "RMATS_EXECUTOR=docker，但未找到 docker"
            RMATS_CMD=(
                docker run --rm
                -u "$(id -u):$(id -g)"
                -v "${RMATS_BIND_ROOT}:${RMATS_BIND_ROOT}"
                -w "${ROOT_DIR}"
                "${RMATS_IMAGE}"
                rmats.py
            )
            return
            ;;
        *)
            die "不支持的 RMATS_EXECUTOR: ${RMATS_EXECUTOR}"
            ;;
    esac

    die "未检测到 rmats.py、singularity 或 docker，无法运行 rMATS"
}

run_rmats() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        log "跳过后转录层分析。"
        return
    fi
    if ! method_enabled "rmats"; then
        die "当前第三阶段要求完整后转录流程，但 RNASEQ_THIRD_METHODS=${RNASEQ_THIRD_METHODS} 未包含 rmats"
    fi

    discover_gtf
    collect_rmats_inputs

    if [[ -z "${RMATS_READ_LENGTH}" ]]; then
        if [[ "${DRY_RUN}" == "1" ]]; then
            RMATS_READ_LENGTH="150"
        else
            RMATS_READ_LENGTH="$(infer_read_length)"
        fi
    fi

    if [[ "${FORCE}" != "1" ]] && [[ "${RMATS_TASK}" == "both" ]] && rmats_outputs_complete && rmats_metadata_matches; then
        log "PostTx-1: 检测到匹配当前 samples.txt 的 rMATS 输出，跳过重跑。"
        return
    fi

    mkdir -p "${RMATS_OUTPUT_DIR}" "${RMATS_TMP_DIR}"
    case "${RMATS_TASK,,}" in
        both|prep)
            RMATS_RUN_TMP_DIR="$(mktemp -d "${RMATS_TMP_DIR}/run_$(date '+%Y%m%d_%H%M%S').XXXXXX")"
            ;;
        *)
            RMATS_RUN_TMP_DIR="${RMATS_TMP_DIR}"
            ;;
    esac
    RMATS_B1_FILE="$(mktemp "${RMATS_RUN_TMP_DIR}/rmats_b1.XXXXXX.txt")"
    RMATS_B2_FILE="$(mktemp "${RMATS_RUN_TMP_DIR}/rmats_b2.XXXXXX.txt")"
    printf '%s' "$(join_by_comma "${RMATS_GROUP1_BAMS[@]}")" > "${RMATS_B1_FILE}"
    printf '%s' "$(join_by_comma "${RMATS_GROUP2_BAMS[@]}")" > "${RMATS_B2_FILE}"
    trap 'rm -f "${RMATS_B1_FILE:-}" "${RMATS_B2_FILE:-}"' EXIT

    build_rmats_command

    local rmats_args=(
        --b1 "${RMATS_B1_FILE}"
        --b2 "${RMATS_B2_FILE}"
        --gtf "${GTF_FILE}"
        --od "${RMATS_OUTPUT_DIR}"
        --tmp "${RMATS_RUN_TMP_DIR}"
        --readLength "${RMATS_READ_LENGTH}"
        --nthread "${RMATS_THREADS}"
        --tstat "${RMATS_THREADS}"
        --libType "${RMATS_LIBTYPE}"
        --task "${RMATS_TASK}"
        --cstat "${RMATS_CSTAT}"
        -t paired
    )

    if is_true "${RMATS_VARIABLE_READ_LENGTH}"; then
        rmats_args+=(--variable-read-length)
    fi
    if is_true "${RMATS_ALLOW_CLIPPING}"; then
        rmats_args+=(--allow-clipping)
    fi
    if is_true "${RMATS_NOVELSS}"; then
        rmats_args+=(--novelSS)
    fi
    if is_true "${RMATS_INDIVIDUAL_COUNTS}"; then
        rmats_args+=(--individual-counts)
    fi

    log "PostTx-1: 运行 rMATS"
    log "  GTF          : ${GTF_FILE}"
    log "  read length  : ${RMATS_READ_LENGTH}"
    log "  ${RMATS_GROUP1} samples: ${RMATS_GROUP1_SAMPLES[*]}"
    log "  ${RMATS_GROUP2} samples: ${RMATS_GROUP2_SAMPLES[*]}"
    log "  output dir   : ${RMATS_OUTPUT_DIR}"
    log "  tmp dir      : ${RMATS_RUN_TMP_DIR}"
    log "  executor     : ${RMATS_CMD[*]}"
    log "  heartbeat    : every ${THIRD_HEARTBEAT_INTERVAL}s"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ $(printf '%q ' "${RMATS_CMD[@]}" "${rmats_args[@]}")"
        return
    fi

    run_with_heartbeat "rMATS" "rmats_heartbeat_status" "${RMATS_CMD[@]}" "${rmats_args[@]}"
    write_rmats_metadata
}

run_exon_usage() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        return
    fi
    if ! method_enabled "exon_usage"; then
        log "跳过 exon-usage，当前方法集: ${RNASEQ_THIRD_METHODS}"
        return
    fi

    discover_gtf
    collect_rmats_inputs

    mkdir -p "${EXON_USAGE_DIR}" "${AS_INTEGRATED_DIR}"

    if [[ "${FORCE}" != "1" ]] && exon_usage_outputs_complete && exon_usage_metadata_matches; then
        log "PostTx-2: 检测到匹配当前 samples.txt 的 exon-usage 输出，跳过重跑。"
        return
    fi

    log "PostTx-2: 运行 exon-usage 差异剪接验证"
    log "  GTF        : ${GTF_FILE}"
    log "  output dir : ${EXON_USAGE_DIR}"
    log "  methods    : ${RNASEQ_THIRD_METHODS}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline Rsubread::featureCounts exon-level counting -> ${EXON_USAGE_DIR}"
        log "+ inline edgeR::diffSpliceDGE -> ${EXON_USAGE_DIR}/diffsplice_*.csv"
        return
    fi

    Rscript --vanilla - "${ROOT_DIR}" "${SAMPLES_FILE}" "${GTF_FILE}" "${EXON_USAGE_DIR}" "${RMATS_GROUP1}" "${RMATS_GROUP2}" "${THREADS}" "${FIG_PNG_DPI}" "${COLOR_SIG}" "${COLOR_NS}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[[1]]
samples_file <- args[[2]]
gtf_file <- args[[3]]
out_dir <- args[[4]]
group1 <- args[[5]]
group2 <- args[[6]]
threads <- as.integer(args[[7]])
fig_png_dpi <- suppressWarnings(as.numeric(args[[8]]))
color_sig <- args[[9]]
color_ns <- args[[10]]

if (!is.finite(threads) || threads <= 0) {
  threads <- 8L
}
if (!is.finite(fig_png_dpi) || fig_png_dpi <= 0) {
  fig_png_dpi <- 300
}

options(scipen = 999)
suppressPackageStartupMessages({
  library(Rsubread)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(ggplot2)
})

pdf_device_fn <- function() {
  if (requireNamespace("Cairo", quietly = TRUE) && "CairoPDF" %in% getNamespaceExports("Cairo")) {
    Cairo::CairoPDF
  } else {
    grDevices::pdf
  }
}

png_device_fn <- function() {
  if (requireNamespace("Cairo", quietly = TRUE) && "CairoPNG" %in% getNamespaceExports("Cairo")) {
    Cairo::CairoPNG
  } else {
    grDevices::png
  }
}

png_path_from_pdf <- function(path) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    sub("\\.pdf$", ".png", path, ignore.case = TRUE)
  } else {
    paste0(path, ".png")
  }
}

save_plot_versions <- function(plot, primary, width = 8, height = 6) {
  dir.create(dirname(primary), recursive = TRUE, showWarnings = FALSE)
  ggsave(primary, plot = plot, width = width, height = height, device = pdf_device_fn())
  ggsave(
    png_path_from_pdf(primary),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = fig_png_dpi,
    device = png_device_fn()
  )
}

drop_na_seqnames <- function(gr) {
  if (is.null(gr) || length(gr) == 0) {
    return(gr)
  }
  keep <- !is.na(as.character(GenomicRanges::seqnames(gr))) &
    nzchar(as.character(GenomicRanges::seqnames(gr)))
  keep[is.na(keep)] <- FALSE
  gr[keep]
}

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
sample_tbl$group <- factor(sample_tbl$group, levels = c(group1, group2))

if (any(is.na(sample_tbl$group))) {
  stop("samples.txt 中包含不属于当前对比的分组")
}

bam_files <- file.path(root_dir, "1.Mapping", paste0(sample_tbl$sample, ".bam"))
if (!all(file.exists(bam_files))) {
  missing_bams <- bam_files[!file.exists(bam_files)]
  stop(sprintf("缺少 BAM 文件: %s", paste(missing_bams, collapse = ", ")))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fc <- Rsubread::featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = FALSE,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  strandSpecific = 2,
  checkFragLength = FALSE,
  countMultiMappingReads = FALSE,
  allowMultiOverlap = FALSE,
  nthreads = threads
)

counts <- as.matrix(fc$counts)
colnames(counts) <- sample_tbl$sample
anno <- as.data.frame(fc$annotation, stringsAsFactors = FALSE, check.names = FALSE)

gene_id_col <- dplyr::coalesce(
  if ("GeneID" %in% colnames(anno)) "GeneID" else NA_character_,
  if ("Geneid" %in% colnames(anno)) "Geneid" else NA_character_
)
if (is.na(gene_id_col)) {
  stop("featureCounts 注释结果中未找到 GeneID/Geneid 列")
}

anno$gene_id <- sub("\\.[0-9]+$", "", anno[[gene_id_col]])
anno$exon_id <- make.unique(paste(anno$Chr, anno$Start, anno$End, anno$Strand, anno$gene_id, sep = ":"))
rownames(counts) <- anno$exon_id

write.table(
  counts,
  file = file.path(out_dir, "exon_count_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.csv(anno, file.path(out_dir, "exon_annotation.tsv"), row.names = FALSE)
write.table(
  edgeR::cpm(counts, log = TRUE, prior.count = 1),
  file = file.path(out_dir, "exon_logCPM.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  data.frame(sample = sample_tbl$sample, bam = bam_files, group = sample_tbl$group, stringsAsFactors = FALSE),
  file = file.path(out_dir, "sample_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

y <- edgeR::DGEList(counts = counts, genes = anno, group = sample_tbl$group)
keep <- edgeR::filterByExpr(y, group = sample_tbl$group)
if (sum(keep) < 2) {
  stop("exon-usage 分析通过过滤后的 exon 数过少")
}
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~0 + group, data = sample_tbl)
colnames(design) <- levels(sample_tbl$group)
y <- edgeR::estimateDisp(y, design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design, robust = TRUE)
contrast <- limma::makeContrasts(contrasts = paste0(group2, "-", group1), levels = design)
qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
ds <- edgeR::diffSpliceDGE(fit, contrast = contrast, geneid = "gene_id", exonid = "exon_id", verbose = TRUE)

gene_simes <- edgeR::topSpliceDGE(ds, test = "Simes", number = Inf, FDR = 1)
gene_f <- edgeR::topSpliceDGE(ds, test = "gene", number = Inf, FDR = 1)
exon_res <- edgeR::topSpliceDGE(ds, test = "exon", number = Inf, FDR = 1)

gene_simes <- gene_simes %>%
  mutate(
    gene_id_clean = sub("\\.[0-9]+$", "", gene_id),
    significant_deu = !is.na(P.Value) & P.Value < 0.05
  ) %>%
  arrange(P.Value, FDR)

gene_f <- gene_f %>%
  mutate(
    gene_id_clean = sub("\\.[0-9]+$", "", gene_id),
    significant_gene_test = !is.na(P.Value) & P.Value < 0.05
  ) %>%
  arrange(P.Value, FDR)

exon_res <- exon_res %>%
  mutate(
    gene_id_clean = sub("\\.[0-9]+$", "", gene_id),
    significant_exon = !is.na(P.Value) & P.Value < 0.05,
    abs_logFC = abs(logFC)
  ) %>%
  arrange(P.Value, desc(abs_logFC))

sig_genes <- gene_simes %>% filter(significant_deu)
sig_exons <- exon_res %>% filter(significant_exon)

gene_burden <- sig_exons %>%
  group_by(gene_id_clean) %>%
  summarise(
    significant_exon_count = n(),
    best_exon_pvalue = min(P.Value, na.rm = TRUE),
    max_abs_exon_logFC = max(abs_logFC, na.rm = TRUE),
    representative_exon = first(exon_id),
    .groups = "drop"
  ) %>%
  left_join(
    gene_simes %>% select(gene_id_clean, P.Value, FDR, NExons),
    by = "gene_id_clean"
  ) %>%
  rename(gene_simes_pvalue = P.Value, gene_simes_fdr = FDR, total_exons = NExons) %>%
  arrange(gene_simes_pvalue, best_exon_pvalue)

write.csv(gene_simes, file.path(out_dir, "diffsplice_gene_simes.csv"), row.names = FALSE)
write.csv(gene_f, file.path(out_dir, "diffsplice_gene_f.csv"), row.names = FALSE)
write.csv(exon_res, file.path(out_dir, "diffsplice_exon.csv"), row.names = FALSE)
write.csv(sig_genes, file.path(out_dir, "significant_genes.csv"), row.names = FALSE)
write.csv(sig_exons, file.path(out_dir, "significant_exons.csv"), row.names = FALSE)
write.csv(gene_burden, file.path(out_dir, "gene_exon_burden.csv"), row.names = FALSE)

summary_tbl <- data.frame(
  category = c("tested_genes", "significant_genes", "tested_exons", "significant_exons"),
  count = c(nrow(gene_simes), nrow(sig_genes), nrow(exon_res), nrow(sig_exons)),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(out_dir, "exon_usage_summary.csv"), row.names = FALSE)

gene_plot <- gene_simes %>%
  filter(is.finite(P.Value)) %>%
  mutate(gene_label = gene_id_clean) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(gene_label, -log10(P.Value + 1e-300)), y = -log10(P.Value + 1e-300), fill = significant_deu)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = color_sig, "FALSE" = color_ns)) +
  labs(title = "Top genes from exon-usage analysis", x = NULL, y = "-log10(PValue)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
save_plot_versions(gene_plot, file.path(out_dir, "exon_usage_top_genes.pdf"), width = 8.5, height = 6)

exon_plot <- exon_res %>%
  filter(is.finite(P.Value), is.finite(logFC)) %>%
  mutate(sig = ifelse(significant_exon, "significant", "ns")) %>%
  ggplot(aes(x = logFC, y = -log10(P.Value + 1e-300), color = sig)) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("significant" = color_sig, "ns" = color_ns)) +
  labs(title = "Exon-level differential usage", x = "Exon logFC", y = "-log10(PValue)") +
  theme_bw(base_size = 11)
save_plot_versions(exon_plot, file.path(out_dir, "exon_usage_exon_scatter.pdf"), width = 8.5, height = 6)
EOF

    write_exon_usage_metadata
}

run_isoform_switch() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        return
    fi
    if ! method_enabled "isoform_switch"; then
        log "跳过 isoform switch，当前方法集: ${RNASEQ_THIRD_METHODS}"
        return
    fi

    discover_gtf
    discover_genome_fasta
    collect_rmats_inputs

    mkdir -p "${ISOFORM_REF_DIR}" "${SALMON_QUANT_DIR}" "${ISOFORM_RESULTS_DIR}" "${AS_INTEGRATED_DIR}"

    if [[ "${FORCE}" != "1" ]] && isoform_outputs_complete && isoform_metadata_matches; then
        log "PostTx-3: 检测到匹配当前 samples.txt 的 isoform switch 输出，跳过重跑。"
        return
    fi

    log "PostTx-3: 运行 isoform switch 分析"
    log "  fastq dir    : ${FASTQC_DIR}"
    log "  genome fasta : ${GENOME_FASTA}"
    log "  transcript   : ${TRANSCRIPT_FASTA}"
    SALMON_BIN="$(resolve_bin "${SALMON_BIN}" salmon)" || die "未找到命令: salmon。请设置 SALMON_BIN 或确认 /home/h1028/miniconda3/bin/salmon 可用。"
    log "  salmon bin   : ${SALMON_BIN}"
    log "  salmon index : ${SALMON_INDEX}"
    log "  results dir  : ${ISOFORM_RESULTS_DIR}"
    log "  libType      : ${SALMON_LIBTYPE}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        while read -r sample_name; do
            [[ -n "${sample_name}" ]] || continue
            log "+ ${SALMON_BIN} quant -i ${SALMON_INDEX} -l ${SALMON_LIBTYPE} -1 $(resolve_clean_fastq "${sample_name}" 1) -2 $(resolve_clean_fastq "${sample_name}" 2) -p ${THREADS} --validateMappings -o ${SALMON_QUANT_DIR}/${sample_name}"
        done < <(read_selected_sample_names)
        log "+ inline R transcript FASTA/tx2gene construction -> ${ISOFORM_REF_DIR}"
        log "+ inline R tximport + IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq -> ${ISOFORM_RESULTS_DIR}"
        return
    fi

    require_cmd Rscript
    require_r_packages "isoform switch" tximport IsoformSwitchAnalyzeR DEXSeq GenomicFeatures Rsamtools Biostrings rtracklayer dplyr readr tibble tidyr

    local sample_name fq1 fq2 quant_dir
    while read -r sample_name; do
        [[ -n "${sample_name}" ]] || continue
        fq1="$(resolve_clean_fastq "${sample_name}" 1)"
        fq2="$(resolve_clean_fastq "${sample_name}" 2)"
        [[ -f "${fq1}" ]] || die "isoform switch 缺少 clean FASTQ: ${fq1}"
        [[ -f "${fq2}" ]] || die "isoform switch 缺少 clean FASTQ: ${fq2}"
    done < <(read_selected_sample_names)

    if [[ "${FORCE}" == "1" || ! -s "${TRANSCRIPT_FASTA}" || ! -s "${TRANSCRIPT_TX2GENE}" ]]; then
        log "  构建 transcript FASTA 与 tx2gene 映射"
        Rscript --vanilla - "${GTF_FILE}" "${GENOME_FASTA}" "${TRANSCRIPT_FASTA}" "${TRANSCRIPT_TX2GENE}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[[1]]
genome_fasta <- args[[2]]
transcript_fasta <- args[[3]]
tx2gene_file <- args[[4]]

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicFeatures)
  library(Rsamtools)
  library(Biostrings)
})

dir.create(dirname(transcript_fasta), recursive = TRUE, showWarnings = FALSE)

gtf_gr <- rtracklayer::import(gtf_file)
exon_gr <- gtf_gr[gtf_gr$type == "exon"]
if (!"transcript_id" %in% colnames(S4Vectors::mcols(exon_gr))) {
  stop("GTF 缺少 transcript_id 注释，无法构建 transcript FASTA")
}
if (!"gene_id" %in% colnames(S4Vectors::mcols(exon_gr))) {
  stop("GTF 缺少 gene_id 注释，无法构建 tx2gene")
}

exon_gr <- exon_gr[!is.na(exon_gr$transcript_id) & exon_gr$transcript_id != "" & !is.na(exon_gr$gene_id) & exon_gr$gene_id != ""]
if (length(exon_gr) == 0) {
  stop("GTF 中没有可用 exon 注释")
}

exon_by_tx <- split(exon_gr, exon_gr$transcript_id)
exon_by_tx <- exon_by_tx[lengths(exon_by_tx) > 0]
exon_by_tx <- endoapply(exon_by_tx, function(gr) {
  if ("exon_number" %in% colnames(S4Vectors::mcols(gr)) && any(!is.na(gr$exon_number))) {
    ord <- order(suppressWarnings(as.numeric(gr$exon_number)), start(gr), end(gr), na.last = TRUE)
  } else {
    ord <- order(start(gr), end(gr))
    if (length(unique(as.character(strand(gr)))) == 1 && unique(as.character(strand(gr))) == "-") {
      ord <- rev(ord)
    }
  }
  gr[ord]
})

if (!file.exists(paste0(genome_fasta, ".fai"))) {
  Rsamtools::indexFa(genome_fasta)
}
dna <- Rsamtools::FaFile(genome_fasta)
open(dna)
on.exit(close(dna), add = TRUE)
tx_seqs <- GenomicFeatures::extractTranscriptSeqs(dna, exon_by_tx)
names(tx_seqs) <- names(exon_by_tx)
Biostrings::writeXStringSet(tx_seqs, filepath = transcript_fasta, format = "fasta")

tx2gene <- do.call(rbind, lapply(seq_along(exon_by_tx), function(i) {
  gr <- exon_by_tx[[i]]
  data.frame(
    isoform_id = names(exon_by_tx)[i],
    gene_id = as.character(gr$gene_id[[1]]),
    gene_name = if ("gene_name" %in% colnames(S4Vectors::mcols(gr))) as.character(gr$gene_name[[1]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}))
tx2gene <- tx2gene[!duplicated(tx2gene$isoform_id), , drop = FALSE]
write.table(tx2gene, tx2gene_file, sep = "\t", quote = FALSE, row.names = FALSE)
EOF
    fi

    if [[ "${FORCE}" == "1" || ! -d "${SALMON_INDEX}" || ! -e "${SALMON_INDEX}/versionInfo.json" ]]; then
        log "  构建 Salmon transcript index"
        mkdir -p "$(dirname "${SALMON_INDEX}")"
        "${SALMON_BIN}" index -t "${TRANSCRIPT_FASTA}" -i "${SALMON_INDEX}" -p "${THREADS}"
    fi

    while read -r sample_name; do
        [[ -n "${sample_name}" ]] || continue
        fq1="$(resolve_clean_fastq "${sample_name}" 1)"
        fq2="$(resolve_clean_fastq "${sample_name}" 2)"
        quant_dir="${SALMON_QUANT_DIR}/${sample_name}"
        if [[ "${FORCE}" != "1" && -s "${quant_dir}/quant.sf" ]]; then
            log "  复用现有 Salmon 定量: ${sample_name}"
            continue
        fi
        rm -rf "${quant_dir}"
        mkdir -p "${quant_dir}"
        run_with_heartbeat \
            "salmon quant ${sample_name}" \
            "" \
            "${SALMON_BIN}" quant \
            -i "${SALMON_INDEX}" \
            -l "${SALMON_LIBTYPE}" \
            -1 "${fq1}" \
            -2 "${fq2}" \
            -p "${THREADS}" \
            --validateMappings \
            -o "${quant_dir}"
    done < <(read_selected_sample_names)

    Rscript --vanilla - "${SAMPLES_FILE}" "${SALMON_QUANT_DIR}" "${TRANSCRIPT_TX2GENE}" "${GTF_FILE}" "${TRANSCRIPT_FASTA}" "${ISOFORM_RESULTS_DIR}" "${RMATS_GROUP1}" "${RMATS_GROUP2}" "${ISOFORM_SWITCH_ALPHA}" "${ISOFORM_SWITCH_DIF_CUTOFF}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
samples_file <- args[[1]]
quant_dir <- args[[2]]
tx2gene_file <- args[[3]]
gtf_file <- args[[4]]
transcript_fasta <- args[[5]]
out_dir <- args[[6]]
group1 <- args[[7]]
group2 <- args[[8]]
alpha <- as.numeric(args[[9]])
dif_cutoff <- as.numeric(args[[10]])

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(tximport)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

quiet_known_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("S4Vectors:::anyMissing", msg, fixed = TRUE) ||
          grepl("Use 'anyNA()' instead.", msg, fixed = TRUE) ||
          grepl("makeTxDbFromGRanges", msg, fixed = TRUE) ||
          grepl("txdbmaker", msg, fixed = TRUE) ||
          grepl("The \"phase\" metadata column contains non-NA values for features of type", msg, fixed = TRUE) ||
          grepl("stop_codon", msg, fixed = TRUE) ||
          grepl("genome version information is not available for this TxDb object", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

sample_tbl <- read.table(
  samples_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
) |>
  dplyr::filter(nzchar(V1), nzchar(V2)) |>
  dplyr::transmute(group = V1, sample = V2)

sample_tbl$group <- factor(sample_tbl$group, levels = c(group1, group2))
if (any(is.na(sample_tbl$group))) {
  stop("samples.txt 中存在不属于当前对比的分组")
}

quant_files <- setNames(file.path(quant_dir, sample_tbl$sample, "quant.sf"), sample_tbl$sample)
if (!all(file.exists(quant_files))) {
  missing_quant <- quant_files[!file.exists(quant_files)]
  stop(sprintf("缺少 Salmon quant.sf: %s", paste(missing_quant, collapse = ", ")))
}

tx2gene <- read.delim(tx2gene_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("isoform_id", "gene_id") %in% colnames(tx2gene))) {
  stop("tx2gene.tsv 缺少 isoform_id/gene_id 列")
}

iso_expr <- IsoformSwitchAnalyzeR::importIsoformExpression(sampleVector = quant_files)
count_df <- iso_expr$counts
abundance_df <- iso_expr$abundance

if (!"isoform_id" %in% colnames(count_df)) {
  count_df <- tibble::rownames_to_column(as.data.frame(count_df), var = "isoform_id")
}
if (!"isoform_id" %in% colnames(abundance_df)) {
  abundance_df <- tibble::rownames_to_column(as.data.frame(abundance_df), var = "isoform_id")
}

design_matrix <- data.frame(
  sampleID = sample_tbl$sample,
  condition = as.character(sample_tbl$group),
  stringsAsFactors = FALSE
)
comparisons_to_make <- data.frame(
  condition_1 = group1,
  condition_2 = group2,
  stringsAsFactors = FALSE
)

switch_list <- quiet_known_warnings(IsoformSwitchAnalyzeR::importRdata(
  isoformCountMatrix = count_df,
  isoformRepExpression = abundance_df,
  designMatrix = design_matrix,
  isoformExonAnnoation = gtf_file,
  isoformNtFasta = transcript_fasta,
  comparisonsToMake = comparisons_to_make,
  addAnnotatedORFs = TRUE,
  ignoreAfterPeriod = TRUE,
  showProgress = FALSE,
  quiet = TRUE
))
switch_list <- quiet_known_warnings(IsoformSwitchAnalyzeR::preFilter(switch_list))
switch_list <- quiet_known_warnings(IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(
  switch_list,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  reduceToSwitchingGenes = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE
))

# Populate alternative splicing classifications before consequence analysis.
switch_list <- quiet_known_warnings(IsoformSwitchAnalyzeR::analyzeAlternativeSplicing(
  switch_list,
  onlySwitchingGenes = FALSE,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  showProgress = FALSE,
  quiet = TRUE
))
switch_list <- quiet_known_warnings(IsoformSwitchAnalyzeR::analyzeIntronRetention(
  switch_list,
  onlySwitchingGenes = FALSE,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  showProgress = FALSE,
  quiet = TRUE
))

consequences_to_use <- c("intron_retention", "ORF_length", "5_utr_length", "3_utr_length", "NMD_status")
switch_list_conseq <- quiet_known_warnings(IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
  switch_list,
  consequencesToAnalyze = consequences_to_use,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  removeNonConseqSwitches = FALSE,
  showProgress = FALSE,
  quiet = TRUE
))

switch_summary <- suppressWarnings(as.data.frame(IsoformSwitchAnalyzeR::extractSwitchSummary(
  switch_list_conseq,
  filterForConsequences = FALSE,
  alpha = alpha,
  dIFcutoff = dif_cutoff
)))
consequence_summary <- suppressWarnings(as.data.frame(IsoformSwitchAnalyzeR::extractConsequenceSummary(
  switch_list_conseq,
  consequencesToAnalyze = consequences_to_use,
  includeCombined = TRUE,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  plot = FALSE,
  returnResult = TRUE
)))
top_switches <- suppressWarnings(as.data.frame(IsoformSwitchAnalyzeR::extractTopSwitches(
  switch_list_conseq,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  filterForConsequences = FALSE,
  sortByQvals = TRUE,
  n = Inf
)))
top_switches_conseq <- suppressWarnings(as.data.frame(IsoformSwitchAnalyzeR::extractTopSwitches(
  switch_list_conseq,
  alpha = alpha,
  dIFcutoff = dif_cutoff,
  filterForConsequences = TRUE,
  sortByQvals = TRUE,
  n = Inf
)))

isoform_features <- as.data.frame(switch_list_conseq$isoformFeatures)
if (!"gene_name" %in% colnames(isoform_features)) {
  isoform_features$gene_name <- NA_character_
}
isoform_features$gene_id_clean <- sub("\\.[0-9]+$", "", isoform_features$gene_id)

isoform_fraction_mat <- as.matrix(abundance_df[, sample_tbl$sample, drop = FALSE])
rownames(isoform_fraction_mat) <- abundance_df$isoform_id
isoform_gene <- tx2gene$gene_id[match(rownames(isoform_fraction_mat), tx2gene$isoform_id)]
valid_iso <- !is.na(isoform_gene) & isoform_gene != ""
isoform_fraction_mat <- isoform_fraction_mat[valid_iso, , drop = FALSE]
isoform_gene <- isoform_gene[valid_iso]
gene_totals <- rowsum(isoform_fraction_mat, group = isoform_gene, reorder = FALSE)
for (i in seq_len(nrow(isoform_fraction_mat))) {
  isoform_fraction_mat[i, ] <- ifelse(gene_totals[isoform_gene[i], ] > 0, isoform_fraction_mat[i, ] / gene_totals[isoform_gene[i], ], 0)
}

isoform_fraction_df <- data.frame(
  isoform_id = rownames(isoform_fraction_mat),
  gene_id = isoform_gene,
  gene_name = tx2gene$gene_name[match(rownames(isoform_fraction_mat), tx2gene$isoform_id)],
  as.data.frame(isoform_fraction_mat, check.names = FALSE),
  check.names = FALSE
)

iso_q_col <- intersect(c("isoform_switch_q_value", "isoform_switch_qval", "iso_q_value"), colnames(isoform_features))
gene_q_col <- intersect(c("gene_switch_q_value", "gene_switch_qval", "gene_q_value"), colnames(isoform_features))
iso_dif_col <- intersect(c("dIF", "dIF_overall"), colnames(isoform_features))
if (length(iso_q_col) == 0 || length(gene_q_col) == 0 || length(iso_dif_col) == 0) {
  stop("isoformFeatures 缺少 isoform switch 关键列，无法整理标准输出")
}
iso_q_col <- iso_q_col[[1]]
gene_q_col <- gene_q_col[[1]]
iso_dif_col <- iso_dif_col[[1]]

significant_switch_isoforms <- isoform_features %>%
  mutate(
    isoform_switch_q_value = .data[[iso_q_col]],
    gene_switch_q_value = .data[[gene_q_col]],
    dIF_value = .data[[iso_dif_col]],
    significant_isoform_switch = is.finite(isoform_switch_q_value) & isoform_switch_q_value < alpha & is.finite(dIF_value) & abs(dIF_value) >= dif_cutoff
  ) %>%
  filter(significant_isoform_switch) %>%
  arrange(isoform_switch_q_value, desc(abs(dIF_value)))

significant_switch_genes <- isoform_features %>%
  mutate(
    isoform_switch_q_value = .data[[iso_q_col]],
    gene_switch_q_value = .data[[gene_q_col]],
    dIF_value = .data[[iso_dif_col]]
  ) %>%
  group_by(gene_id_clean, gene_name) %>%
  summarise(
    gene_switch_q_value = if (all(!is.finite(gene_switch_q_value))) NA_real_ else min(gene_switch_q_value, na.rm = TRUE),
    switching_isoform_count = sum(is.finite(isoform_switch_q_value) & isoform_switch_q_value < alpha & is.finite(dIF_value) & abs(dIF_value) >= dif_cutoff, na.rm = TRUE),
    max_abs_dIF = if (all(!is.finite(dIF_value))) NA_real_ else max(abs(dIF_value), na.rm = TRUE),
    representative_isoform = isoform_id[order(replace(isoform_switch_q_value, !is.finite(isoform_switch_q_value), Inf), -abs(dIF_value))][1],
    .groups = "drop"
  ) %>%
  mutate(significant_switch_gene = is.finite(gene_switch_q_value) & gene_switch_q_value < alpha & is.finite(max_abs_dIF) & max_abs_dIF >= dif_cutoff) %>%
  filter(significant_switch_gene) %>%
  arrange(gene_switch_q_value, desc(max_abs_dIF))

switch_consequence_tbl <- as.data.frame(switch_list_conseq$switchConsequence)
if (nrow(switch_consequence_tbl) > 0 && "gene_id" %in% colnames(switch_consequence_tbl)) {
  switch_consequence_tbl$gene_id_clean <- sub("\\.[0-9]+$", "", switch_consequence_tbl$gene_id)
}
representative_switch_pairs <- if (nrow(switch_consequence_tbl) > 0) {
  pick_first_present <- function(df, candidates, default = NA_character_) {
    existing <- intersect(candidates, colnames(df))
    if (length(existing) == 0) {
      return(rep(default, nrow(df)))
    }
    out <- df[[existing[1]]]
    if (length(existing) > 1) {
      for (nm in existing[-1]) {
        out <- dplyr::coalesce(out, df[[nm]])
      }
    }
    out
  }

  if (!"gene_id_clean" %in% colnames(switch_consequence_tbl)) {
    if ("gene_id" %in% colnames(switch_consequence_tbl)) {
      switch_consequence_tbl$gene_id_clean <- sub("\\.[0-9]+$", "", switch_consequence_tbl$gene_id)
    } else if ("geneID" %in% colnames(switch_consequence_tbl)) {
      switch_consequence_tbl$gene_id_clean <- sub("\\.[0-9]+$", "", switch_consequence_tbl$geneID)
    } else {
      switch_consequence_tbl$gene_id_clean <- NA_character_
    }
  }

  if (!"gene_name" %in% colnames(switch_consequence_tbl)) {
    switch_consequence_tbl$gene_name <- switch_consequence_tbl$gene_id_clean
  }

  switch_consequence_tbl$isoform_up_candidate <- pick_first_present(
    switch_consequence_tbl,
    c("isoformUpregulated", "isoform_1", "isoform_id_1", "upregulatedIsoform", "up_isoform")
  )
  switch_consequence_tbl$isoform_down_candidate <- pick_first_present(
    switch_consequence_tbl,
    c("isoformDownregulated", "isoform_2", "isoform_id_2", "downregulatedIsoform", "down_isoform")
  )
  switch_consequence_tbl$consequence_type_candidate <- pick_first_present(
    switch_consequence_tbl,
    c("switchConsequence", "consequence", "consequenceType"),
    default = ""
  )

  switch_consequence_tbl %>%
    group_by(gene_id_clean, gene_name) %>%
    summarise(
      isoform_up = dplyr::first(isoform_up_candidate),
      isoform_down = dplyr::first(isoform_down_candidate),
      consequence_count = dplyr::n(),
      consequence_types = paste(sort(unique(stats::na.omit(consequence_type_candidate[consequence_type_candidate != ""]))), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(consequence_count), gene_name)
} else {
  data.frame()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(design_matrix, file.path(out_dir, "sample_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = names(quant_files), quant_sf = unname(quant_files), stringsAsFactors = FALSE), file.path(out_dir, "salmon_quant_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(count_df, file.path(out_dir, "tximport_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(abundance_df, file.path(out_dir, "tximport_abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(isoform_fraction_df, file.path(out_dir, "isoform_fraction_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(isoform_features, file.path(out_dir, "isoform_feature_table.csv"), row.names = FALSE)
write.csv(switch_summary, file.path(out_dir, "switch_summary.csv"), row.names = FALSE)
write.csv(top_switches, file.path(out_dir, "top_switches.csv"), row.names = FALSE)
write.csv(top_switches_conseq, file.path(out_dir, "top_switches_with_consequences.csv"), row.names = FALSE)
write.csv(consequence_summary, file.path(out_dir, "consequence_summary.csv"), row.names = FALSE)
write.csv(significant_switch_isoforms, file.path(out_dir, "significant_switch_isoforms.csv"), row.names = FALSE)
write.csv(significant_switch_genes, file.path(out_dir, "significant_switch_genes.csv"), row.names = FALSE)
write.csv(switch_consequence_tbl, file.path(out_dir, "switch_consequence_table.csv"), row.names = FALSE)
write.csv(representative_switch_pairs, file.path(out_dir, "representative_switch_pairs.csv"), row.names = FALSE)
write.csv(as.data.frame(switch_list_conseq$orfAnalysis), file.path(out_dir, "orf_analysis.csv"), row.names = FALSE)
saveRDS(switch_list_conseq, file.path(out_dir, "switchAnalyzeRlist.rds"))
EOF

    write_isoform_metadata
}

run_apa() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        return
    fi
    if ! extra_method_enabled "apa"; then
        log "跳过 APA，当前额外方法集: ${RNASEQ_THIRD_EXTRA_METHODS}"
        return
    fi

    discover_gtf
    discover_genome_fasta
    collect_rmats_inputs

    mkdir -p "${APA_DIR}" "${APA_REF_DIR}" "${AS_INTEGRATED_DIR}"

    if [[ "${FORCE}" != "1" ]] && apa_outputs_complete && apa_metadata_matches; then
        log "PostTx-4: 检测到匹配当前 samples.txt 的 APA 输出，跳过重跑。"
        return
    fi

    local apa_strandtype
    apa_strandtype="$(resolve_apa_strandtype)"

    log "PostTx-4: 运行 Differential APA"
    log "  GTF          : ${GTF_FILE}"
    log "  genome fasta : ${GENOME_FASTA}"
    log "  output dir   : ${APA_DIR}"
    log "  strand type  : ${apa_strandtype}"
    log "  seq type     : ${APA_SEQTYPE}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline APAlyzer::PAS2GEF -> ${APA_REF_DIR}"
        log "+ inline APAlyzer::PASEXP_3UTR/PASEXP_IPA/APAdiff -> ${APA_DIR}"
        return
    fi

    require_r_packages "APA" APAlyzer Rsamtools GenomicAlignments SummarizedExperiment AnnotationDbi org.Mm.eg.db dplyr readr tibble tidyr

    Rscript --vanilla - "${ROOT_DIR}" "${SAMPLES_FILE}" "${GTF_FILE}" "${GENOME_FASTA}" "${APA_DIR}" "${APA_REF_DIR}" "${RMATS_GROUP1}" "${RMATS_GROUP2}" "${apa_strandtype}" "${APA_SEQTYPE}" "${THREADS}" "${APA_PVALUE}" "${APA_DELTA_PAU}" "${APA_TEST_METHOD}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[[1]]
samples_file <- args[[2]]
gtf_file <- args[[3]]
genome_fasta <- args[[4]]
out_dir <- args[[5]]
ref_dir <- args[[6]]
group1 <- args[[7]]
group2 <- args[[8]]
apa_strandtype <- args[[9]]
apa_seqtype <- args[[10]]
threads <- as.integer(args[[11]])
apa_pvalue <- as.numeric(args[[12]])
apa_delta_pau <- as.numeric(args[[13]])
apa_test_method <- args[[14]]

if (!is.finite(threads) || threads <= 0) {
  threads <- 8L
}
if (!is.finite(apa_pvalue) || apa_pvalue <= 0) {
  apa_pvalue <- 0.05
}
if (!is.finite(apa_delta_pau) || apa_delta_pau <= 0) {
  apa_delta_pau <- 0.1
}

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(APAlyzer)
  library(Rsamtools)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

quiet_known_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("S4Vectors:::anyMissing", msg, fixed = TRUE) ||
          grepl("Use 'anyNA()' instead.", msg, fixed = TRUE) ||
          grepl("makeTxDbFromGRanges", msg, fixed = TRUE) ||
          grepl("txdbmaker", msg, fixed = TRUE) ||
          grepl("stop_codon", msg, fixed = TRUE) ||
          grepl("genome version information is not available", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

coerce_coord_cols <- function(df, cols) {
  for (nm in intersect(cols, colnames(df))) {
    df[[nm]] <- suppressWarnings(as.numeric(as.character(df[[nm]])))
  }
  df
}

sample_tbl <- read.table(
  samples_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
) |>
  dplyr::filter(nzchar(V1), nzchar(V2)) |>
  dplyr::transmute(group = V1, sample = V2)

sample_tbl$group <- factor(sample_tbl$group, levels = c(group1, group2))
if (any(is.na(sample_tbl$group))) {
  stop("samples.txt 中存在不属于当前对比的分组")
}

bam_files <- file.path(root_dir, "1.Mapping", paste0(sample_tbl$sample, ".bam"))
if (!all(file.exists(bam_files))) {
  missing_bams <- bam_files[!file.exists(bam_files)]
  stop(sprintf("APA 缺少 BAM 文件: %s", paste(missing_bams, collapse = ", ")))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

sample_manifest <- data.frame(
  sample = sample_tbl$sample,
  group = as.character(sample_tbl$group),
  bam = bam_files,
  stringsAsFactors = FALSE
)
write.table(sample_manifest, file.path(out_dir, "sample_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

bam_inputs <- bam_files
names(bam_inputs) <- sample_tbl$sample
sample_table_apa <- data.frame(
  samplename = sample_tbl$sample,
  condition = as.character(sample_tbl$group),
  stringsAsFactors = FALSE
)

pasref <- quiet_known_warnings(APAlyzer::PAS2GEF(gtf_file, AnnoMethod = "V2"))
if ("dfIPA" %in% names(pasref)) {
  pasref$dfIPA <- coerce_coord_cols(pasref$dfIPA, c("Pos", "upstreamSS", "downstreamSS"))
}
if ("dfLE" %in% names(pasref)) {
  pasref$dfLE <- coerce_coord_cols(pasref$dfLE, c("LEstart", "TES"))
}
if ("refUTRraw" %in% names(pasref)) {
  pasref$refUTRraw <- coerce_coord_cols(pasref$refUTRraw, c("First", "Last", "cdsend"))
}
saveRDS(pasref, file.path(out_dir, "pas_reference.rds"))

utrdb <- quiet_known_warnings(APAlyzer::REF3UTR(pasref$refUTRraw))
utr_raw <- quiet_known_warnings(APAlyzer::PASEXP_3UTR(utrdb, bam_inputs, Strandtype = apa_strandtype))
ipa_raw <- quiet_known_warnings(APAlyzer::PASEXP_IPA(
  pasref$dfIPA,
  pasref$dfLE,
  bam_inputs,
  Strandtype = apa_strandtype,
  nts = threads,
  SeqType = apa_seqtype
))

map_symbol_tbl <- function(symbols) {
  if (length(symbols) == 0) {
    return(data.frame())
  }
  anno <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
  ))
  anno |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::mutate(
      ENSEMBL = dplyr::na_if(ENSEMBL, ""),
      ENTREZID = dplyr::na_if(ENTREZID, "")
    ) |>
    dplyr::arrange(SYMBOL, is.na(ENSEMBL), is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

empty_apa_diff <- function() {
  data.frame(
    gene_symbol = character(),
    RED = numeric(),
    pvalue = numeric(),
    p_adj = numeric(),
    APAreg = character(),
    stringsAsFactors = FALSE
  )
}

empty_standardized_apa <- function() {
  tibble::tibble(
    gene_symbol = character(),
    event_type = character(),
    feature_id = character(),
    ENSEMBL = character(),
    ENTREZID = character(),
    gene_id_clean = character(),
    gene_label = character(),
    event_key = character(),
    delta_PAU = numeric(),
    apa_pvalue = numeric(),
    apa_padj = numeric(),
    significant_apa = logical(),
    apa_direction = character()
  )
}

safe_apadiff <- function(mutiraw, pas_type) {
  mutiraw_df <- as.data.frame(mutiraw, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(mutiraw_df) == 0) {
    message(sprintf("APAdiff 输入为空，%s 将写出空结果表。", pas_type))
    return(empty_apa_diff())
  }
  tryCatch(
    quiet_known_warnings(APAlyzer::APAdiff(
      sampleTable = sample_table_apa,
      mutiraw = mutiraw_df,
      conKET = group1,
      trtKEY = group2,
      PAS = pas_type,
      CUTreads = 0,
      p_adjust_methods = "BH",
      MultiTest = apa_test_method
    )),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("replacement has 1 row, data has 0", msg, fixed = TRUE)) {
        message(sprintf("APAdiff 在 %s 上内部过滤后无可分析事件，写出空结果表。", pas_type))
        return(empty_apa_diff())
      }
      stop(e)
    }
  )
}

utr_diff <- safe_apadiff(utr_raw, "3UTR")
ipa_diff <- safe_apadiff(ipa_raw, "IPA")

utr_gene_map <- map_symbol_tbl(if ("gene_symbol" %in% colnames(utr_diff)) unique(utr_diff$gene_symbol) else character(0))
ipa_gene_map <- map_symbol_tbl(if ("gene_symbol" %in% colnames(ipa_diff)) unique(ipa_diff$gene_symbol) else character(0))
gene_map <- dplyr::bind_rows(utr_gene_map, ipa_gene_map) |>
  dplyr::distinct(SYMBOL, .keep_all = TRUE)
write.table(gene_map, file.path(out_dir, "reference_gene_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

standardize_apa <- function(df, event_type) {
  if (nrow(df) == 0) {
    out <- empty_standardized_apa()
    out$event_type <- event_type
    return(out)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_symbol" %in% colnames(df)) {
    stop(sprintf("%s 结果缺少 gene_symbol 列", event_type))
  }
  if (!"pvalue" %in% colnames(df)) {
    if ("Min_pv" %in% colnames(df)) {
      df$pvalue <- df$Min_pv
    } else {
      stop(sprintf("%s 结果缺少 pvalue/Min_pv 列", event_type))
    }
  }
  if (!"p_adj" %in% colnames(df)) {
    df$p_adj <- p.adjust(df$pvalue, method = "BH")
  }
  if (!"RED" %in% colnames(df)) {
    stop(sprintf("%s 结果缺少 RED 列", event_type))
  }
  df$event_type <- event_type
  df$feature_id <- dplyr::coalesce(
    if ("IPAID" %in% colnames(df)) as.character(df$IPAID) else NA_character_,
    if ("PASid" %in% colnames(df)) as.character(df$PASid) else NA_character_,
    df$gene_symbol
  )
  df <- df |>
    dplyr::left_join(gene_map, by = c("gene_symbol" = "SYMBOL")) |>
    dplyr::mutate(
      gene_id_clean = dplyr::coalesce(ENSEMBL, gene_symbol),
      gene_label = gene_symbol,
      event_key = paste(event_type, gene_id_clean, feature_id, sep = "|"),
      delta_PAU = as.numeric(RED),
      apa_pvalue = as.numeric(pvalue),
      apa_padj = as.numeric(p_adj),
      # 保留当前项目 APA nominal pvalue 口径；更严格汇报时建议同步查看 apa_padj / BH FDR。
      significant_apa = is.finite(apa_pvalue) & apa_pvalue < apa_pvalue_cutoff & is.finite(delta_PAU) & abs(delta_PAU) >= apa_delta_cutoff,
      apa_direction = dplyr::case_when(
        !significant_apa ~ "ns",
        delta_PAU > 0 ~ paste0("higher_in_", group2),
        TRUE ~ paste0("higher_in_", group1)
      )
    ) |>
    dplyr::arrange(apa_pvalue, dplyr::desc(abs(delta_PAU)))
  df
}

apa_pvalue_cutoff <- apa_pvalue
apa_delta_cutoff <- apa_delta_pau
utr_std <- standardize_apa(utr_diff, "3UTR_APA")
ipa_std <- standardize_apa(ipa_diff, "IPA")
apa_all <- dplyr::bind_rows(utr_std, ipa_std) |>
  dplyr::arrange(apa_pvalue, dplyr::desc(abs(delta_PAU)))
apa_sig <- apa_all |>
  dplyr::filter(significant_apa)

if (nrow(apa_all) == 0) {
  apa_gene_summary <- tibble::tibble(
    gene_id_clean = character(),
    gene_label = character(),
    apa_event_count = integer(),
    significant_apa_event_count = integer(),
    apa_event_types = character(),
    best_apa_pvalue = numeric(),
    best_apa_padj = numeric(),
    max_abs_delta_PAU = numeric(),
    representative_feature = character(),
    representative_event_type = character()
  )
  top_apa <- apa_gene_summary
} else {
  apa_gene_summary <- apa_all |>
    dplyr::group_by(gene_id_clean, gene_label) |>
    dplyr::summarise(
      apa_event_count = dplyr::n(),
      significant_apa_event_count = sum(significant_apa, na.rm = TRUE),
      apa_event_types = paste(sort(unique(event_type)), collapse = ";"),
      best_apa_pvalue = if (all(!is.finite(apa_pvalue))) NA_real_ else min(apa_pvalue, na.rm = TRUE),
      best_apa_padj = if (all(!is.finite(apa_padj))) NA_real_ else min(apa_padj, na.rm = TRUE),
      max_abs_delta_PAU = if (all(!is.finite(delta_PAU))) NA_real_ else max(abs(delta_PAU), na.rm = TRUE),
      representative_feature = feature_id[order(replace(apa_pvalue, !is.finite(apa_pvalue), Inf), -abs(delta_PAU))][1],
      representative_event_type = event_type[order(replace(apa_pvalue, !is.finite(apa_pvalue), Inf), -abs(delta_PAU))][1],
      .groups = "drop"
    ) |>
    dplyr::arrange(best_apa_pvalue, dplyr::desc(max_abs_delta_PAU))

  top_apa <- apa_gene_summary |>
    dplyr::slice_head(n = 50)
}

write.table(utr_raw, file.path(out_dir, "3utr_expression_raw.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ipa_raw, file.path(out_dir, "ipa_expression_raw.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(utr_std, file.path(out_dir, "apa_3utr_events.csv"), row.names = FALSE)
write.csv(ipa_std, file.path(out_dir, "apa_ipa_events.csv"), row.names = FALSE)
write.csv(apa_all, file.path(out_dir, "apa_all_events.csv"), row.names = FALSE)
write.csv(apa_sig, file.path(out_dir, "significant_apa_events.csv"), row.names = FALSE)
write.csv(apa_gene_summary, file.path(out_dir, "apa_gene_summary.csv"), row.names = FALSE)
write.csv(top_apa, file.path(out_dir, "top_apa_candidates.csv"), row.names = FALSE)
EOF

    write_apa_metadata
}

run_der() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        return
    fi
    if ! extra_method_enabled "der"; then
        log "跳过 DER，当前额外方法集: ${RNASEQ_THIRD_EXTRA_METHODS}"
        return
    fi

    discover_gtf
    discover_genome_fasta
    collect_rmats_inputs

    if [[ -z "${RMATS_READ_LENGTH}" ]]; then
        if [[ "${DRY_RUN}" == "1" ]]; then
            RMATS_READ_LENGTH="150"
        else
            RMATS_READ_LENGTH="$(infer_read_length)"
        fi
    fi

    mkdir -p "${DER_DIR}" "${AS_INTEGRATED_DIR}"

    if [[ "${FORCE}" != "1" ]] && der_outputs_complete && der_metadata_matches; then
        log "PostTx-5: 检测到匹配当前 samples.txt 的 DER 输出，跳过重跑。"
        return
    fi

    log "PostTx-5: 运行 DER 分析"
    log "  GTF          : ${GTF_FILE}"
    log "  genome fasta : ${GENOME_FASTA}"
    log "  read length  : ${RMATS_READ_LENGTH}"
    log "  output dir   : ${DER_DIR}"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline derfinder::fullCoverage/regionMatrix + limma -> ${DER_DIR}"
        return
    fi

    require_r_packages "DER" derfinder GenomicFeatures GenomicRanges GenomeInfoDb Rsamtools rtracklayer limma AnnotationDbi org.Mm.eg.db dplyr readr tibble

    Rscript --vanilla - "${ROOT_DIR}" "${SAMPLES_FILE}" "${GTF_FILE}" "${DER_DIR}" "${RMATS_GROUP1}" "${RMATS_GROUP2}" "${THREADS}" "${DER_PVALUE}" "${DER_CUTOFF}" "${RMATS_READ_LENGTH}" "${DER_GENOME_STYLE}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[[1]]
samples_file <- args[[2]]
gtf_file <- args[[3]]
out_dir <- args[[4]]
group1 <- args[[5]]
group2 <- args[[6]]
threads <- as.integer(args[[7]])
der_pvalue <- as.numeric(args[[8]])
der_cutoff <- as.numeric(args[[9]])
read_length <- as.numeric(args[[10]])
genome_style <- args[[11]]

if (!is.finite(threads) || threads <= 0) {
  threads <- 8L
}
if (!is.finite(der_pvalue) || der_pvalue <= 0) {
  der_pvalue <- 0.05
}
if (!is.finite(der_cutoff) || der_cutoff <= 0) {
  der_cutoff <- 5
}
if (!is.finite(read_length) || read_length <= 0) {
  read_length <- 150
}

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(derfinder)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Rsamtools)
  library(rtracklayer)
  library(limma)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(tibble)
})

quiet_known_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("S4Vectors:::anyMissing", msg, fixed = TRUE) ||
          grepl("Use 'anyNA()' instead.", msg, fixed = TRUE) ||
          grepl("makeTxDbFromGRanges", msg, fixed = TRUE) ||
          grepl("txdbmaker", msg, fixed = TRUE) ||
          grepl("The \"phase\" metadata column contains non-NA values for features of type", msg, fixed = TRUE) ||
          grepl("stop_codon", msg, fixed = TRUE) ||
          grepl("genome version information is not available for this TxDb object", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

drop_na_seqnames <- function(gr) {
  if (is.null(gr) || length(gr) == 0) {
    return(gr)
  }
  keep <- !is.na(as.character(GenomicRanges::seqnames(gr))) &
    nzchar(as.character(GenomicRanges::seqnames(gr)))
  keep[is.na(keep)] <- FALSE
  gr[keep]
}

sample_tbl <- read.table(
  samples_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
) |>
  dplyr::filter(nzchar(V1), nzchar(V2)) |>
  dplyr::transmute(group = V1, sample = V2)

sample_tbl$group <- factor(sample_tbl$group, levels = c(group1, group2))
if (any(is.na(sample_tbl$group))) {
  stop("samples.txt 中存在不属于当前对比的分组")
}

bam_files <- file.path(root_dir, "1.Mapping", paste0(sample_tbl$sample, ".bam"))
if (!all(file.exists(bam_files))) {
  missing_bams <- bam_files[!file.exists(bam_files)]
  stop(sprintf("DER 缺少 BAM 文件: %s", paste(missing_bams, collapse = ", ")))
}
names(bam_files) <- sample_tbl$sample
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  data.frame(sample = sample_tbl$sample, group = as.character(sample_tbl$group), bam = bam_files, stringsAsFactors = FALSE),
  file.path(out_dir, "sample_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

txdb <- quiet_known_warnings(GenomicFeatures::makeTxDbFromGFF(gtf_file))
gtf_gr <- quiet_known_warnings(rtracklayer::import(gtf_file))
gtf_gr <- drop_na_seqnames(gtf_gr)
gene_annot_tbl <- as.data.frame(S4Vectors::mcols(gtf_gr), stringsAsFactors = FALSE)
if (!"gene_id" %in% colnames(gene_annot_tbl)) {
  stop("GTF 缺少 gene_id 注释，DER 无法做基因注释")
}
if (!"gene_name" %in% colnames(gene_annot_tbl)) {
  gene_annot_tbl$gene_name <- NA_character_
}
gene_annot_tbl <- gene_annot_tbl |>
  dplyr::transmute(gene_id = gene_id, gene_name = gene_name) |>
  dplyr::filter(!is.na(gene_id), gene_id != "") |>
  dplyr::distinct(gene_id, .keep_all = TRUE)

bam_targets <- quiet_known_warnings(names(Rsamtools::scanBamHeader(bam_files[[1]])[[1]]$targets))
bam_targets <- bam_targets[!is.na(bam_targets) & nzchar(bam_targets) & bam_targets != "*"]
txdb_chrs <- GenomeInfoDb::seqlevels(txdb)
txdb_chrs <- txdb_chrs[!is.na(txdb_chrs) & nzchar(txdb_chrs)]
candidate_chrs <- intersect(txdb_chrs, bam_targets)
candidate_chrs <- candidate_chrs[!grepl("(_|random|Un|alt|fix|patch|HLA)", candidate_chrs, ignore.case = TRUE)]
if (length(candidate_chrs) == 0) {
  stop("DER 无法从 GTF 与 BAM 解析公共染色体")
}
if (identical(tolower(genome_style), "auto")) {
  chrs <- candidate_chrs
} else if (tolower(genome_style) == "ucsc") {
  chrs <- candidate_chrs[grepl("^chr", candidate_chrs)]
} else if (tolower(genome_style) == "ncbi") {
  chrs <- candidate_chrs[!grepl("^chr", candidate_chrs)]
} else {
  chrs <- candidate_chrs
}
if (length(chrs) == 0) {
  chrs <- candidate_chrs
}
chrs <- unique(chrs[!is.na(chrs) & nzchar(chrs)])
if (length(chrs) == 0) {
  stop("DER 染色体集合为空")
}

total_mapped <- quiet_known_warnings(vapply(bam_files, derfinder::getTotalMapped, numeric(1)))
target_size <- stats::median(total_mapped, na.rm = TRUE)
mc_cores <- 1L

full_cov <- quiet_known_warnings(derfinder::fullCoverage(
  files = bam_files,
  chrs = as.character(chrs),
  verbose = FALSE,
  mc.cores = mc_cores
))
full_cov <- full_cov[!is.na(names(full_cov)) & nzchar(names(full_cov))]
if (length(full_cov) == 0) {
  stop("DER fullCoverage 未返回有效染色体覆盖结果")
}
region_mat <- quiet_known_warnings(derfinder::regionMatrix(
  fullCov = full_cov,
  cutoff = der_cutoff,
  L = read_length,
  totalMapped = total_mapped,
  targetSize = target_size,
  returnBP = FALSE,
  verbose = FALSE,
  mc.cores = mc_cores
))

coverage_matrix <- do.call(rbind, lapply(region_mat, `[[`, "coverageMatrix"))
regions <- unlist(GenomicRanges::GRangesList(lapply(region_mat, `[[`, "regions")))
valid_region_idx <- !is.na(as.character(GenomicRanges::seqnames(regions))) &
  nzchar(as.character(GenomicRanges::seqnames(regions)))
valid_region_idx[is.na(valid_region_idx)] <- FALSE
regions <- regions[valid_region_idx]
coverage_matrix <- coverage_matrix[valid_region_idx, , drop = FALSE]
if (length(regions) == 0 || is.null(coverage_matrix) || nrow(coverage_matrix) == 0) {
  stop("DER 未产生可分析的 expressed regions")
}

region_ids <- paste0("DER_", seq_along(regions))
names(regions) <- region_ids
rownames(coverage_matrix) <- region_ids
colnames(coverage_matrix) <- sample_tbl$sample

log_cov <- log2(coverage_matrix + 1)
design <- model.matrix(~0 + group, data = sample_tbl)
colnames(design) <- levels(sample_tbl$group)
contrast <- limma::makeContrasts(contrasts = paste0(group2, "-", group1), levels = design)
fit <- limma::lmFit(log_cov, design)
fit <- limma::contrasts.fit(fit, contrast)
fit <- limma::eBayes(fit)
der_tbl <- limma::topTable(fit, number = Inf, sort.by = "P")
der_tbl$region_id <- rownames(der_tbl)

regions_df <- data.frame(
  region_id = region_ids,
  seqnames = as.character(GenomicRanges::seqnames(regions)),
  start = GenomicRanges::start(regions),
  end = GenomicRanges::end(regions),
  width = GenomicRanges::width(regions),
  stringsAsFactors = FALSE
)

genes_gr <- quiet_known_warnings(GenomicFeatures::genes(txdb))
genes_gr <- drop_na_seqnames(genes_gr)
gene_ids <- if ("gene_id" %in% names(S4Vectors::mcols(genes_gr))) as.character(genes_gr$gene_id) else names(genes_gr)
gene_ids[is.na(gene_ids) | gene_ids == ""] <- names(genes_gr)[is.na(gene_ids) | gene_ids == ""]
gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
names(genes_gr) <- gene_ids

exons_gr <- quiet_known_warnings(GenomicFeatures::exons(txdb))
exons_gr <- drop_na_seqnames(exons_gr)
introns_gr <- quiet_known_warnings(unlist(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE), use.names = FALSE))
introns_gr <- drop_na_seqnames(introns_gr)

nearest_idx <- GenomicRanges::distanceToNearest(regions, genes_gr, ignore.strand = TRUE)
nearest_gene <- rep(NA_character_, length(regions))
if (length(nearest_idx) > 0) {
  nearest_gene[S4Vectors::queryHits(nearest_idx)] <- names(genes_gr)[S4Vectors::subjectHits(nearest_idx)]
}

region_class <- rep("intergenic", length(regions))
region_class[countOverlaps(regions, genes_gr, ignore.strand = TRUE) > 0] <- "genic"
region_class[countOverlaps(regions, introns_gr, ignore.strand = TRUE) > 0] <- "intron"
region_class[countOverlaps(regions, exons_gr, ignore.strand = TRUE) > 0] <- "exon"

annot_df <- regions_df |>
  dplyr::mutate(
    gene_id_clean = nearest_gene,
    gene_name = gene_annot_tbl$gene_name[match(nearest_gene, sub("\\.[0-9]+$", "", gene_annot_tbl$gene_id))],
    region_class = region_class
  )

der_results <- der_tbl |>
  dplyr::left_join(annot_df, by = "region_id") |>
  dplyr::mutate(
    significant_der = !is.na(P.Value) & P.Value < der_pvalue,
    der_effect = logFC
  ) |>
  dplyr::arrange(P.Value, dplyr::desc(abs(logFC)))

sig_ders <- der_results |>
  dplyr::filter(significant_der)

der_gene_summary <- der_results |>
  dplyr::filter(!is.na(gene_id_clean), gene_id_clean != "") |>
  dplyr::group_by(gene_id_clean, gene_name) |>
  dplyr::summarise(
    der_region_count = dplyr::n(),
    significant_der_region_count = sum(significant_der, na.rm = TRUE),
    best_der_pvalue = if (all(!is.finite(P.Value))) NA_real_ else min(P.Value, na.rm = TRUE),
    best_der_padj = if (all(!is.finite(adj.P.Val))) NA_real_ else min(adj.P.Val, na.rm = TRUE),
    max_abs_der_effect = if (all(!is.finite(logFC))) NA_real_ else max(abs(logFC), na.rm = TRUE),
    representative_region = region_id[order(replace(P.Value, !is.finite(P.Value), Inf), -abs(logFC))][1],
    major_region_class = region_class[order(replace(P.Value, !is.finite(P.Value), Inf), -abs(logFC))][1],
    .groups = "drop"
  ) |>
  dplyr::arrange(best_der_pvalue, dplyr::desc(max_abs_der_effect))

write.table(coverage_matrix, file.path(out_dir, "coverage_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.csv(der_results, file.path(out_dir, "all_ders.csv"), row.names = FALSE)
write.csv(sig_ders, file.path(out_dir, "significant_ders.csv"), row.names = FALSE)
write.csv(der_results, file.path(out_dir, "annotated_ders.csv"), row.names = FALSE)
write.csv(der_gene_summary, file.path(out_dir, "der_gene_summary.csv"), row.names = FALSE)
write.table(data.frame(chr = chrs, stringsAsFactors = FALSE), file.path(out_dir, "chromosome_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
EOF

    write_der_metadata
}

run_post_tx_summary() {
    if is_true "${RNASEQ_THIRD_SKIP_POST}"; then
        return
    fi

    require_cmd Rscript
    discover_gtf

    local groups=()
    mapfile -t groups < <(get_group_order)
    [[ ${#groups[@]} -eq 2 ]] || die "后转录层整合需要两个分组"

    log "PostTx-6: 汇总与整合后转录层结果"
    log "  rMATS dir: ${RMATS_OUTPUT_DIR}"
    log "  exon dir : ${EXON_USAGE_DIR}"
    log "  isoform  : ${ISOFORM_RESULTS_DIR}"
    log "  APA dir  : ${APA_DIR}"
    log "  DER dir  : ${DER_DIR}"
    log "  integ dir: ${AS_INTEGRATED_DIR}"
    log "  step10 dir: 保留给 rnapost.Rmd 统一生成，rnaseq_third.sh 不直接写图表"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline R post-transcriptional integration from ${RMATS_OUTPUT_DIR} + ${EXON_USAGE_DIR} + ${ISOFORM_RESULTS_DIR} + ${APA_DIR} + ${DER_DIR}"
        return
    fi

    mkdir -p "${AS_INTEGRATED_DIR}"

    Rscript --vanilla - "${RMATS_OUTPUT_DIR}" "${EXON_USAGE_DIR}" "${ISOFORM_RESULTS_DIR}" "${APA_DIR}" "${DER_DIR}" "${AS_INTEGRATED_DIR}" "${DE_RESULTS_FILE}" "${groups[0]}" "${groups[1]}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
rmats_dir <- args[[1]]
exon_usage_dir <- args[[2]]
isoform_dir <- args[[3]]
apa_dir <- args[[4]]
der_dir <- args[[5]]
integrated_dir <- args[[6]]
de_file <- args[[7]]
group1 <- args[[8]]
group2 <- args[[9]]

options(scipen = 999)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

event_types <- c("SE", "A5SS", "A3SS", "RI", "MXE")
count_modes <- c("JC", "JCEC")

split_numeric_string <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") {
    return(numeric(0))
  }
  vals <- trimws(unlist(strsplit(as.character(x), ",")))
  suppressWarnings(vals_num <- as.numeric(vals))
  vals_num[is.finite(vals_num)]
}

mean_numeric_string <- function(x) {
  vals <- split_numeric_string(x)
  if (length(vals) == 0) {
    return(NA_real_)
  }
  mean(vals)
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep]))
}

coalesce_existing_chr <- function(df, candidates) {
  existing <- intersect(candidates, colnames(df))
  if (length(existing) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  out <- as.character(df[[existing[[1]]]])
  if (length(existing) >= 2) {
    for (col_name in existing[-1]) {
      next_vals <- as.character(df[[col_name]])
      missing_idx <- is.na(out) | out == ""
      out[missing_idx] <- next_vals[missing_idx]
    }
  }
  out
}

derive_gene_id_clean <- function(df, candidates) {
  sub("\\.[0-9]+$", "", coalesce_existing_chr(df, candidates))
}

ensure_join_columns <- function(df, cols) {
  if (is.null(df) || !is.data.frame(df)) {
    df <- data.frame(stringsAsFactors = FALSE)
  }
  for (col_name in cols) {
    if (!col_name %in% colnames(df)) {
      fill_value <- if (identical(col_name, "gene_id_clean")) NA_character_ else NA
      df[[col_name]] <- rep(fill_value, nrow(df))
    }
  }
  if ("gene_id_clean" %in% colnames(df)) {
    df$gene_id_clean <- as.character(df$gene_id_clean)
  }
  df[, cols, drop = FALSE]
}

find_rmats_file <- function(base_dir, event_type, count_mode) {
  candidates <- c(
    file.path(base_dir, paste0(event_type, ".MATS.", count_mode, ".txt")),
    file.path(base_dir, paste0("fromGTF.", event_type, ".MATS.", count_mode, ".txt"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    return(NA_character_)
  }
  found[[1]]
}

load_rmats_table <- function(target_file, event_type, count_mode) {
  rmats_df <- read.delim(target_file, check.names = FALSE, stringsAsFactors = FALSE)
  rmats_df$event_type <- event_type
  rmats_df$count_mode <- count_mode
  rmats_df$source_file <- normalizePath(target_file, mustWork = FALSE)
  rmats_df$gene_id_clean <- dplyr::coalesce(
    if ("GeneID" %in% colnames(rmats_df)) sub("\\.[0-9]+$", "", rmats_df$GeneID) else NA_character_,
    if ("geneSymbol" %in% colnames(rmats_df)) sub("\\.[0-9]+$", "", rmats_df$geneSymbol) else NA_character_,
    if ("ID" %in% colnames(rmats_df)) sub("\\.[0-9]+$", "", as.character(rmats_df$ID)) else NA_character_
  )
  rmats_df$gene_label <- dplyr::coalesce(
    if ("geneSymbol" %in% colnames(rmats_df)) rmats_df$geneSymbol else NA_character_,
    if ("GeneID" %in% colnames(rmats_df)) sub("\\.[0-9]+$", "", rmats_df$GeneID) else NA_character_,
    if ("ID" %in% colnames(rmats_df)) as.character(rmats_df$ID) else NA_character_
  )
  rmats_df$event_pvalue <- dplyr::coalesce(
    if ("PValue" %in% colnames(rmats_df)) suppressWarnings(as.numeric(rmats_df$PValue)) else NA_real_,
    if ("pvalue" %in% colnames(rmats_df)) suppressWarnings(as.numeric(rmats_df$pvalue)) else NA_real_
  )
  rmats_df$event_fdr <- dplyr::coalesce(
    if ("FDR" %in% colnames(rmats_df)) suppressWarnings(as.numeric(rmats_df$FDR)) else NA_real_,
    if ("fdr" %in% colnames(rmats_df)) suppressWarnings(as.numeric(rmats_df$fdr)) else NA_real_,
    if ("adj.P.Val" %in% colnames(rmats_df)) suppressWarnings(as.numeric(rmats_df$adj.P.Val)) else NA_real_
  )
  rmats_df$event_pvalue[is.na(rmats_df$event_pvalue)] <- 1
  rmats_df$IncLevelDifference[is.na(rmats_df$IncLevelDifference)] <- 0
  rmats_df$psi_group1_mean <- if ("IncLevel1" %in% colnames(rmats_df)) vapply(rmats_df$IncLevel1, mean_numeric_string, numeric(1)) else NA_real_
  rmats_df$psi_group2_mean <- if ("IncLevel2" %in% colnames(rmats_df)) vapply(rmats_df$IncLevel2, mean_numeric_string, numeric(1)) else NA_real_
  rmats_df$abs_delta_psi <- abs(rmats_df$IncLevelDifference)
  # 保留当前项目 rMATS nominal pvalue < 0.05 + |delta PSI| >= 0.1 口径；更严格汇报时建议同时筛选 rMATS FDR < 0.05。
  rmats_df$significant_as <- rmats_df$event_pvalue < 0.05 & rmats_df$abs_delta_psi >= 0.1
  rmats_df$splice_direction <- dplyr::case_when(
    !rmats_df$significant_as ~ "ns",
    rmats_df$IncLevelDifference > 0 ~ paste0("higher_in_", group2),
    TRUE ~ paste0("higher_in_", group1)
  )

  coord_cols <- intersect(
    c("chr", "strand", "exonStart_0base", "exonEnd", "upstreamEE", "downstreamES",
      "upstreamES", "downstreamEE", "longExonStart_0base", "longExonEnd",
      "shortES", "shortEE", "riExonStart_0base", "riExonEnd",
      "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd"),
    colnames(rmats_df)
  )
  if (length(coord_cols) > 0) {
    rmats_df$event_label <- apply(rmats_df[, coord_cols, drop = FALSE], 1, function(x) {
      paste(x[!is.na(x) & x != ""], collapse = ":")
    })
  } else {
    rmats_df$event_label <- if ("ID" %in% colnames(rmats_df)) as.character(rmats_df$ID) else paste0(event_type, "_", seq_len(nrow(rmats_df)))
  }
  rmats_df$event_key <- paste(rmats_df$event_type, rmats_df$gene_id_clean, rmats_df$event_label, sep = "|")
  rmats_df$event_priority <- -log10(rmats_df$event_pvalue + 1e-300) * pmax(rmats_df$abs_delta_psi, 0.001)
  rmats_df
}

rmats_file_index <- expand.grid(event_type = event_types, count_mode = count_modes, stringsAsFactors = FALSE) %>%
  mutate(
    file_path = mapply(function(one_event, one_mode) find_rmats_file(rmats_dir, one_event, one_mode), event_type, count_mode),
    exists = !is.na(file_path)
  ) %>%
  filter(exists)
if (nrow(rmats_file_index) == 0) {
  stop(sprintf("未在 %s 中检测到 rMATS 结果文件", rmats_dir))
}

rmats_tables <- lapply(seq_len(nrow(rmats_file_index)), function(i) {
  load_rmats_table(rmats_file_index$file_path[i], rmats_file_index$event_type[i], rmats_file_index$count_mode[i])
})
rmats_combined <- bind_rows(rmats_tables)

de_map <- data.frame()
de_sig_ids <- character(0)
if (file.exists(de_file)) {
  de_raw <- read.delim(de_file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("gene_id" %in% colnames(de_raw)) {
    de_raw$gene_id_clean <- sub("\\.[0-9]+$", "", de_raw$gene_id)
    # 保留当前项目 DEG nominal pvalue < 0.05 + |log2FC| >= 0.58 口径；更严格汇报时建议同步查看 padj/FDR。
    de_raw$de_significant <- !is.na(de_raw$pvalue) & de_raw$pvalue < 0.05 & !is.na(de_raw$log2FoldChange) & abs(de_raw$log2FoldChange) >= 0.58
    de_map <- de_raw %>%
      arrange(pvalue, desc(abs(log2FoldChange))) %>%
      distinct(gene_id_clean, .keep_all = TRUE) %>%
      select(gene_id_clean, log2FoldChange, pvalue, padj, de_significant)
    de_sig_ids <- unique(de_map$gene_id_clean[de_map$de_significant])
  }
}

if (nrow(de_map) > 0) {
  rmats_combined <- rmats_combined %>%
    left_join(de_map, by = "gene_id_clean") %>%
    rename(de_log2FoldChange = log2FoldChange, de_pvalue = pvalue, de_padj = padj)
}

exon_gene_tbl <- data.frame()
exon_gene_f_tbl <- data.frame()
exon_exon_tbl <- data.frame()
deu_sig_ids <- character(0)
if (dir.exists(exon_usage_dir)) {
  gene_simes_file <- file.path(exon_usage_dir, "diffsplice_gene_simes.csv")
  gene_f_file <- file.path(exon_usage_dir, "diffsplice_gene_f.csv")
  exon_file <- file.path(exon_usage_dir, "diffsplice_exon.csv")
  if (file.exists(gene_simes_file)) {
    exon_gene_tbl <- read.csv(gene_simes_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "GeneID")), significant_deu = !is.na(P.Value) & P.Value < 0.05)
      }
    deu_sig_ids <- unique(exon_gene_tbl$gene_id_clean[exon_gene_tbl$significant_deu])
  }
  if (file.exists(gene_f_file)) {
    exon_gene_f_tbl <- read.csv(gene_f_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "GeneID")), significant_gene_test = !is.na(P.Value) & P.Value < 0.05)
      }
  }
  if (file.exists(exon_file)) {
    exon_exon_tbl <- read.csv(exon_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "GeneID")), significant_exon = !is.na(P.Value) & P.Value < 0.05)
      }
  }
}

is_gene_tbl <- data.frame()
is_iso_tbl <- data.frame()
is_consequence_tbl <- data.frame()
is_pairs_tbl <- data.frame()
is_ids <- character(0)
if (dir.exists(isoform_dir)) {
  gene_switch_file <- file.path(isoform_dir, "significant_switch_genes.csv")
  iso_switch_file <- file.path(isoform_dir, "significant_switch_isoforms.csv")
  consequence_file <- file.path(isoform_dir, "consequence_summary.csv")
  pairs_file <- file.path(isoform_dir, "representative_switch_pairs.csv")
  if (file.exists(gene_switch_file)) {
    is_gene_tbl <- read.csv(gene_switch_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "gene_name")))
      }
    is_ids <- unique(is_gene_tbl$gene_id_clean[!is.na(is_gene_tbl$gene_id_clean) & is_gene_tbl$gene_id_clean != ""])
  }
  if (file.exists(iso_switch_file)) {
    is_iso_tbl <- read.csv(iso_switch_file, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (file.exists(consequence_file)) {
    is_consequence_tbl <- read.csv(consequence_file, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (file.exists(pairs_file)) {
    is_pairs_tbl <- read.csv(pairs_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "gene_name")))
      }
  }
}

apa_all_tbl <- data.frame()
apa_sig_tbl <- data.frame()
apa_gene_tbl <- data.frame()
apa_ids <- character(0)
if (dir.exists(apa_dir)) {
  apa_all_file <- file.path(apa_dir, "apa_all_events.csv")
  apa_sig_file <- file.path(apa_dir, "significant_apa_events.csv")
  apa_gene_file <- file.path(apa_dir, "apa_gene_summary.csv")
  if (file.exists(apa_all_file)) {
    apa_all_tbl <- read.csv(apa_all_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "ENSEMBL", "gene_id", "gene_symbol", "gene_label")))
      }
  }
  if (file.exists(apa_sig_file)) {
    apa_sig_tbl <- read.csv(apa_sig_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "ENSEMBL", "gene_id", "gene_symbol", "gene_label")))
      }
    apa_ids <- unique(apa_sig_tbl$gene_id_clean[!is.na(apa_sig_tbl$gene_id_clean) & apa_sig_tbl$gene_id_clean != ""])
  } else if (nrow(apa_all_tbl) > 0) {
    apa_sig_tbl <- apa_all_tbl %>%
      # 保留当前项目 APA nominal pvalue < 0.05 + |delta PAU| >= 0.1 口径；更严格汇报时建议同步查看 apa_padj / BH FDR。
      mutate(significant_apa = !is.na(apa_pvalue) & apa_pvalue < 0.05 & !is.na(delta_PAU) & abs(delta_PAU) >= 0.1) %>%
      filter(significant_apa)
    apa_ids <- unique(apa_sig_tbl$gene_id_clean[!is.na(apa_sig_tbl$gene_id_clean) & apa_sig_tbl$gene_id_clean != ""])
  }
  if (file.exists(apa_gene_file)) {
    apa_gene_tbl <- read.csv(apa_gene_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "ENSEMBL", "gene_id", "gene_label", "gene_symbol")))
      }
  }
}

der_all_tbl <- data.frame()
der_sig_tbl <- data.frame()
der_gene_tbl <- data.frame()
der_ids <- character(0)
if (dir.exists(der_dir)) {
  der_all_file <- file.path(der_dir, "all_ders.csv")
  der_sig_file <- file.path(der_dir, "significant_ders.csv")
  der_gene_file <- file.path(der_dir, "der_gene_summary.csv")
  if (file.exists(der_all_file)) {
    der_all_tbl <- read.csv(der_all_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "gene_name")))
      }
  }
  if (file.exists(der_sig_file)) {
    der_sig_tbl <- read.csv(der_sig_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "gene_name")))
      }
    der_ids <- unique(der_sig_tbl$gene_id_clean[!is.na(der_sig_tbl$gene_id_clean) & der_sig_tbl$gene_id_clean != ""])
  } else if (nrow(der_all_tbl) > 0) {
    der_sig_tbl <- der_all_tbl %>%
      mutate(significant_der = !is.na(P.Value) & P.Value < 0.05) %>%
      filter(significant_der)
    der_ids <- unique(der_sig_tbl$gene_id_clean[!is.na(der_sig_tbl$gene_id_clean) & der_sig_tbl$gene_id_clean != ""])
  }
  if (file.exists(der_gene_file)) {
    der_gene_tbl <- read.csv(der_gene_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      {
        mutate(., gene_id_clean = derive_gene_id_clean(., c("gene_id_clean", "gene_id", "gene_name")))
      }
  }
}

preferred_count_mode <- if ("JCEC" %in% unique(rmats_combined$count_mode)) "JCEC" else sort(unique(rmats_combined$count_mode))[1]
rmats_sig <- rmats_combined %>% filter(significant_as) %>% arrange(event_pvalue, desc(abs_delta_psi))

write.csv(rmats_file_index, file.path(integrated_dir, "rmats_file_index.csv"), row.names = FALSE)
write.csv(rmats_combined, file.path(integrated_dir, "rmats_all_events.csv"), row.names = FALSE)
write.csv(rmats_sig, file.path(integrated_dir, "significant_rmats_events.csv"), row.names = FALSE)

rmats_summary <- rmats_combined %>%
  group_by(event_type, count_mode) %>%
  summarise(
    total_events = n(),
    significant_events = sum(significant_as, na.rm = TRUE),
    higher_in_group2 = sum(splice_direction == paste0("higher_in_", group2), na.rm = TRUE),
    higher_in_group1 = sum(splice_direction == paste0("higher_in_", group1), na.rm = TRUE),
    median_abs_delta_psi = ifelse(sum(significant_as, na.rm = TRUE) > 0, median(abs_delta_psi[significant_as], na.rm = TRUE), 0),
    .groups = "drop"
  )
write.csv(rmats_summary, file.path(integrated_dir, "rmats_event_summary.csv"), row.names = FALSE)

direction_summary <- rmats_sig %>%
  group_by(count_mode, event_type, splice_direction) %>%
  summarise(event_count = n(), .groups = "drop")
write.csv(direction_summary, file.path(integrated_dir, "rmats_direction_summary.csv"), row.names = FALSE)

top_events <- rmats_sig %>%
  transmute(event_type, count_mode, gene_id_clean, gene_label, event_label, PValue = event_pvalue, FDR = event_fdr, IncLevelDifference, abs_delta_psi, event_priority, splice_direction, de_log2FoldChange, de_pvalue) %>%
  head(50)
write.csv(top_events, file.path(integrated_dir, "top_rmats_events.csv"), row.names = FALSE)

high_confidence_events <- rmats_sig %>% filter(count_mode == preferred_count_mode)
if (all(c("JC", "JCEC") %in% unique(rmats_combined$count_mode))) {
  jc_events <- rmats_combined %>%
    filter(count_mode == "JC") %>%
    transmute(event_key, event_type, gene_id_clean, gene_label, event_label, jc_pvalue = event_pvalue, jc_fdr = event_fdr, jc_delta_psi = IncLevelDifference, jc_sig = significant_as, jc_direction = splice_direction) %>%
    arrange(jc_pvalue, desc(abs(jc_delta_psi))) %>%
    distinct(event_key, event_type, gene_id_clean, gene_label, event_label, .keep_all = TRUE)
  jcec_events <- rmats_combined %>%
    filter(count_mode == "JCEC") %>%
    transmute(event_key, event_type, gene_id_clean, gene_label, event_label, jcec_pvalue = event_pvalue, jcec_fdr = event_fdr, jcec_delta_psi = IncLevelDifference, jcec_sig = significant_as, jcec_direction = splice_direction) %>%
    arrange(jcec_pvalue, desc(abs(jcec_delta_psi))) %>%
    distinct(event_key, event_type, gene_id_clean, gene_label, event_label, .keep_all = TRUE)
  concordance <- inner_join(jc_events, jcec_events, by = c("event_key", "event_type", "gene_id_clean", "gene_label", "event_label"), relationship = "one-to-one") %>%
    mutate(concordant_direction = jc_direction == jcec_direction, high_confidence = jc_sig & jcec_sig & concordant_direction, delta_psi_abs_diff = abs(jc_delta_psi - jcec_delta_psi))
  concordance_summary <- concordance %>%
    group_by(event_type) %>%
    summarise(shared_events = n(), both_significant = sum(jc_sig & jcec_sig, na.rm = TRUE), high_confidence_events = sum(high_confidence, na.rm = TRUE), delta_psi_cor = safe_cor(jc_delta_psi, jcec_delta_psi), pvalue_cor = safe_cor(-log10(jc_pvalue + 1e-300), -log10(jcec_pvalue + 1e-300)), .groups = "drop")
  write.csv(concordance, file.path(integrated_dir, "jc_jcec_event_concordance.csv"), row.names = FALSE)
  write.csv(concordance_summary, file.path(integrated_dir, "jc_jcec_concordance_summary.csv"), row.names = FALSE)
  high_confidence_events <- concordance %>%
    filter(high_confidence) %>%
    mutate(
      count_mode = "JC_and_JCEC",
      PValue = pmin(jc_pvalue, jcec_pvalue),
      IncLevelDifference = ifelse(abs(jcec_delta_psi) >= abs(jc_delta_psi), jcec_delta_psi, jc_delta_psi),
      abs_delta_psi = abs(IncLevelDifference),
      splice_direction = ifelse(abs(jcec_delta_psi) >= abs(jc_delta_psi), jcec_direction, jc_direction)
    ) %>%
    arrange(PValue, desc(abs_delta_psi))
}
write.csv(high_confidence_events, file.path(integrated_dir, "high_confidence_as_events.csv"), row.names = FALSE)

das_ids <- unique(rmats_sig$gene_id_clean[!is.na(rmats_sig$gene_id_clean) & rmats_sig$gene_id_clean != ""])
overlap_summary <- data.frame(
  category = c("DAS_only", "DE_only", "DE_and_DAS"),
  count = c(length(setdiff(das_ids, de_sig_ids)), length(setdiff(de_sig_ids, das_ids)), length(intersect(das_ids, de_sig_ids))),
  stringsAsFactors = FALSE
)
write.csv(overlap_summary, file.path(integrated_dir, "de_das_overlap_summary.csv"), row.names = FALSE)

deu_summary <- data.frame(
  category = c("DEU_only", "DE_only", "DE_and_DEU"),
  count = c(length(setdiff(deu_sig_ids, de_sig_ids)), length(setdiff(de_sig_ids, deu_sig_ids)), length(intersect(deu_sig_ids, de_sig_ids))),
  stringsAsFactors = FALSE
)
write.csv(deu_summary, file.path(integrated_dir, "de_deu_overlap_summary.csv"), row.names = FALSE)

de_is_summary <- data.frame(
  category = c("IS_only", "DE_only", "DE_and_IS"),
  count = c(length(setdiff(is_ids, de_sig_ids)), length(setdiff(de_sig_ids, is_ids)), length(intersect(is_ids, de_sig_ids))),
  stringsAsFactors = FALSE
)
write.csv(de_is_summary, file.path(integrated_dir, "de_is_overlap_summary.csv"), row.names = FALSE)

de_apa_summary <- data.frame(
  category = c("APA_only", "DE_only", "DE_and_APA"),
  count = c(length(setdiff(apa_ids, de_sig_ids)), length(setdiff(de_sig_ids, apa_ids)), length(intersect(apa_ids, de_sig_ids))),
  stringsAsFactors = FALSE
)
write.csv(de_apa_summary, file.path(integrated_dir, "de_apa_overlap_summary.csv"), row.names = FALSE)

de_der_summary <- data.frame(
  category = c("DER_only", "DE_only", "DE_and_DER"),
  count = c(length(setdiff(der_ids, de_sig_ids)), length(setdiff(de_sig_ids, der_ids)), length(intersect(der_ids, de_sig_ids))),
  stringsAsFactors = FALSE
)
write.csv(de_der_summary, file.path(integrated_dir, "de_der_overlap_summary.csv"), row.names = FALSE)

if (nrow(exon_gene_tbl) > 0) {
  write.csv(exon_gene_tbl, file.path(integrated_dir, "deu_gene_priority_table.csv"), row.names = FALSE)
}
if (nrow(exon_exon_tbl) > 0) {
  write.csv(exon_exon_tbl, file.path(integrated_dir, "deu_exon_table.csv"), row.names = FALSE)
}
if (nrow(apa_all_tbl) > 0) {
  write.csv(apa_all_tbl, file.path(integrated_dir, "apa_event_table.csv"), row.names = FALSE)
}
if (nrow(apa_gene_tbl) > 0) {
  write.csv(apa_gene_tbl, file.path(integrated_dir, "apa_gene_priority_table.csv"), row.names = FALSE)
}
if (nrow(der_all_tbl) > 0) {
  write.csv(der_all_tbl, file.path(integrated_dir, "der_region_table.csv"), row.names = FALSE)
}
if (nrow(der_gene_tbl) > 0) {
  write.csv(der_gene_tbl, file.path(integrated_dir, "der_gene_priority_table.csv"), row.names = FALSE)
}

gene_membership <- data.frame(
  gene_id_clean = sort(unique(c(de_sig_ids, das_ids, deu_sig_ids, is_ids, apa_ids, der_ids))),
  stringsAsFactors = FALSE
) %>%
  mutate(
    DE = gene_id_clean %in% de_sig_ids,
    DAS = gene_id_clean %in% das_ids,
    DEU = gene_id_clean %in% deu_sig_ids,
    IS = gene_id_clean %in% is_ids,
    APA = gene_id_clean %in% apa_ids,
    DER = gene_id_clean %in% der_ids
  )
write.csv(gene_membership, file.path(integrated_dir, "gene_method_membership.csv"), row.names = FALSE)

method_overlap_summary <- gene_membership %>%
  mutate(membership = paste0(ifelse(DE, "DE", ""), ifelse(DAS, "_DAS", ""), ifelse(DEU, "_DEU", ""), ifelse(IS, "_IS", ""), ifelse(APA, "_APA", ""), ifelse(DER, "_DER", ""))) %>%
  mutate(membership = gsub("^_", "", membership), membership = ifelse(membership == "", "none", membership)) %>%
  count(membership, name = "gene_count") %>%
  arrange(desc(gene_count), membership)
write.csv(method_overlap_summary, file.path(integrated_dir, "method_overlap_summary.csv"), row.names = FALSE)

rmats_gene_summary <- rmats_sig %>%
  group_by(gene_id_clean) %>%
  summarise(
    significant_event_count = n(),
    max_abs_delta_psi = max(abs_delta_psi, na.rm = TRUE),
    best_event_pvalue = min(event_pvalue, na.rm = TRUE),
    best_event_fdr = if (all(!is.finite(event_fdr))) NA_real_ else min(event_fdr, na.rm = TRUE),
    representative_gene_label = first(gene_label),
    .groups = "drop"
  )

exon_gene_summary <- exon_exon_tbl %>%
  filter(significant_exon) %>%
  group_by(gene_id_clean) %>%
  summarise(
    significant_exon_count = n(),
    best_exon_pvalue = min(P.Value, na.rm = TRUE),
    max_abs_exon_logFC = max(abs(logFC), na.rm = TRUE),
    representative_exon = first(exon_id),
    .groups = "drop"
  )

isoform_gene_summary <- data.frame()
if (nrow(is_gene_tbl) > 0) {
  isoform_gene_summary <- is_gene_tbl %>%
    transmute(
      gene_id_clean,
      isoform_switch_q_value = dplyr::coalesce(gene_switch_q_value, NA_real_),
      switching_isoform_count = dplyr::coalesce(switching_isoform_count, 0L),
      max_abs_dIF = dplyr::coalesce(max_abs_dIF, 0),
      representative_isoform = dplyr::coalesce(representative_isoform, NA_character_)
    )
}
isoform_gene_summary <- ensure_join_columns(
  isoform_gene_summary,
  c("gene_id_clean", "isoform_switch_q_value", "switching_isoform_count", "max_abs_dIF", "representative_isoform")
)

apa_gene_summary <- data.frame()
if (nrow(apa_gene_tbl) > 0) {
  apa_gene_summary <- apa_gene_tbl %>%
    transmute(
      gene_id_clean,
      apa_event_count = dplyr::coalesce(significant_apa_event_count, apa_event_count, 0L),
      best_apa_pvalue = dplyr::coalesce(best_apa_pvalue, NA_real_),
      best_apa_padj = dplyr::coalesce(best_apa_padj, NA_real_),
      max_abs_delta_PAU = dplyr::coalesce(max_abs_delta_PAU, NA_real_),
      representative_apa_feature = dplyr::coalesce(representative_feature, NA_character_),
      representative_apa_type = dplyr::coalesce(representative_event_type, apa_event_types, NA_character_)
    )
}
apa_gene_summary <- ensure_join_columns(
  apa_gene_summary,
  c("gene_id_clean", "apa_event_count", "best_apa_pvalue", "best_apa_padj", "max_abs_delta_PAU", "representative_apa_feature", "representative_apa_type")
)

der_gene_summary <- data.frame()
if (nrow(der_gene_tbl) > 0) {
  der_gene_summary <- der_gene_tbl %>%
    transmute(
      gene_id_clean,
      der_region_count = dplyr::coalesce(significant_der_region_count, der_region_count, 0L),
      best_der_pvalue = dplyr::coalesce(best_der_pvalue, NA_real_),
      best_der_padj = dplyr::coalesce(best_der_padj, NA_real_),
      max_abs_der_effect = dplyr::coalesce(max_abs_der_effect, NA_real_),
      representative_der_region = dplyr::coalesce(representative_region, NA_character_),
      major_der_region_class = dplyr::coalesce(major_region_class, NA_character_)
    )
}
der_gene_summary <- ensure_join_columns(
  der_gene_summary,
  c("gene_id_clean", "der_region_count", "best_der_pvalue", "best_der_padj", "max_abs_der_effect", "representative_der_region", "major_der_region_class")
)

is_pairs_join <- ensure_join_columns(
  is_pairs_tbl,
  c("gene_id_clean", "isoform_up", "isoform_down", "consequence_count", "consequence_types")
)
exon_gene_join <- ensure_join_columns(
  exon_gene_tbl,
  c("gene_id_clean", "P.Value", "FDR", "NExons")
)
de_map_join <- ensure_join_columns(
  de_map,
  c("gene_id_clean", "log2FoldChange", "pvalue", "padj")
)

gene_evidence <- gene_membership %>%
  left_join(rmats_gene_summary, by = "gene_id_clean") %>%
  left_join(exon_gene_summary, by = "gene_id_clean") %>%
  left_join(isoform_gene_summary, by = "gene_id_clean") %>%
  left_join(apa_gene_summary, by = "gene_id_clean") %>%
  left_join(der_gene_summary, by = "gene_id_clean") %>%
  left_join(is_pairs_join, by = "gene_id_clean") %>%
  left_join(exon_gene_join %>% rename(gene_simes_pvalue = P.Value, gene_simes_fdr = FDR), by = "gene_id_clean") %>%
  left_join(de_map_join %>% rename(de_log2FoldChange = log2FoldChange, de_pvalue = pvalue, de_padj = padj), by = "gene_id_clean") %>%
  mutate(
    method_count = ifelse(DE, 1L, 0L) + ifelse(DAS, 1L, 0L) + ifelse(DEU, 1L, 0L) + ifelse(IS, 1L, 0L) + ifelse(APA, 1L, 0L) + ifelse(DER, 1L, 0L),
    representative_label = dplyr::coalesce(representative_gene_label, gene_id_clean),
    evidence_score = method_count +
      dplyr::coalesce(significant_event_count, 0) / 10 +
      dplyr::coalesce(significant_exon_count, 0) / 10 +
      dplyr::coalesce(switching_isoform_count, 0) / 10 +
      dplyr::coalesce(apa_event_count, 0) / 10 +
      dplyr::coalesce(der_region_count, 0) / 10 +
      dplyr::coalesce(max_abs_dIF, 0) +
      dplyr::coalesce(max_abs_delta_PAU, 0) +
      dplyr::coalesce(max_abs_der_effect, 0)
  ) %>%
  arrange(desc(evidence_score), isoform_switch_q_value, best_apa_pvalue, best_der_pvalue, gene_simes_pvalue, best_event_pvalue, de_pvalue)
write.csv(gene_evidence, file.path(integrated_dir, "as_gene_evidence_table.csv"), row.names = FALSE)
write.csv(gene_evidence %>% slice_head(n = 50), file.path(integrated_dir, "top_as_candidates.csv"), row.names = FALSE)

write.csv(is_gene_tbl, file.path(integrated_dir, "isoform_switch_gene_table.csv"), row.names = FALSE)
write.csv(is_iso_tbl, file.path(integrated_dir, "isoform_switch_isoform_table.csv"), row.names = FALSE)
write.csv(is_consequence_tbl, file.path(integrated_dir, "isoform_switch_consequence_summary.csv"), row.names = FALSE)
write.csv(is_pairs_tbl, file.path(integrated_dir, "isoform_switch_representative_pairs.csv"), row.names = FALSE)
EOF
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
                RMATS_THREADS="$1"
                shift
                ;;
            --skip-de)
                RNASEQ_THIRD_SKIP_DE=true
                shift
                ;;
            --run-de)
                RNASEQ_THIRD_SKIP_DE=false
                shift
                ;;
            --skip-post)
                RNASEQ_THIRD_SKIP_POST=true
                shift
                ;;
            --run-post)
                RNASEQ_THIRD_SKIP_POST=false
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
    log "skip DESeq2     : ${RNASEQ_THIRD_SKIP_DE}"
    log "skip post-tx    : ${RNASEQ_THIRD_SKIP_POST}"

    local step
    for step in "${RUN_ORDER[@]}"; do
        run_step "${step}"
    done

    log "完成: ${RUN_ORDER_LABEL}"
}

main "$@"
