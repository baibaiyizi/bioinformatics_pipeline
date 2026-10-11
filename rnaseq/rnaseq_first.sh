#!/usr/bin/env bash

# Example usage:
#   nohup bash rnaseq_first.sh --steps 2,3,4,5,6,7 > 1.log 2>&1 &
#   bash rnaseq_first.sh --list
#   bash rnaseq_first.sh --dry-run --steps 4,5,6
# Note: use comma/space-separated step ids; ranges like 2-7 are not supported.

set -euo pipefail

# ==============================================================================
# 用户配置区（可用同名环境变量临时覆盖）
# ==============================================================================

# ---------- A. 项目目录、输入与输出 ----------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT_DIR}/scripts/config.sh"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
READS_INFO="${RNASEQ_SAMPLES_TSV}"
FASTQC_DIR="${ROOT_DIR}/fastqc"
MAPPING_DIR="${ROOT_DIR}/1.Mapping"
QUANT_DIR="${ROOT_DIR}/2.Quantification"
MERGE_DIR="${ROOT_DIR}/3.Merge_result"
DE_DIR="${ROOT_DIR}/4.DE_analysis"
RESULT_PCA_DIR="${ROOT_DIR}/result_pca"
RESULT_PCA_SNAPSHOT_DIR="${RESULT_PCA_DIR}/full_sample_inputs"
AUDIT_DIR="${ROOT_DIR}/result/00_audit"
AUDIT_TABLE_DIR="${AUDIT_DIR}/tables"
UPSTREAM_QC_DIR="${AUDIT_DIR}/upstream_qc"
MULTIQC_OUT_DIR="${AUDIT_DIR}/multiqc"

# ---------- B. 参考文件与工具 ----------
REFERENCE_DIR="${REFERENCE_DIR:-${REFERENCE_DIR}}"
HISAT2_INDEX="${HISAT2_INDEX:-}"
GTF_FILE="${GTF_FILE:-}"
FASTP_BIN="${FASTP_BIN:-fastp}"
HISAT2_BIN="${HISAT2_BIN:-hisat2}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-/usr/bin/samtools}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
MULTIQC_BIN="${MULTIQC_BIN:-multiqc}"
RSEQC_INFER_BIN="${RSEQC_INFER_BIN:-infer_experiment.py}"
RSEQC_GENEBODY_BIN="${RSEQC_GENEBODY_BIN:-geneBody_coverage.py}"
RSEQC_REF_BED="${RSEQC_REF_BED:-}"

# ---------- C. 文库方向开关（上游/下游共用） ----------
# 建库方向必须依据实验记录配置。
# 可选值：unstranded、fr-firststrand、fr-secondstrand。
RNASEQ_LIBRARY_STRANDEDNESS="${RNASEQ_LIBRARY_STRANDEDNESS:-unstranded}"

# ---------- D. 运行、资源与缓存开关 ----------
THREADS="${THREADS:-8}"                    # 并行线程数。
FORCE="${FORCE:-0}"                         # 1：忽略已有输出并重跑；链设定变更后建议设为 1。
DRY_RUN="${DRY_RUN:-0}"                     # 1：仅打印计划命令，不执行任何分析。
RNASEQ_SAMPLES_FILE="${RNASEQ_SAMPLES_FILE:-${DE_DIR}/samples.txt}"  # 指定下游选中样本清单。

# ---------- E. 步骤调度 ----------
TOTAL_STEPS=7
# 默认全流程按上游分析的科学顺序执行，并在预先合格样本记录后重建选中样本矩阵。
DEFAULT_STEPS=(1 2 3 4 5 6 7)

declare -A STEP_NAMES=(
    [1]="fastp"
    [2]="mapping"
    [3]="quant"
    [4]="merge_all"
    [5]="eligible_samples"
    [6]="merge_selected"
    [7]="audit"
)

declare -A STEP_TITLES=(
    [1]="fastp QC/trimming"
    [2]="HISAT2 mapping"
    [3]="featureCounts quantification"
    [4]="merge all-sample expression matrices"
    [5]="retain all prespecified eligible samples"
    [6]="merge selected-sample matrices"
    [7]="upstream QC/audit"
)

declare -a SAMPLE_NAMES=()
declare -a SAMPLE_GROUPS=()
declare -a RAW_R1S=()
declare -a RAW_R2S=()
declare -a CLEAN_R1S=()
declare -a CLEAN_R2S=()
declare -a REQUESTED_STEPS=()
declare -a RUN_ORDER=()

SAMPLES_LOADED=0
REQUESTED_STEPS_LABEL=""
RUN_ORDER_LABEL=""

log() {
    echo "[first] $*"
}

die() {
    echo "[first][error] $*" >&2
    exit 1
}

declare -a HISAT2_STRANDNESS_ARGS=()
HISAT2_RNA_STRANDNESS=""
FEATURECOUNTS_STRAND_SPECIFIC=0
RNASEQ_LIBRARY_STRANDEDNESS_LABEL=""
RNASEQ_RMATS_LIBTYPE=""
RNASEQ_APA_STRANDTYPE=""

configure_library_strandedness() {
    local mode="${RNASEQ_LIBRARY_STRANDEDNESS,,}"
    case "${mode}" in
        unstranded|none|non-stranded|non-strand-specific)
            RNASEQ_LIBRARY_STRANDEDNESS="unstranded"
            RNASEQ_LIBRARY_STRANDEDNESS_LABEL="unstranded（链非特异性）"
            HISAT2_STRANDNESS_ARGS=()
            HISAT2_RNA_STRANDNESS="none"
            FEATURECOUNTS_STRAND_SPECIFIC=0
            RNASEQ_RMATS_LIBTYPE="fr-unstranded"
            RNASEQ_APA_STRANDTYPE="NONE"
            ;;
        fr-firststrand|fr-first|rf|firststrand|first-strand)
            RNASEQ_LIBRARY_STRANDEDNESS="fr-firststrand"
            RNASEQ_LIBRARY_STRANDEDNESS_LABEL="fr-firststrand（链特异性，RF）"
            HISAT2_STRANDNESS_ARGS=(--rna-strandness RF)
            HISAT2_RNA_STRANDNESS="RF"
            FEATURECOUNTS_STRAND_SPECIFIC=2
            RNASEQ_RMATS_LIBTYPE="fr-firststrand"
            RNASEQ_APA_STRANDTYPE="invert"
            ;;
        fr-secondstrand|fr-second|fr|secondstrand|second-strand)
            RNASEQ_LIBRARY_STRANDEDNESS="fr-secondstrand"
            RNASEQ_LIBRARY_STRANDEDNESS_LABEL="fr-secondstrand（链特异性，FR）"
            HISAT2_STRANDNESS_ARGS=(--rna-strandness FR)
            HISAT2_RNA_STRANDNESS="FR"
            FEATURECOUNTS_STRAND_SPECIFIC=1
            RNASEQ_RMATS_LIBTYPE="fr-secondstrand"
            RNASEQ_APA_STRANDTYPE="forward"
            ;;
        *)
            die "RNASEQ_LIBRARY_STRANDEDNESS 只支持 unstranded、fr-firststrand 或 fr-secondstrand，收到: ${RNASEQ_LIBRARY_STRANDEDNESS}"
            ;;
    esac
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
        printf '  %d: %-14s %s\n' "${step}" "${STEP_NAMES[${step}]}" "${STEP_TITLES[${step}]}"
    done
}

usage() {
    cat <<EOF
Examples:
  nohup bash ${SCRIPT_NAME} --steps 2,3,4,5,6,7 > 1.log 2>&1 &
  bash ${SCRIPT_NAME} --list
  bash ${SCRIPT_NAME} --dry-run --steps 4,5,6

Usage:
  ${SCRIPT_NAME}
  ${SCRIPT_NAME} --all
  ${SCRIPT_NAME} --list
  ${SCRIPT_NAME} --steps 1,3,7
  ${SCRIPT_NAME} --steps 1 3 7
  ${SCRIPT_NAME} 1 3 7
  ${SCRIPT_NAME} 1,3,7

Options:
  --all              Run the full upstream workflow: $(print_default_order)
  --steps LIST       Run selected steps. LIST can be comma- or space-separated.
  --list             Show available steps.
  --dry-run          Print planned commands without running heavy tools.
  --preflight        只读检查配置、样本表与必要输入。
  --force            Re-run completed sample-level outputs.
  --threads N        Override THREADS for this run.
  -h, --help         Show this help.

Library strandness:
  RNASEQ_LIBRARY_STRANDEDNESS=unstranded|fr-firststrand|fr-secondstrand
  该开关统一派生 HISAT2、featureCounts、rMATS 和 APA 的链方向参数。

Note:
  Use comma/space-separated step ids; ranges like 2-7 are not supported.

Behavior:
  1) Steps are normalized to the scientific order 1-7, not the input order.
  2) Missing hard prerequisites are added automatically. For example, step 5
     needs the full-sample matrix snapshot, and step 6 needs samples.txt plus
     count files.
  3) Full-sample matrices are preserved in:
     ${RESULT_PCA_SNAPSHOT_DIR}
  4) Step 6 rebuilds ${MERGE_DIR} for selected samples and is included in --all.

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
        fastp) printf '1\n' ;;
        mapping) printf '2\n' ;;
        quant|featurecounts) printf '3\n' ;;
        merge_all) printf '4\n' ;;
        eligible_samples) printf '5\n' ;;
        merge_selected) printf '6\n' ;;
        audit) printf '7\n' ;;
        *) return 1 ;;
    esac
}

dedupe_preserve_order() {
    local input=("$@")
    local -a out=()
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

discover_bin() {
    local current="${1:-}"
    shift
    local cmd_name="$1"
    shift || true

    if [[ -n "${current}" && -x "${current}" ]]; then
        printf '%s' "${current}"
        return 0
    fi

    if command -v "${cmd_name}" >/dev/null 2>&1; then
        command -v "${cmd_name}"
        return 0
    fi


    return 1
}

normalize_group() { printf '%s' "$1"; }

resolve_read_path() {
    local hint="$1"
    if [[ "${hint}" = /* ]]; then
        printf '%s\n' "${hint}"
    else
        printf '%s/rawdata/%s\n' "${ROOT_DIR}" "$(basename "${hint}")"
    fi
}

ensure_output_dirs() {
    if [[ "${DRY_RUN}" != "1" ]]; then
        mkdir -p \
            "${FASTQC_DIR}" \
            "${MAPPING_DIR}" \
            "${QUANT_DIR}" \
            "${MERGE_DIR}" \
            "${DE_DIR}" \
            "${RESULT_PCA_SNAPSHOT_DIR}" \
            "${AUDIT_TABLE_DIR}" \
            "${UPSTREAM_QC_DIR}" \
            "${MULTIQC_OUT_DIR}"
    fi
}

discover_reference_assets() {
    [[ -s "${GTF_FILE}" ]] || die "MISSING GTF_FILE: ${GTF_FILE}"
    compgen -G "${HISAT2_INDEX}.1.ht2*" >/dev/null || die "MISSING HISAT2_INDEX: ${HISAT2_INDEX}"
}

ensure_prerequisites() {
    FASTP_BIN="$(discover_bin "${FASTP_BIN}" fastp)" || {
        [[ "${DRY_RUN}" == "1" ]] || die "未找到命令: fastp。可用 FASTP_BIN=/path/to/fastp 覆盖。"
        FASTP_BIN="${FASTP_BIN:-fastp}"
    }
    HISAT2_BIN="$(discover_bin "${HISAT2_BIN}" hisat2)" || {
        [[ "${DRY_RUN}" == "1" ]] || die "未找到命令: hisat2。可用 HISAT2_BIN=/path/to/hisat2 覆盖。"
        HISAT2_BIN="${HISAT2_BIN:-hisat2}"
    }
    SAMTOOLS_BIN="$(discover_bin "${SAMTOOLS_BIN}" samtools)" || {
        [[ "${DRY_RUN}" == "1" ]] || die "未找到命令: samtools。可用 SAMTOOLS_BIN=/path/to/samtools 覆盖。"
        SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
    }
    RSCRIPT_BIN="$(discover_bin "${RSCRIPT_BIN}" Rscript)" || {
        [[ "${DRY_RUN}" == "1" ]] || die "未找到命令: Rscript。可用 RSCRIPT_BIN=/path/to/Rscript 覆盖。"
        RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
    }

    if [[ "${DRY_RUN}" != "1" ]]; then
      "${RSCRIPT_BIN}" --vanilla -e 'suppressPackageStartupMessages(library(DESeq2))' >/dev/null 2>&1 || {
        die "未找到 R 包: DESeq2。请安装 DESeq2 或切换到包含该包的 R 环境。"
      }
    fi

    [[ -f "${READS_INFO}" ]] || die "未找到样本清单: ${READS_INFO}"

    discover_reference_assets
}

load_samples() {
    [[ "${SAMPLES_LOADED}" == "1" ]] && return
    local sample group r1 r2
    local rows
    rows="$("${PYTHON_BIN}" "${ROOT_DIR}/scripts/input_contract.py" rows)" || die "样本表验证失败"
    while IFS='|' read -r sample group r1 r2; do
        SAMPLE_NAMES+=("${sample}"); SAMPLE_GROUPS+=("${group}")
        RAW_R1S+=("${r1}"); RAW_R2S+=("${r2}")
        CLEAN_R1S+=("${FASTQC_DIR}/${sample}_1.fastq.gz")
        CLEAN_R2S+=("${FASTQC_DIR}/${sample}_2.fastq.gz")
    done <<< "${rows}"
    SAMPLES_LOADED=1
}

fastp_complete() {
    local sample="$1"
    [[ -s "${FASTQC_DIR}/${sample}_1.fastq.gz" \
        && -s "${FASTQC_DIR}/${sample}_2.fastq.gz" \
        && -s "${FASTQC_DIR}/${sample}_fastp.html" \
        && -s "${FASTQC_DIR}/${sample}_fastp.json" ]]
}

mapping_complete() {
    local sample="$1"
    [[ -s "${MAPPING_DIR}/${sample}.bam" \
        && -s "${MAPPING_DIR}/${sample}.bam.bai" \
        && -s "${MAPPING_DIR}/${sample}.log" ]]
}

quant_complete() {
    local sample="$1"
    [[ -s "${QUANT_DIR}/${sample}.count" && -s "${QUANT_DIR}/${sample}.log" ]]
}

all_fastp_complete() {
    local sample
    for sample in "${SAMPLE_NAMES[@]}"; do
        fastp_complete "${sample}" || return 1
    done
    return 0
}

all_mapping_complete() {
    local sample
    for sample in "${SAMPLE_NAMES[@]}"; do
        mapping_complete "${sample}" || return 1
    done
    return 0
}

all_quant_complete() {
    local sample
    for sample in "${SAMPLE_NAMES[@]}"; do
        quant_complete "${sample}" || return 1
    done
    return 0
}

full_snapshot_complete() {
    [[ -s "${RESULT_PCA_SNAPSHOT_DIR}/genes.counts.matrix" \
        && ( -s "${RESULT_PCA_SNAPSHOT_DIR}/genes.counts.logCPM.matrix" \
             || -s "${RESULT_PCA_SNAPSHOT_DIR}/genes.DESeq2.normalized_counts.matrix" ) ]]
}

selected_samples_complete() {
    [[ -s "${RNASEQ_SAMPLES_FILE}" ]]
}

add_step_with_deps() {
    local step="$1"
    case "${step}" in
        1)
            add_once 1
            ;;
        2)
            if ! all_fastp_complete || step_scheduled 1; then
                add_step_with_deps 1
            fi
            add_once 2
            ;;
        3)
            if ! all_mapping_complete || step_scheduled 1 || step_scheduled 2; then
                add_step_with_deps 2
            fi
            add_once 3
            ;;
        4)
            if ! all_quant_complete || step_scheduled 1 || step_scheduled 2 || step_scheduled 3; then
                add_step_with_deps 3
            fi
            add_once 4
            ;;
        5)
            if ! full_snapshot_complete || step_scheduled 1 || step_scheduled 2 || step_scheduled 3 || step_scheduled 4; then
                add_step_with_deps 4
            fi
            add_once 5
            ;;
        6)
            if ! selected_samples_complete || step_scheduled 1 || step_scheduled 2 || step_scheduled 3 || step_scheduled 4 || step_scheduled 5; then
                add_step_with_deps 5
            fi
            if ! all_quant_complete || step_scheduled 1 || step_scheduled 2 || step_scheduled 3; then
                add_step_with_deps 3
            fi
            add_once 6
            ;;
        7)
            add_once 7
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
        1) run_fastp_stage ;;
        2) run_mapping_stage ;;
        3) run_quant_stage ;;
        4) run_merge_all_stage ;;
        5) run_eligible_samples_stage ;;
        6) run_merge_selected_stage ;;
        7) run_audit_stage ;;
        *) die "内部错误: 不支持的步骤 ID ${step}" ;;
    esac
}

run_fastp_stage() {
    log "Step 1/7: fastp 质控与修剪"

    for idx in "${!SAMPLE_NAMES[@]}"; do
        local sample="${SAMPLE_NAMES[$idx]}"
        local raw_r1="${RAW_R1S[$idx]}"
        local raw_r2="${RAW_R2S[$idx]}"
        local clean_r1="${CLEAN_R1S[$idx]}"
        local clean_r2="${CLEAN_R2S[$idx]}"
        local html="${FASTQC_DIR}/${sample}_fastp.html"
        local json="${FASTQC_DIR}/${sample}_fastp.json"

        if [[ "${FORCE}" != "1" ]] && fastp_complete "${sample}"; then
            log "skip fastp: ${sample}"
            continue
        fi

        log "fastp: ${sample}"

        if [[ ! -f "${raw_r1}" || ! -f "${raw_r2}" ]]; then
            if [[ "${DRY_RUN}" == "1" ]]; then
                log "warn: 原始 FASTQ 缺失，dry-run 继续: ${sample}"
            else
                [[ -f "${raw_r1}" ]] || die "未找到原始 read1: ${raw_r1}"
                [[ -f "${raw_r2}" ]] || die "未找到原始 read2: ${raw_r2}"
            fi
        fi

        if [[ "${DRY_RUN}" == "1" ]]; then
            log "+ fastp -i ${raw_r1} -I ${raw_r2} -o ${clean_r1} -O ${clean_r2} -h ${html} -j ${json}"
            continue
        fi

        "${FASTP_BIN}" \
            -i "${raw_r1}" \
            -I "${raw_r2}" \
            -o "${clean_r1}" \
            -O "${clean_r2}" \
            -h "${html}" \
            -j "${json}" \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --length_required 36 \
            --thread "${THREADS}"
    done
}

run_mapping_stage() {
    log "Step 2/7: hisat2 比对 + BAM 排序/建索引"
    log "mapping mode  : serial (one sample at a time)"

    for idx in "${!SAMPLE_NAMES[@]}"; do
        local sample="${SAMPLE_NAMES[$idx]}"
        local clean_r1="${CLEAN_R1S[$idx]}"
        local clean_r2="${CLEAN_R2S[$idx]}"
        local bam="${MAPPING_DIR}/${sample}.bam"
        local log_file="${MAPPING_DIR}/${sample}.log"

        if [[ "${FORCE}" != "1" ]] && mapping_complete "${sample}"; then
            log "skip mapping: ${sample}"
            continue
        fi

        log "mapping: ${sample}"

        if [[ "${DRY_RUN}" == "1" ]]; then
            log "+ ${HISAT2_BIN} --new-summary ${HISAT2_STRANDNESS_ARGS[*]} -p ${THREADS} -x ${HISAT2_INDEX} -1 ${clean_r1} -2 ${clean_r2} | ${SAMTOOLS_BIN} sort -@ ${THREADS} -o ${bam} -"
            log "+ ${SAMTOOLS_BIN} index -@ ${THREADS} ${bam}"
            continue
        fi

        [[ -s "${clean_r1}" && -s "${clean_r2}" ]] || die "缺少清洗后的 FASTQ，请先完成 fastp: ${sample}"

        "${HISAT2_BIN}" \
            --new-summary \
            "${HISAT2_STRANDNESS_ARGS[@]}" \
            -p "${THREADS}" \
            -x "${HISAT2_INDEX}" \
            -1 "${clean_r1}" \
            -2 "${clean_r2}" \
            2> "${log_file}" \
            | "${SAMTOOLS_BIN}" sort -@ "${THREADS}" -o "${bam}" -

        "${SAMTOOLS_BIN}" index -@ "${THREADS}" "${bam}"
        [[ -s "${bam}.bai" ]] || die "BAM 建索引失败: ${bam}"
    done
}

run_quant_stage() {
    log "Step 3/7: featureCounts 定量"

    for idx in "${!SAMPLE_NAMES[@]}"; do
        local sample="${SAMPLE_NAMES[$idx]}"
        local bam="${MAPPING_DIR}/${sample}.bam"
        local prefix="${QUANT_DIR}/${sample}"

        if [[ "${FORCE}" != "1" ]] && quant_complete "${sample}"; then
            log "skip quantification: ${sample}"
            continue
        fi

        log "quantification: ${sample}"

        if [[ "${DRY_RUN}" == "1" ]]; then
            log "+ inline Rsubread::featureCounts strandSpecific=${FEATURECOUNTS_STRAND_SPECIFIC} -> ${prefix}.count / ${prefix}.log"
            continue
        fi

        [[ -s "${bam}" ]] || die "缺少 BAM 文件，请先完成比对: ${bam}"

        "${RSCRIPT_BIN}" --vanilla - "${bam}" "${GTF_FILE}" "${prefix}" "${THREADS}" "${FEATURECOUNTS_STRAND_SPECIFIC}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
bam <- args[[1]]
gtf <- args[[2]]
prefix <- args[[3]]
threads <- as.integer(args[[4]])
strand_specific <- as.integer(args[[5]])

options(scipen = 999)
suppressPackageStartupMessages(library(Rsubread))

fc <- featureCounts(
  files = bam,
  annot.ext = gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  strandSpecific = strand_specific,
  nthreads = threads
)

counts <- as.numeric(fc$counts[, 1])
lengths <- as.numeric(fc$annotation$Length)
lengths[lengths <= 0 | is.na(lengths)] <- 1
lib_size <- sum(counts, na.rm = TRUE)

if (lib_size > 0) {
  cpm <- counts / lib_size * 1e6
  fpkm <- counts * 1e9 / (lengths * lib_size)
  rpk <- counts / (lengths / 1000)
  tpm_scale <- sum(rpk, na.rm = TRUE)
  if (tpm_scale > 0) {
    tpm <- rpk / tpm_scale * 1e6
  } else {
    tpm <- rep(0, length(counts))
  }
} else {
  cpm <- rep(0, length(counts))
  fpkm <- rep(0, length(counts))
  tpm <- rep(0, length(counts))
}

stats_df <- as.data.frame(fc$stat, stringsAsFactors = FALSE)
write.table(
  stats_df,
  file = paste0(prefix, ".log"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

result_df <- data.frame(
  gene_id = fc$annotation[, 1],
  counts = counts,
  fpkm = fpkm,
  tpm = tpm,
  cpm = cpm,
  check.names = FALSE
)

write.table(
  result_df,
  file = paste0(prefix, ".count"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
EOF
    done
}

load_selected_samples() {
    local samples_file="$1"
    local -n groups_ref="$2"
    local -n samples_ref="$3"

    groups_ref=()
    samples_ref=()

    [[ -f "${samples_file}" ]] || die "未找到样本文件: ${samples_file}"

    while read -r group sample; do
        [[ -z "${group:-}" || -z "${sample:-}" ]] && continue
        [[ "${group}" =~ ^# ]] && continue
        groups_ref+=("${group}")
        samples_ref+=("${sample}")
    done < "${samples_file}"

    [[ ${#samples_ref[@]} -gt 0 ]] || die "样本文件为空: ${samples_file}"
}

write_quant_file_list_all() {
    local out_file="$1"
    : > "${out_file}"

    for sample in "${SAMPLE_NAMES[@]}"; do
        local count_file="${QUANT_DIR}/${sample}.count"
        [[ -f "${count_file}" ]] || die "缺少 count 文件: ${count_file}"
        echo "${count_file}" >> "${out_file}"
    done
}

write_quant_file_list_from_samples_file() {
    local samples_file="$1"
    local out_file="$2"
    local selected_groups=()
    local selected_samples=()

    load_selected_samples "${samples_file}" selected_groups selected_samples

    : > "${out_file}"
    for sample in "${selected_samples[@]}"; do
        local count_file="${QUANT_DIR}/${sample}.count"
        [[ -f "${count_file}" ]] || die "缺少 count 文件: ${count_file}"
        echo "${count_file}" >> "${out_file}"
    done
}

run_merge_from_quant_list() {
    local quant_list_file="$1"
    local merge_label="$2"

    log "${merge_label}: 合并表达矩阵"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline merge/DESeq2 from quant file list -> ${MERGE_DIR}"
        return
    fi

    [[ -s "${quant_list_file}" ]] || die "quant file list 为空: ${quant_list_file}"

    "${RSCRIPT_BIN}" --vanilla - "${quant_list_file}" "${MERGE_DIR}" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
quant_list_file <- args[[1]]
merge_dir <- args[[2]]

options(scipen = 999)
suppressPackageStartupMessages(library(DESeq2))

quant_files <- readLines(quant_list_file, warn = FALSE)
quant_files <- quant_files[nzchar(quant_files)]
if (length(quant_files) == 0) {
  stop("No quant files listed")
}

quant_tables <- lapply(quant_files, function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
})

sample_names <- sub("\\.count$", "", basename(quant_files))
gene_ids <- quant_tables[[1]]$gene_id

for (tbl in quant_tables) {
  if (!identical(tbl$gene_id, gene_ids)) {
    stop("Gene order is inconsistent across .count files")
  }
}

extract_metric <- function(metric_name) {
  mat <- do.call(cbind, lapply(quant_tables, function(tbl) tbl[[metric_name]]))
  rownames(mat) <- gene_ids
  colnames(mat) <- sample_names
  mat
}

counts_mat <- extract_metric("counts")
tpm_mat <- extract_metric("tpm")

write.table(
  counts_mat,
  file = file.path(merge_dir, "genes.counts.matrix"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  tpm_mat,
  file = file.path(merge_dir, "genes.TPM.matrix"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

sample_tbl <- data.frame(sample = colnames(counts_mat), row.names = colnames(counts_mat), check.names = FALSE)
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData = sample_tbl,
  design = ~ 1
)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

count_lib_size <- colSums(counts_mat)
count_lib_size[count_lib_size <= 0 | is.na(count_lib_size)] <- 1
count_cpm <- sweep(counts_mat, 2, count_lib_size, "/") * 1e6
count_cpm[!is.finite(count_cpm)] <- 0
count_logcpm <- log2(count_cpm + 1)

write.table(
  data.frame(sample = names(size_factors), size_factor = size_factors, check.names = FALSE),
  file = file.path(merge_dir, "genes.DESeq2.size_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  count_cpm,
  file = file.path(merge_dir, "genes.counts.CPM.matrix"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  count_logcpm,
  file = file.path(merge_dir, "genes.counts.logCPM.matrix"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# DESeq2 size-factor normalized counts for downstream expression/PCA modules.
norm_expr_fmt <- apply(norm_counts, 2, function(x) sprintf("%.3f", x))
rownames(norm_expr_fmt) <- rownames(norm_counts)
colnames(norm_expr_fmt) <- colnames(norm_counts)

write.table(
  norm_expr_fmt,
  file = file.path(merge_dir, "genes.DESeq2.normalized_counts.matrix"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
EOF

    log "合并完成:"
    log "  - ${MERGE_DIR}/genes.counts.matrix"
    log "  - ${MERGE_DIR}/genes.TPM.matrix"
    log "  - ${MERGE_DIR}/genes.counts.CPM.matrix"
    log "  - ${MERGE_DIR}/genes.counts.logCPM.matrix"
    log "  - ${MERGE_DIR}/genes.DESeq2.normalized_counts.matrix"
}

refresh_full_sample_snapshot() {
    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ 保存全样本矩阵快照到 ${RESULT_PCA_SNAPSHOT_DIR}"
        return
    fi

    mkdir -p "${RESULT_PCA_SNAPSHOT_DIR}"
    cp -f "${MERGE_DIR}/genes.counts.matrix" "${RESULT_PCA_SNAPSHOT_DIR}/genes.counts.matrix"
    cp -f "${MERGE_DIR}/genes.TPM.matrix" "${RESULT_PCA_SNAPSHOT_DIR}/genes.TPM.matrix"
    cp -f "${MERGE_DIR}/genes.DESeq2.normalized_counts.matrix" "${RESULT_PCA_SNAPSHOT_DIR}/genes.DESeq2.normalized_counts.matrix"
    [[ -f "${MERGE_DIR}/genes.counts.CPM.matrix" ]] && cp -f "${MERGE_DIR}/genes.counts.CPM.matrix" "${RESULT_PCA_SNAPSHOT_DIR}/genes.counts.CPM.matrix"
    [[ -f "${MERGE_DIR}/genes.counts.logCPM.matrix" ]] && cp -f "${MERGE_DIR}/genes.counts.logCPM.matrix" "${RESULT_PCA_SNAPSHOT_DIR}/genes.counts.logCPM.matrix"
}

run_eligible_samples_stage() {
    log "Step 5/7: 保留样本表全部预先合格样本；PCA 仅用于诊断"
    [[ "${DRY_RUN}" == "1" ]] && return 0
    "${PYTHON_BIN}" "${ROOT_DIR}/scripts/input_contract.py" prepare
}

run_merge_all_stage() {
    local quant_list_file="${MERGE_DIR}/genes.quant_files.txt"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "Step 4/7: 合并全样本表达矩阵"
        log "+ inline merge from all sample count files -> ${MERGE_DIR}"
        refresh_full_sample_snapshot
        return
    fi

    mkdir -p "${MERGE_DIR}"
    write_quant_file_list_all "${quant_list_file}"
    run_merge_from_quant_list "${quant_list_file}" "Step 4/7"
    refresh_full_sample_snapshot
}

run_merge_selected_stage() {
    local quant_list_file="${MERGE_DIR}/genes.quant_files.txt"
    local samples_file="${RNASEQ_SAMPLES_FILE}"

    log "Step 6/7: 根据 ${samples_file} 重建选中样本矩阵"

    if [[ "${DRY_RUN}" == "1" ]]; then
        log "+ inline merge from selected sample count files listed in ${samples_file}"
        return
    fi

    mkdir -p "${MERGE_DIR}"
    write_quant_file_list_from_samples_file "${samples_file}" "${quant_list_file}"
    run_merge_from_quant_list "${quant_list_file}" "Step 6/7"
}

write_run_provenance() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${AUDIT_TABLE_DIR}"

    local software_file="${AUDIT_TABLE_DIR}/software_versions.tsv"
    {
        printf 'tool\tpath\tversion\n'
      for tool in fastp hisat2 samtools Rscript DESeq2 multiqc infer_experiment.py geneBody_coverage.py salmon rmats.py; do
        local tool_path="NA"
        local tool_version="NA"
        case "${tool}" in
          fastp) tool_path="${FASTP_BIN}" ;;
          hisat2) tool_path="${HISAT2_BIN}" ;;
          samtools) tool_path="${SAMTOOLS_BIN}" ;;
          Rscript|DESeq2) tool_path="${RSCRIPT_BIN}" ;;
          multiqc) tool_path="$(command -v "${MULTIQC_BIN}" 2>/dev/null || true)" ;;
          infer_experiment.py) tool_path="$(command -v "${RSEQC_INFER_BIN}" 2>/dev/null || true)" ;;
          geneBody_coverage.py) tool_path="$(command -v "${RSEQC_GENEBODY_BIN}" 2>/dev/null || true)" ;;
          salmon) tool_path="$(command -v salmon 2>/dev/null || true)" ;;
          rmats.py) tool_path="$(command -v rmats.py 2>/dev/null || true)" ;;
        esac
        [[ -n "${tool_path}" ]] || tool_path="NA"
        if [[ "${tool_path}" != "NA" && -x "${tool_path}" ]]; then
          case "${tool}" in
            fastp|hisat2|samtools|Rscript|multiqc|infer_experiment.py|geneBody_coverage.py|salmon|rmats.py)
              tool_version="$(${tool_path} --version 2>&1 | head -1 || true)"
              ;;
            DESeq2)
              tool_version="$(${tool_path} --vanilla -e 'cat(as.character(packageVersion("DESeq2")))' 2>&1 | head -1 || true)"
              ;;
          esac
        fi
        printf '%s\t%s\t%s\n' "${tool}" "${tool_path}" "${tool_version:-NA}"
      done
    } > "${software_file}"

    {
        printf 'parameter\tvalue\n'
        printf 'root_dir\t%s\n' "${ROOT_DIR}"
        printf 'reads_info\t%s\n' "${READS_INFO}"
        printf 'reference_dir\t%s\n' "${REFERENCE_DIR}"
        printf 'hisat2_index\t%s\n' "${HISAT2_INDEX}"
        printf 'gtf_file\t%s\n' "${GTF_FILE}"
        printf 'threads\t%s\n' "${THREADS}"
        printf 'force\t%s\n' "${FORCE}"
        printf 'requested_steps\t%s\n' "${REQUESTED_STEPS_LABEL:-NA}"
        printf 'run_order\t%s\n' "${RUN_ORDER_LABEL:-NA}"
        printf 'library_strandedness\t%s\n' "${RNASEQ_LIBRARY_STRANDEDNESS}"
        printf 'library_strandedness_label\t%s\n' "${RNASEQ_LIBRARY_STRANDEDNESS_LABEL}"
        printf 'hisat2_strandness\t%s\n' "${HISAT2_RNA_STRANDNESS:-none}"
        printf 'featurecounts_strandSpecific\t%s\n' "${FEATURECOUNTS_STRAND_SPECIFIC}"
        printf 'rmats_libtype\t%s\n' "${RNASEQ_RMATS_LIBTYPE}"
        printf 'apa_strandtype\t%s\n' "${RNASEQ_APA_STRANDTYPE}"
        printf 'rseqc_ref_bed\t%s\n' "${RSEQC_REF_BED:-NA}"
    } > "${AUDIT_TABLE_DIR}/run_params.tsv"

    {
        printf 'asset\tpath\tsha256\n'
        for asset in "${GTF_FILE}" "${GENOME_FASTA:-}" "${HISAT2_INDEX}.1.ht2" "${HISAT2_INDEX}.1.ht2l"; do
            [[ -n "${asset}" && -f "${asset}" ]] || continue
            printf '%s\t%s\t%s\n' "$(basename "${asset}")" "${asset}" "$(sha256sum "${asset}" | awk '{print $1}')"
        done
    } > "${AUDIT_TABLE_DIR}/reference_manifest.tsv"

    {
        printf 'component\tsetting\tnote\n'
        printf 'library\t%s\t%s\n' "${RNASEQ_LIBRARY_STRANDEDNESS}" "${RNASEQ_LIBRARY_STRANDEDNESS_LABEL}"
        printf 'HISAT2\t%s\tderived from RNASEQ_LIBRARY_STRANDEDNESS\n' "${HISAT2_RNA_STRANDNESS:-no --rna-strandness}"
        printf 'featureCounts\tstrandSpecific=%s\tderived from RNASEQ_LIBRARY_STRANDEDNESS\n' "${FEATURECOUNTS_STRAND_SPECIFIC}"
        printf 'rMATS\t%s\tderived downstream setting\n' "${RNASEQ_RMATS_LIBTYPE}"
        printf 'APA\t%s\tderived downstream setting\n' "${RNASEQ_APA_STRANDTYPE}"
        printf 'validation\tRSeQC infer_experiment.py\tset RSEQC_REF_BED to record empirical strandness evidence\n'
    } > "${AUDIT_TABLE_DIR}/strandness_assumptions.tsv"

    "${RSCRIPT_BIN}" --vanilla -e "writeLines(capture.output(sessionInfo()), '${AUDIT_TABLE_DIR}/R_sessionInfo_from_first_stage.txt')" >/dev/null 2>&1 || true
}

run_alignment_qc_tools() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${UPSTREAM_QC_DIR}"

    local sample bam sample_qc_dir
    for sample in "${SAMPLE_NAMES[@]}"; do
        bam="${MAPPING_DIR}/${sample}.bam"
        sample_qc_dir="${UPSTREAM_QC_DIR}/${sample}"
        mkdir -p "${sample_qc_dir}"
        [[ -s "${bam}" ]] || continue

        if [[ "${FORCE}" == "1" || ! -s "${sample_qc_dir}/${sample}.flagstat.txt" ]]; then
            "${SAMTOOLS_BIN}" flagstat -@ "${THREADS}" "${bam}" > "${sample_qc_dir}/${sample}.flagstat.txt" || true
        fi
        if [[ "${FORCE}" == "1" || ! -s "${sample_qc_dir}/${sample}.idxstats.txt" ]]; then
            "${SAMTOOLS_BIN}" idxstats "${bam}" > "${sample_qc_dir}/${sample}.idxstats.txt" || true
        fi
        if [[ "${FORCE}" == "1" || ! -s "${sample_qc_dir}/${sample}.samtools_stats.txt" ]]; then
            "${SAMTOOLS_BIN}" stats -@ "${THREADS}" "${bam}" > "${sample_qc_dir}/${sample}.samtools_stats.txt" || true
        fi

        if [[ -n "${RSEQC_REF_BED}" && -f "${RSEQC_REF_BED}" ]] && command -v "${RSEQC_INFER_BIN}" >/dev/null 2>&1; then
            if [[ "${FORCE}" == "1" || ! -s "${sample_qc_dir}/${sample}.infer_experiment.txt" ]]; then
                "${RSEQC_INFER_BIN}" -r "${RSEQC_REF_BED}" -i "${bam}" > "${sample_qc_dir}/${sample}.infer_experiment.txt" 2>&1 || true
            fi
        fi
    done

    if [[ -n "${RSEQC_REF_BED}" && -f "${RSEQC_REF_BED}" ]] && command -v "${RSEQC_GENEBODY_BIN}" >/dev/null 2>&1; then
        local bam_list=()
        for sample in "${SAMPLE_NAMES[@]}"; do
            bam="${MAPPING_DIR}/${sample}.bam"
            [[ -s "${bam}" ]] && bam_list+=("${bam}")
        done
        if [[ ${#bam_list[@]} -gt 0 && ( "${FORCE}" == "1" || ! -s "${UPSTREAM_QC_DIR}/geneBody_coverage.geneBodyCoverage.txt" ) ]]; then
            "${RSEQC_GENEBODY_BIN}" -r "${RSEQC_REF_BED}" -i "$(IFS=,; echo "${bam_list[*]}")" -o "${UPSTREAM_QC_DIR}/geneBody_coverage" >/dev/null 2>&1 || true
        fi
    fi
}

write_upstream_qc_summary() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${AUDIT_TABLE_DIR}"
    if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
        log "未找到 ${PYTHON_BIN}，跳过 upstream_qc_summary.tsv。"
        return
    fi

    "${PYTHON_BIN}" - "${READS_INFO}" "${FASTQC_DIR}" "${MAPPING_DIR}" "${QUANT_DIR}" "${UPSTREAM_QC_DIR}" "${AUDIT_TABLE_DIR}/upstream_qc_summary.tsv" <<'PY'
import csv
import json
import re
import sys
from pathlib import Path

reads_info, fastqc_dir, mapping_dir, quant_dir, upstream_qc_dir, out_path = map(Path, sys.argv[1:])

with reads_info.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)
if not rows or "sample" not in rows[0]:
    raise SystemExit("reads_info.txt must have a header with a sample column")

def read_featurecounts(path):
    stats = {}
    if path.exists():
        with path.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    try:
                        stats[parts[0]] = int(float(parts[1]))
                    except ValueError:
                        pass
    total = sum(stats.values()) if stats else None
    assigned = stats.get("Assigned")
    ratio = assigned / total if total else None
    return assigned, total, ratio

def read_hisat2(path):
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    match = re.search(r"([0-9.]+)% overall alignment rate", text)
    return float(match.group(1)) if match else None

def read_fastp(path):
    if not path.exists():
        return {}
    with path.open() as handle:
        obj = json.load(handle)
    before = obj.get("summary", {}).get("before_filtering", {})
    after = obj.get("summary", {}).get("after_filtering", {})
    return {
        "fastp_before_reads": before.get("total_reads"),
        "fastp_after_reads": after.get("total_reads"),
        "fastp_after_q30_rate": after.get("q30_rate"),
        "fastp_after_gc_content": after.get("gc_content"),
    }

def read_flagstat(path):
    if not path.exists():
        return None, None
    total = mapped = None
    for line in path.read_text(errors="replace").splitlines():
        if " in total " in line:
            total = int(line.split()[0])
        elif " mapped (" in line and "mate mapped" not in line:
            mapped = int(line.split()[0])
    return total, mapped

fieldnames = [
    "sample", "group", "stage",
    "fastp_before_reads", "fastp_after_reads", "fastp_retained_ratio",
    "fastp_after_q30_rate", "fastp_after_gc_content",
    "hisat2_overall_alignment_rate",
    "featurecounts_assigned", "featurecounts_total", "featurecounts_assigned_ratio",
    "flagstat_total", "flagstat_mapped", "flagstat_mapped_ratio",
    "has_rseqc_infer_experiment",
]

with out_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        sample = row["sample"]
        group = row["group"]
        out = {
            "sample": sample,
            "group": group,
            "stage": row.get("stage", ""),
        }
        out.update(read_fastp(fastqc_dir / f"{sample}_fastp.json"))
        before = out.get("fastp_before_reads")
        after = out.get("fastp_after_reads")
        out["fastp_retained_ratio"] = after / before if before else None
        out["hisat2_overall_alignment_rate"] = read_hisat2(mapping_dir / f"{sample}.log")
        assigned, total, ratio = read_featurecounts(quant_dir / f"{sample}.log")
        out["featurecounts_assigned"] = assigned
        out["featurecounts_total"] = total
        out["featurecounts_assigned_ratio"] = ratio
        flag_total, flag_mapped = read_flagstat(upstream_qc_dir / sample / f"{sample}.flagstat.txt")
        out["flagstat_total"] = flag_total
        out["flagstat_mapped"] = flag_mapped
        out["flagstat_mapped_ratio"] = flag_mapped / flag_total if flag_total else None
        out["has_rseqc_infer_experiment"] = (upstream_qc_dir / sample / f"{sample}.infer_experiment.txt").exists()
        writer.writerow(out)
PY
}

run_multiqc_summary() {
    [[ "${DRY_RUN}" == "1" ]] && return
    mkdir -p "${MULTIQC_OUT_DIR}"
    if command -v "${MULTIQC_BIN}" >/dev/null 2>&1; then
        "${MULTIQC_BIN}" "${ROOT_DIR}" --outdir "${MULTIQC_OUT_DIR}" --filename multiqc_report.html --force >/dev/null 2>&1 || {
            log "MultiQC 运行失败，已保留其他 upstream QC 表。"
        }
    else
        printf 'MultiQC not found in PATH. Install multiqc or set MULTIQC_BIN to generate this report.\n' > "${MULTIQC_OUT_DIR}/SKIPPED.txt"
    fi
}

run_audit_stage() {
    log "Step 7/7: 上游 QC、MultiQC 与可复现性清单"
    write_run_provenance
    run_alignment_qc_tools
    write_upstream_qc_summary
    run_multiqc_summary
    log "audit 输出:"
    log "  - ${AUDIT_TABLE_DIR}/upstream_qc_summary.tsv"
    log "  - ${AUDIT_TABLE_DIR}/software_versions.tsv"
    log "  - ${AUDIT_TABLE_DIR}/reference_manifest.tsv"
    log "  - ${MULTIQC_OUT_DIR}/multiqc_report.html 或 SKIPPED.txt"
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

    "${PYTHON_BIN:-python3}" "${ROOT_DIR}/scripts/input_contract.py" validate
    configure_library_strandedness
    ensure_prerequisites
    load_samples
    [[ "${DRY_RUN}" == 1 ]] || "${PYTHON_BIN}" "${ROOT_DIR}/scripts/input_contract.py" prepare
    ensure_output_dirs
    resolve_run_order "${REQUESTED_STEPS[@]}"

    REQUESTED_STEPS_LABEL="$(format_steps "${REQUESTED_STEPS[@]}")"
    RUN_ORDER_LABEL="$(format_steps "${RUN_ORDER[@]}")"

    log "root dir        : ${ROOT_DIR}"
    log "reads info      : ${READS_INFO}"
    log "sample number   : ${#SAMPLE_NAMES[@]}"
    log "threads         : ${THREADS}"
    log "force           : ${FORCE}"
    log "dry run         : ${DRY_RUN}"
    log "requested steps : ${REQUESTED_STEPS_LABEL}"
    log "run order       : ${RUN_ORDER_LABEL}"
    log "library strand  : ${RNASEQ_LIBRARY_STRANDEDNESS_LABEL}"
    log "HISAT2 strand   : ${HISAT2_RNA_STRANDNESS:-none}"
    log "featureCounts   : strandSpecific=${FEATURECOUNTS_STRAND_SPECIFIC}"
    log "rMATS libType   : ${RNASEQ_RMATS_LIBTYPE}"
    log "APA strand type : ${RNASEQ_APA_STRANDTYPE}"
    log "reference dir   : ${REFERENCE_DIR}"
    log "hisat2 index    : ${HISAT2_INDEX}"
    log "gtf file        : ${GTF_FILE}"
    log "fastp bin       : ${FASTP_BIN}"
    log "hisat2 bin      : ${HISAT2_BIN}"
    log "samtools bin    : ${SAMTOOLS_BIN}"
    log "Rscript bin     : ${RSCRIPT_BIN}"

    local step
    for step in "${RUN_ORDER[@]}"; do
        run_step "${step}"
    done

    log "完成: ${RUN_ORDER_LABEL}"
}

main "$@"
