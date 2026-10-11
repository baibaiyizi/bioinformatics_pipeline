#!/usr/bin/env bash
set -euo pipefail

# 小鼠 GRCm39 paired-end ATAC 上游模板。
# 先编辑 config/project.env 和 config/samples.tsv，运行 --check 验证设计。
# 上游算法来自 liver 项目，分组全部读取显式样本表。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "$SCRIPT_DIR/scripts/load_config.sh"

# ==============================================================================
# 以下为流程内部状态，一般不修改
# ==============================================================================
READS_INFO_EXPLICIT=0
RENDER_RMD=0
TOTAL_STEPS=15
RUN_ALL_STEPS=0
DRY_RUN=0
FORCE=0
RANGE_EXPLICIT=0
UPSTREAM_WAS_STALE=0
UPSTREAM_WARNING_COUNT=0
UPSTREAM_RUN_START_EPOCH=0
declare -a STEP_ARGS=()
declare -a SELECTED_STEPS=()
declare -A STEP_SELECTED=()

# Prefer a tool from a known conda env, otherwise fall back to PATH.
resolve_cmd() {
    local preferred="$1"
    local fallback="$2"
    if [[ -n "$preferred" && -x "$preferred" ]]; then
        printf '%s\n' "$preferred"
    elif command -v "$fallback" >/dev/null 2>&1; then
        command -v "$fallback"
    else
        printf '%s\n' ""
    fi
}

cmd_available() {
    local cmd="$1"
    if [[ -z "$cmd" ]]; then
        return 1
    fi
    if [[ "$cmd" == */* ]]; then
        [[ -x "$cmd" ]]
    else
        command -v "$cmd" >/dev/null 2>&1
    fi
}

print_user_config_summary() {
    log "Config | fastp: min_length=$FASTP_MIN_LENGTH, window=$FASTP_CUT_WINDOW_SIZE, mean_quality=$FASTP_CUT_MEAN_QUALITY"
    log "Config | mapping: max_fragment=$BOWTIE2_MAX_FRAGMENT, MAPQ=$ATAC_MAPQ, exclude_flags=$ATAC_SAMTOOLS_EXCLUDE_FLAGS"
    log "Config | peaks/RPGC/TSS: MACS2_q=$ATAC_MACS2_QVALUE, genome_size=$ATAC_EFFECTIVE_GENOME_SIZE, bin=$ATAC_BIGWIG_BIN_SIZE, TSS_flank=$ATAC_TSS_FLANK, heatmap_colors=$ATAC_DEEPTOOLS_HEATMAP_COLORS"
    log "Config | IDR: relaxed_p=$ATAC_IDR_MACS2_PVALUE, global_IDR=$ATAC_IDR_THRESHOLD, seed=$ATAC_IDR_SEED, keep_intermediates=$ATAC_KEEP_IDR_INTERMEDIATES"
    log "Config | optional: raw_MACS2_bigWig=$RAW_MACS2_BIGWIG, ataqv=$ATAQV_CMD, TOBIAS=$TOBIAS_CMD"
}

MACS2_CMD="${MACS2_CMD:-$(resolve_cmd "$MACS2_ENV/macs2" "macs2")}"
BAMCOVERAGE_CMD="${BAMCOVERAGE_CMD:-$(resolve_cmd "$MACS2_ENV/bamCoverage" "bamCoverage")}"
COMPUTEMATRIX_CMD="${COMPUTEMATRIX_CMD:-$(resolve_cmd "$MACS2_ENV/computeMatrix" "computeMatrix")}"
PLOTPROFILE_CMD="${PLOTPROFILE_CMD:-$(resolve_cmd "$MACS2_ENV/plotProfile" "plotProfile")}"
PLOTHEATMAP_CMD="${PLOTHEATMAP_CMD:-$(resolve_cmd "$MACS2_ENV/plotHeatmap" "plotHeatmap")}"
IDR_CMD="${IDR_CMD:-$(resolve_cmd "$IDR_ENV/idr" "idr")}"
TOBIAS_CMD="${TOBIAS_CMD:-$(resolve_cmd "$TOBIAS_ENV/TOBIAS" "TOBIAS")}"
BEDGRAPHTOBIGWIG_CMD="${BEDGRAPHTOBIGWIG_CMD:-$(resolve_cmd "bedGraphToBigWig" "bedGraphToBigWig")}"
FEATURECOUNTS_CMD="${FEATURECOUNTS_CMD:-$(resolve_cmd "${FEATURECOUNTS_ENV:+$FEATURECOUNTS_ENV/featureCounts}" "featureCounts")}"
ATAQV_CMD="${ATAQV_CMD:-$(resolve_cmd "$ATAQV_ENV/ataqv" "ataqv")}"
PICARD_CMD="${PICARD_CMD:-picard}"

declare -a SAMPLES=()
declare -a REQUESTED_SAMPLES=()
declare -a PARALLEL_PIDS=()
declare -a GROUP_LEVELS=()
declare -A SAMPLE_GROUP_BY_NAME=()
declare -A SAMPLE_STAGE_BY_NAME=()
declare -A SAMPLE_GROUP_STAGE_BY_NAME=()
declare -A SAMPLE_RAW_R1_BY_NAME=()
declare -A SAMPLE_RAW_R2_BY_NAME=()

# ------------------------------------------------------------------
# 1. CLI and logging helpers
# ------------------------------------------------------------------

usage() {
    cat <<EOF
用法: bash run_upstream.sh [--steps 1,4,08-09 | --all | --list | --check | --dry-run]
配置: config/project.env；样本表由 ATAC_SAMPLES 指定。
默认运行 1-11 和 14；footprint/BINDetect 须先在配置中启用。
--threads N / --parallel-samples N：工具线程与样本并行。
--samples a,b：只执行指定上游样本，不改变下游实验设计。
--force：允许执行完整重算后清除源流程的 STALE 状态。
检查与 dry-run 不创建日志或分析结果；不会证明输入 FASTQ 完整性或工具可用。
EOF
}

step_name() {
    case "$1" in
        1) echo "raw_fastqc" ;;
        2) echo "fastp" ;;
        3) echo "clean_fastqc" ;;
        4) echo "mapping" ;;
        5) echo "peaks" ;;
        6) echo "blacklist" ;;
        7) echo "bigwig" ;;
        8) echo "frip" ;;
        9) echo "tss" ;;
        10) echo "idr" ;;
        11) echo "multiqc" ;;
        12) echo "footprint" ;;
        13) echo "bindetect" ;;
        14) echo "consensus_counts" ;;
        15) echo "ataqv_qc" ;;
    esac
}

step_title_text() {
    case "$1" in
        1) echo "Raw FastQC" ;;
        2) echo "fastp trimming" ;;
        3) echo "Clean FastQC" ;;
        4) echo "Mapping + filtering + dedup + library complexity QC" ;;
        5) echo "Peak calling" ;;
        6) echo "Blacklist filtering" ;;
        7) echo "BigWig generation" ;;
        8) echo "FRiP quantification" ;;
        9) echo "TSS enrichment" ;;
        10) echo "IDR reproducibility" ;;
        11) echo "MultiQC + software/QC aggregation" ;;
        12) echo "Optional TOBIAS footprint tracks" ;;
        13) echo "Optional TOBIAS BINDetect motif-wide analysis" ;;
        14) echo "Consensus peak sets and count matrices" ;;
        15) echo "Optional ataqv standardized ATAC QC" ;;
    esac
}

print_step_table() {
    local i
    for i in $(seq 1 "$TOTAL_STEPS"); do
        printf '  %02d  %-14s  %s\n' "$i" "$(step_name "$i")" "$(step_title_text "$i")"
    done
}

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

warn() {
    UPSTREAM_WARNING_COUNT=$((UPSTREAM_WARNING_COUNT + 1))
    printf '[%s] WARNING: %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
    printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2
    exit 1
}

step_log() {
    local step="$1"
    local title="$2"
    log "Step ${step}/${TOTAL_STEPS}: ${title}"
}

step_enabled() {
    local step="$1"
    [[ -n "${STEP_SELECTED[$step]:-}" ]]
}

step_id_from_token() {
    local token="$1"
    token="$(printf '%s' "$token" | tr '[:upper:]' '[:lower:]')"
    token="${token//_/-}"

    if [[ "$token" =~ ^0*[0-9]+$ ]]; then
        token="$((10#$token))"
        (( token >= 1 && token <= TOTAL_STEPS )) || return 1
        printf '%s\n' "$token"
        return 0
    fi

    case "$token" in
        raw|raw-fastqc|rawqc|fastqc-raw) echo 1 ;;
        fastp|trim|trimming|clean) echo 2 ;;
        clean-fastqc|cleanqc|fastqc-clean) echo 3 ;;
        map|mapping|align|alignment|bam|dedup|filter) echo 4 ;;
        peak|peaks|macs2|callpeak|peak-calling) echo 5 ;;
        blacklist|clean-peaks|filter-peaks) echo 6 ;;
        bigwig|bw|tracks) echo 7 ;;
        frip) echo 8 ;;
        tss|tsse|tss-enrichment) echo 9 ;;
        idr|reproducibility) echo 10 ;;
        multiqc|aggregate) echo 11 ;;
        footprint|footprinting|tobias|tracks-footprint) echo 12 ;;
        bindetect|bin-detect|motif-wide) echo 13 ;;
        consensus|consensus-counts|consensus-count|count-matrix|counts|featurecounts|peak-counts) echo 14 ;;
        ataqv|ataqv-qc|standard-qc|standardized-qc|tsse-standard) echo 15 ;;
        *) return 1 ;;
    esac
}

parse_step_tokens() {
    local tokens=("$@")
    local token part start end i

    for token in "${tokens[@]}"; do
        token="${token//[[:space:]]/}"
        [[ -z "$token" ]] && continue
        IFS=',' read -r -a raw_parts <<<"$token"
        for part in "${raw_parts[@]}"; do
            part="${part//[[:space:]]/}"
            [[ -z "$part" ]] && continue
            if [[ "$part" =~ ^([0-9]+)(-|\\.\\.)([0-9]+)$ ]]; then
                start="$(step_id_from_token "${BASH_REMATCH[1]}")" || die "Invalid step range start: ${BASH_REMATCH[1]}"
                end="$(step_id_from_token "${BASH_REMATCH[3]}")" || die "Invalid step range end: ${BASH_REMATCH[3]}"
                if (( start <= end )); then
                    for i in $(seq "$start" "$end"); do printf '%s\n' "$i"; done
                else
                    for i in $(seq "$start" -1 "$end"); do printf '%s\n' "$i"; done
                fi
            else
                step_id_from_token "$part" || die "Invalid step token: $part"
            fi
        done
    done
}

resolve_selected_steps() {
    local raw_steps=()
    local i step

    if (( RUN_ALL_STEPS == 1 )); then
        raw_steps=($(seq 1 "$TOTAL_STEPS"))
    elif (( ${#STEP_ARGS[@]} > 0 )); then
        raw_steps=($(parse_step_tokens "${STEP_ARGS[@]}"))
    else
        # Standard run: keep BINDetect and ataqv optional, but make the
        # auditable consensus count matrix part of the default upstream output.
        raw_steps=($(seq 1 11) 14)
    fi

    SELECTED_STEPS=()
    STEP_SELECTED=()
    for i in $(seq 1 "$TOTAL_STEPS"); do
        for step in "${raw_steps[@]}"; do
            if [[ "$step" == "$i" && -z "${STEP_SELECTED[$i]:-}" ]]; then
                if [[ "$i" == 12 || "$i" == 13 ]]; then
                    [[ "$ATAC_ENABLE_FOOTPRINT" == 1 ]] || { echo "[SKIPPED] footprint 未启用: $i" >&2; continue; }
                fi
                SELECTED_STEPS+=("$i")
                STEP_SELECTED["$i"]=1
            fi
        done
    done
    (( ${#SELECTED_STEPS[@]} > 0 )) || die "No valid steps selected."
}

print_selected_steps() {
    local step
    for step in "${SELECTED_STEPS[@]}"; do
        printf '  %02d  %-14s  %s\n' "$step" "$(step_name "$step")" "$(step_title_text "$step")"
    done
}

require_cmd() {
    local cmd="$1"
    local label="${2:-$1}"
    cmd_available "$cmd" || die "Missing required command: $label"
}

prepare_logging() {
    mkdir -p "$WORK_DIR/logs"
    local log_file="$WORK_DIR/logs/atac_pipeline_$(date '+%Y%m%d_%H%M%S').log"
    exec > >(tee -i "$log_file") 2>&1
    log "Pipeline log: $log_file"
}

# Create all directories up front so each step can focus on its own work.
ensure_dirs() {
    mkdir -p \
        "$WORK_DIR/metadata" \
        "$WORK_DIR/01rawdata/qc_results" \
        "$WORK_DIR/01rawdata/cleaned_data" \
        "$WORK_DIR/02mapping/bam" \
        "$WORK_DIR/02mapping/qc" \
        "$WORK_DIR/03peaks" \
        "$WORK_DIR/04counts" \
        "$WORK_DIR/05bigwig" \
        "$WORK_DIR/06qc_report" \
        "$WORK_DIR/07footprinting"
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --steps)
                shift
                [[ $# -gt 0 ]] || die "--steps requires a step list"
                while [[ $# -gt 0 && "$1" != --* ]]; do
                    STEP_ARGS+=("$1")
                    shift
                done
                ;;
            --all)
                RUN_ALL_STEPS=1
                shift
                ;;
            --list)
                print_step_table
                exit 0
                ;;
            --dry-run|--check)
                DRY_RUN=1
                shift
                ;;
            --force)
                FORCE=1
                shift
                ;;
            --threads)
                THREADS="$2"
                shift 2
                ;;
            --parallel-samples)
                PARALLEL_SAMPLES="$2"
                shift 2
                ;;
            --samples)
                IFS=',' read -r -a REQUESTED_SAMPLES <<<"$2"
                shift 2
                ;;
            --help|-h)
                usage
                exit 0
                ;;
            *)
                die "Unknown argument: $1"
                ;;
        esac
    done
}

# Basic sanity checks for user-facing numeric arguments.
validate_numeric_args() {
    [[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
    [[ "$PARALLEL_SAMPLES" =~ ^[0-9]+$ ]] || die "--parallel-samples must be an integer"
    [[ "$ATAC_IDR_SEED" =~ ^[0-9]+$ ]] || die "ATAC_IDR_SEED must be a non-negative integer"
    local integer_name integer_value
    for integer_name in FASTP_MIN_LENGTH FASTP_CUT_WINDOW_SIZE FASTP_CUT_MEAN_QUALITY \
        BOWTIE2_MAX_FRAGMENT ATAC_MAPQ ATAC_SAMTOOLS_EXCLUDE_FLAGS \
        ATAC_EFFECTIVE_GENOME_SIZE ATAC_BIGWIG_BIN_SIZE ATAC_TSS_FLANK; do
        integer_value="${!integer_name}"
        [[ "$integer_value" =~ ^[0-9]+$ ]] || die "$integer_name must be a non-negative integer"
    done
    (( FASTP_MIN_LENGTH >= 1 && FASTP_CUT_WINDOW_SIZE >= 1 && BOWTIE2_MAX_FRAGMENT >= 1 )) || die "FASTP_MIN_LENGTH, FASTP_CUT_WINDOW_SIZE and BOWTIE2_MAX_FRAGMENT must be >= 1"
    (( ATAC_EFFECTIVE_GENOME_SIZE >= 1 && ATAC_BIGWIG_BIN_SIZE >= 1 && ATAC_TSS_FLANK >= 1 )) || die "RPGC genome size, bin size and TSS flank must be >= 1"
    (( FASTP_CUT_MEAN_QUALITY <= 40 )) || die "FASTP_CUT_MEAN_QUALITY must be between 0 and 40"
    [[ "$RAW_MACS2_BIGWIG" =~ ^[01]$ ]] || die "RAW_MACS2_BIGWIG must be 0 or 1"
    [[ "$ATAC_KEEP_IDR_INTERMEDIATES" =~ ^[01]$ ]] || die "ATAC_KEEP_IDR_INTERMEDIATES must be 0 or 1"
    awk -v v="$ATAC_MACS2_QVALUE" 'BEGIN{exit !(v>0 && v<=1)}' || die "ATAC_MACS2_QVALUE must be in (0,1]"
    awk -v v="$ATAC_IDR_MACS2_PVALUE" 'BEGIN{exit !(v>0 && v<=1)}' || die "ATAC_IDR_MACS2_PVALUE must be in (0,1]"
    awk -v v="$ATAC_IDR_THRESHOLD" 'BEGIN{exit !(v>0 && v<=1)}' || die "ATAC_IDR_THRESHOLD must be in (0,1]"
    THREADS=$((10#$THREADS))
    PARALLEL_SAMPLES=$((10#$PARALLEL_SAMPLES))
    (( THREADS >= 1 )) || die "--threads must be >= 1"
    (( PARALLEL_SAMPLES >= 1 )) || die "--parallel-samples must be >= 1"
    resolve_selected_steps
}

trim_cr() {
    local value="$1"
    value="${value%$'\r'}"
    printf '%s' "$value"
}

contains_value() {
    local needle="$1"
    shift || true
    local item
    for item in "$@"; do
        [[ "$item" == "$needle" ]] && return 0
    done
    return 1
}

resolve_reads_path() {
    local path_value="$1"
    if [[ "$path_value" = /* ]]; then
        printf '%s\n' "$path_value"
    else
        printf '%s/%s\n' "$WORK_DIR" "$path_value"
    fi
}

sample_requested() {
    local sample="$1"
    if [[ "${#REQUESTED_SAMPLES[@]}" -eq 0 ]]; then
        return 0
    fi
    contains_value "$sample" "${REQUESTED_SAMPLES[@]}"
}

# Load sample metadata from reads_info.txt. This is the source of truth for sample order,
# group labels and raw FASTQ inputs.
discover_samples() {
    [[ -f "$READS_INFO" ]] || die "Sample manifest not found: $READS_INFO"

    local header_line
    IFS= read -r header_line <"$READS_INFO" || true
    [[ -n "$header_line" ]] || die "reads_info.txt is empty: $READS_INFO"

    local -a header_cols=()
    local header_encoded="${header_line//$'\t'/$'\x1f'}"
    IFS=$'\x1f' read -r -a header_cols <<<"$header_encoded"
    declare -A col_idx=()
    local i col
    for i in "${!header_cols[@]}"; do
        col="$(trim_cr "${header_cols[$i]}")"
        [[ -n "$col" ]] && col_idx["$col"]="$i"
    done

    [[ -n "${col_idx[sample]:-}" ]] || die "reads_info.txt missing required column: sample"
    [[ -n "${col_idx[group]:-}" ]] || die "reads_info.txt missing required column: group"
    local raw_fq1_col="raw_fq1"
    local raw_fq2_col="raw_fq2"
    [[ -n "${col_idx[$raw_fq1_col]:-}" ]] || die "reads_info.txt missing required column: raw_fq1/fq1"
    [[ -n "${col_idx[$raw_fq2_col]:-}" ]] || die "reads_info.txt missing required column: raw_fq2/fq2"

    SAMPLES=()
    SAMPLE_GROUP_BY_NAME=()
    SAMPLE_STAGE_BY_NAME=()
    SAMPLE_GROUP_STAGE_BY_NAME=()
    SAMPLE_RAW_R1_BY_NAME=()
    SAMPLE_RAW_R2_BY_NAME=()
    GROUP_LEVELS=()

    local -a seen_samples=()
    local -a group_levels=()
    local line line_no=1
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_no=$((line_no + 1))
        [[ -n "${line//[$'\t\r ']/}" ]] || continue

        local -a fields=()
        local line_encoded="${line//$'\t'/$'\x1f'}"
        IFS=$'\x1f' read -r -a fields <<<"$line_encoded"
        local sample group stage group_stage raw_fq1 raw_fq2 resolved_r1 resolved_r2
        sample="$(trim_cr "${fields[${col_idx[sample]}]:-}")"
        group="$(trim_cr "${fields[${col_idx[group]}]:-}")"
        raw_fq1="$(trim_cr "${fields[${col_idx[$raw_fq1_col]}]:-}")"
        raw_fq2="$(trim_cr "${fields[${col_idx[$raw_fq2_col]}]:-}")"
        stage=""
        group_stage=""
        if [[ -n "${col_idx[stage]:-}" ]]; then
            stage="$(trim_cr "${fields[${col_idx[stage]}]:-}")"
        fi
        if [[ -n "${col_idx[group_stage]:-}" ]]; then
            group_stage="$(trim_cr "${fields[${col_idx[group_stage]}]:-}")"
        fi
        [[ -n "$sample" ]] || die "reads_info.txt line $line_no has empty sample"
        [[ -n "$group" ]] || die "reads_info.txt line $line_no has empty group for sample $sample"
        [[ -n "$raw_fq1" && -n "$raw_fq2" ]] || die "reads_info.txt line $line_no missing raw_fq1/raw_fq2 for sample $sample"
        contains_value "$sample" "${seen_samples[@]}" && die "reads_info.txt has duplicated sample: $sample"
        seen_samples+=("$sample")

        [[ -n "$stage" ]] || die "stage 必须显式提供"
        [[ -n "$group_stage" ]] || group_stage="${group}${stage}"
        contains_value "$group" "${group_levels[@]}" || group_levels+=("$group")

        resolved_r1="$(resolve_reads_path "$raw_fq1")"
        resolved_r2="$(resolve_reads_path "$raw_fq2")"
        if (( DRY_RUN == 0 )) && { step_enabled 1 || step_enabled 2; }; then
            [[ -f "$resolved_r1" ]] || die "FASTQ not found: $resolved_r1"
            [[ -f "$resolved_r2" ]] || die "FASTQ not found: $resolved_r2"
        fi

        SAMPLE_GROUP_BY_NAME["$sample"]="$group"
        SAMPLE_STAGE_BY_NAME["$sample"]="$stage"
        SAMPLE_GROUP_STAGE_BY_NAME["$sample"]="$group_stage"
        SAMPLE_RAW_R1_BY_NAME["$sample"]="$resolved_r1"
        SAMPLE_RAW_R2_BY_NAME["$sample"]="$resolved_r2"

        if sample_requested "$sample"; then
            SAMPLES+=("$sample")
        fi
    done < <(tail -n +2 "$READS_INFO")

    [[ "${#group_levels[@]}" -eq 2 ]] || die "reads_info.txt must contain exactly two groups; got: ${group_levels[*]}"
    GROUP_LEVELS=("${group_levels[@]}")
    [[ "${#SAMPLES[@]}" -gt 0 ]] || die "No selected samples found in reads_info.txt"
    local requested
    for requested in "${REQUESTED_SAMPLES[@]}"; do
        [[ -n "${SAMPLE_GROUP_BY_NAME[$requested]:-}" ]] || die "--samples requested sample not found in reads_info.txt: $requested"
    done
}

# Infer biological group from reads_info, falling back to the sample name prefix.
sample_group() {
    local sample="$1"
    if [[ -n "${SAMPLE_GROUP_BY_NAME[$sample]:-}" ]]; then
        printf '%s\n' "${SAMPLE_GROUP_BY_NAME[$sample]}"
    else
        die "样本未在表中声明: $sample"
    fi
}

# A compact manifest makes the downstream Rmd independent from hard-coded filenames.
write_sample_manifest() {
    local manifest="$WORK_DIR/metadata/sample_manifest.tsv"
    {
        printf 'sample\tgroup\traw_r1\traw_r2\tclean_r1\tclean_r2\tmarked_bam\tbam\tpeak\tbigwig\tbigwig_normalization\n'
        local sample
        for sample in "${SAMPLES[@]}"; do
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$sample" \
                "$(sample_group "$sample")" \
                "${SAMPLE_RAW_R1_BY_NAME[$sample]}" \
                "${SAMPLE_RAW_R2_BY_NAME[$sample]}" \
                "01rawdata/cleaned_data/${sample}_1.fastq.gz" \
                "01rawdata/cleaned_data/${sample}_2.fastq.gz" \
                "02mapping/bam/${sample}.marked.bam" \
                "02mapping/bam/${sample}.final.bam" \
                "03peaks/${sample}_peaks.clean.narrowPeak" \
                "05bigwig/${sample}_RPGC.bw" \
                "RPGC"
        done
    } >"$manifest"
    log "Sample manifest written: $manifest"
}

write_upstream_freshness_signature() {
    local phase="${1:-start}"
    local signature_file="$WORK_DIR/metadata/upstream_config_signature.tsv"
    local payload_file previous_signature previous_status current_signature status sample path size mtime
    payload_file="$(mktemp "$WORK_DIR/metadata/.atac_signature.XXXXXX")"
    {
        printf 'threads=%s\nparallel_samples=%s\nreference=%s\nblacklist=%s\ntss=%s\n' \
            "$THREADS" "$PARALLEL_SAMPLES" "$REF_GENOME_INDEX" "$BLACKLIST_BED" "$TSS_BED"
        printf 'idr_macs2_p=%s\nidr_threshold=%s\nidr_seed=%s\nataqv_organism=%s\n' \
            "$ATAC_IDR_MACS2_PVALUE" "$ATAC_IDR_THRESHOLD" "$ATAC_IDR_SEED" "$ATAQV_ORGANISM"
        printf 'fastp=%s|%s|%s\nbowtie_max_fragment=%s\nmapq=%s\nexclude_flags=%s\n' \
            "$FASTP_MIN_LENGTH" "$FASTP_CUT_WINDOW_SIZE" "$FASTP_CUT_MEAN_QUALITY" \
            "$BOWTIE2_MAX_FRAGMENT" "$ATAC_MAPQ" "$ATAC_SAMTOOLS_EXCLUDE_FLAGS"
        printf 'macs2_q=%s\nrpgc=%s|%s\ntss_flank=%s\nraw_macs2_bigwig=%s\n' \
            "$ATAC_MACS2_QVALUE" "$ATAC_EFFECTIVE_GENOME_SIZE" "$ATAC_BIGWIG_BIN_SIZE" \
            "$ATAC_TSS_FLANK" "$RAW_MACS2_BIGWIG"
        for path in "$WORK_DIR/run_upstream.sh" "$READS_INFO" "$BLACKLIST_BED" "$TSS_BED"; do
            if [[ -e "$path" ]]; then
                size="$(stat -c '%s' "$path")"; mtime="$(stat -c '%Y' "$path")"
                printf 'file=%s|%s|%s\n' "$path" "$size" "$mtime"
            else
                printf 'file=%s|MISSING\n' "$path"
            fi
        done
        for sample in "${SAMPLES[@]}"; do
            for path in "${SAMPLE_RAW_R1_BY_NAME[$sample]}" "${SAMPLE_RAW_R2_BY_NAME[$sample]}"; do
                if [[ -e "$path" ]]; then
                    size="$(stat -c '%s' "$path")"; mtime="$(stat -c '%Y' "$path")"
                    printf 'file=%s|%s|%s\n' "$path" "$size" "$mtime"
                else
                    printf 'file=%s|MISSING\n' "$path"
                fi
            done
            printf 'sample=%s|group=%s\n' "$sample" "$(sample_group "$sample")"
        done
    } >"$payload_file"
    if command -v sha256sum >/dev/null 2>&1; then
        current_signature="$(sha256sum "$payload_file" | awk '{print $1}')"
    else
        current_signature="$(cksum "$payload_file" | awk '{print $1 "-" $2}')"
    fi
    rm -f "$payload_file"
    previous_signature=""
    previous_status="MISSING_SIGNATURE"
    if [[ -f "$signature_file" ]]; then
        previous_signature="$(awk -F '\t' '$1=="Signature" {print $2; exit}' "$signature_file")"
        previous_status="$(awk -F '\t' '$1=="Status" {print $2; exit}' "$signature_file")"
    fi
    status="CURRENT"
    if [[ -z "$previous_signature" ]] && compgen -G "$WORK_DIR/02mapping/bam/*.final.bam" >/dev/null; then
        status="STALE_MISSING_SIGNATURE"
        UPSTREAM_WAS_STALE=1
    elif [[ -n "$previous_signature" && "$previous_signature" != "$current_signature" ]]; then
        status="STALE_PREVIOUS_RESULTS"
        UPSTREAM_WAS_STALE=1
    elif [[ "$previous_status" != "CURRENT" || "$UPSTREAM_WAS_STALE" == "1" ]]; then
        UPSTREAM_WAS_STALE=1
        status="STALE_PREVIOUS_RESULTS"
    fi
    # A repeated signature does not make previously stale files current. Only a
    # complete forced regeneration without warnings and with fresh core outputs
    # can clear the stale state. Optional TOBIAS/ataqv have separate status files.
    if [[ "$phase" == "start" ]]; then
        UPSTREAM_RUN_START_EPOCH="$(date +%s)"
        status="STALE_RUNNING"
    else
        local required_step all_required_selected=1 outputs_fresh=1
        for required_step in 1 2 3 4 5 6 7 8 9 10 11 14; do
            if ! step_enabled "$required_step"; then all_required_selected=0; break; fi
        done
        for sample in "${SAMPLES[@]}"; do
            for path in \
                "${SAMPLE_RAW_R1_BY_NAME[$sample]}" "${SAMPLE_RAW_R2_BY_NAME[$sample]}" \
                "$WORK_DIR/01rawdata/cleaned_data/${sample}_1.fastq.gz" \
                "$WORK_DIR/01rawdata/cleaned_data/${sample}_2.fastq.gz"; do
                [[ -s "$path" ]] || outputs_fresh=0
            done
            for path in "$WORK_DIR/02mapping/bam/${sample}.final.bam" \
                        "$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak" \
                        "$WORK_DIR/05bigwig/${sample}_RPGC.bw"; do
                if [[ ! -s "$path" ]] || [[ "$(stat -c '%Y' "$path" 2>/dev/null || echo 0)" -lt "$UPSTREAM_RUN_START_EPOCH" ]]; then outputs_fresh=0; fi
            done
        done
        for path in "$WORK_DIR/04counts/consensus_peak_counts.raw.tsv" \
                    "$WORK_DIR/04counts/FRiP_scores.txt" \
                    "$WORK_DIR/03peaks/IDR_Reproducibility_Summary.tsv" \
                    "$WORK_DIR/06qc_report/atac_multiqc_report.html"; do
            if [[ ! -s "$path" ]] || [[ "$(stat -c '%Y' "$path" 2>/dev/null || echo 0)" -lt "$UPSTREAM_RUN_START_EPOCH" ]]; then outputs_fresh=0; fi
        done
        if [[ "$all_required_selected" == "1" && "$FORCE" == "1" && "$outputs_fresh" == "1" && "$UPSTREAM_WARNING_COUNT" == "0" ]]; then
            status="CURRENT"
        else
            status="STALE_PARTIAL_RERUN"
        fi
    fi

    {
        printf 'Field\tValue\n'
        printf 'Signature\t%s\n' "$current_signature"
        printf 'PreviousSignature\t%s\n' "${previous_signature:-NA}"
        printf 'Status\t%s\n' "$status"
        printf 'GeneratedAt\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    } >"$signature_file"
    [[ "$status" == "CURRENT" ]] || log "Upstream status: $status; current source execution alone does not verify historical outputs."
}

# ------------------------------------------------------------------
# 2. Pre-flight checks
# ------------------------------------------------------------------

check_reference_files() {
    [[ -f "${REF_GENOME_INDEX}.1.bt2" || -f "${REF_GENOME_INDEX}.1.bt2l" ]] || die "Bowtie2 index not found: $REF_GENOME_INDEX"
    [[ -f "$BLACKLIST_BED" ]] || warn "Blacklist bed not found: $BLACKLIST_BED"
    [[ -f "$TSS_BED" ]] || warn "TSS bed not found: $TSS_BED"
    [[ -f "$CHROM_SIZE" ]] || warn "Chrom size file not found: $CHROM_SIZE"
    [[ -f "$GENOME_FASTA" ]] || warn "Genome FASTA not found: $GENOME_FASTA"
    [[ -f "$GENE_GTF" ]] || warn "Gene GTF not found: $GENE_GTF"
}

check_required_tools() {
    local required=()
    step_enabled 1 && required+=(fastqc)
    step_enabled 2 && required+=(fastp)
    step_enabled 3 && required+=(fastqc)
    if step_enabled 4; then
        required+=(bowtie2 samtools awk "$PICARD_CMD")
    fi
    if step_enabled 5; then
        required+=("$MACS2_CMD")
    fi
    if step_enabled 6; then
        command -v bedtools >/dev/null 2>&1 || require_cmd python3 "python3 (bedtools fallback for blacklist filtering)"
    fi
    if step_enabled 8; then
        required+=(samtools awk)
    fi
    if step_enabled 9; then
        if [[ ! -x "$COMPUTEMATRIX_CMD" || ! -x "$PLOTPROFILE_CMD" ]]; then
            warn "deepTools not available; Step 9 will be skipped."
        fi
    fi
    if step_enabled 10 && [[ ! -x "$IDR_CMD" ]]; then
        warn "IDR not available; Step 10 will be skipped."
    fi
    if step_enabled 11 && ! command -v multiqc >/dev/null 2>&1; then
        warn "multiqc not available; Step 11 will be skipped."
    fi
    if step_enabled 7; then
        if [[ ! -x "$BAMCOVERAGE_CMD" ]]; then
            warn "bamCoverage not available; normalized bigWig files will be skipped."
        fi
        if [[ "$RAW_MACS2_BIGWIG" == "1" ]] && ! cmd_available "$BEDGRAPHTOBIGWIG_CMD"; then
            warn "bedGraphToBigWig not available; raw MACS2 bigWig files will be skipped."
        fi
    fi
    if step_enabled 12; then
        if ! cmd_available "$TOBIAS_CMD"; then
            warn "TOBIAS not available; Step 12 will only prepare consensus peaks and merged BAMs."
        fi
        if [[ ! -f "$GENOME_FASTA" ]]; then
            warn "Genome FASTA missing; TOBIAS Step 12 will be skipped."
        fi
    fi
    if step_enabled 13; then
        if ! cmd_available "$TOBIAS_CMD"; then
            warn "TOBIAS not available; Step 13 BINDetect will be skipped."
        fi
        if [[ ! -f "$GENOME_FASTA" ]]; then
            warn "Genome FASTA missing; TOBIAS Step 13 will be skipped."
        fi
        if [[ -z "$TOBIAS_MOTIFS" || ! -f "$TOBIAS_MOTIFS" ]]; then
            warn "TOBIAS motif database not configured; Step 13 BINDetect will be skipped."
        fi
    fi
    if step_enabled 14; then
        if ! cmd_available "$FEATURECOUNTS_CMD"; then
            warn "featureCounts not available; Step 14 will write SKIPPED rather than substitute alignment counts."
        fi
    fi
    if step_enabled 15; then
        if ! cmd_available "$ATAQV_CMD"; then
            warn "ataqv not available; Step 15 will write skipped status."
        fi
        if [[ -z "$ATAQV_ORGANISM" ]]; then
            warn "ATAQV_ORGANISM is not set; Step 15 will write skipped status unless configured."
        fi
    fi

    local cmd
    for cmd in "${required[@]}"; do
        require_cmd "$cmd"
    done
}

# ------------------------------------------------------------------
# 3. Parallel execution helpers
# ------------------------------------------------------------------

# Track only worker PIDs started by this script.
# This avoids waiting on the log tee process created by prepare_logging().
run_with_parallel_limit() {
    local limit="$1"
    shift
    trim_parallel_pids
    "$@" &
    PARALLEL_PIDS+=("$!")
    while (( "${#PARALLEL_PIDS[@]}" >= limit )); do
        wait -n "${PARALLEL_PIDS[@]}"
        trim_parallel_pids
    done
}

trim_parallel_pids() {
    local -a alive_pids=()
    local pid
    for pid in "${PARALLEL_PIDS[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            alive_pids+=("$pid")
        fi
    done
    PARALLEL_PIDS=("${alive_pids[@]}")
}

# Wait for the sample jobs of the current step and propagate the first error code.
wait_for_parallel_jobs() {
    local pid
    local exit_code=0
    local status=0
    for pid in "${PARALLEL_PIDS[@]}"; do
        if wait "$pid"; then
            continue
        else
            status=$?
            if (( exit_code == 0 )); then
                exit_code="$status"
            fi
        fi
    done
    PARALLEL_PIDS=()
    return "$exit_code"
}

# Run one worker function over all selected samples using the shared parallel limit.
run_parallel_for_samples() {
    local worker="$1"
    local sample
    for sample in "${SAMPLES[@]}"; do
        run_with_parallel_limit "$PARALLEL_SAMPLES" "$worker" "$sample"
    done
    wait_for_parallel_jobs
}

# Write ENCODE-style fragment complexity metrics from the pre-dedup filtered BAM.
write_library_complexity_metrics() {
    local sample="$1"
    local bam_file="$2"
    local out_file="$WORK_DIR/02mapping/qc/${sample}.library_complexity.tsv"

    printf 'Sample\tTotalFragments\tDistinctFragments\tOneReadFragments\tTwoReadFragments\tNRF\tPBC1\tPBC2\tMetricSource\n' >"$out_file"
    if [[ ! -f "$bam_file" ]]; then
        warn "Missing filtered BAM for library complexity: $sample"
        return 0
    fi

    if ! samtools view -f 64 -F "$ATAC_SAMTOOLS_EXCLUDE_FLAGS" "$bam_file" | awk -v sample="$sample" '
        BEGIN {FS="\t"; OFS="\t"}
        {
            tlen = $9 + 0
            if(tlen == 0 || $3 == "*") next
            if(tlen > 0) {
                start = $4
                end = $4 + tlen - 1
            } else {
                start = $4 + tlen
                end = $4 - 1
            }
            if(start < 1) start = 1
            key = $3 ":" start ":" end
            total++
            count[key]++
        }
        END {
            for(key in count) {
                distinct++
                if(count[key] == 1) one++
                else if(count[key] == 2) two++
            }
            nrf = (total > 0 ? distinct / total : "NA")
            pbc1 = (distinct > 0 ? one / distinct : "NA")
            pbc2 = (two > 0 ? one / two : "NA")
            printf "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\tpre_dedup_filtered_bam_first_mates\n", sample, total, distinct, one, two, nrf, pbc1, pbc2
        }' >>"$out_file"; then
        warn "Failed to compute library complexity metrics for $sample"
    fi
}

# Summarize ATAC insert-size classes used by the Rmd QC layer.
write_fragment_class_metrics() {
    local sample="$1"
    local insert_file="$2"
    local out_file="$WORK_DIR/02mapping/qc/${sample}.fragment_classes.tsv"

    printf 'Sample\tClass\tMinSize\tMaxSize\tFragmentCount\tFraction\tMetricSource\n' >"$out_file"
    if [[ ! -f "$insert_file" ]]; then
        warn "Missing insert-size file for fragment classes: $sample"
        return 0
    fi

    awk -v sample="$sample" '
        function emit(class_name, min_size, max_size, count) {
            frac = (total > 0 ? count / total : 0)
            printf "%s\t%s\t%s\t%s\t%d\t%.6f\tinsert_sizes_txt\n", sample, class_name, min_size, max_size, count, frac
        }
        {
            size = $1 + 0
            if(size <= 0) next
            total++
            if(size < 100) nfr++
            else if(size < 180) subnuc++
            else if(size <= 247) mono++
            else if(size >= 315 && size <= 473) di++
            else other++
        }
        END {
            emit("NFR_lt100bp", 1, 99, nfr)
            emit("Subnucleosomal_100_179bp", 100, 179, subnuc)
            emit("Mononucleosome_180_247bp", 180, 247, mono)
            emit("Dinucleosome_315_473bp", 315, 473, di)
            emit("Other", "NA", "NA", other)
        }' "$insert_file" >>"$out_file"
}

aggregate_sample_qc_tsvs() {
    local suffix="$1"
    local out_file="$2"
    local sample file wrote_header=0
    : >"$out_file"
    for sample in "${SAMPLES[@]}"; do
        file="$WORK_DIR/02mapping/qc/${sample}.${suffix}.tsv"
        [[ -f "$file" ]] || continue
        if (( wrote_header == 0 )); then
            awk '1' "$file" >>"$out_file"
            wrote_header=1
        else
            awk 'FNR > 1 {print}' "$file" >>"$out_file"
        fi
    done
    if (( wrote_header == 0 )); then
        warn "No per-sample ${suffix} metrics found for aggregation."
        rm -f "$out_file"
    fi
}

write_optional_tool_status() {
    local out_file="$WORK_DIR/06qc_report/Optional_Tool_Status.tsv"
    local tool cmd status note
    printf 'Tool\tResolvedCommand\tStatus\tNote\n' >"$out_file"
    for tool in preseq ataqv featureCounts bedGraphToBigWig; do
        case "$tool" in
            bedGraphToBigWig) cmd="$BEDGRAPHTOBIGWIG_CMD" ;;
            featureCounts) cmd="$FEATURECOUNTS_CMD" ;;
            *) cmd="$(command -v "$tool" 2>/dev/null || true)" ;;
        esac
        if cmd_available "$cmd"; then
            status="available"
            note="optional QC/tool can be run"
        else
            status="missing"
            note="optional P0/P1 QC recorded as skipped"
        fi
        printf '%s\t%s\t%s\t%s\n' "$tool" "${cmd:-NA}" "$status" "$note" >>"$out_file"
    done
}

tool_version_line() {
    local label="$1"
    local cmd="$2"
    local version="NA"
    if cmd_available "$cmd"; then
        version="$("$cmd" --version 2>&1 | head -n 1 || true)"
        [[ -n "$version" ]] || version="$("$cmd" version 2>&1 | head -n 1 || true)"
        [[ -n "$version" ]] || version="available"
    else
        cmd="NA"
        version="missing"
    fi
    printf '%s\t%s\t%s\n' "$label" "$cmd" "$version"
}

write_software_versions() {
    local out_file="$WORK_DIR/06qc_report/software_versions.tsv"
    {
        printf 'Tool\tResolvedCommand\tVersion\n'
        tool_version_line "fastqc" "$(command -v fastqc 2>/dev/null || true)"
        tool_version_line "fastp" "$(command -v fastp 2>/dev/null || true)"
        tool_version_line "bowtie2" "$(command -v bowtie2 2>/dev/null || true)"
        tool_version_line "samtools" "$(command -v samtools 2>/dev/null || true)"
        tool_version_line "picard" "$(command -v "$PICARD_CMD" 2>/dev/null || printf '%s' "$PICARD_CMD")"
        tool_version_line "macs2" "$MACS2_CMD"
        tool_version_line "bamCoverage" "$BAMCOVERAGE_CMD"
        tool_version_line "computeMatrix" "$COMPUTEMATRIX_CMD"
        tool_version_line "plotProfile" "$PLOTPROFILE_CMD"
        tool_version_line "plotHeatmap" "$PLOTHEATMAP_CMD"
        tool_version_line "idr" "$IDR_CMD"
        tool_version_line "multiqc" "$(command -v multiqc 2>/dev/null || true)"
        tool_version_line "TOBIAS" "$TOBIAS_CMD"
        tool_version_line "bedGraphToBigWig" "$BEDGRAPHTOBIGWIG_CMD"
        tool_version_line "featureCounts" "$FEATURECOUNTS_CMD"
        tool_version_line "preseq" "$(command -v preseq 2>/dev/null || true)"
        tool_version_line "ataqv" "$ATAQV_CMD"
    } >"$out_file"
}

write_tss_method_manifest() {
    local out_file="$WORK_DIR/06qc_report/TSS_Enrichment_Method.tsv"
    local sample matrix profile heatmap
    printf 'Sample\tTSS_Method\tTSS_Reference\tTSS_Threshold_System\tMatrixFile\tProfilePlot\tHeatmapPlot\tInterpretationNote\n' >"$out_file"
    for sample in "${SAMPLES[@]}"; do
        matrix="$WORK_DIR/06qc_report/${sample}_TSS_matrix.gz"
        profile="$WORK_DIR/06qc_report/${sample}_TSS_profile.pdf"
        heatmap="$WORK_DIR/06qc_report/${sample}_TSS_heatmap.pdf"
        printf '%s\tdeepTools computeMatrix reference-point TSS signal-over-flank approximation\t%s\tProject Rmd approximate TSS: >=5 PASS; ENCODE-style mouse TSSE often expects higher values and should be interpreted with method-specific thresholds\t%s\t%s\t%s\tUse this as a within-project QC trend unless ATACseqQC/ENCODE TSSE is run with matching thresholds.\n' \
            "$sample" "$TSS_BED" "$matrix" "$profile" "$heatmap" >>"$out_file"
    done
}

# ------------------------------------------------------------------
# 4. Small Python fallbacks
# ------------------------------------------------------------------

# Fallback for blacklist filtering when bedtools is unavailable.
filter_blacklist_python() {
    local input_peak="$1"
    local blacklist_bed="$2"
    local output_peak="$3"
    python3 - "$input_peak" "$blacklist_bed" "$output_peak" <<'PY'
import bisect
import sys
from collections import defaultdict

peak_file, blacklist_file, output_file = sys.argv[1:4]
intervals = defaultdict(list)
starts = {}

with open(blacklist_file) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, *_ = line.rstrip("\n").split("\t")
        intervals[chrom].append((int(start), int(end)))

for chrom in intervals:
    intervals[chrom].sort()
    starts[chrom] = [item[0] for item in intervals[chrom]]

def overlaps(chrom, start, end):
    if chrom not in intervals:
        return False
    idx = bisect.bisect_left(starts[chrom], end)
    lo = max(0, idx - 4)
    hi = min(len(intervals[chrom]), idx + 4)
    for i in range(lo, hi):
        a, b = intervals[chrom][i]
        if a < end and start < b:
            return True
    return False

with open(peak_file) as src, open(output_file, "w") as dst:
    for line in src:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        if not overlaps(chrom, start, end):
            dst.write(line)
PY
}

# Fallback for bedtools merge when preparing consensus peaks.
merge_intervals_python() {
    local input_bed="$1"
    local output_bed="$2"
    python3 - "$input_bed" "$output_bed" <<'PY'
import sys

input_bed, output_bed = sys.argv[1:3]
records = []
with open(input_bed) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end, *_ = line.rstrip("\n").split("\t")
        records.append((chrom, int(start), int(end)))

def chrom_key(chrom):
    core = chrom[3:] if chrom.startswith("chr") else chrom
    if core.isdigit():
        return (0, int(core))
    if core == "X":
        return (1, 23)
    if core == "Y":
        return (1, 24)
    if core in {"M", "MT"}:
        return (1, 25)
    return (2, core)

records.sort(key=lambda x: (chrom_key(x[0]), x[1], x[2]))
merged = []
for chrom, start, end in records:
    if not merged or merged[-1][0] != chrom or start > merged[-1][2]:
        merged.append([chrom, start, end])
    else:
        merged[-1][2] = max(merged[-1][2], end)

with open(output_bed, "w") as out:
    for chrom, start, end in merged:
        out.write(f"{chrom}\t{start}\t{end}\n")
PY
}

bed_file_has_chr_prefix() {
    local bed_file="$1"
    [[ -f "$bed_file" ]] || return 1
    awk 'NF >= 1 && $1 !~ /^#/ {print ($1 ~ /^chr/ ? "yes" : "no"); exit}' "$bed_file"
}

normalize_bed_chrom_style() {
    local input_bed="$1"
    local output_bed="$2"
    local target_has_chr="$3"
    awk -v target_has_chr="$target_has_chr" 'BEGIN{OFS="\t"} {
        if (NF == 0 || $0 ~ /^#/) {
            print
            next
        }
        chrom = $1
        sub(/^chr/, "", chrom)
        if (chrom == "M") chrom = "MT"
        if (target_has_chr == "yes") {
            if (chrom == "MT") chrom = "M"
            $1 = "chr" chrom
        } else {
            $1 = chrom
        }
        print
    }' "$input_bed" >"$output_bed"
}

# ------------------------------------------------------------------
# 5. Footprinting helpers
# ------------------------------------------------------------------

# Merge all clean peak intervals into one consensus peak set for TOBIAS.
build_consensus_peaks() {
    local raw_bed="$TOBIAS_DIR/consensus_peaks.raw.bed"
    local consensus_bed="$TOBIAS_DIR/consensus_peaks.bed"
    : >"$raw_bed"

    local sample peak_file
    for sample in "${SAMPLES[@]}"; do
        peak_file="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        [[ -f "$peak_file" ]] || continue
        awk 'BEGIN{OFS="\t"} NF >= 3 {print $1, $2, $3}' "$peak_file" >>"$raw_bed"
    done

    if [[ ! -s "$raw_bed" ]]; then
        warn "No clean peak files found; cannot build consensus peaks."
        return 1
    fi

    if command -v bedtools >/dev/null 2>&1; then
        sort -k1,1 -k2,2n "$raw_bed" | bedtools merge -i - >"$consensus_bed"
    else
        merge_intervals_python "$raw_bed" "$consensus_bed"
    fi

    log "Consensus peaks prepared: $consensus_bed ($(wc -l <"$consensus_bed") regions)"
}

# Merge replicate BAMs into one per-condition BAM for footprinting.
merge_group_bam() {
    local group="$1"
    local merged_bam="$TOBIAS_DIR/${group}.merged.bam"
    local -a bam_inputs=()
    local sample bam_file

    for sample in "${SAMPLES[@]}"; do
        if [[ "$(sample_group "$sample")" == "$group" ]]; then
            bam_file="$WORK_DIR/02mapping/bam/${sample}.final.bam"
            [[ -f "$bam_file" ]] && bam_inputs+=("$bam_file")
        fi
    done

    if [[ "${#bam_inputs[@]}" -eq 0 ]]; then
        warn "No BAM files found for group $group."
        return 1
    fi

    rm -f "$merged_bam" "$merged_bam.bai"
    if [[ "${#bam_inputs[@]}" -eq 1 ]]; then
        ln -sf "${bam_inputs[0]}" "$merged_bam"
    else
        samtools merge -@ "$THREADS" -f "$merged_bam" "${bam_inputs[@]}"
    fi
    samtools index -@ "$THREADS" "$merged_bam"
    log "Merged BAM prepared for $group: $merged_bam"
}

# Run TOBIAS bias correction and footprint scoring for one condition.
run_tobias_condition() {
    local group="$1"
    local merged_bam="$TOBIAS_DIR/${group}.merged.bam"
    local prefix="$2"
    local corrected_bw="$TOBIAS_DIR/${prefix}_corrected.bw"
    local footprints_bw="$TOBIAS_DIR/${prefix}_footprints.bw"
    local consensus_bed="$TOBIAS_DIR/consensus_peaks.bed"

    "$TOBIAS_CMD" ATACorrect \
        --bam "$merged_bam" \
        --genome "$GENOME_FASTA" \
        --peaks "$consensus_bed" \
        --outdir "$TOBIAS_DIR" \
        --prefix "$prefix" \
        --cores "$THREADS" \
        >"$TOBIAS_DIR/${prefix}_ATACorrect.log" 2>&1

    "$TOBIAS_CMD" ScoreBigwig \
        --signal "$corrected_bw" \
        --regions "$consensus_bed" \
        --output "$footprints_bw" \
        --cores "$THREADS" \
        >"$TOBIAS_DIR/${prefix}_ScoreBigwig.log" 2>&1
}

# ------------------------------------------------------------------
# 6. Step implementations
# ------------------------------------------------------------------

step1_raw_fastqc() {
    step_log 1 "Raw FastQC"
    local -a raw_fastqs=()
    local sample
    for sample in "${SAMPLES[@]}"; do
        raw_fastqs+=("${SAMPLE_RAW_R1_BY_NAME[$sample]}" "${SAMPLE_RAW_R2_BY_NAME[$sample]}")
    done
    fastqc -t "$THREADS" -o "$WORK_DIR/01rawdata/qc_results" "${raw_fastqs[@]}"
}

# Clean one sample with fastp and write both JSON/HTML QC outputs.
run_fastp_sample() {
    local sample="$1"
    local raw_r1="${SAMPLE_RAW_R1_BY_NAME[$sample]:-$WORK_DIR/01rawdata/${sample}_1.fastq.gz}"
    local raw_r2="${SAMPLE_RAW_R2_BY_NAME[$sample]:-$WORK_DIR/01rawdata/${sample}_2.fastq.gz}"
    [[ -f "$raw_r1" && -f "$raw_r2" ]] || {
        warn "Missing paired FASTQ for $sample, skipping fastp."
        return 0
    }

    fastp \
        -i "$raw_r1" \
        -I "$raw_r2" \
        -o "$WORK_DIR/01rawdata/cleaned_data/${sample}_1.fastq.gz" \
        -O "$WORK_DIR/01rawdata/cleaned_data/${sample}_2.fastq.gz" \
        -h "$WORK_DIR/01rawdata/qc_results/${sample}_fastp.html" \
        -j "$WORK_DIR/01rawdata/qc_results/${sample}_fastp.json" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --length_required "$FASTP_MIN_LENGTH" \
        --cut_front --cut_tail \
        --cut_window_size "$FASTP_CUT_WINDOW_SIZE" \
        --cut_mean_quality "$FASTP_CUT_MEAN_QUALITY"
}

step2_fastp() {
    step_log 2 "fastp trimming"
    run_parallel_for_samples run_fastp_sample
}

step3_clean_fastqc() {
    step_log 3 "Clean FastQC"
    local -a clean_fastqs=()
    local sample
    for sample in "${SAMPLES[@]}"; do
        clean_fastqs+=(
            "$WORK_DIR/01rawdata/cleaned_data/${sample}_1.fastq.gz"
            "$WORK_DIR/01rawdata/cleaned_data/${sample}_2.fastq.gz"
        )
    done
    fastqc -t "$THREADS" -o "$WORK_DIR/01rawdata/qc_results" "${clean_fastqs[@]}"
}

# For one sample: align, remove low-quality/mitochondrial reads, deduplicate, and write QC.
process_sample_bam() {
    local sample="$1"
    local clean_r1="$WORK_DIR/01rawdata/cleaned_data/${sample}_1.fastq.gz"
    local clean_r2="$WORK_DIR/01rawdata/cleaned_data/${sample}_2.fastq.gz"
    local raw_bam="$WORK_DIR/02mapping/bam/${sample}.raw.bam"
    local sorted_bam="$WORK_DIR/02mapping/bam/${sample}.sorted.bam"
    local filtered_bam="$WORK_DIR/02mapping/bam/${sample}.filtered.bam"
    local marked_bam="$WORK_DIR/02mapping/bam/${sample}.marked.bam"
    local final_bam="$WORK_DIR/02mapping/bam/${sample}.final.bam"
    local mito_stats="$WORK_DIR/02mapping/qc/${sample}.mito_stats.txt"

    [[ -f "$clean_r1" && -f "$clean_r2" ]] || {
        warn "Missing cleaned FASTQ for $sample, skipping mapping."
        return 0
    }

    log "Mapping and filtering $sample"

    # 1) Align paired-end reads and keep confidently mapped fragments.
    bowtie2 --very-sensitive -X "$BOWTIE2_MAX_FRAGMENT" -p "$THREADS" -x "$REF_GENOME_INDEX" \
        -1 "$clean_r1" \
        -2 "$clean_r2" \
        2>"$WORK_DIR/02mapping/qc/${sample}.bowtie2.log" | \
        samtools view -bS -q "$ATAC_MAPQ" - >"$raw_bam"

    samtools sort -@ "$THREADS" -o "$sorted_bam" "$raw_bam"
    samtools index "$sorted_bam"

    # 2) Summarize mitochondrial contamination before filtering.
    local total_reads mt_reads mt_pct filtered_reads dedup_reads dup_rate
    total_reads="$(samtools view -c "$sorted_bam")"
    mt_reads="$(samtools idxstats "$sorted_bam" | awk '$1=="MT" || $1=="chrM" {sum += $3} END {print sum+0}')"
    mt_pct="$(awk -v mt="$mt_reads" -v total="$total_reads" 'BEGIN {if (total == 0) print "0.00"; else printf "%.2f", (mt/total)*100}')"
    printf '%s: Total=%s, MT=%s (%s%%)\n' "$sample" "$total_reads" "$mt_reads" "$mt_pct" >"$mito_stats"

    # 3) Keep proper paired fragments on primary nuclear chromosomes only.
    samtools idxstats "$sorted_bam" | \
        awk '$1 ~ /^(chr)?([1-9]|1[0-9]|X|Y)$/ {print $1}' | \
        xargs -r samtools view -b -q "$ATAC_MAPQ" -f 2 -F "$ATAC_SAMTOOLS_EXCLUDE_FLAGS" "$sorted_bam" >"$filtered_bam"

    filtered_reads="$(samtools view -c "$filtered_bam")"
    printf '%s: After filtering=%s\n' "$sample" "$filtered_reads" >>"$mito_stats"
    write_library_complexity_metrics "$sample" "$filtered_bam"

    "$PICARD_CMD" MarkDuplicates \
        I="$filtered_bam" \
        O="$marked_bam" \
        M="$WORK_DIR/02mapping/qc/${sample}.picard_dup_metrics.txt" \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT

    samtools index "$marked_bam"
    local duplicate_reads
    duplicate_reads="$(samtools view -c -f 1024 "$marked_bam")"
    dedup_reads="$((filtered_reads - duplicate_reads))"
    dup_rate="$(awk -v dup="$duplicate_reads" -v filtered="$filtered_reads" 'BEGIN {if (filtered == 0) print "0.00"; else printf "%.2f", (dup/filtered)*100}')"
    printf '%s: Marked=%s, duplicates=%s, final_dedup=%s (DupRate=%s%%)\n' "$sample" "$filtered_reads" "$duplicate_reads" "$dedup_reads" "$dup_rate" >>"$mito_stats"

    # 4) Keep the marked BAM for duplicate-aware QC and remove duplicate-flagged reads only from final BAM.
    samtools view -@ "$THREADS" -b -F 1024 "$marked_bam" >"$final_bam"
    samtools index "$final_bam"
    [[ "$(samtools view -c -f 1024 "$marked_bam")" -gt 0 ]] || warn "$sample marked BAM contains no duplicate flags."
    [[ "$(samtools view -c -f 1024 "$final_bam")" -eq 0 ]] || die "$sample final BAM still contains duplicate flags."
    samtools flagstat "$final_bam" >"$WORK_DIR/02mapping/qc/${sample}.flagstat"
    samtools view "$final_bam" | awk '{if ($9 > 0 && $9 < 1000) print $9}' >"$WORK_DIR/02mapping/qc/${sample}.insert_sizes.txt"
    write_fragment_class_metrics "$sample" "$WORK_DIR/02mapping/qc/${sample}.insert_sizes.txt"

    rm -f "$raw_bam" "$sorted_bam" "$sorted_bam.bai" "$filtered_bam"
}

step4_mapping() {
    step_log 4 "Mapping, filtering, deduplication"
    run_parallel_for_samples process_sample_bam
    aggregate_sample_qc_tsvs "library_complexity" "$WORK_DIR/06qc_report/Library_Complexity_Metrics.tsv"
    aggregate_sample_qc_tsvs "fragment_classes" "$WORK_DIR/06qc_report/Fragment_Size_Classes.tsv"
}

# Call peaks for one sample using the final filtered BAM.
call_peaks_for_sample() {
    local sample="$1"
    [[ -f "$WORK_DIR/02mapping/bam/${sample}.final.bam" ]] || {
        warn "Final BAM missing for $sample, skipping peak calling."
        return 0
    }

    # BAMPE mode uses the observed paired-end fragment extents; keep local lambda correction.
    "$MACS2_CMD" callpeak \
        -t "$WORK_DIR/02mapping/bam/${sample}.final.bam" \
        -f BAMPE \
        -n "$sample" \
        -g mm \
        --outdir "$WORK_DIR/03peaks" \
        --keep-dup all \
        -q "$ATAC_MACS2_QVALUE" \
        --bdg \
        2>"$WORK_DIR/03peaks/${sample}.macs2.log"
}

step5_call_peaks() {
    step_log 5 "MACS2 peak calling"
    run_parallel_for_samples call_peaks_for_sample
}

step6_blacklist() {
    step_log 6 "Blacklist filtering"
    : >"$WORK_DIR/03peaks/blacklist_filter_stats.txt"
    local sample input_peak output_peak before after blacklist_for_sample peak_has_chr tmp_blacklist

    if [[ ! -f "$BLACKLIST_BED" ]]; then
        warn "Blacklist file missing; clean peak files will mirror raw peaks."
        for sample in "${SAMPLES[@]}"; do
            input_peak="$WORK_DIR/03peaks/${sample}_peaks.narrowPeak"
            output_peak="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
            [[ -f "$input_peak" ]] && cp "$input_peak" "$output_peak"
        done
        return 0
    fi

    for sample in "${SAMPLES[@]}"; do
        input_peak="$WORK_DIR/03peaks/${sample}_peaks.narrowPeak"
        output_peak="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        [[ -f "$input_peak" ]] || {
            warn "Peak file missing for $sample, skipping blacklist filtering."
            continue
        }

        peak_has_chr="$(bed_file_has_chr_prefix "$input_peak")"
        blacklist_for_sample="$BLACKLIST_BED"
        tmp_blacklist=""
        if [[ -n "$peak_has_chr" ]]; then
            tmp_blacklist="$(mktemp)"
            normalize_bed_chrom_style "$BLACKLIST_BED" "$tmp_blacklist" "$peak_has_chr"
            blacklist_for_sample="$tmp_blacklist"
        fi

        if command -v bedtools >/dev/null 2>&1; then
            bedtools intersect -v -a "$input_peak" -b "$blacklist_for_sample" >"$output_peak"
        else
            filter_blacklist_python "$input_peak" "$blacklist_for_sample" "$output_peak"
        fi

        [[ -n "$tmp_blacklist" ]] && rm -f "$tmp_blacklist"

        before="$(wc -l <"$input_peak")"
        after="$(wc -l <"$output_peak")"
        printf '%s: Peaks before blacklist filter=%s, after=%s\n' "$sample" "$before" "$after" >>"$WORK_DIR/03peaks/blacklist_filter_stats.txt"
    done
}

step7_bigwig() {
    step_log 7 "BigWig generation"
    local sample bdg_file sorted_bdg bam_file raw_status_file raw_status
    raw_status_file="$WORK_DIR/05bigwig/raw_macs2_bigwig_status.tsv"
    printf 'Sample\tRawBigWig\tStatus\tReason\n' >"$raw_status_file"
    for sample in "${SAMPLES[@]}"; do
        raw_status="not_requested"
        if [[ "$RAW_MACS2_BIGWIG" == "1" ]]; then
            bdg_file="$WORK_DIR/03peaks/${sample}_treat_pileup.bdg"
            sorted_bdg="$WORK_DIR/03peaks/${sample}_treat_pileup.sorted.bdg"
            if [[ -f "$bdg_file" ]]; then
                sort -k1,1 -k2,2n "$bdg_file" | \
                    awk 'BEGIN{OFS="\t"} {
                        if ($1 ~ /^[0-9]+$/ && $1 >= 1 && $1 <= 19) {$1="chr"$1; print}
                        else if ($1 ~ /^[XY]$/) {$1="chr"$1; print}
                        else if ($1 ~ /^chr([0-9]+|X|Y|M)$/) {print}
                    }' >"$sorted_bdg"

                if cmd_available "$BEDGRAPHTOBIGWIG_CMD" && [[ -f "$CHROM_SIZE" ]]; then
                    if ! "$BEDGRAPHTOBIGWIG_CMD" "$sorted_bdg" "$CHROM_SIZE" "$WORK_DIR/05bigwig/${sample}.bw"; then
                        warn "bedGraphToBigWig failed for $sample; skipping raw MACS2 bigWig and continuing with bamCoverage."
                        raw_status="skipped_bedGraphToBigWig_failed"
                    else
                        raw_status="ok"
                    fi
                else
                    warn "bedGraphToBigWig or chrom.sizes missing; skipping ${sample}.bw"
                    raw_status="skipped_missing_bedGraphToBigWig_or_chrom_sizes"
                fi
            else
                warn "Missing bedGraph for $sample, skipping raw MACS2 signal track."
                raw_status="skipped_missing_macs2_bedGraph"
            fi
        fi
        printf '%s\t%s\t%s\t%s\n' "$sample" "$WORK_DIR/05bigwig/${sample}.bw" "$raw_status" "normalized bamCoverage tracks are generated independently" >>"$raw_status_file"

        # Normalized coverage track for downstream visualization/QC.
        bam_file="$WORK_DIR/02mapping/bam/${sample}.final.bam"
        [[ -f "$bam_file" ]] || {
            warn "Missing final BAM for $sample, skipping normalized BigWig."
            continue
        }
        if [[ -x "$BAMCOVERAGE_CMD" ]]; then
            "$BAMCOVERAGE_CMD" \
                -b "$bam_file" \
                -o "$WORK_DIR/05bigwig/${sample}_RPGC.bw" \
                --binSize "$ATAC_BIGWIG_BIN_SIZE" \
                --normalizeUsing RPGC \
                --effectiveGenomeSize "$ATAC_EFFECTIVE_GENOME_SIZE" \
                --extendReads \
                -p "$THREADS" || warn "bamCoverage failed for $sample"
        else
            warn "bamCoverage unavailable; skipping normalized BigWig for $sample"
        fi
    done
}

step8_frip() {
    step_log 8 "FRiP quantification"
    local out_file="$WORK_DIR/04counts/FRiP_scores.txt"
    local method_file="$WORK_DIR/04counts/FRiP_Method.tsv"
    printf 'Sample\tTotal_Reads\tReads_in_Peaks\tFRiP\n' >"$out_file"
    printf 'Sample\tBamFile\tPeakFile\tMethod\tRecommendedThreshold\tInterpretationNote\n' >"$method_file"

    local sample bam_file peak_file total_reads reads_in_peaks frip
    for sample in "${SAMPLES[@]}"; do
        bam_file="$WORK_DIR/02mapping/bam/${sample}.final.bam"
        peak_file="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        [[ -f "$bam_file" && -f "$peak_file" ]] || {
            warn "Missing BAM or peak file for $sample, skipping FRiP."
            continue
        }
        total_reads="$(samtools view -c "$bam_file")"
        reads_in_peaks="$(samtools view -c -M -L "$peak_file" "$bam_file")"
        frip="$(awk -v in_peaks="$reads_in_peaks" -v total="$total_reads" 'BEGIN {if (total == 0) print "NA"; else printf "%.4f", in_peaks/total}')"
        printf '%s\t%s\t%s\t%s\n' "$sample" "$total_reads" "$reads_in_peaks" "$frip" >>"$out_file"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$bam_file" "$peak_file" "read-level alignment overlaps: samtools view -M -L final dedup BAM against blacklist-filtered sample peaks" "NA; no validated read-level acceptance threshold configured" "Each mate is an alignment; not fragment FRiP. Descriptive QC with TSS enrichment and complexity." >>"$method_file"
    done
    cat "$out_file"
}

step9_tss() {
    step_log 9 "TSS enrichment"
    if [[ ! -x "$COMPUTEMATRIX_CMD" || ! -x "$PLOTPROFILE_CMD" || ! -f "$TSS_BED" ]]; then
        warn "deepTools or TSS BED missing; skipping Step 9."
        write_tss_method_manifest
        return 0
    fi

    local sample bw_file
    for sample in "${SAMPLES[@]}"; do
        bw_file="$WORK_DIR/05bigwig/${sample}_RPGC.bw"
        [[ -f "$bw_file" ]] || {
            warn "Missing normalized BigWig for $sample, skipping TSS enrichment."
            continue
        }

        "$COMPUTEMATRIX_CMD" reference-point \
            --referencePoint TSS \
            -b "$ATAC_TSS_FLANK" -a "$ATAC_TSS_FLANK" \
            -R "$TSS_BED" \
            -S "$bw_file" \
            -o "$WORK_DIR/06qc_report/${sample}_TSS_matrix.gz" \
            --missingDataAsZero \
            -p "$THREADS"

        "$PLOTPROFILE_CMD" \
            -m "$WORK_DIR/06qc_report/${sample}_TSS_matrix.gz" \
            -out "$WORK_DIR/06qc_report/${sample}_TSS_profile.pdf" \
            --perGroup \
            --plotTitle "TSS Enrichment - $sample"

        if [[ -x "$PLOTHEATMAP_CMD" ]]; then
            "$PLOTHEATMAP_CMD" \
                -m "$WORK_DIR/06qc_report/${sample}_TSS_matrix.gz" \
                -out "$WORK_DIR/06qc_report/${sample}_TSS_heatmap.pdf" \
                --plotTitle "TSS Enrichment Heatmap - $sample" \
                --refPointLabel TSS \
                --xAxisLabel "gene distance (bp)" \
                --colorList "$ATAC_DEEPTOOLS_HEATMAP_COLORS" \
                --whatToShow "plot, heatmap and colorbar" \
                --heatmapHeight 18 \
                --heatmapWidth 6 \
                --sortRegions descend \
                --sortUsing mean
        else
            warn "plotHeatmap unavailable; skipped TSS heatmap for $sample."
        fi
    done
    write_tss_method_manifest
}

call_relaxed_idr_peaks() {
    local bam_file="$1"
    local prefix="$2"
    local tmp_dir="$3"
    local raw_peak="$tmp_dir/${prefix}_peaks.narrowPeak"
    local clean_peak="$tmp_dir/${prefix}.relaxed.clean.narrowPeak"
    "$MACS2_CMD" callpeak -t "$bam_file" -f BAMPE -n "$prefix" -g mm \
        --outdir "$tmp_dir" --keep-dup all --pvalue "$ATAC_IDR_MACS2_PVALUE" \
        2>"$tmp_dir/${prefix}.macs2.log"
    if [[ -f "$BLACKLIST_BED" ]] && command -v bedtools >/dev/null 2>&1; then
        bedtools intersect -v -a "$raw_peak" -b "$BLACKLIST_BED" >"$clean_peak"
    elif [[ -f "$BLACKLIST_BED" ]]; then
        filter_blacklist_python "$raw_peak" "$BLACKLIST_BED" "$clean_peak"
    else
        cp "$raw_peak" "$clean_peak"
    fi
    printf '%s\n' "$clean_peak"
}

split_bam_pseudoreplicates() {
    local input_bam="$1"
    local output1="$2"
    local output2="$3"
    local seed="$4"
    local tmp1="${output1}.sam"
    local tmp2="${output2}.sam"
    : >"$tmp1"
    : >"$tmp2"
    samtools collate -@ "$THREADS" -u -O "$input_bam" | samtools view -h - | \
        awk -v seed="$seed" -v out1="$tmp1" -v out2="$tmp2" '
          BEGIN{srand(seed); previous=""; target=1}
          /^@/ {print > out1; print > out2; next}
          {if($1 != previous){target=(rand()<0.5 ? 1 : 2); previous=$1}
           if(target==1) print > out1; else print > out2}'
    samtools view -@ "$THREADS" -b "$tmp1" | samtools sort -@ "$THREADS" -o "$output1" -
    samtools view -@ "$THREADS" -b "$tmp2" | samtools sort -@ "$THREADS" -o "$output2" -
    samtools index "$output1"
    samtools index "$output2"
    rm -f "$tmp1" "$tmp2"
}

threshold_idr_output() {
    local unfiltered="$1"
    local thresholded="$2"
    awk -v threshold="$ATAC_IDR_THRESHOLD" 'BEGIN{min_score=-log(threshold)/log(10)} NF>=12 && $12>=min_score' "$unfiltered" >"$thresholded"
}

run_idr_pair() {
    local peak1="$1"
    local peak2="$2"
    local prefix="$3"
    local unfiltered="${prefix}.unfiltered.narrowPeak"
    local thresholded="${prefix}.thresholded.narrowPeak"
    "$IDR_CMD" --samples "$peak1" "$peak2" --input-file-type narrowPeak --rank p.value \
        --output-file "$unfiltered" --plot --log-output-file "${prefix}.log"
    threshold_idr_output "$unfiltered" "$thresholded"
    printf '%s\n' "$thresholded"
}

step10_idr() {
    step_log 10 "ENCODE-style IDR reproducibility"
    local status_file="$WORK_DIR/03peaks/IDR_Reproducibility_Summary.tsv"
    printf 'Group\tStatus\tNt\tNp\tN1\tN2\tRescueRatio\tSelfConsistencyRatio\tConservativeFile\tOptimalFile\tNote\n' >"$status_file"
    if [[ ! -x "$IDR_CMD" || ! -x "$MACS2_CMD" ]]; then
        printf 'ALL\tSKIPPED\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tMissing IDR or MACS2 executable\n' >>"$status_file"
        warn "IDR or MACS2 unavailable; wrote skipped status."
        return 0
    fi

    local group safe_group sample tmp_dir pooled_bam pooled_pr1 pooled_pr2
    local true_threshold pooled_threshold self1_threshold self2_threshold
    local nt np n1 n2 rescue self_ratio status optimal_source conservative_file optimal_file
    local -a group_samples=()

    for group in "${GROUP_LEVELS[@]}"; do
        group_samples=()
        for sample in "${SAMPLES[@]}"; do
            [[ "$(sample_group "$sample")" == "$group" ]] && group_samples+=("$sample")
        done
        safe_group="$(printf '%s' "$group" | tr -c '[:alnum:]_.-' '_')"
        if [[ "${#group_samples[@]}" -ne 2 ]]; then
            printf '%s\tSKIPPED\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tExactly two biological replicates are required; got %s\n' "$group" "${#group_samples[@]}" >>"$status_file"
            warn "Full IDR requires exactly two biological replicates for $group; skipping without selecting samples."
            continue
        fi

        tmp_dir="$WORK_DIR/03peaks/.idr_${safe_group}"
        rm -rf "$tmp_dir"
        mkdir -p "$tmp_dir"
        local bam1="$WORK_DIR/02mapping/bam/${group_samples[0]}.final.bam"
        local bam2="$WORK_DIR/02mapping/bam/${group_samples[1]}.final.bam"
        if [[ ! -f "$bam1" || ! -f "$bam2" ]]; then
            printf '%s\tSKIPPED\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tMissing final BAM for one or both biological replicates\n' "$group" >>"$status_file"
            continue
        fi

        local rep1_peak rep2_peak
        rep1_peak="$(call_relaxed_idr_peaks "$bam1" "${safe_group}_rep1" "$tmp_dir")"
        rep2_peak="$(call_relaxed_idr_peaks "$bam2" "${safe_group}_rep2" "$tmp_dir")"

        local rep1_pr1="$tmp_dir/${safe_group}_rep1_pr1.bam" rep1_pr2="$tmp_dir/${safe_group}_rep1_pr2.bam"
        local rep2_pr1="$tmp_dir/${safe_group}_rep2_pr1.bam" rep2_pr2="$tmp_dir/${safe_group}_rep2_pr2.bam"
        split_bam_pseudoreplicates "$bam1" "$rep1_pr1" "$rep1_pr2" "$ATAC_IDR_SEED"
        split_bam_pseudoreplicates "$bam2" "$rep2_pr1" "$rep2_pr2" "$((ATAC_IDR_SEED + 1))"

        pooled_bam="$tmp_dir/${safe_group}_pooled.bam"
        samtools merge -@ "$THREADS" -f "$pooled_bam" "$bam1" "$bam2"
        samtools index "$pooled_bam"
        pooled_pr1="$tmp_dir/${safe_group}_pooled_pr1.bam"
        pooled_pr2="$tmp_dir/${safe_group}_pooled_pr2.bam"
        split_bam_pseudoreplicates "$pooled_bam" "$pooled_pr1" "$pooled_pr2" "$((ATAC_IDR_SEED + 2))"

        local rep1_pr1_peak rep1_pr2_peak rep2_pr1_peak rep2_pr2_peak pooled_pr1_peak pooled_pr2_peak
        rep1_pr1_peak="$(call_relaxed_idr_peaks "$rep1_pr1" "${safe_group}_rep1_pr1" "$tmp_dir")"
        rep1_pr2_peak="$(call_relaxed_idr_peaks "$rep1_pr2" "${safe_group}_rep1_pr2" "$tmp_dir")"
        rep2_pr1_peak="$(call_relaxed_idr_peaks "$rep2_pr1" "${safe_group}_rep2_pr1" "$tmp_dir")"
        rep2_pr2_peak="$(call_relaxed_idr_peaks "$rep2_pr2" "${safe_group}_rep2_pr2" "$tmp_dir")"
        pooled_pr1_peak="$(call_relaxed_idr_peaks "$pooled_pr1" "${safe_group}_pooled_pr1" "$tmp_dir")"
        pooled_pr2_peak="$(call_relaxed_idr_peaks "$pooled_pr2" "${safe_group}_pooled_pr2" "$tmp_dir")"

        true_threshold="$(run_idr_pair "$rep1_peak" "$rep2_peak" "$WORK_DIR/03peaks/${safe_group}_true_replicates_IDR")"
        self1_threshold="$(run_idr_pair "$rep1_pr1_peak" "$rep1_pr2_peak" "$WORK_DIR/03peaks/${safe_group}_self_rep1_IDR")"
        self2_threshold="$(run_idr_pair "$rep2_pr1_peak" "$rep2_pr2_peak" "$WORK_DIR/03peaks/${safe_group}_self_rep2_IDR")"
        pooled_threshold="$(run_idr_pair "$pooled_pr1_peak" "$pooled_pr2_peak" "$WORK_DIR/03peaks/${safe_group}_pooled_pseudoreplicates_IDR")"

        nt="$(wc -l <"$true_threshold")"; np="$(wc -l <"$pooled_threshold")"
        n1="$(wc -l <"$self1_threshold")"; n2="$(wc -l <"$self2_threshold")"
        rescue="$(awk -v a="$nt" -v b="$np" 'BEGIN{if(a==0||b==0) print "NA"; else if(a>b) printf "%.4f",a/b; else printf "%.4f",b/a}')"
        self_ratio="$(awk -v a="$n1" -v b="$n2" 'BEGIN{if(a==0||b==0) print "NA"; else if(a>b) printf "%.4f",a/b; else printf "%.4f",b/a}')"
        status="FAIL"
        if [[ "$rescue" != "NA" && "$self_ratio" != "NA" ]] && awk -v r="$rescue" -v s="$self_ratio" 'BEGIN{exit !(r<2 && s<2)}'; then status="PASS"; fi

        conservative_file="$WORK_DIR/03peaks/${safe_group}_IDR_conservative.narrowPeak"
        optimal_file="$WORK_DIR/03peaks/${safe_group}_IDR_optimal.narrowPeak"
        cp "$true_threshold" "$conservative_file"
        optimal_source="$true_threshold"; [[ "$np" -gt "$nt" ]] && optimal_source="$pooled_threshold"
        cp "$optimal_source" "$optimal_file"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tglobal IDR<=%s; ratios must both be <2; n=2/group evidence\n' \
            "$group" "$status" "$nt" "$np" "$n1" "$n2" "$rescue" "$self_ratio" "$conservative_file" "$optimal_file" "$ATAC_IDR_THRESHOLD" >>"$status_file"

        if [[ "$ATAC_KEEP_IDR_INTERMEDIATES" != "1" ]]; then rm -rf "$tmp_dir"; fi
    done
}

step11_multiqc() {
    step_log 11 "MultiQC aggregation"
    write_optional_tool_status
    write_software_versions
    if ! command -v multiqc >/dev/null 2>&1; then
        warn "multiqc unavailable; skipping Step 11."
        return 0
    fi
    local -a multiqc_inputs=(
        "$WORK_DIR/01rawdata/qc_results"
        "$WORK_DIR/02mapping/qc"
        "$WORK_DIR/03peaks"
    )
    multiqc "${multiqc_inputs[@]}" \
        -o "$WORK_DIR/06qc_report" \
        --title "ATAC-seq QC Report" \
        --filename atac_multiqc_report.html \
        --force
}

step12_footprinting() {
    step_log 12 "Optional TOBIAS footprint tracks"
    mkdir -p "$TOBIAS_DIR"

    # Always prepare the footprinting inputs so downstream can inspect them even if
    # TOBIAS itself is skipped.
    build_consensus_peaks || return 0
    merge_group_bam "$ATAC_REFERENCE_GROUP" || true
    merge_group_bam "$ATAC_COMPARISON_GROUP" || true

    cat >"$TOBIAS_DIR/Footprinting_Input_Manifest.tsv" <<EOF
type	group	path
consensus_peaks	ALL	$TOBIAS_DIR/consensus_peaks.bed
merged_bam	$ATAC_REFERENCE_GROUP	$TOBIAS_DIR/${ATAC_REFERENCE_GROUP}.merged.bam
merged_bam	$ATAC_COMPARISON_GROUP	$TOBIAS_DIR/${ATAC_COMPARISON_GROUP}.merged.bam
motif_db	ALL	${TOBIAS_MOTIFS:-NA}
genome_fasta	ALL	$GENOME_FASTA
EOF

    if ! cmd_available "$TOBIAS_CMD"; then
        warn "TOBIAS command not found. Prepared consensus peaks and merged BAMs only."
        return 0
    fi
    if [[ ! -f "$GENOME_FASTA" ]]; then
        warn "Genome FASTA missing. Skipping TOBIAS execution."
        return 0
    fi
    if [[ ! -f "$TOBIAS_DIR/${ATAC_REFERENCE_GROUP}.merged.bam" || ! -f "$TOBIAS_DIR/${ATAC_COMPARISON_GROUP}.merged.bam" ]]; then
        warn "Merged BAMs incomplete. Skipping TOBIAS execution."
        return 0
    fi

    # Correct Tn5 bias and derive footprint scores for each group.
    run_tobias_condition "$ATAC_REFERENCE_GROUP" "$TOBIAS_PREFIX_REFERENCE" || warn "TOBIAS processing failed for reference."
    run_tobias_condition "$ATAC_COMPARISON_GROUP" "$TOBIAS_PREFIX_COMPARISON" || warn "TOBIAS processing failed for comparison."

    if [[ ! -f "$TOBIAS_DIR/${TOBIAS_PREFIX_COMPARISON}_footprints.bw" || ! -f "$TOBIAS_DIR/${TOBIAS_PREFIX_REFERENCE}_footprints.bw" ]]; then
        warn "Corrected/footprint bigWig files are incomplete. Track manifest not written."
        return 0
    fi

    cat >"$TOBIAS_DIR/Footprinting_Track_Manifest.tsv" <<EOF
group	track_type	path
$ATAC_COMPARISON_GROUP	corrected_bw	$TOBIAS_DIR/${TOBIAS_PREFIX_COMPARISON}_corrected.bw
$ATAC_COMPARISON_GROUP	footprints_bw	$TOBIAS_DIR/${TOBIAS_PREFIX_COMPARISON}_footprints.bw
$ATAC_REFERENCE_GROUP	corrected_bw	$TOBIAS_DIR/${TOBIAS_PREFIX_REFERENCE}_corrected.bw
$ATAC_REFERENCE_GROUP	footprints_bw	$TOBIAS_DIR/${TOBIAS_PREFIX_REFERENCE}_footprints.bw
EOF
}

write_bindetect_status() {
    printf 'Status\tReason\tGeneratedAt\n%s\t%s\t%s\n' "$1" "$2" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >"$TOBIAS_DIR/BINDetect_Status.tsv"
}

step13_bindetect() {
    step_log 13 "Optional TOBIAS BINDetect motif-wide analysis"
    mkdir -p "$TOBIAS_DIR"
    local provenance_file="$TOBIAS_DIR/BINDetect_Provenance.tsv"
    local signal_comparison="$TOBIAS_DIR/${TOBIAS_PREFIX_COMPARISON}_footprints.bw"
    local signal_reference="$TOBIAS_DIR/${TOBIAS_PREFIX_REFERENCE}_footprints.bw"
    local peaks="$TOBIAS_DIR/consensus_peaks.bed"
    local output="$TOBIAS_DIR/bindetect_results.txt"
    write_bindetect_status "RUNNING" "Current attempt has not completed; historical outputs are not admitted"
    if ! cmd_available "$TOBIAS_CMD"; then
        write_bindetect_status "SKIPPED" "TOBIAS executable missing; historical results are not current"
        warn "TOBIAS command not found. Skipping BINDetect."
        return 0
    fi
    if [[ ! -s "$GENOME_FASTA" || -z "$TOBIAS_MOTIFS" || ! -s "$TOBIAS_MOTIFS" || ! -s "$peaks" || ! -s "$signal_comparison" || ! -s "$signal_reference" ]]; then
        write_bindetect_status "SKIPPED" "Required genome, motif, consensus peak or footprint signal input is missing"
        warn "BINDetect inputs incomplete; run and verify Step 12 before Step 13."
        return 0
    fi
    local temp_manifest role path expected_hash
    temp_manifest="$(mktemp "$TOBIAS_DIR/.bindetect_provenance.XXXXXX")"
    printf 'Role\tPath\tSHA256\n' >"$temp_manifest"
    while IFS=$'\t' read -r role path; do
        printf '%s\t%s\t%s\n' "$role" "$(realpath "$path")" "$(sha256sum "$path" | awk '{print $1}')" >>"$temp_manifest"
    done <<EOF
signal_comparison	$signal_comparison
signal_reference	$signal_reference
genome	$GENOME_FASTA
motifs	$TOBIAS_MOTIFS
peaks	$peaks
EOF
    if ! "$TOBIAS_CMD" BINDetect \
        --motifs "$TOBIAS_MOTIFS" --signals "$signal_comparison" "$signal_reference" \
        --genome "$GENOME_FASTA" --peaks "$peaks" --outdir "$TOBIAS_DIR" \
        --cond_names comparison reference --cores "$THREADS" >"$TOBIAS_DIR/BINDetect.log" 2>&1; then
        rm -f "$temp_manifest"
        write_bindetect_status "FAILED" "BINDetect failed; existing result files are not admitted"
        warn "TOBIAS BINDetect failed."
        return 1
    fi
    if [[ ! -s "$output" || ! "$output" -nt "$temp_manifest" ]]; then
        rm -f "$temp_manifest"
        write_bindetect_status "FAILED" "BINDetect did not produce a newly written nonempty result table"
        return 1
    fi
    # Check inputs again after execution so mutation during a run cannot be hidden.
    while IFS=$'\t' read -r role path expected_hash; do
        [[ "$role" == "Role" ]] && continue
        if [[ "$(sha256sum "$path" | awk '{print $1}')" != "$expected_hash" ]]; then
            rm -f "$temp_manifest"
            write_bindetect_status "FAILED" "An input changed during BINDetect execution"
            return 1
        fi
    done <"$temp_manifest"
    printf 'results\t%s\t%s\n' "$(realpath "$output")" "$(sha256sum "$output" | awk '{print $1}')" >>"$temp_manifest"
    mv "$temp_manifest" "$provenance_file"
    write_bindetect_status "COMPLETED" "SHA256 recorded for all actual inputs and result; comparison versus reference"
}

write_consensus_counts_skipped() {
    local reason="$1"
    local status_file="$WORK_DIR/04counts/Consensus_Count_Matrix_Status.tsv"
    mkdir -p "$WORK_DIR/04counts"
    {
        printf 'Item\tPath\tStatus\tNote\n'
        printf 'consensus_peak_counts\t%s\tSKIPPED\t%s\n' "$WORK_DIR/04counts/consensus_peak_counts.raw.tsv" "$reason"
    } >"$status_file"
    warn "$reason"
}

write_consensus_count_method() {
    local count_tool="$1"
    local peak_set="$2"
    local method_file="$WORK_DIR/04counts/Consensus_Count_Matrix_Method.tsv"
    {
        printf 'Field\tValue\n'
        printf 'ConsensusPeakSet\t%s\n' "$peak_set"
        printf 'PeakInputs\t%s\n' "blacklist-filtered per-sample narrowPeak files in 03peaks"
        printf 'CountingTool\t%s\n' "$count_tool"
        printf 'ReadInputs\t%s\n' "02mapping/bam/*.final.bam"
        printf 'CountingUnit\t%s\n' "paired-end fragments; featureCounts -p --countReadPairs -B -C --primary"
        printf 'MAPQ\t%s\n' ">=$ATAC_MAPQ for featureCounts; BAMs were already MAPQ-filtered upstream"
        printf 'DuplicateHandling\t%s\n' "final BAMs are duplicate-removed upstream"
        printf 'Normalization\t%s\n' "CPM = raw count / sample assigned counts across consensus peaks * 1e6"
        printf 'CoordinateSystem\t%s\n' "BED outputs use 0-based half-open coordinates; SAF conversion uses 1-based inclusive start for featureCounts"
        printf 'GeneratedAt\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    } >"$method_file"
}

build_count_peak_sets() {
    local count_dir="$WORK_DIR/04counts"
    local raw_bed="$count_dir/consensus_peaks_union.raw.bed"
    local union_bed="$count_dir/consensus_peaks_union.bed"
    local min2_bed="$count_dir/consensus_peaks_reproducible_min2.bed"
    local idr_raw="$count_dir/consensus_peaks_idr_based.raw.bed"
    local idr_bed="$count_dir/consensus_peaks_idr_based.bed"
    local presence_out="$count_dir/consensus_peak_presence_absence.tsv"
    local tmp_dir="$count_dir/.consensus_tmp"
    local sample peak_file count_file idr_file

    mkdir -p "$count_dir" "$tmp_dir"
    rm -f "$raw_bed" "$union_bed" "$min2_bed" "$idr_raw" "$idr_bed" "$presence_out"
    rm -f "$tmp_dir"/*
    : >"$raw_bed"

    for sample in "${SAMPLES[@]}"; do
        peak_file="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        if [[ -f "$peak_file" ]]; then
            awk 'BEGIN{OFS="\t"} NF >= 3 {print $1, $2, $3}' "$peak_file" >>"$raw_bed"
        else
            warn "Missing clean peak file for consensus count matrix: $peak_file"
        fi
    done

    if [[ ! -s "$raw_bed" ]]; then
        write_consensus_counts_skipped "No clean peak intervals were available for consensus peak construction."
        return 1
    fi

    if command -v bedtools >/dev/null 2>&1; then
        sort -k1,1 -k2,2n "$raw_bed" | bedtools merge -i - >"$union_bed"
    else
        merge_intervals_python "$raw_bed" "$union_bed"
    fi

    local -a count_files=()
    for sample in "${SAMPLES[@]}"; do
        peak_file="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        count_file="$tmp_dir/${sample}.presence_count.txt"
        count_files+=("$count_file")
        if [[ -f "$peak_file" && -s "$union_bed" && "$(command -v bedtools 2>/dev/null || true)" != "" ]]; then
            bedtools intersect -c -a "$union_bed" -b "$peak_file" | awk '{print $NF}' >"$count_file"
        else
            awk 'BEGIN{OFS="\t"} {print 0}' "$union_bed" >"$count_file"
        fi
    done

    python3 - "$union_bed" "$presence_out" "$min2_bed" "$(IFS=,; echo "${SAMPLES[*]}")" "$(IFS=,; for s in "${SAMPLES[@]}"; do printf '%s,' "$(sample_group "$s")"; done | sed 's/,$//')" "${count_files[@]}" <<'PY'
import csv
import sys

union_bed, presence_out, min2_bed, sample_arg, group_arg, *count_files = sys.argv[1:]
samples = sample_arg.split(",") if sample_arg else []
groups = group_arg.split(",") if group_arg else [""] * len(samples)
group_by_sample = dict(zip(samples, groups))

with open(union_bed) as bed:
    peaks = [line.rstrip("\n").split("\t")[:3] for line in bed if line.strip() and not line.startswith("#")]

counts = []
for path in count_files:
    with open(path) as handle:
        values = [int(float(line.strip() or 0)) for line in handle]
    counts.append(values)

with open(presence_out, "w", newline="") as out, open(min2_bed, "w") as min2:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["peak_id", "CHR", "START", "END", *samples, "PeakSupport", "SupportGroups"])
    for i, (chrom, start, end) in enumerate(peaks):
        flags = []
        present_groups = []
        for sample_idx, sample in enumerate(samples):
            present = 1 if sample_idx < len(counts) and i < len(counts[sample_idx]) and counts[sample_idx][i] > 0 else 0
            flags.append(present)
            if present:
                present_groups.append(group_by_sample.get(sample, ""))
        support = sum(flags)
        support_groups = ",".join(sorted({g for g in present_groups if g}))
        peak_id = f"{chrom}:{start}:{end}"
        writer.writerow([peak_id, chrom, start, end, *flags, support, support_groups])
        if support >= 2:
            min2.write(f"{chrom}\t{start}\t{end}\n")
PY

    : >"$idr_raw"
    for idr_file in "$WORK_DIR/03peaks"/*_IDR_optimal.narrowPeak; do
        [[ -f "$idr_file" ]] || continue
        awk 'BEGIN{OFS="\t"} NF >= 3 {print $1, $2, $3}' "$idr_file" >>"$idr_raw"
    done
    if [[ -s "$idr_raw" ]]; then
        if command -v bedtools >/dev/null 2>&1; then
            sort -k1,1 -k2,2n "$idr_raw" | bedtools merge -i - >"$idr_bed"
        else
            merge_intervals_python "$idr_raw" "$idr_bed"
        fi
    else
        : >"$idr_bed"
    fi

    rm -rf "$tmp_dir"
    log "Consensus peak sets written: $union_bed ($(wc -l <"$union_bed") union peaks), $min2_bed ($(wc -l <"$min2_bed") min2 peaks), $idr_bed ($(wc -l <"$idr_bed") IDR peaks)"
}

write_peak_saf_from_bed() {
    local bed_file="$1"
    local saf_file="$2"
    awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} NF >= 3 {
        peak_id = $1 ":" $2 ":" $3
        print peak_id, $1, $2 + 1, $3, "."
    }' "$bed_file" >"$saf_file"
}

standardize_featurecounts_output() {
    local featurecounts_out="$1"
    local raw_out="$2"
    local cpm_out="$3"
    local library_out="$4"
    local sample_arg="$5"
    local group_arg="$6"
    shift 6
    python3 - "$featurecounts_out" "$raw_out" "$cpm_out" "$library_out" "$sample_arg" "$group_arg" "$@" <<'PY'
import csv
import os
import sys

fc_out, raw_out, cpm_out, library_out, sample_arg, group_arg, *bam_paths = sys.argv[1:]
samples = sample_arg.split(",") if sample_arg else []
groups = group_arg.split(",") if group_arg else [""] * len(samples)

with open(fc_out) as handle:
    reader = csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")
    rows = list(reader)
if not rows:
    raise SystemExit("featureCounts output is empty")

header = rows[0]
data = rows[1:]
if len(header) < 7:
    raise SystemExit("featureCounts output does not contain count columns")

count_col_count = len(header) - 6
if len(samples) != count_col_count:
    samples = [os.path.basename(p).replace(".final.bam", "") for p in header[6:]]
if len(groups) != len(samples):
    groups = [""] * len(samples)

raw_header = ["peak_id", "CHR", "START", "END", "PeakWidth", *samples]
raw_records = []
library_sizes = [0.0] * len(samples)
for row in data:
    if len(row) < 6 + len(samples):
        continue
    peak_id = row[0]
    chrom = row[1]
    start_0 = int(float(row[2])) - 1
    end = int(float(row[3]))
    width = max(0, end - start_0)
    vals = [float(x) for x in row[6:6 + len(samples)]]
    library_sizes = [a + b for a, b in zip(library_sizes, vals)]
    raw_records.append([peak_id, chrom, start_0, end, width, *vals])

with open(raw_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(raw_header)
    for rec in raw_records:
        writer.writerow([int(x) if isinstance(x, float) and x.is_integer() else x for x in rec])

with open(cpm_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(raw_header)
    for rec in raw_records:
        vals = rec[5:]
        cpm = [(v / library_sizes[i] * 1_000_000.0) if library_sizes[i] > 0 else 0.0 for i, v in enumerate(vals)]
        writer.writerow([*rec[:5], *[f"{v:.6f}" for v in cpm]])

with open(library_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Sample", "Group", "BamFile", "AssignedFragmentsInConsensusPeaks"])
    for sample, group, bam, lib in zip(samples, groups, bam_paths, library_sizes):
        writer.writerow([sample, group, bam, int(lib) if float(lib).is_integer() else lib])
PY
}

step14_consensus_counts() {
    step_log 14 "Consensus peak sets and count matrices"
    local count_dir="$WORK_DIR/04counts"
    local union_bed="$count_dir/consensus_peaks_union.bed"
    local saf_file="$count_dir/consensus_peaks_union.saf"
    local fc_out="$count_dir/featureCounts.consensus_peak_counts.txt"
    local raw_out="$count_dir/consensus_peak_counts.raw.tsv"
    local cpm_out="$count_dir/consensus_peak_counts.cpm.tsv"
    local library_out="$count_dir/consensus_peak_counts.library_size.tsv"
    local status_file="$count_dir/Consensus_Count_Matrix_Status.tsv"
    local sample_arg group_arg sample bam_file pending_dir matrix_file
    local -a bam_inputs=()

    mkdir -p "$count_dir" || return 1
    printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tRUNNING\tCurrent attempt is incomplete; existing matrices are not admitted\n' "$raw_out" >"$status_file" || return 1
    if ! build_count_peak_sets; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tConsensus peak construction failed; existing matrices are not current\n' "$raw_out" >"$status_file"
        return 1
    fi
    [[ -s "$union_bed" ]] || {
        write_consensus_counts_skipped "Consensus union peak BED is empty."
        return 0
    }

    for sample in "${SAMPLES[@]}"; do
        bam_file="$WORK_DIR/02mapping/bam/${sample}.final.bam"
        if [[ -f "$bam_file" ]]; then
            bam_inputs+=("$bam_file")
        else
            warn "Missing final BAM for consensus counting: $bam_file"
        fi
    done
    if [[ "${#bam_inputs[@]}" -ne "${#SAMPLES[@]}" ]]; then
        write_consensus_counts_skipped "One or more final BAM files were missing; consensus count matrix was not generated."
        return 0
    fi

    sample_arg="$(IFS=,; echo "${SAMPLES[*]}")"
    group_arg="$(IFS=,; for s in "${SAMPLES[@]}"; do printf '%s,' "$(sample_group "$s")"; done | sed 's/,$//')"

    if ! cmd_available "$FEATURECOUNTS_CMD"; then
        write_consensus_counts_skipped "featureCounts is required for paired-fragment counts; no alignment-count fallback is permitted."
        return 0
    fi
    if ! write_peak_saf_from_bed "$union_bed" "$saf_file" || ! awk 'END {exit(NR < 2)}' "$saf_file"; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tSAF generation failed or produced no feature rows\n' "$raw_out" >"$status_file"
        return 1
    fi
    if ! "$FEATURECOUNTS_CMD" \
        -T "$THREADS" -F SAF -a "$saf_file" -o "$fc_out" \
        -p --countReadPairs -B -C -Q "$ATAC_MAPQ" --primary \
        "${bam_inputs[@]}" >"$count_dir/featureCounts.consensus_peak_counts.log" 2>&1; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tfeatureCounts failed; existing counts are not current\n' "$raw_out" >"$status_file"
        warn "featureCounts failed for consensus peaks; see $count_dir/featureCounts.consensus_peak_counts.log"
        return 1
    fi
    if [[ ! -s "$fc_out" || ! "$fc_out" -nt "$status_file" ]]; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tfeatureCounts did not produce a new nonempty output table\n' "$raw_out" >"$status_file"
        return 1
    fi
    if ! pending_dir="$(mktemp -d "$count_dir/.counts_pending.XXXXXX")"; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tCould not create temporary matrix output directory\n' "$raw_out" >"$status_file"
        return 1
    fi
    if ! standardize_featurecounts_output "$fc_out" "$pending_dir/raw.tsv" "$pending_dir/cpm.tsv" "$pending_dir/library.tsv" "$sample_arg" "$group_arg" "${bam_inputs[@]}"; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tfeatureCounts matrix standardization failed; existing matrices are not current\n' "$raw_out" >"$status_file"
        return 1
    fi
    for matrix_file in "$pending_dir/raw.tsv" "$pending_dir/cpm.tsv" "$pending_dir/library.tsv"; do
        if [[ ! -s "$matrix_file" ]] || ! awk 'END {exit(NR < 2)}' "$matrix_file"; then
            printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tA required new matrix is missing, empty, or header-only\n' "$raw_out" >"$status_file"
            return 1
        fi
    done
    if ! write_consensus_count_method "featureCounts" "$union_bed"; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tCount-method manifest could not be written\n' "$raw_out" >"$status_file"
        return 1
    fi
    if ! mv "$pending_dir/raw.tsv" "$raw_out" ||
       ! mv "$pending_dir/cpm.tsv" "$cpm_out" ||
       ! mv "$pending_dir/library.tsv" "$library_out"; then
        printf 'Item\tPath\tStatus\tNote\nconsensus_peak_counts\t%s\tFAILED\tCould not publish all validated matrices\n' "$raw_out" >"$status_file"
        return 1
    fi
    rmdir "$pending_dir" || true

    {
        printf 'Item\tPath\tStatus\tNote\n'
        printf 'union_peaks\t%s\t%s\t%s\n' "$union_bed" "$(if [[ -s "$union_bed" ]]; then echo PASS; else echo MISSING; fi)" "merged clean peaks from all selected samples"
        printf 'min2_reproducible_peaks\t%s\t%s\t%s\n' "$count_dir/consensus_peaks_reproducible_min2.bed" "$(if [[ -s "$count_dir/consensus_peaks_reproducible_min2.bed" ]]; then echo PASS; else echo MISSING; fi)" "union peaks overlapping peaks in at least two selected samples"
        printf 'idr_based_peaks\t%s\t%s\t%s\n' "$count_dir/consensus_peaks_idr_based.bed" "$(if [[ -s "$count_dir/consensus_peaks_idr_based.bed" ]]; then echo PASS; else echo EMPTY_OR_MISSING; fi)" "merged IDR peaks when IDR outputs exist"
        printf 'presence_absence\t%s\t%s\t%s\n' "$count_dir/consensus_peak_presence_absence.tsv" "$(if [[ -s "$count_dir/consensus_peak_presence_absence.tsv" ]]; then echo PASS; else echo MISSING; fi)" "per-sample peak overlap support for each union region"
        printf 'raw_count_matrix\t%s\t%s\t%s\n' "$raw_out" "$(if [[ -s "$raw_out" ]]; then echo PASS; else echo MISSING; fi)" "raw consensus peak counts by sample"
        printf 'cpm_count_matrix\t%s\t%s\t%s\n' "$cpm_out" "$(if [[ -s "$cpm_out" ]]; then echo PASS; else echo MISSING; fi)" "CPM-normalized consensus peak counts by sample"
        printf 'library_size\t%s\t%s\t%s\n' "$library_out" "$(if [[ -s "$library_out" ]]; then echo PASS; else echo MISSING; fi)" "assigned paired fragments used for CPM normalization"
    } >"$status_file"

    log "Consensus count matrix written: $raw_out ($(($(wc -l <"$raw_out") - 1)) peaks)"
}

step15_ataqv_qc() {
    step_log 15 "Optional ataqv standardized ATAC QC"
    local ataqv_dir="$WORK_DIR/06qc_report/ataqv"
    local status_file="$ataqv_dir/Ataqv_QC_Status.tsv"
    local sample bam_file peak_file out_prefix autosome_file
    mkdir -p "$ataqv_dir"
    printf 'Sample\tStatus\tCommand\tMetricsFile\tNote\n' >"$status_file"

    if ! cmd_available "$ATAQV_CMD"; then
        printf 'ALL\tSKIPPED\tNA\tNA\tataqv executable was not found in PATH/ATAQV_CMD\n' >>"$status_file"
        warn "ataqv executable was not found; wrote skipped status."
        return 0
    fi
    for sample in "${SAMPLES[@]}"; do
        bam_file="$WORK_DIR/02mapping/bam/${sample}.marked.bam"
        peak_file="$WORK_DIR/03peaks/${sample}_peaks.clean.narrowPeak"
        out_prefix="$ataqv_dir/${sample}"
        autosome_file="${out_prefix}.autosomes.txt"
        if [[ ! -f "$bam_file" || ! -f "$peak_file" ]]; then
            printf '%s\tSKIPPED\tNA\tNA\tMissing duplicate-marked BAM or clean peak file\n' "$sample" >>"$status_file"
            continue
        fi
        samtools idxstats "$bam_file" | awk '$1 ~ /^(chr)?([1-9]|1[0-9])$/ {print $1}' >"$autosome_file"
        local -a cmd_args=(
            "--peak-file" "$peak_file"
            "--metrics-file" "${out_prefix}.ataqv.json.gz"
            "--autosomal-reference-file" "$autosome_file"
        )
        [[ -f "$BLACKLIST_BED" ]] && cmd_args+=("--excluded-region-file" "$BLACKLIST_BED")
        [[ -f "$TSS_BED" ]] && cmd_args+=("--tss-file" "$TSS_BED")
        if [[ -n "$ATAQV_EXTRA_ARGS" ]]; then
            # shellcheck disable=SC2206
            local -a extra_args=( $ATAQV_EXTRA_ARGS )
            cmd_args+=("${extra_args[@]}")
        fi
        cmd_args+=("$ATAQV_ORGANISM" "$bam_file")
        if "$ATAQV_CMD" "${cmd_args[@]}" >"${out_prefix}.ataqv.stdout.txt" 2>"${out_prefix}.ataqv.stderr.txt"; then
            printf '%s\tPASS\t%s\t%s\tataqv completed\n' "$sample" "$ATAQV_CMD ${cmd_args[*]}" "${out_prefix}.ataqv.json.gz" >>"$status_file"
        else
            printf '%s\tFAILED\t%s\t%s\tSee per-sample stderr/stdout files\n' "$sample" "$ATAQV_CMD ${cmd_args[*]}" "${out_prefix}.ataqv.json.gz" >>"$status_file"
            warn "ataqv failed for $sample; see ${out_prefix}.ataqv.stderr.txt"
        fi
    done

    python3 - "$ataqv_dir" "$WORK_DIR/06qc_report/Ataqv_QC_Metrics_Long.tsv" <<'PY'
import gzip
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])

def scalar_items(value, prefix=""):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from scalar_items(child, f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(value, list):
        return
    elif isinstance(value, (str, int, float, bool)) or value is None:
        yield prefix, value

with output.open("w") as out:
    out.write("Sample\tMetric\tValue\tSource\n")
    for path in sorted(root.glob("*.ataqv.json.gz")):
        sample = path.name.removesuffix(".ataqv.json.gz")
        try:
            with gzip.open(path, "rt") as handle:
                payload = json.load(handle)
            for metric, value in scalar_items(payload):
                out.write(f"{sample}\t{metric}\t{value}\tataqv_marked_bam\n")
        except Exception as exc:
            out.write(f"{sample}\tparse_error\t{str(exc).replace(chr(9), ' ')}\tataqv_marked_bam\n")
PY
}

# ------------------------------------------------------------------
# 7. Main entry
# ------------------------------------------------------------------

# Keep the main flow intentionally linear so the execution order is easy to follow.
main() {
    parse_args "$@"
    validate_numeric_args
    python3 "$SCRIPT_DIR/scripts/input_audit.py" --root "$WORK_DIR" --prefix ATAC
    [[ -f "$READS_INFO" ]] || die "reads_info.txt not found: $READS_INFO"

    if (( DRY_RUN == 1 )); then
        discover_samples
        log "ATAC-seq upstream pipeline dry run"
        log "Work directory: $WORK_DIR"
        log "Threads: $THREADS | Parallel samples: $PARALLEL_SAMPLES"
        log "Force: $FORCE"
        log "Samples: ${SAMPLES[*]}"
        print_user_config_summary
        log "Selected steps:"
        print_selected_steps
        exit 0
    fi

    python3 "$SCRIPT_DIR/scripts/input_audit.py" --root "$WORK_DIR" --prefix ATAC --record upstream
    ensure_dirs
    prepare_logging

    log "=============================================="
    log "ATAC-seq upstream pipeline"
    log "Work directory: $WORK_DIR"
    log "Selected steps: ${SELECTED_STEPS[*]}"
    log "Threads: $THREADS | Parallel samples: $PARALLEL_SAMPLES"
    log "Force: $FORCE"
    print_user_config_summary
    log "=============================================="

    check_reference_files
    discover_samples
    check_required_tools
    write_sample_manifest
    write_upstream_freshness_signature start
    log "Samples: ${SAMPLES[*]}"

    step_enabled 1 && step1_raw_fastqc
    step_enabled 2 && step2_fastp
    step_enabled 3 && step3_clean_fastqc
    step_enabled 4 && step4_mapping
    step_enabled 5 && step5_call_peaks
    step_enabled 6 && step6_blacklist
    step_enabled 7 && step7_bigwig
    step_enabled 8 && step8_frip
    step_enabled 9 && step9_tss
    step_enabled 10 && step10_idr
    step_enabled 11 && step11_multiqc
    step_enabled 12 && step12_footprinting
    step_enabled 13 && step13_bindetect
    step_enabled 14 && step14_consensus_counts
    step_enabled 15 && step15_ataqv_qc
    write_upstream_freshness_signature finish

    log "=============================================="
    log "Upstream pipeline finished."
    log "Key outputs:"
    log "  BAM:      02mapping/bam/*.final.bam"
    log "  Peaks:    03peaks/*_peaks.clean.narrowPeak"
    log "  BigWig:   05bigwig/*_RPGC.bw"
    log "  QC:       06qc_report/"
    log "  FRiP:     04counts/FRiP_scores.txt"
    log "  Counts:   04counts/consensus_peak_counts.raw.tsv"
    log "  Footprint:07footprinting/"
    log "  BINDetect:07footprinting/bindetect_results.txt (Step 13, optional)"
    log "  ataqv:    06qc_report/ataqv/Ataqv_QC_Status.tsv (Step 15, optional)"
    log "  Manifest: metadata/sample_manifest.tsv"
    log "=============================================="
}

main "$@"
