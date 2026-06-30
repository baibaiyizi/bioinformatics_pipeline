#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# WGBS upstream pipeline for sperm samples (GRCm39)
#
# 用法:
#   bash sperm.sh
#     自动选择 cached/raw/clean 运行模式并运行全部上游 step。
#
#   bash sperm.sh --steps 3,5-8
#     只运行指定 step，不自动前置其它 step。这是唯一的选步骤入口。
#
#   bash sperm.sh --list
#     列出所有可运行 step。
#
# 常用重跑示例:
#   cd /home/h1028/workspace/wgbs/wgbs_sperm
#   nohup bash sperm.sh --steps 3 > nohup.out 2>&1 &
#
# 输入约定:
#   01rawdata/reads_info.txt
#     必需。至少包含 sample、group、raw_fq1/raw_fq2 或 fq1/fq2。
#
#   reference/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa
#     step 3/4/5/6 需要。Bismark 还需要该目录已建好 bisulfite index。
#
# Step 说明:
#   1 raw_fastqc
#     作用: 对 reads_info.txt 中 raw_fq1/raw_fq2 指向的原始 FASTQ 运行 FastQC。
#     输入: 01rawdata/reads_info.txt；原始 FASTQ。
#     输出: 01rawdata/qc_results_GRCm39/*fastqc.html 和 *fastqc.zip。
#     工具: fastqc。
#
#   2 fastp_trim
#     作用: 对原始双端 FASTQ 做质控、接头检测和修剪。
#     输入: 01rawdata/reads_info.txt；原始 FASTQ。
#     输出:
#       01rawdata/cleaned_data_GRCm39/sample_1.clean.fastq.gz
#       01rawdata/cleaned_data_GRCm39/sample_2.clean.fastq.gz
#       01rawdata/qc_results_GRCm39/sample_fastp.html
#       01rawdata/qc_results_GRCm39/sample_fastp.json
#     工具: fastp。
#
#   3 bismark_align_dedup_sort
#     作用: Bismark 比对、去重、排序、建索引，并记录读数指标。
#     输入:
#       01rawdata/cleaned_data_GRCm39/sample_1.clean.fastq.gz
#       01rawdata/cleaned_data_GRCm39/sample_2.clean.fastq.gz
#       reference/GRCm39/ Bismark index
#     输出:
#       02bismark_results_GRCm39/sample.sorted.bam
#       02bismark_results_GRCm39/sample.sorted.bam.bai
#       02bismark_results_GRCm39/sample_1.clean_bismark_bt2_PE_report.txt
#       02bismark_results_GRCm39/sample_1.clean_bismark_bt2_pe.deduplication_report.txt
#       01rawdata/qc_results_GRCm39/sample_alignment_metrics.tsv
#       设置 WGBS_RUN_MBIAS=1 时生成 03methylation_results_GRCm39/bismark_mbias/*M-bias.txt(.gz)
#     工具: bismark、samtools、deduplicate_bismark、可选 bismark_methylation_extractor。
#
#   4 bam2nuc_qc
#     作用: 统计 BAM 核苷酸覆盖组成。该步骤较重。
#     输入: 02bismark_results_GRCm39/sample.sorted.bam；reference/GRCm39。
#     输出:
#       01rawdata/qc_results_GRCm39/sample.bam2nuc.done
#       bam2nuc 在 01rawdata/qc_results_GRCm39/ 下生成的报告文件。
#     工具: bam2nuc。
#
#   5 methylkit_cpg
#     作用: 用 MethylDackel 导出 CpG methylKit 文本。
#     输入: 02bismark_results_GRCm39/sample.sorted.bam；GRCm39 FASTA。
#     输出: 04methykit/sample_CpG.methylKit。
#     工具: MethylDackel。
#
#   6 methylkit_noncg
#     作用: 可选导出 CHG/CHH methylKit 文本。
#     输入: 02bismark_results_GRCm39/sample.sorted.bam；GRCm39 FASTA。
#     输出: 04methykit/sample_CHG* 和 04methykit/sample_CHH* methylKit 文件。
#     默认: 跳过；设置 WGBS_FULL_NONCG=1 后运行实际提取。
#     工具: MethylDackel。
#
#   7 mbias_qc
#     作用: 检查或补生成 Bismark M-bias QC 报告。
#     输入:
#       优先使用 02bismark_results_GRCm39/sample_1.clean_bismark_bt2_pe.deduplicated.bam
#       若只保留 sorted BAM，可设置 WGBS_MBIAS_NAME_SORT_EXISTING_BAM=1 临时 name-sort 后补做。
#     输出: 03methylation_results_GRCm39/bismark_mbias/*M-bias.txt(.gz)。
#     工具: bismark_methylation_extractor、可选 samtools sort -n。
#
#   8 mbias_index
#     作用: 生成 M-bias 报告索引表，便于人工检查是否需要 ignore/trim。
#     输入: 03methylation_results_GRCm39/bismark_mbias/*M-bias.txt*。
#     输出: 01rawdata/qc_results_GRCm39/mbias_ignore_advice/mbias_report_index.tsv。
#
#   9 bs_directionality
#     作用: 从 Bismark PE report 中提取 OT/OB/CTOT/CTOB 方向性计数。
#     输入: 02bismark_results_GRCm39/*_PE_report.txt。
#     输出: 01rawdata/qc_results_GRCm39/sample_bs_directionality.tsv。
#     工具: awk。
#
#   10 insert_size
#     作用: 从 sorted BAM 提取 insert size 分布。该步骤会扫描 BAM。
#     输入: 02bismark_results_GRCm39/sample.sorted.bam。
#     输出: 01rawdata/qc_results_GRCm39/sample_insert_size.tsv。
#     工具: samtools stats。
#
#   11 bismark_report
#     作用: 生成 Bismark 单样本 HTML 报告和跨样本汇总。
#     输入: 02bismark_results_GRCm39/ 下的 Bismark report/dedup report。
#     输出:
#       01rawdata/qc_results_GRCm39/bismark_reports/
#       01rawdata/qc_results_GRCm39/bismark_summary*
#     工具: bismark2report、bismark2summary。
#
#   12 methylkit_rdata
#     作用: 用 R/methylKit 读取 CpG methylKit 文本并保存 RData 缓存。
#     输入: 04methykit/sample_CpG.methylKit；01rawdata/reads_info.txt。
#     输出:
#       04methykit/CpG.RData
#       04methykit/cpg_coverage_threshold_summary.csv
#       可选 04methykit/CpG_mincovN.RData
#     工具: Rscript；R 包 methylKit。
#
#   13 qc_summary
#     作用: 汇总 fastp、Bismark、去重、M-bias 等 QC 指标。
#     输入:
#       01rawdata/qc_results_GRCm39/*_fastp.json
#       01rawdata/qc_results_GRCm39/*_alignment_metrics.tsv
#       02bismark_results_GRCm39/*_PE_report.txt
#       02bismark_results_GRCm39/*deduplication_report.txt
#       03methylation_results_GRCm39/bismark_mbias/*M-bias.txt*
#     输出: 01rawdata/qc_results_GRCm39/qc_read_alignment_summary.csv。
#     工具: python3。
#
#   14 snp_bed
#     作用: 将 dbSNP VCF 转成 BED，供后续 SNP 过滤使用。
#     输入: WGBS_DBSNP_VCF 指向的 VCF/VCF.GZ。
#     输出: 04methykit/snp_filtered/snps.bed。
#     工具: bcftools、awk。
#
#   15 multiqc
#     作用: 汇总 FastQC/fastp/Bismark/M-bias 等 QC 结果。
#     输入:
#       01rawdata/qc_results_GRCm39/
#       02bismark_results_GRCm39/
#       03methylation_results_GRCm39/
#     输出: 01rawdata/qc_results_GRCm39/multiqc/。
#     工具: multiqc。
###############################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${ROOT_DIR}/../.." && pwd)"
READS_INFO="${ROOT_DIR}/01rawdata/reads_info.txt"

REF_BASE_DIR="${REPO_ROOT}/reference/GRCm39"
REF_GENOME_DIR="${REF_BASE_DIR}"
REF_GENOME_FA="${REF_BASE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# 可选参数集中配置。运行前可用同名环境变量覆盖，例如: WGBS_BISMARK_PARALLEL=4 bash sperm.sh --steps 3
WGBS_ENV_BIN="${WGBS_ENV_BIN:-/home/h1028/miniconda3/envs/wgbs/bin}" # 若目录存在，优先使用该 conda 环境里的工具。
WGBS_BISMARK_PARALLEL="${WGBS_BISMARK_PARALLEL:-10}" # Bismark 拆分并行数；设为 1 时改用 Bowtie2 -p，降低 temp BAM 合并风险。
WGBS_BOWTIE2_THREADS="${WGBS_BOWTIE2_THREADS:-10}" # WGBS_BISMARK_PARALLEL=1 时传给 Bowtie2 的线程数。
WGBS_BISMARK_USE_TEMP_DIR="${WGBS_BISMARK_USE_TEMP_DIR:-0}" # Bismark --parallel 模式是否显式指定 temp_${sample} 目录。
WGBS_WRITE_UNMAPPED="${WGBS_WRITE_UNMAPPED:-0}" # 是否输出 unmapped FASTQ；默认关闭以减少大文件和合并风险。
WGBS_SAMTOOLS_SORT_THREADS="${WGBS_SAMTOOLS_SORT_THREADS:-8}" # samtools sort 使用的线程数。
WGBS_REUSE_SORTED_BAM="${WGBS_REUSE_SORTED_BAM:-1}" # 已有 sample.sorted.bam 且 quickcheck 通过时是否复用。
WGBS_REUSE_METHYLKIT="${WGBS_REUSE_METHYLKIT:-1}" # 已有 sample_CpG.methylKit 时是否复用。
WGBS_RUN_BAM2NUC="${WGBS_RUN_BAM2NUC:-0}" # 是否运行 bam2nuc；该步骤会扫描 BAM，默认跳过。
WGBS_REUSE_BAM2NUC="${WGBS_REUSE_BAM2NUC:-1}" # 已有 bam2nuc done 标记时是否复用。
WGBS_RUN_MBIAS="${WGBS_RUN_MBIAS:-0}" # 是否生成 M-bias QC；该步骤会额外扫描 BAM，默认跳过。
WGBS_REUSE_MBIAS="${WGBS_REUSE_MBIAS:-1}" # 已有 M-bias 报告时是否复用。
WGBS_MBIAS_NAME_SORT_EXISTING_BAM="${WGBS_MBIAS_NAME_SORT_EXISTING_BAM:-0}" # 仅有 sorted BAM 时，是否临时 name-sort 后补做 M-bias。
WGBS_MBIAS_SORT_THREADS="${WGBS_MBIAS_SORT_THREADS:-4}" # M-bias 临时 name-sort 使用的线程数。
WGBS_RUN_INSERT_SIZE="${WGBS_RUN_INSERT_SIZE:-0}" # 是否运行 insert size 统计；该步骤会扫描 BAM，默认跳过。
WGBS_FULL_NONCG="${WGBS_FULL_NONCG:-0}" # 是否导出完整 CHG/CHH methylKit 文件，默认只保留 CpG。
WGBS_METHREAD_MINCOV="${WGBS_METHREAD_MINCOV:-4}" # methylKit::methRead 读取 CpG methylKit 文本的最小覆盖度。
WGBS_EXTRA_MINCOV_LEVELS="${WGBS_EXTRA_MINCOV_LEVELS:-}" # 逗号分隔的额外 mincov 缓存，例如 1,5,10；为空则不生成。
WGBS_DBSNP_VCF="${WGBS_DBSNP_VCF:-}" # step 14 使用的 dbSNP VCF/VCF.GZ；为空则跳过 SNP BED。

if [[ -d "${WGBS_ENV_BIN}" ]]; then
  export PATH="${WGBS_ENV_BIN}:${PATH}"
fi

STEPS=(
  "raw_fastqc"
  "fastp_trim"
  "bismark_align_dedup_sort"
  "bam2nuc_qc"
  "methylkit_cpg"
  "methylkit_noncg"
  "mbias_qc"
  "mbias_index"
  "bs_directionality"
  "insert_size"
  "bismark_report"
  "methylkit_rdata"
  "qc_summary"
  "snp_bed"
  "multiqc"
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

log() {
  echo "[INFO] $*"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || reads_info_die "缺少命令: $1"
}

have_command() {
  command -v "$1" >/dev/null 2>&1
}

usage() {
  cat <<EOF
用法:
  bash sperm.sh                 # 自动选择 cached/raw/clean 输入运行全部上游步骤
  bash sperm.sh --list          # 列出上游步骤
  bash sperm.sh --steps 3,5-8   # 只运行指定步骤

步骤:
EOF
  list_steps
  cat <<EOF

说明:
  无参数时，脚本按现有文件自动选择 cached/raw/clean 输入。
  只保留 --steps 这一种选步骤方式；范围和逗号可以组合，例如 --steps 3,5-8。
  例如重新运行 Bismark: bash sperm.sh --steps 3
EOF
}

list_steps() {
  local i
  for i in "${!STEPS[@]}"; do
    printf '%2d  %s\n' "$((i + 1))" "${STEPS[$i]}"
  done
}

split_step_tokens() {
  local token="$1"
  token="${token//,/ }"
  printf '%s\n' ${token}
}

append_step_token() {
  local token="$1"
  local step

  [[ -n "${token}" ]] || return 0
  if [[ "${token}" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    local start="${BASH_REMATCH[1]}"
    local end="${BASH_REMATCH[2]}"
    if (( start < 1 || end > ${#STEPS[@]} || start > end )); then
      reads_info_die "无效步骤范围: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    for step in $(seq "${start}" "${end}"); do
      REQUESTED_STEPS+=("${step}")
    done
  elif [[ "${token}" =~ ^[0-9]+$ ]]; then
    if (( token < 1 || token > ${#STEPS[@]} )); then
      reads_info_die "无效步骤: ${token}，有效范围为 1-${#STEPS[@]}"
    fi
    REQUESTED_STEPS+=("${token}")
  else
    reads_info_die "无法识别步骤参数: ${token}"
  fi
}

select_all_steps() {
  local i

  STEP_SELECTION_MODE="all"
  REQUESTED_STEPS=()
  for i in $(seq 1 "${#STEPS[@]}"); do
    REQUESTED_STEPS+=("${i}")
  done
}

parse_args() {
  REQUESTED_STEPS=()
  STEP_SELECTION_MODE="all"

  if (( $# == 0 )); then
    select_all_steps
  else
    while (( $# > 0 )); do
      case "$1" in
        --help|-h)
          usage
          exit 0
          ;;
        --list)
          list_steps
          exit 0
          ;;
        --steps)
          shift
          (( $# > 0 )) || reads_info_die "--steps 后需要指定步骤，例如 --steps 3,5-8"
          [[ "${STEP_SELECTION_MODE}" != "explicit" ]] || reads_info_die "只能指定一次 --steps"
          STEP_SELECTION_MODE="explicit"
          local token
          while IFS= read -r token; do
            append_step_token "${token}"
          done < <(split_step_tokens "$1")
          shift
          ;;
        --*)
          reads_info_die "未知参数: $1；选步骤只支持 --steps 3,5-8"
          ;;
        *)
          reads_info_die "未知位置参数: $1；选步骤只支持 --steps 3,5-8"
          ;;
      esac
    done
  fi

  (( ${#REQUESTED_STEPS[@]} > 0 )) || reads_info_die "未指定要运行的步骤"

  local -A seen=()
  local -a sorted=()
  local step
  for step in "${REQUESTED_STEPS[@]}"; do
    seen["${step}"]=1
  done

  for step in $(seq 1 "${#STEPS[@]}"); do
    if [[ -n "${seen[${step}]:-}" ]]; then
      sorted+=("${step}")
    fi
  done

  REQUESTED_STEPS=("${sorted[@]}")
  declare -gA REQUESTED_STEP_SET=()
  for step in "${REQUESTED_STEPS[@]}"; do
    REQUESTED_STEP_SET["${step}"]=1
  done
}

step_selected() {
  local step="$1"
  [[ -n "${REQUESTED_STEP_SET[${step}]:-}" ]]
}

step_explicitly_selected() {
  local step="$1"
  [[ "${STEP_SELECTION_MODE}" == "explicit" && -n "${REQUESTED_STEP_SET[${step}]:-}" ]]
}

step_runs_with_raw_inputs() {
  local step="$1"
  step_selected "${step}" && { [[ "${STEP_SELECTION_MODE}" == "explicit" ]] || [[ "${run_mode}" == "raw" ]]; }
}

step_runs_with_clean_inputs() {
  local step="$1"
  step_selected "${step}" && { [[ "${STEP_SELECTION_MODE}" == "explicit" ]] || [[ "${run_mode}" == "raw" || "${run_mode}" == "clean" ]]; }
}

render_requested_steps() {
  local -a step_names=()
  local step
  for step in "${REQUESTED_STEPS[@]}"; do
    step_names+=("${step}:${STEPS[$((step - 1))]}")
  done
  printf '%s' "${step_names[*]}"
}

positive_integer_or_die() {
  local name="$1"
  local value="$2"
  [[ "${value}" =~ ^[1-9][0-9]*$ ]] || reads_info_die "${name} 必须是正整数，当前为: ${value}"
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
  boolean_or_die WGBS_BISMARK_USE_TEMP_DIR
  boolean_or_die WGBS_WRITE_UNMAPPED
  boolean_or_die WGBS_REUSE_SORTED_BAM
  boolean_or_die WGBS_REUSE_METHYLKIT
  boolean_or_die WGBS_RUN_BAM2NUC
  boolean_or_die WGBS_REUSE_BAM2NUC
  boolean_or_die WGBS_RUN_MBIAS
  boolean_or_die WGBS_REUSE_MBIAS
  boolean_or_die WGBS_MBIAS_NAME_SORT_EXISTING_BAM
  boolean_or_die WGBS_RUN_INSERT_SIZE
  boolean_or_die WGBS_FULL_NONCG
  positive_integer_or_die "WGBS_BISMARK_PARALLEL" "${WGBS_BISMARK_PARALLEL}"
  positive_integer_or_die "WGBS_BOWTIE2_THREADS" "${WGBS_BOWTIE2_THREADS}"
  positive_integer_or_die "WGBS_SAMTOOLS_SORT_THREADS" "${WGBS_SAMTOOLS_SORT_THREADS}"
  positive_integer_or_die "WGBS_MBIAS_SORT_THREADS" "${WGBS_MBIAS_SORT_THREADS}"
  positive_integer_or_die "WGBS_METHREAD_MINCOV" "${WGBS_METHREAD_MINCOV}"
}

validate_bam_file() {
  local bam="$1"
  local label="$2"

  [[ -s "${bam}" ]] || reads_info_die "${label} BAM 不存在或为空: ${bam}"
  samtools quickcheck -v "${bam}" || reads_info_die "${label} BAM 未完整写入或已损坏: ${bam}"
}

clean_r1_path() {
  local sample="$1"
  printf '%s/01rawdata/cleaned_data_GRCm39/%s_1.clean.fastq.gz' "${ROOT_DIR}" "${sample}"
}

clean_r2_path() {
  local sample="$1"
  printf '%s/01rawdata/cleaned_data_GRCm39/%s_2.clean.fastq.gz' "${ROOT_DIR}" "${sample}"
}

run_fastp_for_sample() {
  local sample="$1"
  local input_r1="$2"
  local input_r2="$3"

  log "Step 2: fastp 处理 ${sample}"
  fastp \
    -i "${input_r1}" \
    -I "${input_r2}" \
    -o "$(clean_r1_path "${sample}")" \
    -O "$(clean_r2_path "${sample}")" \
    -q 30 -u 40 -l 50 \
    --cut_tail \
    --cut_tail_window_size 1 \
    --cut_tail_mean_quality 30 \
    --detect_adapter_for_pe \
    -w 8 \
    -h "${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}_fastp.html" \
    -j "${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}_fastp.json"
}

process_sample() {
  local sample="$1"
  local r1="$2"
  local r2="$3"
  local bismark_dir="${ROOT_DIR}/02bismark_results_GRCm39"
  local bismark_bam="${bismark_dir}/${sample}_1.clean_bismark_bt2_pe.bam"
  local dedup_bam="${bismark_dir}/${sample}_1.clean_bismark_bt2_pe.deduplicated.bam"
  local sorted_bam="${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam"
  local sorted_bai="${sorted_bam}.bai"
  local temp_dir="${bismark_dir}/temp_${sample}"
  local -a bismark_args=()

  log "处理样本: ${sample}"
  if [[ "${WGBS_REUSE_SORTED_BAM}" == "1" && -f "${sorted_bam}" ]]; then
    if samtools quickcheck "${sorted_bam}" >/dev/null 2>&1; then
      if [[ ! -f "${sorted_bai}" ]]; then
        log "  检测到已有 sorted BAM 但缺少索引，补建 ${sample}.sorted.bam.bai"
        samtools index "${sorted_bam}"
      fi
      log "  检测到已有 sorted BAM，复用 ${sample}.sorted.bam；如需重新比对设置 WGBS_REUSE_SORTED_BAM=0"
      return 0
    fi
    log "  检测到已有 sorted BAM 但 quickcheck 失败，将重新比对 ${sample}"
  fi

  rm -f "${bismark_bam}" "${dedup_bam}"
  rm -f "${bismark_dir}/${sample}_1.clean_bismark_bt2_PE_report.txt"
  rm -f "${bismark_dir}/${sample}_1.clean.fastq.gz_unmapped_reads_1.fq.gz"
  rm -f "${bismark_dir}/${sample}_2.clean.fastq.gz_unmapped_reads_2.fq.gz"
  rm -f "${bismark_dir}/${sample}_1.clean.fastq.gz.temp."*
  rm -f "${bismark_dir}/${sample}_2.clean.fastq.gz.temp."*
  rm -rf "${temp_dir}"

  bismark_args=(
    --genome "${REF_GENOME_DIR}"
    -1 "${r1}"
    -2 "${r2}"
    --output_dir "${bismark_dir}"
    --bowtie2
    --bam
    -L 30
    -N 1
  )

  if [[ "${WGBS_WRITE_UNMAPPED}" == "1" ]]; then
    bismark_args+=(-un)
  fi

  if (( WGBS_BISMARK_PARALLEL > 1 )); then
    bismark_args+=(--parallel "${WGBS_BISMARK_PARALLEL}")
    if [[ "${WGBS_BISMARK_USE_TEMP_DIR}" == "1" ]]; then
      bismark_args+=(--temp_dir "${temp_dir}")
    fi
  elif (( WGBS_BOWTIE2_THREADS > 1 )); then
    bismark_args+=(-p "${WGBS_BOWTIE2_THREADS}")
  fi

  log "  Bismark 比对 ${sample}"
  bismark "${bismark_args[@]}"
  validate_bam_file "${bismark_bam}" "${sample} Bismark"

  local mapped_primary_reads
  mapped_primary_reads="$(
    samtools view -c -F 2308 \
      "${bismark_bam}"
  )"

  log "  去重 ${sample}"
  deduplicate_bismark \
    --bam \
    --paired \
    "${bismark_bam}" \
    --output_dir "${bismark_dir}"
  validate_bam_file "${dedup_bam}" "${sample} deduplicated"

  if [[ "${WGBS_RUN_MBIAS}" == "1" || ( "${STEP_SELECTION_MODE}" == "explicit" && -n "${REQUESTED_STEP_SET[7]:-}" ) ]]; then
    if have_command bismark_methylation_extractor; then
      log "  生成 ${sample} M-bias QC (去重后、甲基化提取前)"
      mkdir -p "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias"
      run_mbias_for_sample "${sample}"
    else
      log "  未找到 bismark_methylation_extractor，跳过 ${sample} M-bias QC"
    fi
  else
    log "  跳过 ${sample} M-bias QC；如需生成设置 WGBS_RUN_MBIAS=1 或显式指定 --steps 3,7"
  fi

  log "  排序与索引 ${sample}"
  samtools sort -@ "${WGBS_SAMTOOLS_SORT_THREADS}" \
    -o "${sorted_bam}" \
    "${dedup_bam}"
  validate_bam_file "${sorted_bam}" "${sample} sorted"
  samtools index "${sorted_bam}"

  local deduplicated_primary_reads
  deduplicated_primary_reads="$(
    samtools view -c -F 2308 "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam"
  )"

  printf 'sample\tmapped_primary_reads\tdeduplicated_primary_reads\n%s\t%s\t%s\n' \
    "${sample}" \
    "${mapped_primary_reads}" \
    "${deduplicated_primary_reads}" \
    > "${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}_alignment_metrics.tsv"

  log "  清理 ${sample} 中间文件"
  rm -f "${bismark_bam}"
  rm -f "${dedup_bam}"
  rm -rf "${temp_dir}"
}

run_mbias_for_sample() {
  local sample="$1"
  local mbias_input="${ROOT_DIR}/02bismark_results_GRCm39/${sample}_1.clean_bismark_bt2_pe.deduplicated.bam"
  local temp_name_sorted=""
  local existing_mbias=""

  existing_mbias="$(compgen -G "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias/${sample}*M-bias.txt*" 2>/dev/null | head -1 || true)"
  if [[ "${WGBS_REUSE_MBIAS}" == "1" && -n "${existing_mbias}" && -s "${existing_mbias}" ]]; then
    log "  ${sample}: 检测到已有 M-bias 报告，复用；如需重算设置 WGBS_REUSE_MBIAS=0"
    return 0
  fi

  if [[ ! -f "${mbias_input}" ]]; then
    if [[ "${WGBS_MBIAS_NAME_SORT_EXISTING_BAM}" == "1" && -f "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam" ]]; then
      temp_name_sorted="${ROOT_DIR}/02bismark_results_GRCm39/${sample}.name_sorted_for_mbias.bam"
      log "  ${sample}: 缺少未坐标排序 deduplicated BAM，临时按 read name 排序以生成 M-bias"
      samtools sort -n -@ "${WGBS_MBIAS_SORT_THREADS}" \
        -o "${temp_name_sorted}" \
        "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam"
      mbias_input="${temp_name_sorted}"
    else
      log "  ${sample}: 缺少 read-name 顺序的 deduplicated BAM，跳过 M-bias；如确需补做，设置 WGBS_MBIAS_NAME_SORT_EXISTING_BAM=1"
      return 0
    fi
  fi

  bismark_methylation_extractor \
    --paired-end \
    --mbias_only \
    --gzip \
    --output "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias" \
    "${mbias_input}" || \
    log "  ${sample} M-bias 报告生成失败，继续后续步骤"

  if [[ -n "${temp_name_sorted}" ]]; then
    rm -f "${temp_name_sorted}"
  fi
}

parse_args "$@"
validate_config
load_reads_info "${READS_INFO}" "${ROOT_DIR}" 0

have_fastqs=1
have_clean_fastqs=1
have_methylkit_files=1
for idx in "${!READS_INFO_SAMPLES[@]}"; do
  if [[ ! -f "${READS_INFO_RESOLVED_FQ1[$idx]}" || ! -f "${READS_INFO_RESOLVED_FQ2[$idx]}" ]]; then
    have_fastqs=0
  fi
  if [[ ! -f "$(clean_r1_path "${READS_INFO_SAMPLES[$idx]}")" || ! -f "$(clean_r2_path "${READS_INFO_SAMPLES[$idx]}")" ]]; then
    have_clean_fastqs=0
  fi
done

for sample in "${READS_INFO_SAMPLES[@]}"; do
  if [[ ! -f "${ROOT_DIR}/04methykit/${sample}_CpG.methylKit" ]]; then
    have_methylkit_files=0
  fi
done

run_mode="cached"
if [[ "${have_methylkit_files}" == "1" ]]; then
  run_mode="cached"
elif [[ "${have_fastqs}" == "1" ]]; then
  run_mode="raw"
elif [[ "${have_clean_fastqs}" == "1" ]]; then
  run_mode="clean"
else
  reads_info_die "原始 FASTQ 和 clean FASTQ 都不可用，且 04methykit methylKit 文本缓存不完整，无法继续运行"
fi

if step_runs_with_clean_inputs 3 || step_runs_with_clean_inputs 4 || step_runs_with_clean_inputs 5 || step_runs_with_clean_inputs 6; then
  [[ -d "${REF_GENOME_DIR}" ]] || reads_info_die "找不到参考基因组目录: ${REF_GENOME_DIR}"
  [[ -f "${REF_GENOME_FA}" ]] || reads_info_die "找不到参考基因组 FASTA: ${REF_GENOME_FA}"
fi

if step_runs_with_raw_inputs 1 || step_runs_with_raw_inputs 2; then
  [[ "${have_fastqs}" == "1" ]] || reads_info_die "step 1/2 需要 reads_info 指向的原始 FASTQ 均存在"
fi

if step_runs_with_raw_inputs 1; then
  require_command fastqc
fi

if step_runs_with_raw_inputs 2; then
  require_command fastp
fi

if step_runs_with_clean_inputs 3; then
  if [[ "${have_clean_fastqs}" != "1" ]] && ! step_runs_with_raw_inputs 2; then
    reads_info_die "step 3 需要完整 clean FASTQ；可先运行 step 2 生成 01rawdata/cleaned_data_GRCm39/*clean.fastq.gz"
  fi
  require_command bismark
  require_command samtools
  require_command deduplicate_bismark
fi

if step_runs_with_clean_inputs 4 && { step_explicitly_selected 4 || [[ "${WGBS_RUN_BAM2NUC}" == "1" ]]; }; then
  require_command bam2nuc
fi

if step_runs_with_clean_inputs 5 || step_runs_with_clean_inputs 6; then
  require_command MethylDackel
fi

if step_selected 7 && { step_explicitly_selected 7 || [[ "${WGBS_RUN_MBIAS}" == "1" ]]; }; then
  require_command bismark_methylation_extractor
fi

if step_selected 10 && { step_explicitly_selected 10 || [[ "${WGBS_RUN_INSERT_SIZE}" == "1" ]]; }; then
  require_command samtools
fi

if step_selected 12; then
  require_command Rscript
fi

if step_selected 13; then
  require_command python3
fi

log "运行模式: ${run_mode}"
if [[ "${run_mode}" == "raw" || "${run_mode}" == "clean" ]]; then
  log "Bismark 设置: WGBS_BISMARK_PARALLEL=${WGBS_BISMARK_PARALLEL}, WGBS_BOWTIE2_THREADS=${WGBS_BOWTIE2_THREADS}, WGBS_WRITE_UNMAPPED=${WGBS_WRITE_UNMAPPED}"
fi
log "请求步骤: $(render_requested_steps)"

mkdir -p "${ROOT_DIR}/01rawdata/qc_results_GRCm39"
mkdir -p "${ROOT_DIR}/01rawdata/cleaned_data_GRCm39"
mkdir -p "${ROOT_DIR}/02bismark_results_GRCm39"
mkdir -p "${ROOT_DIR}/03methylation_results_GRCm39"
mkdir -p "${ROOT_DIR}/04methykit"

log "读取 reads_info: ${READS_INFO}"
log "分组顺序: ${READS_INFO_GROUP1_NAME} -> ${READS_INFO_GROUP2_NAME}"
log "样本数量: ${#READS_INFO_SAMPLES[@]}"

if step_runs_with_raw_inputs 1; then
  log "Step 1: 原始数据 FastQC"
  fastqc_inputs=()
  for idx in "${!READS_INFO_SAMPLES[@]}"; do
    fastqc_inputs+=("${READS_INFO_RESOLVED_FQ1[$idx]}")
    fastqc_inputs+=("${READS_INFO_RESOLVED_FQ2[$idx]}")
  done
  fastqc -t 20 -o "${ROOT_DIR}/01rawdata/qc_results_GRCm39" "${fastqc_inputs[@]}"
elif step_selected 1; then
  log "Step 1: 当前 run_mode=${run_mode}，跳过原始数据 FastQC"
fi

if step_runs_with_raw_inputs 2; then
  log "Step 2: fastp 质控修剪"
  for current_group in "${READS_INFO_GROUP_LEVELS[@]}"; do
    log "开始处理 fastp 分组: ${current_group}"
    for idx in "${!READS_INFO_SAMPLES[@]}"; do
      [[ "${READS_INFO_GROUPS[$idx]}" == "${current_group}" ]] || continue
      run_fastp_for_sample \
        "${READS_INFO_SAMPLES[$idx]}" \
        "${READS_INFO_RESOLVED_FQ1[$idx]}" \
        "${READS_INFO_RESOLVED_FQ2[$idx]}"
    done
  done
elif step_selected 2; then
  log "Step 2: 当前 run_mode=${run_mode}，跳过 fastp"
fi

if step_runs_with_clean_inputs 3; then
  log "Step 3: Bismark 比对、去重、排序与索引"
  for current_group in "${READS_INFO_GROUP_LEVELS[@]}"; do
    log "开始处理分组: ${current_group}"
    for idx in "${!READS_INFO_SAMPLES[@]}"; do
      [[ "${READS_INFO_GROUPS[$idx]}" == "${current_group}" ]] || continue
      process_sample \
        "${READS_INFO_SAMPLES[$idx]}" \
        "$(clean_r1_path "${READS_INFO_SAMPLES[$idx]}")" \
        "$(clean_r2_path "${READS_INFO_SAMPLES[$idx]}")"
    done
  done
elif step_selected 3; then
  log "Step 3: 当前 run_mode=${run_mode}，跳过 Bismark 比对"
fi

if step_runs_with_clean_inputs 4; then
  if [[ "${WGBS_RUN_BAM2NUC}" == "1" || "${STEP_SELECTION_MODE}" == "explicit" ]]; then
    log "Step 4: bam2nuc 核苷酸覆盖度分析"
    for sample in "${READS_INFO_SAMPLES[@]}"; do
      bam2nuc_done="${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}.bam2nuc.done"
      if [[ "${WGBS_REUSE_BAM2NUC}" == "1" && -s "${bam2nuc_done}" ]]; then
        log "  ${sample}: 检测到已有 bam2nuc 完成标记，复用；如需重跑设置 WGBS_REUSE_BAM2NUC=0"
        continue
      fi
      if bam2nuc \
          --genome_folder "${REF_GENOME_DIR}" \
          --dir "${ROOT_DIR}/01rawdata/qc_results_GRCm39" \
          "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam"; then
        printf 'completed\t%s\n' "$(date -Is)" > "${bam2nuc_done}"
      else
        log "  ${sample} bam2nuc 分析失败，继续后续步骤"
      fi
    done
  else
    log "Step 4: 默认跳过 bam2nuc 重型 QC；如需运行设置 WGBS_RUN_BAM2NUC=1 或显式指定 --steps 4"
  fi
elif step_selected 4; then
  log "Step 4: 当前 run_mode=${run_mode}，跳过 bam2nuc"
fi

if step_runs_with_clean_inputs 5; then
  log "Step 5: 使用 MethylDackel 提取 methylKit 格式 (CpG)"
  for sample in "${READS_INFO_SAMPLES[@]}"; do
    if [[ "${WGBS_REUSE_METHYLKIT}" == "1" && -s "${ROOT_DIR}/04methykit/${sample}_CpG.methylKit" ]]; then
      log "  检测到已有 ${sample}_CpG.methylKit，复用；如需重新提取设置 WGBS_REUSE_METHYLKIT=0"
      continue
    fi
    log "Processing ${sample} for methylKit format (CpG)"
    MethylDackel extract \
      --methylKit \
      -@ 10 \
      "${REF_GENOME_FA}" \
      "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam" \
      -o "${ROOT_DIR}/04methykit/${sample}"
  done
elif step_selected 5; then
  log "Step 5: 当前 run_mode=${run_mode}，跳过 CpG methylKit 提取"
fi

if step_runs_with_clean_inputs 6; then
  # 转换率下界优先从 Bismark report 的 CHG/CHH 甲基化率估算。
  # 如需完整 CHG/CHH methylKit 文件，设置 WGBS_FULL_NONCG=1 后重新运行。
  if [[ "${WGBS_FULL_NONCG}" == "1" ]]; then
    log "Step 6: 提取 CHG/CHH methylKit 格式"
    for sample in "${READS_INFO_SAMPLES[@]}"; do
      log "Processing ${sample} for methylKit format (CHG - full)"
      MethylDackel extract \
        --methylKit \
        --CHG \
        -@ 10 \
        "${REF_GENOME_FA}" \
        "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam" \
        -o "${ROOT_DIR}/04methykit/${sample}_CHG" || \
        log "  ${sample} CHG 提取失败，继续后续步骤"
    done
    for sample in "${READS_INFO_SAMPLES[@]}"; do
      log "Processing ${sample} for methylKit format (CHH - full)"
      MethylDackel extract \
        --methylKit \
        --CHH \
        -@ 10 \
        "${REF_GENOME_FA}" \
        "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam" \
        -o "${ROOT_DIR}/04methykit/${sample}_CHH" || \
        log "  ${sample} CHH 提取失败，继续后续步骤"
    done
  else
    log "Step 6: 跳过逐位点 CHG/CHH 提取；如需完整非 CpG 结果设置 WGBS_FULL_NONCG=1"
  fi
elif step_selected 6; then
  log "Step 6: 当前 run_mode=${run_mode}，跳过 CHG/CHH 提取"
fi

if step_selected 7; then
  log "Step 7: 检查/补生成 Bismark M-bias QC 报告"
  missing_mbias_samples=()
  for sample in "${READS_INFO_SAMPLES[@]}"; do
    mbias_found="$(compgen -G "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias/${sample}*M-bias.txt*" 2>/dev/null | head -1 || true)"
    if [[ -f "${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam" && ( -z "${mbias_found}" || ! -s "${mbias_found}" ) ]]; then
      missing_mbias_samples+=("${sample}")
    fi
  done

  if (( ${#missing_mbias_samples[@]} > 0 )); then
    if [[ "${WGBS_RUN_MBIAS}" != "1" && "${STEP_SELECTION_MODE}" != "explicit" ]]; then
      log "Step 7: 检测到缺失 M-bias 报告；默认跳过，设置 WGBS_RUN_MBIAS=1 或显式运行 --steps 7 可补生成"
    elif have_command bismark_methylation_extractor; then
      log "Step 7: 检测到缺失 M-bias 报告；仅在保留 deduplicated BAM 或 WGBS_MBIAS_NAME_SORT_EXISTING_BAM=1 时补生成"
      mkdir -p "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias"
      for sample in "${missing_mbias_samples[@]}"; do
        run_mbias_for_sample "${sample}"
      done
    else
      log "Step 7: 缺少 bismark_methylation_extractor，无法补生成 M-bias QC 报告"
    fi
  else
    log "Step 7: 未发现缺失的 M-bias 报告"
  fi
fi

if step_selected 8; then
  # M-bias 是否需要 ignore/trim 应结合图形人工判断；默认不自动推断参数。
  log "Step 8: 索引已有 M-bias 报告 (不生成新报告)"
mkdir -p "${ROOT_DIR}/01rawdata/qc_results_GRCm39/mbias_ignore_advice"
echo -e "sample\tmbias_report\tstatus\tnote" > "${ROOT_DIR}/01rawdata/qc_results_GRCm39/mbias_ignore_advice/mbias_report_index.tsv"
for sample in "${READS_INFO_SAMPLES[@]}"; do
  mbias_file="$(compgen -G "${ROOT_DIR}/03methylation_results_GRCm39/bismark_mbias/${sample}*M-bias.txt*" 2>/dev/null | head -1 || true)"
  if [[ -n "${mbias_file}" && -s "${mbias_file}" ]]; then
    echo -e "${sample}\t${mbias_file}\tfound\tinspect plot/table before setting ignore parameters" >> "${ROOT_DIR}/01rawdata/qc_results_GRCm39/mbias_ignore_advice/mbias_report_index.tsv"
  else
    echo -e "${sample}\t\tmissing\tno M-bias report found" >> "${ROOT_DIR}/01rawdata/qc_results_GRCm39/mbias_ignore_advice/mbias_report_index.tsv"
  fi
done
log "  M-bias 索引表: 01rawdata/qc_results_GRCm39/mbias_ignore_advice/mbias_report_index.tsv"
fi

if step_selected 9; then
log "Step 9: Bisulfite 方向性检查"
for sample in "${READS_INFO_SAMPLES[@]}"; do
  report_file="$(compgen -G "${ROOT_DIR}/02bismark_results_GRCm39/${sample}*PE_report.txt" 2>/dev/null | head -1 || true)"
  if [[ -n "${report_file}" && -f "${report_file}" ]]; then
    awk -v sample="${sample}" '
      BEGIN {
        OFS="\t";
        print "sample", "strand_pattern", "read_pairs", "strand_note";
      }
      /^(CT\/GA\/CT|GA\/CT\/CT|GA\/CT\/GA|CT\/GA\/GA):/ {
        pattern=$1;
        sub(":", "", pattern);
        note=$0;
        sub(/^[^\t]+\t[0-9]+\t*/, "", note);
        print sample, pattern, $2, note;
      }
    ' "${report_file}" > "${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}_bs_directionality.tsv"
  fi
done
fi

if step_selected 10; then
if [[ "${WGBS_RUN_INSERT_SIZE}" == "1" || "${STEP_SELECTION_MODE}" == "explicit" ]]; then
  log "Step 10: 提取 Insert size 分布"
  for sample in "${READS_INFO_SAMPLES[@]}"; do
    sorted_bam="${ROOT_DIR}/02bismark_results_GRCm39/${sample}.sorted.bam"
    if [[ -f "${sorted_bam}" ]]; then
      samtools view -c -f 2 "${sorted_bam}" > /dev/null 2>&1 || true
      samtools stats "${sorted_bam}" 2>/dev/null | \
        grep "^IS" | \
        awk -F'\t' '{print $2"\t"$3}' > \
        "${ROOT_DIR}/01rawdata/qc_results_GRCm39/${sample}_insert_size.tsv" || \
        log "  ${sample} insert size 提取失败"
    fi
  done
else
  log "Step 10: 默认跳过 Insert size 全 BAM 扫描；如需运行设置 WGBS_RUN_INSERT_SIZE=1 或显式指定 --steps 10"
fi
fi

if step_selected 11; then
if have_command bismark2report; then
  log "Step 11: 生成 bismark2report HTML 报告"
  mkdir -p "${ROOT_DIR}/01rawdata/qc_results_GRCm39/bismark_reports"
  bismark2report \
    --dir "${ROOT_DIR}/02bismark_results_GRCm39" \
    --output_dir "${ROOT_DIR}/01rawdata/qc_results_GRCm39/bismark_reports" 2>/dev/null || \
    log "  bismark2report 生成失败"
else
  log "Step 11: 未找到 bismark2report，跳过"
fi

if have_command bismark2summary; then
  log "Step 11: 生成 bismark2summary 跨样本汇总"
  bismark2summary \
    "${ROOT_DIR}/02bismark_results_GRCm39" \
    -o "${ROOT_DIR}/01rawdata/qc_results_GRCm39/bismark_summary" 2>/dev/null || \
    log "  bismark2summary 生成失败"
else
  log "Step 11: 未找到 bismark2summary，跳过"
fi
fi

if step_selected 12; then
if [[ -f "${ROOT_DIR}/04methykit/CpG.RData" && -z "${WGBS_EXTRA_MINCOV_LEVELS}" ]]; then
  log "Step 12: 检测到现有 04methykit/CpG.RData，跳过重建"
else
  log "Step 12: 运行 R 读取 methylKit 并保存 CpG.RData / coverage 敏感性缓存"
  tmp_r_script="$(mktemp "${ROOT_DIR}/tmp.methylkit.XXXXXX.R")"
  cat <<'EOF' > "${tmp_r_script}"
library(methylKit)

read_reads_info <- function(path) {
  reads_info <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("sample", "group")
  missing_cols <- setdiff(required_cols, colnames(reads_info))
  if (length(missing_cols) > 0) {
    stop("reads_info.txt 缺少列: ", paste(missing_cols, collapse = ", "))
  }

  if (!"raw_fq1" %in% colnames(reads_info) && "fq1" %in% colnames(reads_info)) {
    reads_info$raw_fq1 <- reads_info$fq1
  }
  if (!"raw_fq2" %in% colnames(reads_info) && "fq2" %in% colnames(reads_info)) {
    reads_info$raw_fq2 <- reads_info$fq2
  }
  if (!all(c("raw_fq1", "raw_fq2") %in% colnames(reads_info))) {
    stop("reads_info.txt 缺少 raw_fq1/raw_fq2 列")
  }

  reads_info <- reads_info[!duplicated(reads_info$sample), , drop = FALSE]
  reads_info
}

root_dir <- normalizePath(Sys.getenv("WGBS_ROOT_DIR"), winslash = "/", mustWork = TRUE)
reads_info_path <- file.path(root_dir, "01rawdata", "reads_info.txt")
reads_info <- read_reads_info(reads_info_path)
group_levels <- unique(reads_info$group)
if (length(group_levels) != 2) {
  stop("reads_info.txt 必须恰好包含两个分组，当前得到: ", paste(group_levels, collapse = ", "))
}

sample_ids <- reads_info$sample
file_list <- file.path(root_dir, "04methykit", paste0(sample_ids, "_CpG.methylKit"))
missing_files <- file_list[!file.exists(file_list)]
if (length(missing_files) > 0) {
  stop("缺少 methylKit 文件: ", paste(missing_files, collapse = ", "))
}

treatment <- ifelse(reads_info$group == group_levels[[1]], 0L, 1L)
primary_mincov <- as.integer(Sys.getenv("WGBS_METHREAD_MINCOV", unset = "4"))
if (is.na(primary_mincov) || primary_mincov < 1) {
  stop("WGBS_METHREAD_MINCOV 必须是 >= 1 的整数")
}
extra_mincov_env <- Sys.getenv("WGBS_EXTRA_MINCOV_LEVELS", unset = "")
extra_mincov_levels <- integer()
if (nzchar(extra_mincov_env)) {
  extra_mincov_levels <- as.integer(strsplit(extra_mincov_env, ",", fixed = TRUE)[[1]])
  extra_mincov_levels <- sort(unique(extra_mincov_levels[!is.na(extra_mincov_levels) & extra_mincov_levels >= 1]))
}

read_methylkit_obj <- function(mincov_value) {
  methRead(
    location = as.list(file_list),
    sample.id = as.list(sample_ids),
    assembly = "GRCm39",
    treatment = treatment,
    context = "CpG",
    mincov = mincov_value,
    pipeline = "amp"
  )
}

cat("正在读取 methylKit 文本文件，主缓存 mincov=", primary_mincov, "...\n", sep = "")
obj.cpg <- read_methylkit_obj(primary_mincov)

if (length(obj.cpg) != length(sample_ids)) {
  stop("CpG 对象样本数与 reads_info 不一致: ", length(obj.cpg), " vs ", length(sample_ids))
}

cat("读取完成，正在保存为 RData...\n")
print(obj.cpg)
save(obj.cpg, file = file.path(root_dir, "04methykit", "CpG.RData"))

for (extra_mincov in setdiff(extra_mincov_levels, primary_mincov)) {
  extra_file <- file.path(root_dir, "04methykit", paste0("CpG_mincov", extra_mincov, ".RData"))
  cat("正在生成额外 coverage 敏感性缓存: ", basename(extra_file), "\n", sep = "")
  obj.extra <- read_methylkit_obj(extra_mincov)
  save(obj.extra, file = extra_file)
  rm(obj.extra)
}

qc_dir <- file.path(root_dir, "04methykit")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

coverage_summary <- do.call(rbind, lapply(seq_along(obj.cpg), function(i) {
  df <- getData(obj.cpg[[i]])
  data.frame(
    sample = sample_ids[[i]],
    group = reads_info$group[[i]],
    primary_mincov = primary_mincov,
    analyzed_cpg_sites_primary_mincov = nrow(df),
    analyzed_cpg_sites_mincov4 = if (primary_mincov == 4) nrow(df) else NA_integer_,
    ge_4_sites = sum(df$coverage >= 4, na.rm = TRUE),
    ge_4_pct_of_analyzed = round(mean(df$coverage >= 4, na.rm = TRUE) * 100, 2),
    ge_10_sites = sum(df$coverage >= 10, na.rm = TRUE),
    ge_10_pct_of_analyzed = round(mean(df$coverage >= 10, na.rm = TRUE) * 100, 2),
    ge_20_sites = sum(df$coverage >= 20, na.rm = TRUE),
    ge_20_pct_of_analyzed = round(mean(df$coverage >= 20, na.rm = TRUE) * 100, 2),
    mean_coverage = round(mean(df$coverage, na.rm = TRUE), 4),
    global_methylation_pct = round(sum(df$numCs, na.rm = TRUE) / sum(df$coverage, na.rm = TRUE) * 100, 4)
  )
}))

write.csv(
  coverage_summary,
  file = file.path(qc_dir, "cpg_coverage_threshold_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("成功保存至 04methykit/CpG.RData，并导出 coverage QC 汇总表\n")
EOF
  WGBS_ROOT_DIR="${ROOT_DIR}" \
  WGBS_METHREAD_MINCOV="${WGBS_METHREAD_MINCOV}" \
  WGBS_EXTRA_MINCOV_LEVELS="${WGBS_EXTRA_MINCOV_LEVELS}" \
  Rscript "${tmp_r_script}"
  rm -f "${tmp_r_script}"
fi

fi

if step_selected 13; then
log "Step 13: 导出 QC 汇总表"
WGBS_ROOT_DIR="${ROOT_DIR}" \
WGBS_SAMPLE_ORDER="$(IFS=';'; echo "${READS_INFO_SAMPLES[*]}")" \
python3 - <<'PY'
import csv
import json
import os
import re
from pathlib import Path


def maybe_pct(numerator, denominator):
    if numerator in (None, "") or denominator in (None, "", 0):
        return ""
    return round(float(numerator) / float(denominator) * 100, 2)


def maybe_pair_count(read_count):
    if read_count in (None, ""):
        return ""
    read_count = int(read_count)
    if read_count % 2 != 0:
        return ""
    return read_count // 2


def read_fastp_json(json_file):
    with json_file.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    raw_reads = before.get("total_reads")
    clean_reads = after.get("total_reads")

    return {
        "sample": json_file.name.replace("_fastp.json", ""),
        "raw_reads": raw_reads,
        "raw_read_pairs_est": maybe_pair_count(raw_reads),
        "clean_reads": clean_reads,
        "clean_read_pairs_est": maybe_pair_count(clean_reads),
        "clean_raw_pct": maybe_pct(clean_reads, raw_reads),
        "q30_before_pct": round(before.get("q30_rate", 0) * 100, 2) if before.get("q30_rate") is not None else "",
        "q30_after_pct": round(after.get("q30_rate", 0) * 100, 2) if after.get("q30_rate") is not None else "",
    }


def read_alignment_metrics(metrics_file):
    with metrics_file.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader)

    mapped_reads = int(row["mapped_primary_reads"]) if row.get("mapped_primary_reads") else None
    deduplicated_reads = int(row["deduplicated_primary_reads"]) if row.get("deduplicated_primary_reads") else None

    return {
        "sample": row["sample"],
        "mapped_primary_reads": mapped_reads,
        "mapped_primary_pairs_est": maybe_pair_count(mapped_reads),
        "deduplicated_primary_reads": deduplicated_reads,
        "deduplicated_primary_pairs_est": maybe_pair_count(deduplicated_reads),
    }


def sample_from_report_name(path, suffix_pattern):
    match = re.match(suffix_pattern, path.name)
    if match:
        return match.group("sample")
    return None


def read_bismark_alignment_report(report_file):
    sample = sample_from_report_name(
        report_file,
        r"(?P<sample>.+)_1\.clean_bismark_bt2_PE_report\.txt$",
    )
    record = {"sample": sample or report_file.stem}

    with report_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("Sequence pairs analysed in total:"):
                record["bismark_sequence_pairs_total"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Number of paired-end alignments with a unique best hit:"):
                record["bismark_unique_best_pairs"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Mapping efficiency:"):
                value = line.split(":", 1)[1].strip().replace("%", "")
                record["bismark_mapping_efficiency_pct"] = float(value)
            elif line.startswith("C methylated in CpG context:"):
                record["cpg_methylation_pct_bismark"] = float(line.split(":", 1)[1].strip().replace("%", ""))
            elif line.startswith("C methylated in CHG context:"):
                record["chg_methylation_pct_bismark"] = float(line.split(":", 1)[1].strip().replace("%", ""))
            elif line.startswith("C methylated in CHH context:"):
                record["chh_methylation_pct_bismark"] = float(line.split(":", 1)[1].strip().replace("%", ""))

    non_cg_values = [
        record.get("chg_methylation_pct_bismark"),
        record.get("chh_methylation_pct_bismark"),
    ]
    non_cg_values = [value for value in non_cg_values if value is not None]
    if non_cg_values:
        record["non_cg_methylation_max_pct_bismark"] = max(non_cg_values)
        record["bisulfite_conversion_lower_bound_pct"] = round(100 - max(non_cg_values), 2)

    return record


def read_bismark_dedup_report(report_file):
    sample = sample_from_report_name(
        report_file,
        r"(?P<sample>.+)_1\.clean_bismark_bt2_pe\.deduplication_report\.txt$",
    )
    record = {"sample": sample or report_file.stem}

    with report_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("Total number duplicated alignments removed:"):
                match = re.search(r"\(([-0-9.]+)%\)", line)
                if match:
                    record["bismark_duplicate_rate_pct"] = float(match.group(1))
            elif line.startswith("Total count of deduplicated leftover sequences:"):
                match = re.search(r"^Total count of deduplicated leftover sequences:\s+([0-9]+)\s+\(([-0-9.]+)%", line)
                if match:
                    record["bismark_dedup_leftover_reads"] = int(match.group(1))
                    record["bismark_dedup_retention_pct"] = float(match.group(2))

    return record


def write_csv(rows, output_file, fieldnames):
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


root_dir = Path(os.environ["WGBS_ROOT_DIR"]).resolve()
sample_order = [sample for sample in os.environ.get("WGBS_SAMPLE_ORDER", "").split(";") if sample]
qc_dir = root_dir / "01rawdata" / "qc_results_GRCm39"
out_dir = root_dir / "01rawdata" / "qc_results_GRCm39"
bismark_dir = root_dir / "02bismark_results_GRCm39"
mbias_dir = root_dir / "03methylation_results_GRCm39" / "bismark_mbias"

fastp_records = {}
for json_file in qc_dir.glob("*_fastp.json"):
    record = read_fastp_json(json_file)
    fastp_records[record["sample"]] = record

align_records = {}
for metrics_file in qc_dir.glob("*_alignment_metrics.tsv"):
    record = read_alignment_metrics(metrics_file)
    align_records[record["sample"]] = record

bismark_records = {}
for report_file in bismark_dir.glob("*_PE_report.txt"):
    record = read_bismark_alignment_report(report_file)
    bismark_records[record["sample"]] = record

dedup_report_records = {}
for report_file in bismark_dir.glob("*deduplication_report.txt"):
    record = read_bismark_dedup_report(report_file)
    dedup_report_records[record["sample"]] = record

mbias_report_samples = {}
for mbias_file in mbias_dir.glob("*M-bias.txt*"):
    sample = re.sub(r"(\.sorted|\.name_sorted_for_mbias|_1\.clean_bismark_bt2_pe\.deduplicated)?\.M-bias\.txt(\.gz)?$", "", mbias_file.name)
    mbias_report_samples[sample] = True

remaining = sorted(
    (set(fastp_records) | set(align_records) | set(bismark_records) | set(dedup_report_records)) -
    set(sample_order)
)
all_samples = sample_order + remaining

read_alignment_rows = []
for sample in all_samples:
    fastp_record = fastp_records.get(sample, {})
    align_record = align_records.get(sample, {})
    bismark_record = bismark_records.get(sample, {})
    dedup_report_record = dedup_report_records.get(sample, {})

    clean_reads = fastp_record.get("clean_reads")
    mapped_reads = align_record.get("mapped_primary_reads")
    dedup_reads = align_record.get("deduplicated_primary_reads")

    read_alignment_rows.append(
        {
            "sample": sample,
            "raw_reads": fastp_record.get("raw_reads", ""),
            "raw_read_pairs_est": fastp_record.get("raw_read_pairs_est", ""),
            "clean_reads": clean_reads if clean_reads is not None else "",
            "clean_read_pairs_est": fastp_record.get("clean_read_pairs_est", ""),
            "clean_raw_pct": fastp_record.get("clean_raw_pct", ""),
            "q30_before_pct": fastp_record.get("q30_before_pct", ""),
            "q30_after_pct": fastp_record.get("q30_after_pct", ""),
            "mapped_primary_reads": mapped_reads if mapped_reads is not None else "",
            "mapped_primary_pairs_est": align_record.get("mapped_primary_pairs_est", ""),
            "alignment_rate_pct": maybe_pct(mapped_reads, clean_reads),
            "bismark_sequence_pairs_total": bismark_record.get("bismark_sequence_pairs_total", ""),
            "bismark_unique_best_pairs": bismark_record.get("bismark_unique_best_pairs", ""),
            "bismark_mapping_efficiency_pct": bismark_record.get("bismark_mapping_efficiency_pct", ""),
            "deduplicated_primary_reads": dedup_reads if dedup_reads is not None else "",
            "deduplicated_primary_pairs_est": align_record.get("deduplicated_primary_pairs_est", ""),
            "dedup_retention_pct": maybe_pct(dedup_reads, mapped_reads),
            "bismark_duplicate_rate_pct": dedup_report_record.get("bismark_duplicate_rate_pct", ""),
            "bismark_dedup_retention_pct": dedup_report_record.get("bismark_dedup_retention_pct", ""),
            "cpg_methylation_pct_bismark": bismark_record.get("cpg_methylation_pct_bismark", ""),
            "chg_methylation_pct_bismark": bismark_record.get("chg_methylation_pct_bismark", ""),
            "chh_methylation_pct_bismark": bismark_record.get("chh_methylation_pct_bismark", ""),
            "non_cg_methylation_max_pct_bismark": bismark_record.get("non_cg_methylation_max_pct_bismark", ""),
            "bisulfite_conversion_lower_bound_pct": bismark_record.get("bisulfite_conversion_lower_bound_pct", ""),
            "mbias_report_available": "yes" if mbias_report_samples.get(sample) else "no",
        }
    )

write_csv(
    read_alignment_rows,
    out_dir / "qc_read_alignment_summary.csv",
    [
        "sample",
        "raw_reads",
        "raw_read_pairs_est",
        "clean_reads",
        "clean_read_pairs_est",
        "clean_raw_pct",
        "q30_before_pct",
        "q30_after_pct",
        "mapped_primary_reads",
        "mapped_primary_pairs_est",
        "alignment_rate_pct",
        "bismark_sequence_pairs_total",
        "bismark_unique_best_pairs",
        "bismark_mapping_efficiency_pct",
        "deduplicated_primary_reads",
        "deduplicated_primary_pairs_est",
        "dedup_retention_pct",
        "bismark_duplicate_rate_pct",
        "bismark_dedup_retention_pct",
        "cpg_methylation_pct_bismark",
        "chg_methylation_pct_bismark",
        "chh_methylation_pct_bismark",
        "non_cg_methylation_max_pct_bismark",
        "bisulfite_conversion_lower_bound_pct",
        "mbias_report_available",
    ],
)

print(f"QC read/alignment summary written to: {out_dir / 'qc_read_alignment_summary.csv'}")
PY
fi

if step_selected 14; then
# C/T 多态性在 WGBS 中会被误判为甲基化差异；本步骤只准备 BED，不直接改写 methylKit/DMR 输入。
if [[ -n "${WGBS_DBSNP_VCF}" && -f "${WGBS_DBSNP_VCF}" ]]; then
  log "Step 14: SNP BED 准备 (使用 ${WGBS_DBSNP_VCF})"
  mkdir -p "${ROOT_DIR}/04methykit/snp_filtered"
  if have_command bcftools; then
    # 将 VCF 转换为 BED 格式供 MethylDackel 使用
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${WGBS_DBSNP_VCF}" | \
      awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
      > "${ROOT_DIR}/04methykit/snp_filtered/snps.bed" || \
      log "  SNP BED 转换失败"
    log "  SNP BED 已保存至 04methykit/snp_filtered/snps.bed"
    log "  注意: 当前步骤未实际过滤 methylKit/DMR 输入；后续可在 MethylDackel 或 R 层面使用此 BED"
  else
    log "Step 14: 未找到 bcftools，跳过 SNP BED 准备"
  fi
else
  log "Step 14: 未设置 WGBS_DBSNP_VCF 环境变量或文件不存在，跳过 SNP BED 准备"
  log "  提示: 设置 WGBS_DBSNP_VCF=/path/to/dbsnp.vcf.gz 后运行可准备 SNP BED"
fi
fi

if step_selected 15; then
if have_command multiqc; then
  log "Step 15: 运行 MultiQC 汇总 QC 报告"
  mkdir -p "${ROOT_DIR}/01rawdata/qc_results_GRCm39/multiqc"
  multiqc \
    --force \
    --outdir "${ROOT_DIR}/01rawdata/qc_results_GRCm39/multiqc" \
    "${ROOT_DIR}/01rawdata/qc_results_GRCm39" \
    "${ROOT_DIR}/02bismark_results_GRCm39" \
    "${ROOT_DIR}/03methylation_results_GRCm39" || \
    log "  MultiQC 汇总失败，保留已有单工具 QC 输出"
else
  log "Step 15: 未找到 multiqc，跳过 MultiQC 汇总"
fi
fi

log "WGBS 上游流程运行完成"
