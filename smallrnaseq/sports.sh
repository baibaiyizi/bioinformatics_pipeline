#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE=run
while (($#)); do
  case "$1" in
    --help|-h) echo 'Usage: sports.sh [--config project.env] [--check|--dry-run]'; exit 0 ;;
    --config) export SMALLRNA_CONFIG="$2"; shift 2 ;;
    --check) MODE=check; shift ;;
    --dry-run) MODE=dry-run; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
source "$ROOT_DIR/scripts/load_config.sh"
RAW_DIR="${ROOT_DIR}"
RESULT_DIR="${SPORTS_RESULT_DIR}"
READS_INFO_FILE="${SAMPLES_TSV}"

METADATA_DIR="${RESULT_DIR}/02_metadata"
QC_DIR="${RESULT_DIR}/03_qc"
SPORTS_ANNOTATION_DIR="${RESULT_DIR}/04_sports_annotation"
SAMPLE_ROOT="${SPORTS_ANNOTATION_DIR}/samples"
ENDOSIRNA_DIR="${RESULT_DIR}/05_endosirna"
ENDOSIRNA_SAMPLE_ROOT="${ENDOSIRNA_DIR}/samples"
ENDOSIRNA_SUMMARY_DIR="${ENDOSIRNA_DIR}/summary"
LOG_DIR="${RESULT_DIR}/06_logs"
SAMPLE_LOG_DIR="${LOG_DIR}/samples"
CUTADAPT_LOG_DIR="${LOG_DIR}/cutadapt"
ENDOSIRNA_LOG_DIR="${LOG_DIR}/endosirna"
WORK_DIR="${RESULT_DIR}/99_work"
TMP_DIR="${TMP_DIR:-${WORK_DIR}/fastq}"
METADATA_TSV="${METADATA_DIR}/sample_metadata.tsv"
NORMALIZED_READS_INFO_TSV="${METADATA_DIR}/reads_info.normalized.tsv"
OUTPUT_LAYOUT_TSV="${METADATA_DIR}/sports_output_layout.tsv"

SPORTS_PL="${SPORTS_PL:?explicit config required}"
REFERENCE_ROOT="${REFERENCE_ROOT:?explicit config required}"
# PANDORA-local small RNA reference root. Override with DB_ROOT=/path/to/Mus_musculus if needed.
DB_ROOT="${DB_ROOT:?explicit config required}"
RMD_RUNNER="${RMD_RUNNER:-${ROOT_DIR}/run_rmd.sh}"

THREADS="${THREADS:-8}"
PARALLEL_JOBS="${PARALLEL_JOBS:-3}"
# 索引已建立时直接并行；索引未建立时，先串行预热若干样本，再切到并行。
INDEX_AWARE_SCHEDULING="${INDEX_AWARE_SCHEDULING:-1}"
INDEX_WARMUP_SAMPLES="${INDEX_WARMUP_SAMPLES:-1}"
# Length range retained after cutadapt and accepted by SPORTS.
MIN_LEN="${MIN_LEN:-15}"
MAX_LEN="${MAX_LEN:-45}"
# Bowtie mismatch tolerance used by SPORTS annotation.
MISMATCH="${MISMATCH:-1}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-1}"
KEEP_FASTQ="${KEEP_FASTQ:-0}"
# Recommended mode: cutadapt pre-trims FASTQ, then SPORTS annotates without
# repeating adapter trimming.
PRETRIM_FASTQ="${PRETRIM_FASTQ:-1}"
TRIM_ADAPTER="${TRIM_ADAPTER:-0}"
# VAHTS Small RNA Library Prep Kit for Illumina V2 read-through 3' adapter.
ADAPTER3="${ADAPTER3:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
# Two-color sequencers often create artificial terminal G runs in raw reads.
CUTADAPT_NEXTSEQ_TRIM="${CUTADAPT_NEXTSEQ_TRIM:-20}"
# Optional ordinary quality trimming, for example CUTADAPT_QUALITY_CUTOFF=20.
CUTADAPT_QUALITY_CUTOFF="${CUTADAPT_QUALITY_CUTOFF:-}"
RUN_ANTISENSE="${RUN_ANTISENSE:-0}"
RUN_DOWNSTREAM="${RUN_DOWNSTREAM:-0}"
RUN_FASTQC="${RUN_FASTQC:-auto}"
RUN_ENDOSIRNA="${RUN_ENDOSIRNA:-0}"
SAMPLES_FILE="${SAMPLES_FILE:-}"
SKIP_EXISTING=1
FASTQC_DIR="${QC_DIR}/fastqc"

case "${RUN_ENDOSIRNA}" in
  true|TRUE|yes|YES) RUN_ENDOSIRNA=1 ;;
  false|FALSE|no|NO) RUN_ENDOSIRNA=0 ;;
esac

GENOME_FASTA="${GENOME_FASTA:?explicit config required}"
GENOME_BASE="${GENOME_BASE:?explicit config required}"
MIRNA_BASE="${MIRNA_BASE:?explicit config required}"
RRNA_BASE="${RRNA_BASE:?explicit config required}"
TRNA_BASE="${TRNA_BASE:?explicit config required}"
TRNA_MITO_BASE="${TRNA_MITO_BASE:?explicit config required}"
PIRNA_BASE="${PIRNA_BASE:?explicit config required}"
ENSEMBL_BASE="${ENSEMBL_BASE:?explicit config required}"
RFAM_BASE="${RFAM_BASE:?explicit config required}"
ENDOSIRNA_REFS="${ENDOSIRNA_REFS:-}"
ENDOSIRNA_MISMATCH="${ENDOSIRNA_MISMATCH:-0}"
ENDOSIRNA_MAX_ALIGNMENTS="${ENDOSIRNA_MAX_ALIGNMENTS:-20}"
ENDOSIRNA_SUMMARY_TSV="${ENDOSIRNA_SUMMARY_DIR}/endo_siRNA_summary.tsv"

declare -a RRNA_SUBTYPES=(
  "4.5S"
  "5S"
  "5.8S"
  "12S"
  "16S"
  "18S"
  "28S"
  "45S"
  "RNY1"
  "RNY3"
)

declare -a DISCOVERED_SAMPLES=()
declare -a DISCOVERED_FASTQS=()
declare -a DISCOVERED_GROUPS=()
declare -a DISCOVERED_REPLICATES=()
declare -a DISCOVERED_CONDITIONS=()
declare -a DISCOVERED_AGE_GROUPS=()
declare -a DISCOVERED_SPECIES=()
declare -a DISCOVERED_STAGES=()
declare -A REQUESTED_SAMPLE_SET=()
declare -a BOWTIE_INDEX_PREFIXES=(
  "${GENOME_BASE}"
  "${MIRNA_BASE}"
  "${TRNA_BASE}"
  "${TRNA_BASE}_CCA"
  "${TRNA_MITO_BASE}"
  "${TRNA_MITO_BASE}_CCA"
  "${PIRNA_BASE}"
  "${ENSEMBL_BASE}"
  "${RFAM_BASE}"
)

for rrna_subtype in "${RRNA_SUBTYPES[@]}"; do
  BOWTIE_INDEX_PREFIXES+=("${RRNA_BASE}_${rrna_subtype}")
done

if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
  IFS=';' read -r -a _endosirna_ref_specs <<< "${ENDOSIRNA_REFS}"
  for _endosirna_spec in "${_endosirna_ref_specs[@]}"; do
    _endosirna_spec="${_endosirna_spec#"${_endosirna_spec%%[![:space:]]*}"}"
    _endosirna_spec="${_endosirna_spec%"${_endosirna_spec##*[![:space:]]}"}"
    [[ -z "${_endosirna_spec}" ]] && continue
    if [[ "${_endosirna_spec}" == *"="* ]]; then
      BOWTIE_INDEX_PREFIXES+=("${_endosirna_spec#*=}")
    else
      BOWTIE_INDEX_PREFIXES+=("${_endosirna_spec}")
    fi
  done
  unset _endosirna_ref_specs _endosirna_spec
fi

log_info() {
  echo "[INFO] $*"
}

log_warn() {
  echo "[WARN] $*" >&2
}

log_error() {
  echo "[ERROR] $*" >&2
}

sample_output_dir() {
  printf '%s\n' "${SAMPLE_ROOT}/$1"
}

endosirna_sample_output_dir() {
  printf '%s\n' "${ENDOSIRNA_SAMPLE_ROOT}/$1"
}

write_output_layout_manifest() {
  mkdir -p "${METADATA_DIR}"
  {
    printf 'section\tpath\tcontents\n'
    printf 'metadata\t%s\tsample metadata, normalized reads_info, output layout manifest\n' "${METADATA_DIR}"
    printf 'qc\t%s\traw and trimmed FastQC reports\n' "${FASTQC_DIR}"
    printf 'sports_annotation\t%s\tprimary SPORTS annotation outputs, one directory per sample\n' "${SAMPLE_ROOT}"
    printf 'endosirna\t%s\toptional RUN_ENDOSIRNA Bowtie evidence layer\n' "${ENDOSIRNA_DIR}"
    printf 'logs\t%s\tpipeline logs grouped by stage\n' "${LOG_DIR}"
    printf 'work\t%s\ttemporary FASTQ files used during processing\n' "${WORK_DIR}"
  } > "${OUTPUT_LAYOUT_TSV}"
}

trim_shell_string() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s\n' "${value}"
}

sanitize_endosirna_tag() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
}

endosirna_ref_specs() {
  local -a specs=()
  local spec tag prefix

  IFS=';' read -r -a specs <<< "${ENDOSIRNA_REFS}"
  for spec in "${specs[@]}"; do
    spec="$(trim_shell_string "${spec}")"
    [[ -z "${spec}" ]] && continue

    if [[ "${spec}" == *"="* ]]; then
      tag="$(trim_shell_string "${spec%%=*}")"
      prefix="$(trim_shell_string "${spec#*=}")"
    else
      prefix="$(trim_shell_string "${spec}")"
      tag="$(basename "${prefix}")"
    fi

    tag="$(sanitize_endosirna_tag "${tag}")"
    [[ -z "${tag}" || -z "${prefix}" ]] && continue
    printf '%s\t%s\n' "${tag}" "${prefix}"
  done
}

reference_fasta_exists() {
  local prefix="$1"
  [[ -f "${prefix}.fa" || -f "${prefix}.fasta" ]]
}

check_endosirna_refs() {
  local ref_count=0
  local ref_tag ref_prefix

  if ! [[ "${ENDOSIRNA_MISMATCH}" =~ ^[0-9]+$ ]]; then
    log_error "ENDOSIRNA_MISMATCH 必须是非负整数，当前为: ${ENDOSIRNA_MISMATCH}"
    exit 1
  fi
  if ! [[ "${ENDOSIRNA_MAX_ALIGNMENTS}" =~ ^[1-9][0-9]*$ ]]; then
    log_error "ENDOSIRNA_MAX_ALIGNMENTS 必须是正整数，当前为: ${ENDOSIRNA_MAX_ALIGNMENTS}"
    exit 1
  fi

  while IFS=$'\t' read -r ref_tag ref_prefix; do
    ref_count=$((ref_count + 1))
    if ! reference_fasta_exists "${ref_prefix}"; then
      log_error "缺少 endo-siRNA FASTA: ${ref_prefix}.fa 或 ${ref_prefix}.fasta"
      exit 1
    fi
    if ! bowtie_index_exists "${ref_prefix}"; then
      log_error "缺少 endo-siRNA Bowtie 索引: ${ref_prefix}.[1-4,rev.1,rev.2].ebwt"
      exit 1
    fi
  done < <(endosirna_ref_specs)

  if [[ "${ref_count}" -eq 0 ]]; then
    log_error "RUN_ENDOSIRNA=1 但 ENDOSIRNA_REFS 为空"
    exit 1
  fi
}

initialize_outputs() {
  mkdir -p \
    "${RESULT_DIR}" \
    "${METADATA_DIR}" \
    "${TMP_DIR}" \
    "${SAMPLE_ROOT}" \
    "${ENDOSIRNA_SAMPLE_ROOT}" \
    "${ENDOSIRNA_SUMMARY_DIR}" \
    "${LOG_DIR}" \
    "${SAMPLE_LOG_DIR}" \
    "${CUTADAPT_LOG_DIR}" \
    "${ENDOSIRNA_LOG_DIR}" \
    "${FASTQC_DIR}"
  write_output_layout_manifest

}

check_prereqs() {

  if [[ ! -x "${SPORTS_PL}" ]]; then
    log_error "sports.pl 不存在或不可执行: ${SPORTS_PL}"
    exit 1
  fi

  local executables=(perl bowtie bowtie-build gzip xargs)
  local exe
  for exe in "${executables[@]}"; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
      log_error "未找到依赖程序: ${exe}"
      exit 1
    fi
  done

  if [[ "${PRETRIM_FASTQ}" == "1" || "${TRIM_ADAPTER}" == "1" ]] && ! command -v cutadapt >/dev/null 2>&1; then
    log_error "PRETRIM_FASTQ=1 或 TRIM_ADAPTER=1 时需要 cutadapt"
    exit 1
  fi

  if [[ "${RUN_DOWNSTREAM}" == "1" ]] && ! command -v Rscript >/dev/null 2>&1; then
    log_error "RUN_DOWNSTREAM=1 时需要 Rscript"
    exit 1
  fi

  if [[ "${RUN_FASTQC}" == "1" ]] && ! command -v fastqc >/dev/null 2>&1; then
    log_error "RUN_FASTQC=1 时需要 fastqc"
    exit 1
  fi
}

check_databases() {
  local required_fa=(
    "${MIRNA_BASE}.fa"
    "${TRNA_BASE}.fa"
    "${TRNA_MITO_BASE}.fa"
    "${PIRNA_BASE}.fa"
    "${ENSEMBL_BASE}.fa"
    "${RFAM_BASE}.fa"
  )

  if [[ ! -s "${GENOME_FASTA}" ]]; then
    log_error "配置的基因组 FASTA 不存在或为空: ${GENOME_FASTA}"
    exit 1
  fi

  if ! bowtie_index_exists "${GENOME_BASE}"; then
    log_error "配置的 Bowtie1 索引不完整: ${GENOME_BASE}"
    exit 1
  fi

  local f
  for f in "${required_fa[@]}"; do
    if [[ ! -f "${f}" ]]; then
      log_error "缺少数据库文件: ${f}"
      exit 1
    fi
  done

  local rrna_subtype
  for rrna_subtype in "${RRNA_SUBTYPES[@]}"; do
    if [[ ! -f "${RRNA_BASE}_${rrna_subtype}.fa" ]]; then
      log_error "缺少 configured rRNA 数据库文件: ${RRNA_BASE}_${rrna_subtype}.fa"
      exit 1
    fi
  done

  if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
    check_endosirna_refs
  fi

  log_info "中央全基因组参考: ${REFERENCE_ROOT}"
  log_info "PANDORA 小 RNA 专用库: ${DB_ROOT}"
  log_info "PRETRIM_FASTQ=${PRETRIM_FASTQ}, TRIM_ADAPTER=${TRIM_ADAPTER}, RUN_ANTISENSE=${RUN_ANTISENSE}, RUN_ENDOSIRNA=${RUN_ENDOSIRNA}"
  log_info "ADAPTER3=${ADAPTER3}, CUTADAPT_NEXTSEQ_TRIM=${CUTADAPT_NEXTSEQ_TRIM}, CUTADAPT_QUALITY_CUTOFF=${CUTADAPT_QUALITY_CUTOFF:-<none>}"
  if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
    log_info "endo-siRNA 独立证据层已启用: ENDOSIRNA_MISMATCH=${ENDOSIRNA_MISMATCH}, ENDOSIRNA_MAX_ALIGNMENTS=${ENDOSIRNA_MAX_ALIGNMENTS}"
    log_info "endo-siRNA 参考库: ${ENDOSIRNA_REFS}"
    if [[ "${PRETRIM_FASTQ}" != "1" ]]; then
      log_warn "RUN_ENDOSIRNA=1 时建议 PRETRIM_FASTQ=1；当前 endo-siRNA 将使用未预剪切的 SPORTS 输入 FASTQ。"
    fi
  else
    log_info "当前流程使用配置指定的参考库；endo-siRNA 独立证据层未启用，仅使用 genome/miRNA/rRNA/tRNA/piRNA/Ensembl/Rfam。"
  fi
}

bowtie_index_exists() {
  local prefix="$1"

  [[ -f "${prefix}.1.ebwt" && -f "${prefix}.2.ebwt" && -f "${prefix}.3.ebwt" && -f "${prefix}.4.ebwt" && -f "${prefix}.rev.1.ebwt" && -f "${prefix}.rev.2.ebwt" ]]
}

missing_bowtie_indexes() {
  local prefix

  for prefix in "${BOWTIE_INDEX_PREFIXES[@]}"; do
    if ! bowtie_index_exists "${prefix}"; then
      printf '%s\n' "${prefix}"
    fi
  done
}

all_bowtie_indexes_ready() {
  local missing=()

  mapfile -t missing < <(missing_bowtie_indexes)
  [[ "${#missing[@]}" -eq 0 ]]
}

load_requested_samples() {
  local entry sample

  REQUESTED_SAMPLE_SET=()
  [[ -z "${SAMPLES_FILE}" ]] && return 0

  if [[ ! -f "${SAMPLES_FILE}" ]]; then
    log_error "SAMPLES_FILE 指定的文件不存在: ${SAMPLES_FILE}"
    exit 1
  fi

  while IFS= read -r entry || [[ -n "${entry}" ]]; do
    entry="${entry%%#*}"
    entry="$(printf '%s' "${entry}" | xargs)"
    [[ -z "${entry}" ]] && continue

    if [[ "${entry}" == *.fq.gz ]]; then
      sample="$(basename "${entry}")"
      sample="${sample%.fq.gz}"
    elif [[ "${entry}" == *.fastq.gz ]]; then
      sample="$(basename "${entry}")"
      sample="${sample%.fastq.gz}"
    else
      sample="${entry}"
    fi

    REQUESTED_SAMPLE_SET["${sample}"]=1
  done < "${SAMPLES_FILE}"
}

sample_requested() {
  local sample="$1"

  if [[ "${#REQUESTED_SAMPLE_SET[@]}" -eq 0 ]]; then
    return 0
  fi

  [[ -n "${REQUESTED_SAMPLE_SET[${sample}]+x}" ]]
}

normalize_fastq_path() {
  if [[ "$1" == /* ]]; then printf '%s\n' "$1"; else printf '%s/%s\n' "$ROOT_DIR" "$1"; fi
}

append_reads_info_record() {
  local group="$1"
  local sample="$2"
  local replicate="$3"
  local fq_gz="$4"
  local condition="${5:-}"
  local age_group="${6:-}"
  local species="${7:-}"
  local stage="${8:-}"

  if ! sample_requested "${sample}"; then
    return 0
  fi

  if [[ -z "${group}" || -z "${sample}" || -z "${replicate}" || -z "${fq_gz}" ]]; then
    log_error "reads_info.txt 中存在缺少必要字段的记录: group=${group:-<empty>}, sample=${sample:-<empty>}"
    exit 1
  fi

  if [[ ! -f "${fq_gz}" ]]; then
    log_error "reads_info.txt 中指定的 fastq 不存在: ${fq_gz}"
    exit 1
  fi

  DISCOVERED_GROUPS+=("${group}")
  DISCOVERED_SAMPLES+=("${sample}")
  DISCOVERED_REPLICATES+=("${replicate}")
  DISCOVERED_FASTQS+=("${fq_gz}")
  DISCOVERED_CONDITIONS+=("${condition:-${group}}")
  DISCOVERED_AGE_GROUPS+=("${age_group}")
  DISCOVERED_SPECIES+=("${species}")
  DISCOVERED_STAGES+=("${stage}")
}

discover_samples_from_reads_info() {
  local group sample replicate raw_fq condition age_group species stage fq_gz
  local records
  records="$(python3 "$ROOT_DIR/scripts/preflight.py" --mode raw --emit)" || return 1
  printf 'group\tsample\treplicate\traw_fq\tcondition\tage_group\tspecies\tstage\n' > "$NORMALIZED_READS_INFO_TSV"
  while IFS=$'\x1f' read -r group sample replicate raw_fq condition age_group species stage; do
    fq_gz="$(normalize_fastq_path "$raw_fq")"
    append_reads_info_record "$group" "$sample" "$replicate" "$fq_gz" "$condition" "$age_group" "$species" "$stage"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$group" "$sample" "$replicate" "$fq_gz" "$condition" "$age_group" "$species" "$stage" >> "$NORMALIZED_READS_INFO_TSV"
  done <<< "$records"
}

discover_samples() {
  DISCOVERED_SAMPLES=()
  DISCOVERED_FASTQS=()
  DISCOVERED_GROUPS=()
  DISCOVERED_REPLICATES=()
  DISCOVERED_CONDITIONS=()
  DISCOVERED_AGE_GROUPS=()
  DISCOVERED_SPECIES=()
  DISCOVERED_STAGES=()

  load_requested_samples
  discover_samples_from_reads_info

  if [[ "${#DISCOVERED_SAMPLES[@]}" -eq 0 ]]; then
    log_error "未从 reads_info.txt 识别到可运行样本，请检查 ${READS_INFO_FILE}"
    exit 1
  fi

  log_info "共识别到 ${#DISCOVERED_SAMPLES[@]} 个样本"
}

write_metadata() {
  local idx
  local sample
  local group
  local replicate
  local fq_gz

  printf 'sample\tgroup\tcondition\tage_group\tspecies\tstage\treplicate\traw_fastq_gz\n' > "${METADATA_TSV}"
  for idx in "${!DISCOVERED_SAMPLES[@]}"; do
    sample="${DISCOVERED_SAMPLES[${idx}]}"
    group="${DISCOVERED_GROUPS[${idx}]}"
    replicate="${DISCOVERED_REPLICATES[${idx}]}"
    fq_gz="${DISCOVERED_FASTQS[${idx}]}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample}" \
      "${group}" \
      "${DISCOVERED_CONDITIONS[${idx}]}" \
      "${DISCOVERED_AGE_GROUPS[${idx}]}" \
      "${DISCOVERED_SPECIES[${idx}]}" \
      "${DISCOVERED_STAGES[${idx}]}" \
      "${replicate}" \
      "${fq_gz}" >> "${METADATA_TSV}"
  done

  log_info "样本 metadata 已写出: ${METADATA_TSV}"
}

should_run_fastqc() {
  case "${RUN_FASTQC}" in
    1|true|TRUE|yes|YES)
      return 0
      ;;
    0|false|FALSE|no|NO)
      return 1
      ;;
    auto|AUTO)
      command -v fastqc >/dev/null 2>&1
      return $?
      ;;
    *)
      log_warn "未知 RUN_FASTQC=${RUN_FASTQC}，按 auto 处理"
      command -v fastqc >/dev/null 2>&1
      return $?
      ;;
  esac
}

run_raw_fastqc_if_requested() {
  local status_file="${FASTQC_DIR}/fastqc_status.tsv"
  local -a missing_reports=()
  local idx sample fastq report_zip

  mkdir -p "${FASTQC_DIR}/raw"

  if ! should_run_fastqc; then
    printf 'stage\tstatus\treason\nraw\tskipped\tfastqc not requested or not available\n' > "${status_file}"
    log_info "RUN_FASTQC=${RUN_FASTQC}; raw FastQC 跳过"
    return 0
  fi

  for idx in "${!DISCOVERED_SAMPLES[@]}"; do
    sample="${DISCOVERED_SAMPLES[${idx}]}"
    report_zip="${FASTQC_DIR}/raw/${sample}_R1_fastqc.zip"
    [[ -s "${report_zip}" ]] || missing_reports+=("${DISCOVERED_FASTQS[${idx}]}")
  done

  if [[ "${#missing_reports[@]}" -eq 0 ]]; then
    printf 'stage\tstatus\treason\nraw\tcompleted\texisting fastqc reports reused\n' > "${status_file}"
    log_info "raw FastQC 结果已存在，跳过重新运行"
    return 0
  fi

  log_info "开始 raw FastQC: ${#missing_reports[@]} 个 fastq"
  fastqc -t "${THREADS}" -o "${FASTQC_DIR}/raw" "${missing_reports[@]}"
  printf 'stage\tstatus\treason\nraw\tcompleted\tfastqc finished\n' > "${status_file}"

  if command -v multiqc >/dev/null 2>&1; then
    multiqc -o "${FASTQC_DIR}/raw" "${FASTQC_DIR}/raw" >/dev/null 2>&1 || log_warn "MultiQC 汇总 raw FastQC 失败"
  fi
}

run_trimmed_fastqc_for_sample() {
  local sample="$1"
  local fastq="$2"
  local report_zip="${FASTQC_DIR}/trimmed/${sample}_fastqc.zip"

  if ! should_run_fastqc; then
    return 0
  fi
  if [[ ! -s "${fastq}" ]]; then
    log_warn "[${sample}] trimmed FastQC 跳过，文件不存在或为空: ${fastq}"
    return 0
  fi
  if [[ -s "${report_zip}" ]]; then
    return 0
  fi

  mkdir -p "${FASTQC_DIR}/trimmed"
  fastqc -t 1 -o "${FASTQC_DIR}/trimmed" "${fastq}" >/dev/null 2>&1 || log_warn "[${sample}] trimmed FastQC 失败"
}

run_trimmed_multiqc_if_requested() {
  local status_file="${FASTQC_DIR}/fastqc_status.tsv"

  if [[ "${PRETRIM_FASTQ}" != "1" ]]; then
    return 0
  fi
  if ! should_run_fastqc; then
    printf 'stage\tstatus\treason\nraw\tskipped\tfastqc not requested or not available\ntrimmed\tskipped\tfastqc not requested or not available\n' > "${status_file}"
    return 0
  fi

  if ! compgen -G "${FASTQC_DIR}/trimmed/*_fastqc.zip" >/dev/null; then
    {
      printf 'stage\tstatus\treason\n'
      printf 'raw\tcompleted\traw FastQC requested; see raw directory\n'
      printf 'trimmed\tskipped\tno trimmed FastQC reports found; rerun upstream samples with SKIP_EXISTING=0 to regenerate\n'
    } > "${status_file}"
    return 0
  fi

  if command -v multiqc >/dev/null 2>&1 && [[ -d "${FASTQC_DIR}/trimmed" ]]; then
    multiqc -o "${FASTQC_DIR}/trimmed" "${FASTQC_DIR}/trimmed" >/dev/null 2>&1 || log_warn "MultiQC 汇总 trimmed FastQC 失败"
  fi

  {
    printf 'stage\tstatus\treason\n'
    printf 'raw\tcompleted\traw FastQC requested; see raw directory\n'
    printf 'trimmed\tcompleted\ttrimmed FastQC finished or existing reports reused\n'
  } > "${status_file}"
}

is_sample_complete() {
  local sample="$1"
  local sample_out
  sample_out="$(sample_output_dir "${sample}")"
  local summary_file="${sample_out}/${sample}_summary.txt"
  local length_file="${sample_out}/${sample}_length_distribution.txt"
  local output_file="${sample_out}/${sample}_output.txt"

  [[ -s "${summary_file}" && -s "${length_file}" && -s "${output_file}" ]]
}

is_endosirna_complete() {
  local sample="$1"
  local ref_count=0
  local ref_tag ref_prefix stats_file sample_endosirna_out

  [[ "${RUN_ENDOSIRNA}" != "1" ]] && return 0
  sample_endosirna_out="$(endosirna_sample_output_dir "${sample}")"

  while IFS=$'\t' read -r ref_tag ref_prefix; do
    ref_count=$((ref_count + 1))
    stats_file="${sample_endosirna_out}/${sample}.${ref_tag}.stats.tsv"
    [[ -s "${stats_file}" ]] || return 1
  done < <(endosirna_ref_specs)

  [[ "${ref_count}" -gt 0 ]]
}

sample_needs_processing() {
  local sample="$1"

  if [[ "${SKIP_EXISTING}" == "1" ]] && is_sample_complete "${sample}"; then
    if [[ "${RUN_ENDOSIRNA}" != "1" ]] || is_endosirna_complete "${sample}"; then
      return 1
    fi
  fi

  return 0
}

write_endosirna_alignment_stats() {
  local sample="$1"
  local ref_tag="$2"
  local ref_prefix="$3"
  local total_reads="$4"
  local alignment_tsv="$5"
  local bowtie_log="$6"
  local stats_tsv="$7"

  {
    printf 'sample\tref_tag\tref_prefix\ttotal_reads\tmapped_reads\tmapped_fraction\ttotal_alignments\tunique_mapped_sequences\tunique_reference_ids\tmodal_length\tmodal_length_alignments\tendosirna_mismatch\tmax_alignments\tbowtie_output\tbowtie_log\n'
    awk -F $'\t' -v OFS=$'\t' \
      -v sample="${sample}" \
      -v ref_tag="${ref_tag}" \
      -v ref_prefix="${ref_prefix}" \
      -v total_reads="${total_reads}" \
      -v mismatch="${ENDOSIRNA_MISMATCH}" \
      -v max_alignments="${ENDOSIRNA_MAX_ALIGNMENTS}" \
      -v bowtie_output="${alignment_tsv}" \
      -v bowtie_log="${bowtie_log}" '
      {
        total_alignments++
        read_id = $1
        ref_id = $3
        seq = $5
        if (read_id != "" && !(read_id in read_seen)) {
          read_seen[read_id] = 1
          mapped_reads++
        }
        if (seq != "" && !(seq in seq_seen)) {
          seq_seen[seq] = 1
          unique_sequences++
        }
        if (ref_id != "" && !(ref_id in ref_seen)) {
          ref_seen[ref_id] = 1
          unique_reference_ids++
        }
        len = length(seq)
        if (len > 0) {
          length_alignments[len]++
          if (length_alignments[len] > modal_length_alignments) {
            modal_length_alignments = length_alignments[len]
            modal_length = len
          }
        }
      }
      END {
        if (total_reads + 0 > 0) {
          mapped_fraction = mapped_reads / total_reads
        } else {
          mapped_fraction = 0
        }
        if (modal_length == "") {
          modal_length = "NA"
          modal_length_alignments = 0
        }
        print sample, ref_tag, ref_prefix, total_reads + 0, mapped_reads + 0, sprintf("%.8f", mapped_fraction),
          total_alignments + 0, unique_sequences + 0, unique_reference_ids + 0,
          modal_length, modal_length_alignments + 0, mismatch, max_alignments,
          bowtie_output, bowtie_log
      }
    ' "${alignment_tsv}"
  } > "${stats_tsv}"
}

write_endosirna_reference_counts() {
  local sample="$1"
  local ref_tag="$2"
  local alignment_tsv="$3"
  local counts_tsv="$4"

  {
    printf 'sample\tref_tag\treference_id\tmapped_reads\ttotal_alignments\n'
    awk -F $'\t' -v OFS=$'\t' -v sample="${sample}" -v ref_tag="${ref_tag}" '
      {
        read_id = $1
        reference_id = $3
        if (reference_id == "") {
          next
        }
        total_alignments[reference_id]++
        read_ref_key = read_id SUBSEP reference_id
        if (read_id != "" && !(read_ref_key in read_ref_seen)) {
          read_ref_seen[read_ref_key] = 1
          mapped_reads[reference_id]++
        }
      }
      END {
        for (reference_id in total_alignments) {
          print sample, ref_tag, reference_id, mapped_reads[reference_id] + 0, total_alignments[reference_id] + 0
        }
      }
    ' "${alignment_tsv}"
  } > "${counts_tsv}"
}

write_endosirna_sequence_counts() {
  local sample="$1"
  local ref_tag="$2"
  local alignment_tsv="$3"
  local counts_tsv="$4"

  {
    printf 'sample\tref_tag\tsequence\tlength\tmapped_reads\ttotal_alignments\n'
    awk -F $'\t' -v OFS=$'\t' -v sample="${sample}" -v ref_tag="${ref_tag}" '
      {
        read_id = $1
        seq = $5
        if (seq == "") {
          next
        }
        total_alignments[seq]++
        read_seq_key = read_id SUBSEP seq
        if (read_id != "" && !(read_seq_key in read_seq_seen)) {
          read_seq_seen[read_seq_key] = 1
          mapped_reads[seq]++
        }
      }
      END {
        for (seq in total_alignments) {
          print sample, ref_tag, seq, length(seq), mapped_reads[seq] + 0, total_alignments[seq] + 0
        }
      }
    ' "${alignment_tsv}"
  } > "${counts_tsv}"
}

run_endosirna_for_sample() {
  local sample="$1"
  local fastq="$2"
  local out_dir
  local total_reads
  local ref_tag ref_prefix alignment_tsv bowtie_log stats_tsv ref_counts_tsv seq_counts_tsv

  [[ "${RUN_ENDOSIRNA}" != "1" ]] && return 0
  out_dir="$(endosirna_sample_output_dir "${sample}")"

  if [[ ! -s "${fastq}" ]]; then
    log_error "[${sample}] endo-siRNA 输入 FASTQ 不存在或为空: ${fastq}"
    return 1
  fi

  mkdir -p "${out_dir}"
  total_reads="$(awk 'END {printf "%d", NR / 4}' "${fastq}")"

  while IFS=$'\t' read -r ref_tag ref_prefix; do
    alignment_tsv="${out_dir}/${sample}.${ref_tag}.bowtie.tsv"
    bowtie_log="${ENDOSIRNA_LOG_DIR}/${sample}.${ref_tag}.log"
    stats_tsv="${out_dir}/${sample}.${ref_tag}.stats.tsv"
    ref_counts_tsv="${out_dir}/${sample}.${ref_tag}.reference_counts.tsv"
    seq_counts_tsv="${out_dir}/${sample}.${ref_tag}.sequence_counts.tsv"

    if [[ "${SKIP_EXISTING}" == "1" && -s "${stats_tsv}" ]]; then
      log_info "[${sample}] endo-siRNA ${ref_tag} 统计已存在，跳过"
      continue
    fi

    if [[ "${total_reads}" -eq 0 ]]; then
      : > "${alignment_tsv}"
      printf '[INFO] No reads available after trimming for %s %s\n' "${sample}" "${ref_tag}" > "${bowtie_log}"
    else
      log_info "[${sample}] endo-siRNA 比对 ${ref_tag}"
      bowtie \
        -q \
        -v "${ENDOSIRNA_MISMATCH}" \
        -k "${ENDOSIRNA_MAX_ALIGNMENTS}" \
        --best \
        --strata \
        -p "${THREADS}" \
        "${ref_prefix}" \
        "${fastq}" \
        > "${alignment_tsv}" \
        2> "${bowtie_log}"
    fi

    write_endosirna_alignment_stats "${sample}" "${ref_tag}" "${ref_prefix}" "${total_reads}" "${alignment_tsv}" "${bowtie_log}" "${stats_tsv}"
    write_endosirna_reference_counts "${sample}" "${ref_tag}" "${alignment_tsv}" "${ref_counts_tsv}"
    write_endosirna_sequence_counts "${sample}" "${ref_tag}" "${alignment_tsv}" "${seq_counts_tsv}"
  done < <(endosirna_ref_specs)
}

summarize_endosirna_results() {
  local sample stats_file
  local wrote_header=0

  [[ "${RUN_ENDOSIRNA}" != "1" ]] && return 0

  mkdir -p "$(dirname "${ENDOSIRNA_SUMMARY_TSV}")"
  : > "${ENDOSIRNA_SUMMARY_TSV}"

  for sample in "${DISCOVERED_SAMPLES[@]}"; do
    for stats_file in "$(endosirna_sample_output_dir "${sample}")/"*.stats.tsv; do
      [[ -s "${stats_file}" ]] || continue
      if [[ "${wrote_header}" -eq 0 ]]; then
        sed -n '1p' "${stats_file}" >> "${ENDOSIRNA_SUMMARY_TSV}"
        wrote_header=1
      fi
      sed -n '2,$p' "${stats_file}" >> "${ENDOSIRNA_SUMMARY_TSV}"
    done
  done

  if [[ "${wrote_header}" -eq 0 ]]; then
    printf 'sample\tref_tag\tref_prefix\ttotal_reads\tmapped_reads\tmapped_fraction\ttotal_alignments\tunique_mapped_sequences\tunique_reference_ids\tmodal_length\tmodal_length_alignments\tendosirna_mismatch\tmax_alignments\tbowtie_output\tbowtie_log\n' > "${ENDOSIRNA_SUMMARY_TSV}"
    log_warn "RUN_ENDOSIRNA=1 但没有找到 endo-siRNA stats 输出"
  else
    log_info "endo-siRNA 汇总表已写出: ${ENDOSIRNA_SUMMARY_TSV}"
  fi
}

flatten_sample_output() {
  local sample="$1"
  local sample_out
  local nested_root nested_prefix result_dir src_dir target_dir path base renamed idx path_dir
  local -a nested_roots=()
  local -a nested_prefixes=()

  sample_out="$(sample_output_dir "${sample}")"
  nested_roots=("${sample_out}/1_${sample}" "${sample_out}/1_${sample}.trimmed")
  nested_prefixes=("${sample}" "${sample}.trimmed")

  for idx in "${!nested_roots[@]}"; do
    nested_root="${nested_roots[${idx}]}"
    nested_prefix="${nested_prefixes[${idx}]}"
    [[ -d "${nested_root}" ]] || continue

    result_dir="${nested_root}/${nested_prefix}_result"
    if [[ -d "${result_dir}" ]]; then
      mv "${result_dir}/"* "${sample_out}/" 2>/dev/null || true
    fi

    if [[ "${KEEP_INTERMEDIATE}" == "1" ]]; then
      for extra_dir in "processed" "cutadapter" "fa"; do
        src_dir="${nested_root}/${nested_prefix}_${extra_dir}"
        target_dir="${sample_out}/${sample}_${extra_dir}"
        if [[ -d "${src_dir}" ]]; then
          rm -rf "${target_dir}"
          cp -a "${src_dir}" "${target_dir}"
        fi
      done
    fi

    rm -rf "${nested_root}"
  done

  while IFS= read -r -d '' path; do
    base="$(basename "${path}")"
    path_dir="$(dirname "${path}")"
    renamed="${path_dir}/${sample}_${base#${sample}.trimmed_}"
    mv "${path}" "${renamed}"
  done < <(find "${sample_out}" -depth -name "${sample}.trimmed_*" -print0)
}

run_sample() {
  local sample="$1"
  local fq_gz="$2"
  local sample_fastq="${TMP_DIR}/${sample}.fastq"
  local sample_out
  local sample_log="${SAMPLE_LOG_DIR}/${sample}.sports.log"
  local cutadapt_log="${CUTADAPT_LOG_DIR}/${sample}.cutadapt.log"
  local -a cutadapt_cmd=()
  local cmd
  local run_sports=1
  sample_out="$(sample_output_dir "${sample}")"

  if [[ "${SKIP_EXISTING}" == "1" ]] && is_sample_complete "${sample}"; then
    run_sports=0
    if [[ "${RUN_ENDOSIRNA}" != "1" ]] || is_endosirna_complete "${sample}"; then
      log_info "[${sample}] 已检测到完整结果，跳过"
      return 0
    fi
    log_info "[${sample}] SPORTS 结果已完整，仅补跑缺失的 endo-siRNA 层"
  fi

  if [[ "${run_sports}" == "1" ]]; then
    rm -rf "${sample_out}"
    mkdir -p "${sample_out}"
    : > "${sample_log}"
  else
    mkdir -p "${sample_out}"
    : >> "${sample_log}"
  fi

  if [[ "${PRETRIM_FASTQ}" == "1" ]]; then
    log_info "[${sample}] cutadapt 预剪切 3' adapter / polyG / N，并过滤 ${MIN_LEN}-${MAX_LEN} nt"
    cutadapt_cmd=(
      cutadapt
      -j "${THREADS}"
      -a "${ADAPTER3}"
      --trim-n
      -m "${MIN_LEN}"
      -M "${MAX_LEN}"
    )
    if [[ -n "${CUTADAPT_NEXTSEQ_TRIM}" && "${CUTADAPT_NEXTSEQ_TRIM}" != "0" ]]; then
      cutadapt_cmd+=(--nextseq-trim "${CUTADAPT_NEXTSEQ_TRIM}")
    fi
    if [[ -n "${CUTADAPT_QUALITY_CUTOFF}" ]]; then
      cutadapt_cmd+=(-q "${CUTADAPT_QUALITY_CUTOFF}")
    fi
    cutadapt_cmd+=(-o "${sample_fastq}" "${fq_gz}")

    printf '%q ' "${cutadapt_cmd[@]}" >> "${sample_log}"
    printf '\n' >> "${sample_log}"
    "${cutadapt_cmd[@]}" > "${cutadapt_log}" 2>&1
    run_trimmed_fastqc_for_sample "${sample}" "${sample_fastq}"
  else
    log_info "[${sample}] 解压 fastq"
    gzip -dc "${fq_gz}" > "${sample_fastq}"
  fi

  if [[ "${run_sports}" == "1" ]]; then
    cmd=(
      perl "${SPORTS_PL}"
      -i "${sample_fastq}"
      -p "${THREADS}"
      -l "${MIN_LEN}"
      -L "${MAX_LEN}"
      -M "${MISMATCH}"
      -g "${GENOME_BASE}"
      -m "${MIRNA_BASE}"
      -r "${RRNA_BASE}"
      -t "${TRNA_BASE}"
      -w "${PIRNA_BASE}"
      -e "${ENSEMBL_BASE}"
      -f "${RFAM_BASE}"
      -o "${sample_out}/"
    )

    if [[ "${RUN_ANTISENSE}" == "1" ]]; then
      cmd+=(-s)
    fi

    if [[ "${PRETRIM_FASTQ}" == "1" && "${TRIM_ADAPTER}" == "1" ]]; then
      log_warn "[${sample}] PRETRIM_FASTQ=1 已完成 adapter 预剪切，忽略 TRIM_ADAPTER=1，避免 SPORTS 内部重复剪切"
    elif [[ "${TRIM_ADAPTER}" == "1" ]]; then
      cmd+=(-a -y "${ADAPTER3}")
    fi

    if [[ "${KEEP_INTERMEDIATE}" == "1" ]]; then
      cmd+=(-k)
    fi

    printf '%q ' "${cmd[@]}" >> "${sample_log}"
    printf '\n' >> "${sample_log}"

    log_info "[${sample}] 开始 SPORTS1.1 分析"
    "${cmd[@]}" 2>&1 | tee -a "${sample_log}"

    flatten_sample_output "${sample}"
  else
    log_info "[${sample}] 跳过 SPORTS1.1，复用已有结果"
  fi

  run_endosirna_for_sample "${sample}" "${sample_fastq}"

  if [[ "${KEEP_FASTQ}" != "1" ]]; then
    rm -f "${sample_fastq}"
  fi

  log_info "[${sample}] 完成"
}

render_downstream() {
  if [[ "${RUN_DOWNSTREAM}" != "1" ]]; then
    return 0
  fi

  if [[ ! -f "${RMD_RUNNER}" ]]; then
    log_error "未找到下游 Rmd 调度脚本: ${RMD_RUNNER}"
    exit 1
  fi

  log_info "开始渲染下游 Rmd 模块: ${RMD_RUNNER} --all"
  bash "${RMD_RUNNER}" --all
}

export_worker_context() {
  export ROOT_DIR RAW_DIR RESULT_DIR READS_INFO_FILE
  export METADATA_DIR QC_DIR SPORTS_ANNOTATION_DIR SAMPLE_ROOT ENDOSIRNA_DIR ENDOSIRNA_SAMPLE_ROOT ENDOSIRNA_SUMMARY_DIR
  export LOG_DIR SAMPLE_LOG_DIR CUTADAPT_LOG_DIR ENDOSIRNA_LOG_DIR WORK_DIR TMP_DIR METADATA_TSV NORMALIZED_READS_INFO_TSV OUTPUT_LAYOUT_TSV
  export SPORTS_PL REFERENCE_ROOT DB_ROOT THREADS PARALLEL_JOBS MIN_LEN MAX_LEN MISMATCH
  export KEEP_INTERMEDIATE KEEP_FASTQ PRETRIM_FASTQ TRIM_ADAPTER ADAPTER3 CUTADAPT_NEXTSEQ_TRIM CUTADAPT_QUALITY_CUTOFF
  export RUN_ANTISENSE RUN_FASTQC RUN_ENDOSIRNA FASTQC_DIR
  export ENDOSIRNA_REFS ENDOSIRNA_MISMATCH ENDOSIRNA_MAX_ALIGNMENTS ENDOSIRNA_SUMMARY_TSV
  export GENOME_FASTA GENOME_BASE MIRNA_BASE RRNA_BASE TRNA_BASE TRNA_MITO_BASE PIRNA_BASE ENSEMBL_BASE RFAM_BASE
  export SKIP_EXISTING
  export -f log_info log_warn log_error sample_output_dir endosirna_sample_output_dir trim_shell_string sanitize_endosirna_tag endosirna_ref_specs reference_fasta_exists
  export -f should_run_fastqc run_trimmed_fastqc_for_sample is_sample_complete is_endosirna_complete sample_needs_processing
  export -f write_endosirna_alignment_stats write_endosirna_reference_counts write_endosirna_sequence_counts run_endosirna_for_sample
  export -f flatten_sample_output run_sample
}

main() {
  local idx
  local sample
  local fq_gz
  local -a pending_samples=()
  local -a pending_fastqs=()
  local -a remaining_pairs=()
  local warmup_count="${INDEX_WARMUP_SAMPLES}"
  local warmup_limit=0
  local -a missing_indexes=()
  local start_index=0

  python3 "$ROOT_DIR/scripts/preflight.py" --mode config
  if [[ "$MODE" == "dry-run" ]]; then
    printf 'SPORTS=%s\nINPUT=%s\nOUTPUT=%s\nRUN_ENDOSIRNA=%s\n' "$SPORTS_PL" "$READS_INFO_FILE" "$RESULT_DIR" "$RUN_ENDOSIRNA"
    return 0
  fi
  python3 "$ROOT_DIR/scripts/preflight.py" --mode raw
  if [[ "$PRETRIM_FASTQ" == "1" && -z "$ADAPTER3" ]]; then echo 'ADAPTER3 must be set for actual library kit' >&2; exit 1; fi
  check_prereqs
  check_databases
  if [[ "$MODE" == "check" ]]; then return 0; fi
  if [[ -e "$RESULT_DIR" ]]; then echo "Output exists: $RESULT_DIR; choose a new SPORTS_RESULT_DIR, no stale resume" >&2; exit 1; fi
  initialize_outputs
  cp "$SMALLRNA_CONFIG" "$METADATA_DIR/project.env.frozen"
  cp "$SAMPLES_TSV" "$METADATA_DIR/samples.tsv.frozen"
  sha256sum "$SMALLRNA_CONFIG" "$SAMPLES_TSV" > "$METADATA_DIR/config_sample_sha256.txt"
  python3 "$ROOT_DIR/scripts/preflight.py" --mode raw --sha256 > "$METADATA_DIR/input_sha256.tsv"
  discover_samples
  write_metadata
  run_raw_fastqc_if_requested

  for idx in "${!DISCOVERED_SAMPLES[@]}"; do
    sample="${DISCOVERED_SAMPLES[${idx}]}"
    fq_gz="${DISCOVERED_FASTQS[${idx}]}"
    if ! sample_needs_processing "${sample}"; then
      if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
        log_info "[${sample}] SPORTS 和 endo-siRNA 结果已完整，纳入 metadata 但不重新运行"
      else
        log_info "[${sample}] 已存在完整结果，纳入 metadata 但不重新运行"
      fi
      continue
    fi
    pending_samples+=("${sample}")
    pending_fastqs+=("${fq_gz}")
  done

  export_worker_context

  if ! [[ "${warmup_count}" =~ ^[0-9]+$ ]]; then
    warmup_count=0
  fi

  if [[ "${INDEX_AWARE_SCHEDULING}" == "1" ]]; then
    if all_bowtie_indexes_ready; then
      warmup_count=0
      log_info "Bowtie 索引已全部建立，直接按并行数 ${PARALLEL_JOBS} 调度待处理样本"
    else
      mapfile -t missing_indexes < <(missing_bowtie_indexes)
      log_info "Bowtie 索引尚未全部建立，缺失前缀数: ${#missing_indexes[@]}"
      log_info "先串行预热 ${warmup_count} 个样本，再切换到并行 ${PARALLEL_JOBS} 个样本"
      if [[ "${#missing_indexes[@]}" -gt 0 ]]; then
        log_info "未就绪索引前缀: ${missing_indexes[*]}"
      fi
    fi
  else
    log_info "已关闭索引自动检测，固定先串行预热 ${warmup_count} 个样本，再并行 ${PARALLEL_JOBS} 个样本"
  fi

  warmup_limit="${warmup_count}"
  if [[ "${warmup_limit}" -gt "${#pending_samples[@]}" ]]; then
    warmup_limit="${#pending_samples[@]}"
  fi

  if [[ "${#pending_samples[@]}" -gt 0 ]]; then
    if [[ "${warmup_limit}" -gt 0 ]]; then
      for idx in $(seq 0 $((warmup_limit - 1))); do
        log_info "串行预热样本 [${pending_samples[${idx}]}]"
        run_sample "${pending_samples[${idx}]}" "${pending_fastqs[${idx}]}"
      done
    fi

    start_index="${warmup_limit}"
    if [[ "${start_index}" -lt "${#pending_samples[@]}" ]]; then
      for idx in "${!pending_samples[@]}"; do
        if [[ "${idx}" -lt "${start_index}" ]]; then
          continue
        fi
        remaining_pairs+=("${pending_samples[${idx}]}"$'\t'"${pending_fastqs[${idx}]}")
      done

      log_info "索引准备完成，启动并行处理其余 ${#remaining_pairs[@]} 个样本（并行数: ${PARALLEL_JOBS}）"
      printf '%s\0' "${remaining_pairs[@]}" | xargs -0 -n 1 -P "${PARALLEL_JOBS}" bash -lc '
        pair="$1"
        IFS=$'\''\t'\'' read -r sample fq_gz <<< "${pair}"
        run_sample "${sample}" "${fq_gz}"
      ' _
    else
      log_info "所有待处理样本已完成串行预热，无需并行"
    fi
  else
    log_info "所有样本结果已完整存在，跳过上游重新运行"
  fi

  run_trimmed_multiqc_if_requested
  summarize_endosirna_results
  render_downstream

  log_info "上游分析完成"
  log_info "结果目录: ${RESULT_DIR}"
  log_info "目录结构说明: ${OUTPUT_LAYOUT_TSV}"
  log_info "SPORTS 样本结果: ${SAMPLE_ROOT}"
  if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
    log_info "endo-siRNA 汇总: ${ENDOSIRNA_SUMMARY_TSV}"
  fi
  log_info "样本信息表: ${READS_INFO_FILE}"
  log_info "元数据文件: ${METADATA_TSV}"
  if [[ "${RUN_DOWNSTREAM}" == "1" ]]; then
    log_info "下游 Rmd 已执行: ${RMD_RUNNER} --all"
  else
    log_info "如需生成下游报告，可执行:"
    log_info "RUN_DOWNSTREAM=1 bash ${ROOT_DIR}/sports.sh"
  fi
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
