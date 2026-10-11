#!/usr/bin/env bash
set -euo pipefail

# SPORTS 公共执行主体；sports_homo.sh 提供人类参考与原有剪切策略。
# 样本及设计只读取 config/samples.tsv；此文件不覆盖任何设计元数据。

ROOT_DIR="${PANDORA_PROJECT_ROOT:?请通过项目 run.sh 启动}"
RAW_DIR="${ROOT_DIR}/01rawdata"
RESULT_DIR="${ROOT_DIR}/sportsresult"
READS_INFO_FILE="${PANDORA_SAMPLES_FILE:-${ROOT_DIR}/config/samples.tsv}"

METADATA_DIR="${RESULT_DIR}/02_metadata"
QC_DIR="${RESULT_DIR}/03_qc"
SPORTS_ANNOTATION_DIR="${RESULT_DIR}/04_sports_annotation"
SAMPLE_ROOT="${PANDORA_SPORTS_SAMPLE_ROOT:-${SPORTS_ANNOTATION_DIR}/samples}"
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

SPORTS_PL="${SPORTS_PL:?缺少 SPORTS_PL 配置}"
REFERENCE_ROOT="${REFERENCE_ROOT:?缺少参考目录配置}"
# PANDORA-local small RNA reference root. Override with DB_ROOT=/path/to/Mus_musculus if needed.
DB_ROOT="${DB_ROOT:?缺少物种数据库配置}"

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
RUN_FASTQC="${RUN_FASTQC:-auto}"
RUN_ENDOSIRNA="${RUN_ENDOSIRNA:?流程模式必须由分支入口提供}"
SKIP_EXISTING="${SKIP_EXISTING:-1}"
FASTQC_DIR="${QC_DIR}/fastqc"

case "${RUN_ENDOSIRNA}" in
  true|TRUE|yes|YES) RUN_ENDOSIRNA=1 ;;
  false|FALSE|no|NO) RUN_ENDOSIRNA=0 ;;
esac

GENOME_FASTA="${GENOME_FASTA:-${REFERENCE_ROOT}/Mus_musculus.GRCm39.dna.primary_assembly.fa}"
GENOME_BASE="${GENOME_BASE:-${REFERENCE_ROOT}/indexes/bowtie1/Mus_musculus.GRCm39.dna.primary_assembly}"
MIRNA_BASE="${MIRNA_BASE:-${DB_ROOT}/miRBase/22.1/miRBase_22.1-mmu}"
RRNA_BASE="${RRNA_BASE:-${DB_ROOT}/rRNAdb/mouse_rRNA}"
TRNA_BASE="${TRNA_BASE:-${DB_ROOT}/GtRNAdb/mm39/mm39-tRNAs}"
TRNA_MITO_BASE="${TRNA_MITO_BASE:-${DB_ROOT}/GtRNAdb/mm39/mm39-mt_tRNAs}"
PIRNA_BASE="${PIRNA_BASE:-${DB_ROOT}/piRNA/piRBase/mmu.v3.0}"
ENSEMBL_BASE="${ENSEMBL_BASE:-${DB_ROOT}/Ensembl/release-115/Mus_musculus.GRCm39.ncrna}"
RFAM_BASE="${RFAM_BASE:-${DB_ROOT}/Rfam/15.0/Rfam-15.0-mouse}"
ENDOSIRNA_REFS="${ENDOSIRNA_REFS:-wu2016=${DB_ROOT}/endo-siRNA/wuall/endo-siRNA-wu2016;mouse_embryos_RNase6=${DB_ROOT}/endo-siRNA/mouse_embryos/endo-siRNA-mouse_embryos_RNase6;mouse_embryos_UltraPandora=${DB_ROOT}/endo-siRNA/mouse_embryos/endo-siRNA-mouse_embryos_UltraPandora}"
ENDOSIRNA_MISMATCH="${ENDOSIRNA_MISMATCH:-0}"
ENDOSIRNA_MAX_ALIGNMENTS="${ENDOSIRNA_MAX_ALIGNMENTS:-20}"
ENDOSIRNA_SUMMARY_TSV="${ENDOSIRNA_SUMMARY_DIR}/endo_siRNA_summary.tsv"

read -r -a RRNA_SUBTYPES <<< "${RRNA_SUBTYPE_LIST:-4.5S 5S 5.8S 12S 16S 18S 28S 45S RNY1 RNY3}"

declare -a DISCOVERED_SAMPLES=()
declare -a DISCOVERED_FASTQS=()
declare -a DISCOVERED_GROUPS=()
declare -a DISCOVERED_REPLICATES=()
declare -a DISCOVERED_CONDITIONS=()
declare -a DISCOVERED_AGE_GROUPS=()
declare -a DISCOVERED_SPECIES=()
declare -a DISCOVERED_STAGES=()
declare -a BOWTIE_INDEX_PREFIXES=(
  "${GENOME_BASE}"
  "${MIRNA_BASE}"
  "${TRNA_BASE}"
  "${TRNA_MITO_BASE}"
  "${PIRNA_BASE}"
  "${ENSEMBL_BASE}"
  "${RFAM_BASE}"
)

if [[ "${PANDORA_SPECIES}" == Mus ]]; then
  BOWTIE_INDEX_PREFIXES+=("${TRNA_BASE}_CCA" "${TRNA_MITO_BASE}_CCA")
fi

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

check_prereqs() {
  mkdir -p \
    "${RESULT_DIR}" \
    "${METADATA_DIR}" \
    "${TMP_DIR}" \
    "${SAMPLE_ROOT}" \
    "${LOG_DIR}" \
    "${SAMPLE_LOG_DIR}" \
    "${CUTADAPT_LOG_DIR}" \
    "${FASTQC_DIR}"
  if [[ "${RUN_ENDOSIRNA}" == 1 ]]; then
    mkdir -p "${ENDOSIRNA_SAMPLE_ROOT}" "${ENDOSIRNA_SUMMARY_DIR}" "${ENDOSIRNA_LOG_DIR}"
  fi
  write_output_layout_manifest

  if [[ ! -x "${SPORTS_PL}" ]]; then
    log_error "sports.pl 不存在或不可执行: ${SPORTS_PL}"
    exit 1
  fi

  local executables=(perl bowtie bowtie-build gzip xargs python3 sha256sum)
  local exe
  for exe in "${executables[@]}"; do
    if ! command -v "${exe}" >/dev/null 2>&1; then
      log_error "未找到依赖程序: ${exe}"
      exit 1
    fi
  done

  if [[ "${PRETRIM_FASTQ}" == "1" || "${TRIM_ADAPTER}" == "1" || "${RUN_ENDOSIRNA}" == "1" ]] && ! command -v cutadapt >/dev/null 2>&1; then
    log_error "PRETRIM_FASTQ=1 或 TRIM_ADAPTER=1 时需要 cutadapt"
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
    log_error "中央物种基因组 FASTA 不存在或为空: ${GENOME_FASTA}"
    exit 1
  fi

  if ! bowtie_index_exists "${GENOME_BASE}"; then
    log_error "中央物种 Bowtie1 索引不完整: ${GENOME_BASE}"
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
      log_error "缺少 物种 rRNA 数据库文件: ${RRNA_BASE}_${rrna_subtype}.fa"
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
      log_warn "RUN_ENDOSIRNA=1 时建议 PRETRIM_FASTQ=1；endo-siRNA 将按 SPORTS 的 adapter/N/长度参数另行准备同口径输入。"
    fi
  else
    log_info "当前流程使用所选物种参考库；endo-siRNA 独立证据层未启用，仅使用 genome/miRNA/rRNA/tRNA/piRNA/Ensembl/Rfam。"
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

# 仅从权威样本表读取精确路径；不再猜文件名或覆盖样本设计。
discover_samples() {
  DISCOVERED_SAMPLES=(); DISCOVERED_FASTQS=(); DISCOVERED_GROUPS=()
  DISCOVERED_REPLICATES=(); DISCOVERED_CONDITIONS=(); DISCOVERED_AGE_GROUPS=()
  DISCOVERED_SPECIES=(); DISCOVERED_STAGES=()
  python3 - "${READS_INFO_FILE}" "${ROOT_DIR}" > "${NORMALIZED_READS_INFO_TSV}" <<'PY_READER'
import csv, re, sys
from pathlib import Path
manifest, root = Path(sys.argv[1]), Path(sys.argv[2])
with manifest.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"sample", "group", "replicate", "raw_fastq_gz"}
    if not required.issubset(reader.fieldnames or []):
        raise SystemExit("config/samples.tsv 缺少必需列: " + ", ".join(sorted(required - set(reader.fieldnames or []))))
    rows = list(reader)
seen = set()
for row in rows:
    sample = row["sample"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample) or sample in {".", ".."} or sample in seen:
        raise SystemExit("无效或重复样本 ID: " + sample)
    seen.add(sample)
    raw = Path(row["raw_fastq_gz"])
    raw = raw if raw.is_absolute() else root / raw
    if not raw.is_file():
        raise SystemExit("上游不能运行：原始 FASTQ 缺失: " + str(raw) + "；已有 SPORTS 输入仍可运行 describe。")
    if not row["group"] or not row["replicate"]:
        raise SystemExit("group/replicate 不能为空: " + sample)
    values = [row["group"], sample, row["replicate"], str(raw), row.get("condition", ""), row.get("age_group", ""), row.get("species", ""), row.get("stage", "")]
    print("\t".join(value or "." for value in values))
if not rows:
    raise SystemExit("config/samples.tsv 没有样本")
PY_READER
  local group sample replicate fq_gz condition age_group species stage
  while IFS=$'\t' read -r group sample replicate fq_gz condition age_group species stage; do
    DISCOVERED_GROUPS+=("${group}"); DISCOVERED_SAMPLES+=("${sample}")
    DISCOVERED_REPLICATES+=("${replicate}"); DISCOVERED_FASTQS+=("${fq_gz}")
    DISCOVERED_CONDITIONS+=("${condition}"); DISCOVERED_AGE_GROUPS+=("${age_group}")
    DISCOVERED_SPECIES+=("${species}"); DISCOVERED_STAGES+=("${stage}")
  done < "${NORMALIZED_READS_INFO_TSV}"
  log_info "从权威配置读取 ${#DISCOVERED_SAMPLES[@]} 个样本: ${READS_INFO_FILE}"
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

# 复用要求样本输入、参考内容与参数一致，且所有派生计数完整。
endosirna_expected_signature() {
  local sample="$1" ref_prefix="$2" raw_path ref_fasta
  raw_path="$(awk -F '\t' -v sample="${sample}" 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next} $(h["sample"])==sample {print $(h["raw_fastq_gz"]); exit}' "${READS_INFO_FILE}")"
  [[ "${raw_path}" == /* ]] || raw_path="${ROOT_DIR}/${raw_path}"
  [[ -f "${raw_path}" ]] || return 1
  if [[ -f "${ref_prefix}.fa" ]]; then ref_fasta="${ref_prefix}.fa"; else ref_fasta="${ref_prefix}.fasta"; fi
  {
    printf '%s\n' 'endosirna_input_v1' "${sample}" "${ADAPTER3}" "${PRETRIM_FASTQ}" "${TRIM_ADAPTER}" "${CUTADAPT_NEXTSEQ_TRIM}" "${CUTADAPT_QUALITY_CUTOFF}" "${MIN_LEN}" "${MAX_LEN}" "${ENDOSIRNA_MISMATCH}" "${ENDOSIRNA_MAX_ALIGNMENTS}"
    stat -Lc '%n %s %Y' "${raw_path}"
    sha256sum "${ref_fasta}"
    bowtie --version | head -n 1
  } | sha256sum | cut -d ' ' -f 1
}

endosirna_reference_complete() {
  local sample="$1" ref_tag="$2" ref_prefix="$3" out_dir signature expected
  out_dir="$(endosirna_sample_output_dir "${sample}")"
  signature="${out_dir}/${sample}.${ref_tag}.run_signature.sha256"
  [[ -s "${signature}" && -s "${out_dir}/${sample}.${ref_tag}.stats.tsv" &&
     -s "${out_dir}/${sample}.${ref_tag}.reference_counts.tsv" &&
     -s "${out_dir}/${sample}.${ref_tag}.sequence_counts.tsv" &&
     -f "${out_dir}/${sample}.${ref_tag}.bowtie.tsv" ]] || return 1
  expected="$(endosirna_expected_signature "${sample}" "${ref_prefix}")" || return 1
  [[ "$(cat "${signature}")" == "${expected}" ]]
}

is_endosirna_complete() {
  local sample="$1" ref_count=0 ref_tag ref_prefix
  [[ "${RUN_ENDOSIRNA}" != 1 ]] && return 0
  while IFS=$'\t' read -r ref_tag ref_prefix; do
    ref_count=$((ref_count + 1))
    endosirna_reference_complete "${sample}" "${ref_tag}" "${ref_prefix}" || return 1
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

  if [[ ! -f "${fastq}" ]]; then
    log_error "[${sample}] endo-siRNA 输入 FASTQ 不存在: ${fastq}"
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

    if [[ "${SKIP_EXISTING}" == "1" ]] && endosirna_reference_complete "${sample}" "${ref_tag}" "${ref_prefix}"; then
      log_info "[${sample}] endo-siRNA ${ref_tag} 统计已存在，跳过"
      continue
    fi

    rm -f "${out_dir}/${sample}.${ref_tag}.run_signature.sha256"

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
    endosirna_expected_signature "${sample}" "${ref_prefix}" > "${out_dir}/${sample}.${ref_tag}.run_signature.sha256.tmp"
    mv "${out_dir}/${sample}.${ref_tag}.run_signature.sha256.tmp" "${out_dir}/${sample}.${ref_tag}.run_signature.sha256"
  done < <(endosirna_ref_specs)
}

summarize_endosirna_results() {
  local sample stats_file ref_tag ref_prefix
  local wrote_header=0

  [[ "${RUN_ENDOSIRNA}" != "1" ]] && return 0

  mkdir -p "$(dirname "${ENDOSIRNA_SUMMARY_TSV}")"
  : > "${ENDOSIRNA_SUMMARY_TSV}"

  for sample in "${DISCOVERED_SAMPLES[@]}"; do
    while IFS=$'\t' read -r ref_tag ref_prefix; do
      stats_file="$(endosirna_sample_output_dir "${sample}")/${sample}.${ref_tag}.stats.tsv"
      [[ -s "${stats_file}" ]] || { log_error "缺少当前参考统计: ${stats_file}"; return 1; }
      if [[ "${wrote_header}" -eq 0 ]]; then
        sed -n '1p' "${stats_file}" >> "${ENDOSIRNA_SUMMARY_TSV}"
        wrote_header=1
      fi
      sed -n '2,$p' "${stats_file}" >> "${ENDOSIRNA_SUMMARY_TSV}"
    done < <(endosirna_ref_specs)
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

  if [[ "${RUN_ENDOSIRNA}" == 1 && "${PRETRIM_FASTQ}" != 1 ]]; then
    # Homo 保留 SPORTS 内部剪切；证据层单独复制同样 adapter/N/长度口径。
    local endo_fastq="${TMP_DIR}/${sample}.endo.fastq"
    local -a endo_trim=(cutadapt -j "${THREADS}" --max-n 0 -m "${MIN_LEN}" -M "${MAX_LEN}")
    if [[ "${TRIM_ADAPTER}" == 1 ]]; then endo_trim+=(-a "${ADAPTER3}"); fi
    "${endo_trim[@]}" -o "${endo_fastq}" "${sample_fastq}" > "${CUTADAPT_LOG_DIR}/${sample}.endo.log" 2>&1
    run_endosirna_for_sample "${sample}" "${endo_fastq}"
    [[ "${KEEP_FASTQ}" == 1 ]] || rm -f "${endo_fastq}"
  else
    run_endosirna_for_sample "${sample}" "${sample_fastq}"
  fi

  if [[ "${KEEP_FASTQ}" != "1" ]]; then
    rm -f "${sample_fastq}"
  fi

  log_info "[${sample}] 完成"
}

export_worker_context() {
  export ROOT_DIR RAW_DIR RESULT_DIR READS_INFO_FILE PANDORA_SPECIES PANDORA_PROJECT_NAME
  export METADATA_DIR QC_DIR SPORTS_ANNOTATION_DIR SAMPLE_ROOT ENDOSIRNA_DIR ENDOSIRNA_SAMPLE_ROOT ENDOSIRNA_SUMMARY_DIR
  export LOG_DIR SAMPLE_LOG_DIR CUTADAPT_LOG_DIR ENDOSIRNA_LOG_DIR WORK_DIR TMP_DIR METADATA_TSV NORMALIZED_READS_INFO_TSV OUTPUT_LAYOUT_TSV
  export SPORTS_PL REFERENCE_ROOT DB_ROOT THREADS PARALLEL_JOBS MIN_LEN MAX_LEN MISMATCH
  export KEEP_INTERMEDIATE KEEP_FASTQ PRETRIM_FASTQ TRIM_ADAPTER ADAPTER3 CUTADAPT_NEXTSEQ_TRIM CUTADAPT_QUALITY_CUTOFF
  export RUN_ANTISENSE RUN_FASTQC RUN_ENDOSIRNA FASTQC_DIR
  export ENDOSIRNA_REFS ENDOSIRNA_MISMATCH ENDOSIRNA_MAX_ALIGNMENTS ENDOSIRNA_SUMMARY_TSV
  export GENOME_FASTA GENOME_BASE MIRNA_BASE RRNA_BASE TRNA_BASE TRNA_MITO_BASE PIRNA_BASE ENSEMBL_BASE RFAM_BASE
  export SKIP_EXISTING
  export -f log_info log_warn log_error sample_output_dir endosirna_sample_output_dir trim_shell_string sanitize_endosirna_tag endosirna_ref_specs reference_fasta_exists endosirna_expected_signature endosirna_reference_complete
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

  check_prereqs
  check_databases
  discover_samples
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

  log_info "上游分析完成"
  log_info "结果目录: ${RESULT_DIR}"
  log_info "目录结构说明: ${OUTPUT_LAYOUT_TSV}"
  log_info "SPORTS 样本结果: ${SAMPLE_ROOT}"
  if [[ "${RUN_ENDOSIRNA}" == "1" ]]; then
    log_info "endo-siRNA 汇总: ${ENDOSIRNA_SUMMARY_TSV}"
  fi
  log_info "样本信息表: ${READS_INFO_FILE}"
  log_info "权威样本表（只读）: ${READS_INFO_FILE}"
  log_info "描述性结果请通过项目 run.sh --project ${PANDORA_PROJECT_NAME} --stage describe 生成"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
