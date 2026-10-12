#!/usr/bin/env bash
set -euo pipefail
# 项目根由入口位置确定，配置中的资源路径由具体数据集根解析。
SMALLRNA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PANDORA_SHARED_ROOT="${SMALLRNA_ROOT}/shared"
project_name=.
analysis_stage=describe
modules=()
list_only=0
check_only=0
usage() {
  cat <<'USAGE'
用法: bash run.sh --project <数据集名|.|绝对数据集路径> [--stage upstream|describe|differential|enrichment] [--modules 名称,...] [--check]
--list 只列出本地数据集。默认 --project . --stage describe。
人类方法比较使用物种目录下的 integration/run.sh；本入口只执行单组学。
--check 只读检查所选阶段的输入，不创建结果或缓存；RUN_ENDOSIRNA 和物种只读取数据集配置。
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project|--stage)
      [[ $# -ge 2 && "$2" != --* ]] || { echo "[ERROR] $1 缺少值" >&2; exit 2; }
      if [[ "$1" == --project ]]; then project_name="$2"; else analysis_stage="$2"; fi
      shift 2 ;;
    --modules)
      shift
      [[ $# -gt 0 && "$1" != --* ]] || { echo '[ERROR] --modules 缺少模块名' >&2; exit 2; }
      while [[ $# -gt 0 && "$1" != --* ]]; do modules+=("$1"); shift; done ;;
    --list) list_only=1; shift ;;
    --check) check_only=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "[ERROR] 未知参数: $1" >&2; usage >&2; exit 2 ;;
  esac
done
case "${analysis_stage}" in upstream|describe|differential|enrichment) ;; *) echo '[ERROR] 无效 stage' >&2; exit 2 ;; esac
if [[ "${analysis_stage}" == upstream && ${#modules[@]} -gt 0 ]]; then
  echo '[ERROR] upstream 不接受 --modules' >&2; exit 2
fi
if [[ "${list_only}" == 1 ]]; then
  [[ ! -f "${SMALLRNA_ROOT}/config/project.sh" ]] || echo .
  for project_config in "${SMALLRNA_ROOT}"/*/config/project.sh; do
    [[ -f "${project_config}" ]] || continue
    basename "$(dirname "$(dirname "${project_config}")")"
  done
  exit 0
fi
if [[ "${project_name}" == /* ]]; then
  project_root="${project_name}"
elif [[ "${project_name}" == . ]]; then
  project_root="${SMALLRNA_ROOT}"
else
  [[ "${project_name}" != */* && "${project_name}" != .. ]] || { echo '[ERROR] 使用数据集名、. 或绝对数据集路径' >&2; exit 2; }
  project_root="${SMALLRNA_ROOT}/${project_name}"
fi
[[ -d "${project_root}" ]] || { echo "[ERROR] 数据集不存在: ${project_root}" >&2; exit 2; }
export PANDORA_PROJECT_ROOT="$(cd "${project_root}" && pwd)"
[[ -f "${PANDORA_PROJECT_ROOT}/config/project.sh" ]] || { echo '[ERROR] 缺少 config/project.sh' >&2; exit 2; }
# 不从调用进程继承项目身份或模式。
unset PROJECT_SPECIES RUN_ENDOSIRNA DB_ROOT REFERENCE_ROOT SPORTS_PL ENDOSIRNA_REFS
set -a
source "${PANDORA_PROJECT_ROOT}/config/project.sh"
set +a
case "${PROJECT_SPECIES:-}" in 'Mus musculus') export PANDORA_SPECIES=Mus ;; 'Homo sapiens') export PANDORA_SPECIES=Homo ;; *) echo '[ERROR] PROJECT_SPECIES 必须显式配置为 Mus musculus 或 Homo sapiens' >&2; exit 2 ;; esac
case "${RUN_ENDOSIRNA:-}" in 1) export PANDORA_MODE=with_endo ;; 0) export PANDORA_MODE=without_endo ;; *) echo '[ERROR] RUN_ENDOSIRNA 必须显式配置为 0 或 1' >&2; exit 2 ;; esac
: "${PANDORA_PROJECT_NAME:?缺少 PANDORA_PROJECT_NAME}"
: "${PANDORA_DESIGN:?缺少 PANDORA_DESIGN}"
: "${SPORTS_LAYOUT:?缺少 SPORTS_LAYOUT}"
resolve_path() {
  python3 - "$PANDORA_PROJECT_ROOT" "$1" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[2])
print((p if p.is_absolute() else Path(sys.argv[1]) / p).resolve())
PY
}
for resource in DB_ROOT REFERENCE_ROOT SPORTS_PL; do
  [[ -n "${!resource:-}" ]] || { echo "[ERROR] 缺少 ${resource} 配置" >&2; exit 2; }
  printf -v "${resource}" '%s' "$(resolve_path "${!resource}")"
  export "${resource}"
done
if [[ -n "${EXOGENOUS_DB_ROOT:-}" ]]; then export EXOGENOUS_DB_ROOT="$(resolve_path "${EXOGENOUS_DB_ROOT}")"; fi
if [[ "${RUN_ENDOSIRNA}" == 1 ]]; then
  [[ -n "${ENDOSIRNA_REFS:-}" ]] || { echo '[ERROR] RUN_ENDOSIRNA=1 时必须配置 ENDOSIRNA_REFS=标签=参考前缀;...' >&2; exit 2; }
  IFS=';' read -r -a reference_specs <<< "${ENDOSIRNA_REFS}"
  resolved_refs=()
  for spec in "${reference_specs[@]}"; do
    [[ "$spec" == *=* && -n "${spec%%=*}" && -n "${spec#*=}" ]] || { echo '[ERROR] 无效 ENDOSIRNA_REFS 项' >&2; exit 2; }
    resolved_refs+=("${spec%%=*}=$(resolve_path "${spec#*=}")")
  done
  ENDOSIRNA_REFS="$(IFS=';'; printf '%s' "${resolved_refs[*]}")"
  export ENDOSIRNA_REFS
else
  export ENDOSIRNA_REFS=''
fi
export PANDORA_STAGE="${analysis_stage}"
export PANDORA_RESULT_DIR="${PANDORA_PROJECT_ROOT}/result/${analysis_stage}"
export PANDORA_SAMPLES_FILE="${PANDORA_PROJECT_ROOT}/config/samples.tsv"
export SAMPLES_FILE="${PANDORA_SAMPLES_FILE}" READS_INFO_FILE="${PANDORA_SAMPLES_FILE}"
export ROOT_DIR="${PANDORA_PROJECT_ROOT}" RESULT_DIR="${PANDORA_RESULT_DIR}"
export SHARED_ROOT="${PANDORA_SHARED_ROOT}" SPECIES="${PANDORA_SPECIES}" STAGE="${analysis_stage}"
export PROJECT_CONTRAST_BASELINE="${CONTRAST_BASELINE:-}" PROJECT_CONTRAST_CASE="${CONTRAST_CASE:-}"
export ENDOSIRNA_MISMATCH=0 ENDOSIRNA_MAX_ALIGNMENTS=20
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}" R_DATATABLE_NUM_THREADS="${R_DATATABLE_NUM_THREADS:-1}"
module_arg="$(IFS=,; printf '%s' "${modules[*]}")"
python3 "${PANDORA_SHARED_ROOT}/preflight.py" --stage "${analysis_stage}" --modules "${module_arg}"
[[ "${check_only}" == 0 ]] || exit 0
echo "[INFO] project=${PANDORA_PROJECT_NAME} species=${PANDORA_SPECIES} mode=${PANDORA_MODE} stage=${analysis_stage}"
cd "${PANDORA_PROJECT_ROOT}"
if [[ "${analysis_stage}" == upstream ]]; then
  bash "${PANDORA_SHARED_ROOT}/upstream/sports_$(printf '%s' "${PANDORA_SPECIES}" | tr '[:upper:]' '[:lower:]').sh"
else
  r_args=(--stage "${analysis_stage}")
  [[ -z "${module_arg}" ]] || r_args+=(--modules "${module_arg}")
  Rscript "${PANDORA_SHARED_ROOT}/R/run_analysis.R" "${r_args[@]}"
fi
