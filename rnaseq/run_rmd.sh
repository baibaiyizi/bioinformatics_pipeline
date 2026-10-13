#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 纯绘图必须在清理、缓存补跑和setup执行之前分支。
if [[ "${1:-}" == "--plots-only" ]]; then
  shift
  exec python3 "${ROOT_DIR}/rmd/render_plots.py" "$@"
fi

source "${ROOT_DIR}/scripts/config.sh"
export RNASEQ_RMD_DIR="${ROOT_DIR}/rmd"
export RNASEQ_CLEAN_MODULE_RESULT=false RNASEQ_CLEAN_MODULE_CACHE=false
files=(unused 01_main_setup.R 02_selection.Rmd 03_qc.Rmd 04_pca.Rmd 05_de.Rmd 06_enrichment.Rmd 07_go.Rmd 08_gsea.Rmd 09_gsva.Rmd 10_tf_regulator.Rmd 11_ppi.Rmd 12_wgcna.Rmd 13_apa.Rmd 14_der.Rmd 15_splicing.Rmd 16_immune.Rmd 17_pathway_activity.Rmd 18_circle.Rmd)
results=(unused 01_setup 02_selection 03_qc 04_pca 05_de 06_enrichment 07_go 08_gsea 09_gsva 10_tf_regulator 11_ppi 12_wgcna 13_apa 14_der 15_splicing 16_immune 17_pathway_activity 18_circle)
usage() {
 cat <<'EOF'
用法：bash run_rmd.sh [--all | --modules 01,02,05] [--clean | --no-clean]
--all / 无参数：运行 config/project.env 中 RNASEQ_ENABLED_MODULES。
--plots-only：从原图来源及冻结统计结果重绘；必须作为第一个参数。
--list：列出模块；--preflight：只读检查输入；--report-only：从现有结果刷新报告。
默认保留结果与缓存。扩展模块必须显式加入配置；停用模块不会隐式运行。
差异筛选默认 nominal P，结果仅作探索性筛选；保留 BH/FDR 列。
EOF
}
requested=(); report_only=false
while [[ $# -gt 0 ]]; do
 case "$1" in
  --help|-h) usage; exit 0 ;;
  --list) for i in {1..18}; do printf '%02d\t%s\n' "$i" "${files[$i]}"; done; exit 0 ;;
  --all) shift ;;
  --modules) shift ;;
  --clean) export RNASEQ_CLEAN_MODULE_RESULT=true RNASEQ_CLEAN_MODULE_CACHE=true; shift ;;
  --no-clean) shift ;;
  --report-only) report_only=true; shift ;;
  --*) echo "未知参数: $1" >&2; exit 1 ;;
  *) for token in ${1//,/ }; do
       [[ "$token" =~ ^[0-9]+$ ]] || { echo "无效模块: $token" >&2; exit 1; }
       id=$((10#$token)); ((id>=1 && id<=18)) || exit 1; requested+=("$id")
     done; shift ;;
 esac
done
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/input_contract.py" validate
if [[ "$report_only" == true ]]; then
 "${RSCRIPT_BIN}" --vanilla "${ROOT_DIR}/rmd/generate_result_report.R"
 exit 0
fi
read -r -a enabled_raw <<< "${RNASEQ_ENABLED_MODULES//,/ }"
enabled=()
for token in "${enabled_raw[@]}"; do
 [[ "$token" =~ ^[0-9]+$ ]] || { echo "无效配置模块: $token" >&2; exit 1; }
 id=$((10#$token)); (( id>=1 && id<=18 )) || exit 1; enabled+=("$id")
done
if [[ ${#requested[@]} -eq 0 ]]; then requested=("${enabled[@]}"); fi
order=(); seen=' '
add() {
 local id="$1"
 [[ "$seen" == *" $id "* ]] && return 0
 [[ " ${enabled[*]} " == *" $id "* ]] || { echo "模块 $id 未启用；如适用请加入 RNASEQ_ENABLED_MODULES（含依赖）。" >&2; exit 1; }
 [[ "$id" == 1 ]] || add 1
 if (( id>=6 )); then add 5; fi
 case "$id" in 7|8|9|10|13|14|15|16|17|18)
   [[ -f "${ROOT_DIR}/.cache/rmd/06_main_enrichment_tables.rds" && "${RNASEQ_CLEAN_MODULE_CACHE}" == false ]] || add 6 ;;
 esac
 if [[ "$id" == 10 ]]; then
   [[ -f "${ROOT_DIR}/.cache/rmd/09_gsva_state.rds" && "${RNASEQ_CLEAN_MODULE_CACHE}" == false ]] || add 9
 fi
 if [[ "$id" == 17 ]]; then add 13; add 14; add 15; add 16; fi
 seen+="$id "; order+=("$id")
}
for id in "${requested[@]}"; do add "$((10#$id))"; done
"${PYTHON_BIN}" "${ROOT_DIR}/scripts/input_contract.py" prepare
mkdir -p "${ROOT_DIR}/result/logs" "${ROOT_DIR}/.cache/rmd"
status_file="${ROOT_DIR}/result/module_status.tsv"
printf 'module\tstatus\tdetail\n' > "$status_file"
for id in {1..18}; do
 [[ "$seen" == *" $id "* ]] && state=not_run || state=skipped
 printf '%s\t%s\t%s\n' "$id" "$state" "${files[$id]}" >> "$status_file"
done
log_file="${ROOT_DIR}/result/logs/run_rmd_$(date +%Y%m%d_%H%M%S)_$$.log"
exec > >(tee -a "$log_file") 2>&1
for id in "${order[@]}"; do
 echo "[run] module $id ${files[$id]}"
 export RNASEQ_ACTIVE_RESULT_STEP="${results[$id]}"
 if "${RSCRIPT_BIN}" --vanilla "${ROOT_DIR}/rmd/run_module.R" "$id" "${files[$id]}"; then
   sed -i "s/^${id}\tnot_run\t/${id}\tcompleted\t/" "$status_file"
 else
   sed -i "s/^${id}\tnot_run\t/${id}\tfailed\t/" "$status_file"
   exit 1
 fi
done
"${RSCRIPT_BIN}" --vanilla "${ROOT_DIR}/rmd/generate_result_report.R"
