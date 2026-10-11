#!/usr/bin/env bash
# 在新临时目录重放合成数据闭环；不运行 FASTQ 上游或联网扩展模块。
set -euo pipefail
VALIDATION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_DIR="$(cd "${VALIDATION_DIR}/.." && pwd)"
R_BIN="$(command -v "${RSCRIPT_BIN:-Rscript}")"
RUN_DIR="$(mktemp -d "${TMPDIR:-/tmp}/bulk replay core.XXXXXX")"
PROJECT_DIR="${RUN_DIR}/project"
mkdir -p "${PROJECT_DIR}/rawdata"
cp -a "${TEMPLATE_DIR}/rmd" "${TEMPLATE_DIR}/scripts" "${TEMPLATE_DIR}/config" "${PROJECT_DIR}/"
for script in rnaseq_first.sh rnaseq_second.sh rnaseq_third.sh run_rmd.sh; do
  cp "${TEMPLATE_DIR}/${script}" "${PROJECT_DIR}/${script}"
done
cat >> "${PROJECT_DIR}/config/project.env" <<'EOF'
RNASEQ_INPUT_MODE=counts
RNASEQ_COUNTS_FILE=rawdata/counts.tsv
RNASEQ_SAMPLES_TSV=config/samples.tsv
RNASEQ_DESIGN=independent_two_group
RNASEQ_GROUP1=control
RNASEQ_GROUP2=treatment
RNASEQ_SPECIES=mouse
RNASEQ_ENABLED_MODULES="1 2 3 4 5"
RNASEQ_IMMUNE_APPLICABLE=false
RNASEQ_DE_USE_PADJ=false
RNASEQ_DE_PVALUE_CUTOFF=0.05
RNASEQ_DE_LOG2FC_CUTOFF=1
EOF
printf 'RSCRIPT_BIN=%q\n' "${R_BIN}" >> "${PROJECT_DIR}/config/project.env"
"${R_BIN}" --vanilla "${VALIDATION_DIR}/generate_counts_fixture.R" "${PROJECT_DIR}"
bash "${PROJECT_DIR}/rnaseq_second.sh" --all > "${RUN_DIR}/core_de.log" 2>&1
bash "${PROJECT_DIR}/run_rmd.sh" --modules 2,3,4,5 > "${RUN_DIR}/core_rmd.log" 2>&1
printf 'PASS core replay: %s\n' "${RUN_DIR}"
cat "${PROJECT_DIR}/result/module_status.tsv"
