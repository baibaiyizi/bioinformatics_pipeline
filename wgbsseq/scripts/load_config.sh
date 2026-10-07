#!/usr/bin/env bash
# 必须由项目入口 source；不创建目录或输出。
TEMPLATE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
set -a
source "${TEMPLATE_ROOT}/config/project.env"
set +a
export WGBS_ROOT_DIR="${TEMPLATE_ROOT}"
# 所有项目输入相对路径均相对复制后的模板根目录。
for path_key in WGBS_SAMPLES WGBS_REF_DIR WGBS_FASTA WGBS_GTF WGBS_RNA_SEQ_FILE WGBS_RNA_KEGG_FILE WGBS_ATAC_PEAK_DIR WGBS_ATAC_DIFF_DIR WGBS_ENV_BIN WGBS_LOG_DIR WGBS_SPERM_ICR_FILE WGBS_SPERMATOGENESIS_GENE_FILE; do
  path_value="${!path_key:-}"
  if [[ -n "$path_value" && "$path_value" != /* ]]; then
    printf -v "$path_key" '%s/%s' "$TEMPLATE_ROOT" "$path_value"
    export "$path_key"
  fi
done
