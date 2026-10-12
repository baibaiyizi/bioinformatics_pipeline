#!/usr/bin/env bash
# 必须由项目入口 source；不创建目录或输出。
TEMPLATE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
set -a
source "${TEMPLATE_ROOT}/config/project.env"
set +a
export ATAC_ROOT_DIR="${TEMPLATE_ROOT}"
# 所有项目输入相对路径均相对复制后的模板根目录。
for path_key in ATAC_COUNTS_FILE ATAC_SAMPLES WORK_DIR READS_INFO TOBIAS_MOTIFS REF_DIR REF_GENOME_INDEX BLACKLIST_BED TSS_BED CHROM_SIZE GENOME_FASTA GENE_GTF TOBIAS_DIR MACS2_ENV IDR_ENV TOBIAS_ENV ATAQV_ENV FEATURECOUNTS_ENV ATAC_MM10_TO_MM39_CHAIN ATAC_ENSEMBL_GTF ATAC_ENSEMBL_TXDB_SQLITE; do
  path_value="${!path_key:-}"
  if [[ -n "$path_value" && "$path_value" != /* ]]; then
    printf -v "$path_key" '%s/%s' "$TEMPLATE_ROOT" "$path_value"
    export "$path_key"
  fi
done
