#!/usr/bin/env bash
set -euo pipefail
: "${PANDORA_PROJECT_ROOT:?请通过项目 run.sh 启动}"
: "${PANDORA_SHARED_ROOT:?缺少共享目录}"
: "${DB_ROOT:?缺少人类数据库目录}"
: "${REFERENCE_ROOT:?缺少 GRCh38 目录}"
export PRETRIM_FASTQ=0 TRIM_ADAPTER=1 CUTADAPT_NEXTSEQ_TRIM=0
export PANDORA_SPORTS_SAMPLE_ROOT="${PANDORA_PROJECT_ROOT}/sportsresult"
export GENOME_FASTA="${REFERENCE_ROOT}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
export GENOME_BASE="${REFERENCE_ROOT}/indexes/bowtie1/Homo_sapiens.GRCh38.dna.primary_assembly"
export MIRNA_BASE="${DB_ROOT}/miRBase/22.1/miRBase_22.1-hsa"
export RRNA_BASE="${DB_ROOT}/rRNAdb/human_rRNA"
export TRNA_BASE="${DB_ROOT}/GtRNAdb/hg38/hg38-tRNAs"
export TRNA_MITO_BASE="${DB_ROOT}/GtRNAdb/hg38/hg38-mt_tRNAs"
export PIRNA_BASE="${DB_ROOT}/piRNA/piR_human"
export ENSEMBL_BASE="${DB_ROOT}/Ensembl/release-115/Homo_sapiens.GRCh38.ncrna"
export RFAM_BASE="${DB_ROOT}/Rfam/15.0/Rfam-15.0-human"
export RRNA_SUBTYPE_LIST="5S 5.8S 12S 16S 18S 28S 45S other RNY1 RNY3 RNY4 RNY5"
export ENDOSIRNA_REFS="${ENDOSIRNA_REFS:-}"
source "${PANDORA_SHARED_ROOT}/upstream/sports_mus.sh"
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
