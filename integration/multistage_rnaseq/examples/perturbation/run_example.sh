#!/usr/bin/env bash
set -euo pipefail

EXAMPLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_DIR="$(cd "$EXAMPLE_DIR/../.." && pwd)"
INPUT_DIR="$EXAMPLE_DIR/input"
OUTPUT_DIR="$EXAMPLE_DIR/output"
mkdir -p "$INPUT_DIR"
rm -rf "$OUTPUT_DIR"

Rscript - "$INPUT_DIR" <<'RS'
args <- commandArgs(TRUE)
out <- args[[1L]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
genes <- paste0("G", sprintf("%02d", 1:40))
ranked <- do.call(rbind, lapply(c("StageA", "StageB"), function(stage) {
  stat <- if (stage == "StageA") seq(4, -4, length.out=40) else seq(-2, 2, length.out=40)
  data.frame(stage=stage, gene_id=paste0("ENSMUSG", sprintf("%011d", 1:40)), gene_symbol=genes,
             stat=stat, log2FC=stat/2, padj=p.adjust(2*pnorm(-abs(stat)), method="BH"), stringsAsFactors=FALSE)
}))
write.table(ranked, file.path(out,"ranked.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
samples <- expand.grid(stage=c("StageA","StageB"), age_group=c("MY","MO"), replicate=1:3, stringsAsFactors=FALSE)
samples$sample_id <- paste(samples$stage,samples$age_group,samples$replicate,sep="_")
write.table(samples[c("sample_id","stage","age_group")], file.path(out,"samples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
vst <- data.frame(gene_id=paste0("ENSMUSG", sprintf("%011d", 1:40)), stringsAsFactors=FALSE)
for (i in seq_len(nrow(samples))) {
  base <- seq(5,9,length.out=40)
  effect <- ifelse(samples$age_group[i]=="MO", ranked$log2FC[ranked$stage==samples$stage[i]], 0)
  vst[[samples$sample_id[i]]] <- base + effect + (samples$replicate[i]-2)*0.05
}
write.table(vst, file.path(out,"vst.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
signatures <- rbind(
  data.frame(signature_id="S_A_up",source_stage="StageA",direction="up_in_mzKO",gene_symbol=genes[1:6]),
  data.frame(signature_id="S_A_down",source_stage="StageA",direction="down_in_mzKO",gene_symbol=genes[35:40]),
  data.frame(signature_id="S_B_up",source_stage="StageB",direction="up_in_mzKO",gene_symbol=genes[35:40]),
  data.frame(signature_id="S_missing",source_stage="StageB",direction="down_in_mzKO",gene_symbol=paste0("X",1:6))
)
write.table(signatures, file.path(out,"obox.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
pabpn <- do.call(rbind,lapply(c("StageA","StageB"), function(stage) data.frame(
  stage=stage,gene_symbol=genes,mean_KD=2^(seq(-1,1,length.out=40)+2),mean_WT=4,
  pass_expression=TRUE,effect_eps1=seq(-1,1,length.out=40),effect_eps0_1=seq(-1,1,length.out=40)*1.02)))
write.table(pabpn, file.path(out,"pabpn.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
RS

Rscript "$TEMPLATE_DIR/scripts/perturbation_signatures.R" \
  --ranked "$INPUT_DIR/ranked.tsv" \
  --vst "$INPUT_DIR/vst.tsv" \
  --samples "$INPUT_DIR/samples.tsv" \
  --obox-signatures "$INPUT_DIR/obox.tsv" \
  --pabpn-effects "$INPUT_DIR/pabpn.tsv" \
  --config "$EXAMPLE_DIR/config.json" \
  --output "$OUTPUT_DIR"

Rscript - "$OUTPUT_DIR" <<'RS'
args <- commandArgs(TRUE)
x <- read.delim(file.path(args[[1L]],"tables","obox_fgsea_results.tsv"), check.names=FALSE)
stopifnot(nrow(x)==8L)
stopifnot(all(x$test_status[x$signature_id=="S_missing"]=="not_tested_below_min_size"))
stopifnot(x$concordance_NES[x$signature_id=="S_A_up" & x$comparison_stage=="StageA"] > 0)
stopifnot(x$concordance_NES[x$signature_id=="S_A_down" & x$comparison_stage=="StageA"] > 0)
m <- read.delim(file.path(args[[1L]],"tables","obox_signature_mapping.tsv"), check.names=FALSE)
stopifnot(any(m$mapped), any(!m$mapped))
cat("Replayable perturbation example passed direction, exact-symbol mapping and missing-set checks.\n")
RS
