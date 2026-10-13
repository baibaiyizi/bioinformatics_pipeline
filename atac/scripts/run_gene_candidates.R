#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1]), fixed=TRUE)
root <- if(length(args)) normalizePath(args[[1]], mustWork=TRUE) else dirname(dirname(normalizePath(script, mustWork=TRUE)))
input <- file.path(root,'result/07_peak_gene_enrichment/tables/Annotated_DiffPeaks.csv')
setup <- file.path(root,'rmd/01_main_setup.R')
stopifnot(file.exists(input),file.exists(setup))
suppressPackageStartupMessages(library(dplyr))
for(e in parse(setup)) if(is.call(e) && identical(e[[1]], as.name('<-')) && is.symbol(e[[2]]) &&
 identical(as.character(e[[2]]),'build_atac_gene_candidates')) eval(e, envir=globalenv())
peaks <- read.csv(input,stringsAsFactors=FALSE)
if(!'Comparison' %in% names(peaks) || length(unique(peaks$Comparison)) != 1L) stop('ATAC 来源缺少唯一比较方向')
genes <- build_atac_gene_candidates(peaks)
output <- file.path(root,'result/13_keygenes_locus_reports')
dir.create(file.path(output,'tables'),recursive=TRUE,showWarnings=FALSE)
write.csv(genes,file.path(output,'tables/ATAC_Gene_Candidates.csv'),row.names=FALSE)
key <- head(genes$geneSymbol,12)
write.csv(peaks[peaks$geneSymbol %in% key,],file.path(output,'tables/key_gene_peaks.csv'),row.names=FALSE)
write.csv(data.frame(geneSymbol=key,PriorityRank=seq_along(key)),file.path(output,'tables/Representative_Loci.csv'),row.names=FALSE)
readr::write_tsv(data.frame(Input='result/07_peak_gene_enrichment/tables/Annotated_DiffPeaks.csv',
 SHA256=digest::digest(file=input,algo='sha256'),SourceStatus='PRESERVED_ANNOTATED_PEAK_OUTPUT',
 Comparison=unique(peaks$Comparison),Rows=nrow(peaks),SourceSetupSHA256=digest::digest(file=setup,algo='sha256')),file.path(output,'input_manifest.tsv'))
readr::write_tsv(data.frame(time=format(Sys.time(),tz='UTC',usetz=TRUE),stage='pure_atac_gene_candidates',status='completed',
 source_status='PRESERVED_ANNOTATED_PEAK_OUTPUT',peak_rows=nrow(peaks),candidate_genes=nrow(genes)),file.path(output,'execution.tsv'))
cat('completed ATAC candidates:',nrow(genes),'genes from',nrow(peaks),'annotated peaks\n')
