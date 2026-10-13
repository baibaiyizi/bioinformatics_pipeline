#!/usr/bin/env Rscript
# 从本组学 DMR 注释重建纯甲基化汇总，不读旧联合缓存或表达表。
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1L) stop('用法：Rscript run_gene_methylation.R 分析模块目录')
root <- normalizePath(args[[1]], mustWork = TRUE)
input <- file.path(root, 'result/03_dmr/02_annotation/tables/dmr_annotated.csv')
setup <- file.path(root, 'rmd/01_setup.Rmd')
output <- file.path(root, 'result/08_gene_methylation')
stopifnot(file.exists(input), file.exists(setup))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(stringr);library(ggplot2)})
source(file.path(root,'scripts/figure_style.R'),local=globalenv())
# 仅载入三个纯函数，不执行 setup 的数据恢复和大对象处理。
purl <- tempfile(fileext = '.R')
knitr::purl(setup, output = purl, quiet = TRUE)
expressions <- parse(purl)
unlink(purl)
needed <- c('normalize_feature_label', 'derive_feature_layer', 'build_gene_methylation_summary')
loaded <- character()
for(e in expressions) if(is.call(e) && identical(e[[1]], as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% needed) {
 eval(e, envir = globalenv()); loaded <- c(loaded, as.character(e[[2]]))
}
if(!setequal(needed, loaded)) stop('纯甲基化汇总函数不完整')
dmr <- read.csv(input, stringsAsFactors = FALSE)
if(!all(c('group1_label','group2_label') %in% names(dmr))) stop('DMR 缺少比较方向')
contrast <- unique(paste(dmr$group2_label, dmr$group1_label, sep = ' - '))
if(length(contrast) != 1L) stop('DMR 存在混合比较方向')
long <- build_gene_methylation_summary(dmr)
wide <- long %>% select(gene_id_clean, SYMBOL, Layer, diffMethy, DMR_id, RegionCount, DirectionConflict) %>%
 pivot_wider(names_from = Layer, values_from = c(diffMethy, DMR_id, RegionCount, DirectionConflict))
tables <- file.path(output,'01_group_comparison/tables');plots <- file.path(output,'01_group_comparison/plots')
dir.create(tables,recursive = TRUE,showWarnings = FALSE);dir.create(plots,recursive = TRUE,showWarnings = FALSE)
write_csv(long,file.path(tables,'gene_methylation_long.csv'))
write_csv(wide,file.path(tables,'gene_methylation_wide.csv'))
p <- ggplot(long,aes(Layer,diffMethy)) + geom_boxplot(outlier.shape = NA) +
 geom_hline(yintercept = 0,linetype='dashed') + pfos_theme() +
 labs(x='Genomic layer',y=paste('Methylation proportion difference:',contrast))
pfos_save(p,file.path(plots,'gene_methylation_layers.pdf'),height_mm=137.25)
manifest <- data.frame(Input='result/03_dmr/02_annotation/tables/dmr_annotated.csv',
 SHA256=digest::digest(file=input,algo='sha256'),SourceStatus='PRESERVED_DMR_OUTPUT',
 Contrast=contrast,Rows=nrow(dmr),SourceSetupSHA256=digest::digest(file=setup,algo='sha256'))
write_tsv(manifest,file.path(output,'input_manifest.tsv'))
write_tsv(data.frame(time=format(Sys.time(),tz='UTC',usetz=TRUE),stage='pure_gene_methylation',status='completed',
 source_status='PRESERVED_DMR_OUTPUT',dmr_rows=nrow(dmr),gene_layer_rows=nrow(long)),file.path(output,'execution.tsv'))
cat('completed',basename(root),'DMRs',nrow(dmr),'gene layers',nrow(long),'\n')
