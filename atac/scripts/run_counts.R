#!/usr/bin/env Rscript
# 显式 consensus fragment-count 输入模式；原区间 DESeq2，与原生 DiffBind 分开记录。
args <- commandArgs(trailingOnly = TRUE)
check_only <- '--check' %in% args
root <- normalizePath(Sys.getenv('ATAC_ROOT_DIR'), mustWork = TRUE)
input <- Sys.getenv('ATAC_COUNTS_FILE')
if(!nzchar(input) || !file.exists(input)) stop('MISSING: ATAC_COUNTS_FILE 必须指向原始整数片段计数 TSV')
meta <- read.delim(Sys.getenv('ATAC_SAMPLES'), check.names = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
meta[] <- lapply(meta, trimws)
reference <- Sys.getenv('ATAC_REFERENCE_GROUP'); comparison <- Sys.getenv('ATAC_COMPARISON_GROUP')
meta <- meta[order(match(meta$group, c(reference,comparison)),meta$sample), , drop=FALSE]
raw <- read.delim(input,check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(peak_id='character',CHR='character'))
coord <- c('peak_id','CHR','START','END')
if(!all(coord %in% names(raw))) stop('计数表必须含 peak_id, CHR, START, END（0-based half-open）')
if(anyDuplicated(names(raw))) stop('计数表列名重复')
count_cols <- setdiff(names(raw),c(coord,'PeakWidth'))
if(!setequal(count_cols,meta$sample)) stop('计数列必须与样本表身份一一对应；禁止按列位置推断')
if(!nrow(raw) || anyNA(raw[,coord]) || anyDuplicated(raw$peak_id) || any(!nzchar(raw$peak_id)) || any(!nzchar(raw$CHR))) stop('peak 身份或坐标缺失/重复')
if(any(!is.finite(raw$START)) || any(!is.finite(raw$END)) || any(raw$START<0) || any(raw$END<=raw$START) || any(raw$START != floor(raw$START)) || any(raw$END != floor(raw$END))) stop('非法 BED 坐标')
if(anyDuplicated(paste(raw$CHR,raw$START,raw$END,sep=':'))) stop('重复检验区间')
if('PeakWidth' %in% names(raw) && anyNA(raw$PeakWidth)) stop('PeakWidth 不能缺失')
if('PeakWidth' %in% names(raw) && any(raw$PeakWidth != raw$END-raw$START)) stop('PeakWidth 与 BED 区间不符')
counts <- as.matrix(raw[,meta$sample,drop=FALSE]); storage.mode(counts) <- 'numeric'
if(any(!is.finite(counts)) || any(counts<0) || any(counts != floor(counts)) || any(counts > .Machine$integer.max)) stop('必须为有限、非负、可用 R integer 表示的原始整数片段计数；CPM/TPM 不可用于拟合')
if(any(colSums(counts)==0)) stop('至少一个样本总计数为零')
min_total <- as.numeric(Sys.getenv('ATAC_COUNTS_MIN_TOTAL','10'))
if(length(min_total)!=1 || !is.finite(min_total) || min_total<0) stop('ATAC_COUNTS_MIN_TOTAL 非法')
keep <- rowSums(counts)>=min_total
if(sum(keep)<2L) stop('过滤后不足两个可检验 peak')
for(key in c('ATAC_PRIMARY_P_CUTOFF','ATAC_DIFF_FDR_CUTOFF')) {
  value <- suppressWarnings(as.numeric(Sys.getenv(key)))
  if(length(value)!=1L || !is.finite(value) || value<=0 || value>1) stop(key, ' 必须为 (0,1] 的有限数值')
}
fc_threshold <- suppressWarnings(as.numeric(Sys.getenv('ATAC_DIFF_LOG2FC_CUTOFF')))
if(length(fc_threshold)!=1L || !is.finite(fc_threshold) || fc_threshold<0) stop('ATAC_DIFF_LOG2FC_CUTOFF 必须非负有限')
if(!Sys.getenv('ATAC_SIGNIFICANCE_MODE') %in% c('pvalue','padj','fdr')) stop('ATAC_SIGNIFICANCE_MODE 非法')
cat('COUNTS_INPUT_PASS:',sum(keep),'/',nrow(raw),'peaks;',comparison,'-',reference,'; integer fragments\n')
if(check_only) quit(status=0)
for(pkg in c('DESeq2','digest')) if(!requireNamespace(pkg,quietly=TRUE)) stop('MISSING R package: ',pkg)
source(file.path(root,'rmd','significance_helpers.R'))
tables <- file.path(root,'result','05_differential_accessibility','tables');plots<-file.path(root,'result','05_differential_accessibility','plots')
dir.create(tables,recursive=TRUE,showWarnings=FALSE);dir.create(plots,recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,name) write.csv(x,file.path(tables,name),row.names=FALSE,na='')
counts <- counts[keep,,drop=FALSE];raw <- raw[keep,,drop=FALSE];storage.mode(counts)<-'integer';rownames(counts)<-raw$peak_id
rownames(meta)<-meta$sample;meta$group<-factor(meta$group,levels=c(reference,comparison))
design<-model.matrix(~group,meta)
if(qr(design)$rank != ncol(design)) stop('设计矩阵不满秩')
set.seed(as.integer(Sys.getenv('ATAC_COUNTS_SEED','20260929')))
dds <- DESeq2::DESeqDataSetFromMatrix(counts,meta,design=~group)
dds <- DESeq2::DESeq(dds,quiet=TRUE,parallel=FALSE)
fdr_cutoff<-as.numeric(Sys.getenv('ATAC_DIFF_FDR_CUTOFF','0.05'))
res<-as.data.frame(DESeq2::results(dds,contrast=c('group',comparison,reference),alpha=fdr_cutoff))
res<-res[raw$peak_id,,drop=FALSE]
out<-data.frame(peak_id=raw$peak_id,seqnames=raw$CHR,start=raw$START+1L,end=raw$END,
                baseMean=res$baseMean,Fold=res$log2FoldChange,lfcSE=res$lfcSE,stat=res$stat,
                pvalue=res$pvalue,padj=res$padj,FDR=res$padj,
                ReferenceGroup=reference,ComparisonGroup=comparison,stringsAsFactors=FALSE)
flags<-compute_atac_diff_peak_flags(out$pvalue,out$FDR,out$Fold,
 significance_mode=Sys.getenv('ATAC_SIGNIFICANCE_MODE','pvalue'),
 primary_p_cutoff=as.numeric(Sys.getenv('ATAC_PRIMARY_P_CUTOFF','0.05')),
 fdr_cutoff=fdr_cutoff,log2fc_cutoff=as.numeric(Sys.getenv('ATAC_DIFF_LOG2FC_CUTOFF','0.585')),
 primary_require_fc=tolower(Sys.getenv('ATAC_PRIMARY_REQUIRE_FC','0')) %in% c('1','true','yes'))
out<-cbind(out,flags)
write_csv(out,'Differential_Peaks_AllSites.csv')
write_csv(out[out$PrimaryDiffPeak,,drop=FALSE],'Differential_Peaks.csv')
write_csv(out[out$PvalueDiffPeak,,drop=FALSE],'Differential_Peaks_Pvalue.csv')
write_csv(out[out$FDRDiffPeak,,drop=FALSE],'Differential_Peaks_FDR.csv')
normalized<-DESeq2::counts(dds,normalized=TRUE)
write_csv(data.frame(peak_id=rownames(normalized),normalized,check.names=FALSE),'Normalized_Consensus_Counts.csv')
write_csv(data.frame(sample=meta$sample,biological_sample_id=meta$biological_sample_id,group=meta$group,
                    size_factor=DESeq2::sizeFactors(dds),total_fragments=colSums(counts)), 'Sample_Count_QC.csv')
write_csv(data.frame(AnalysisSource='consensus_counts_DESeq2',AnalysisMethod='DESeq2 independent-two-group Wald test; original consensus intervals',
 InputFile=normalizePath(input),InputSHA256=digest::digest(file=input,algo='sha256'),
 InputCoordinateSystem='BED 0-based half-open',OutputCoordinateSystem='1-based inclusive',
 ReferenceGroup=reference,ComparisonGroup=comparison,Unit='independent biological sample; technical libraries already aggregated',
 NReference=sum(meta$group==reference),NComparison=sum(meta$group==comparison),MinTotalCount=min_total,
 DESeq2Version=as.character(utils::packageVersion('DESeq2')),EvidenceRole='exploratory raw-p screen with BH/FDR retained',
 UpstreamExecution='not_run_in_counts_mode',DiffBindExecution='not_run_in_counts_mode'), 'Differential_Peak_Source_Audit.csv')
# PCA 使用全部保留样本，仅诊断；没有基于分离程度的选样。
logcounts<-log2(normalized+1);vary<-apply(logcounts,1,var)>0
if(sum(vary)>=2L) {
 pc<-prcomp(t(logcounts[vary,,drop=FALSE]),center=TRUE,scale.=FALSE)
 scores<-data.frame(sample=rownames(pc$x),group=meta[rownames(pc$x),'group'],pc$x[,seq_len(min(2,ncol(pc$x))),drop=FALSE])
 write_csv(scores,'Sample_PCA_Scores.csv')
}
render_pair<-function(stem,draw) {
 grDevices::cairo_pdf(file.path(plots,paste0(stem,'.pdf')),width=183/25.4,height=130/25.4,family='sans');draw();grDevices::dev.off()
 grDevices::png(file.path(plots,paste0(stem,'.png')),width=183,height=130,units='mm',res=150,type='cairo');draw();grDevices::dev.off()
}
if(exists('scores') && 'PC2' %in% names(scores)) render_pair('Sample_PCA',function(){
 plot(scores$PC1,scores$PC2,col=ifelse(scores$group==comparison,'#D95F59','#3974A5'),pch=19,xlab='PC1',ylab='PC2',main='All eligible biological samples');text(scores$PC1,scores$PC2,labels=scores$sample,pos=3,cex=.7)
})
valid<-is.finite(out$Fold)&is.finite(out$SelectedSignificanceValue)
if(any(valid)) render_pair('Volcano_Plot',function(){
 plot(out$Fold[valid],-log10(pmax(out$SelectedSignificanceValue[valid],.Machine$double.xmin)),
  col=ifelse(out$PrimaryDiffPeak[valid],'#D95F59','#B8B8B8'),pch=19,cex=.6,
  xlab=paste('log2 fold change:',comparison,'-',reference),ylab=paste('-log10',unique(out$SelectedSignificanceMetric)),main='Consensus peak DESeq2; exploratory screen')
})
writeLines(capture.output(sessionInfo()),file.path(tables,'sessionInfo.txt'))
cat('COUNTS_DESEQ2_COMPLETED:',nrow(out),'tested peaks;',sum(out$PrimaryDiffPeak),'selected; all eligible samples retained\n')
