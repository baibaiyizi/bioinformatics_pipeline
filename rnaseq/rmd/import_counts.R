# 显式 counts 输入；沿用 DESeq2 size factor 和 counts-derived logCPM。
root <- Sys.getenv('RNASEQ_ROOT_DIR')
input <- Sys.getenv('RNASEQ_COUNTS_FILE')
if (!file.exists(input)) stop('MISSING RNASEQ_COUNTS_FILE: ', input)
if (!requireNamespace('DESeq2', quietly=TRUE)) stop('MISSING R package: DESeq2')
x <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
if (ncol(x)<5L || anyDuplicated(x[[1]]) || anyNA(x[[1]]) || any(!nzchar(x[[1]]))) stop('counts 需要唯一 gene_id 和至少 4 个样本列')
if (anyDuplicated(names(x))) stop('counts 存在重复列名')
samples <- read.delim(file.path(root,'4.DE_analysis/samples.txt'), header=FALSE, col.names=c('group','sample'), stringsAsFactors=FALSE)
if (!setequal(names(x)[-1], samples$sample)) stop('counts 列必须与全部合格样本完全匹配')
mat <- as.matrix(x[,samples$sample,drop=FALSE]); storage.mode(mat) <- 'numeric'; rownames(mat) <- x[[1]]
if (anyNA(mat) || any(!is.finite(mat)) || any(mat<0) || any(mat!=round(mat))) stop('counts 必须是非负整数；不能输入 TPM/FPKM/归一化值')
if (any(colSums(mat)==0)) stop('counts 含空文库')
cd <- data.frame(row.names=colnames(mat), group=samples$group)
dds <- DESeq2::DESeqDataSetFromMatrix(mat, cd, ~1)
dds <- DESeq2::estimateSizeFactors(dds)
norm <- DESeq2::counts(dds,normalized=TRUE)
logcpm <- log2(sweep(mat,2,colSums(mat),'/')*1e6+1)
for (directory in c('3.Merge_result','result_pca/full_sample_inputs')) {
 out <- file.path(root,directory); dir.create(out,recursive=TRUE,showWarnings=FALSE)
 write.table(mat,file.path(out,'genes.counts.matrix'),sep='\t',quote=FALSE,col.names=NA)
 write.table(norm,file.path(out,'genes.DESeq2.normalized_counts.matrix'),sep='\t',quote=FALSE,col.names=NA)
 write.table(logcpm,file.path(out,'genes.counts.logCPM.matrix'),sep='\t',quote=FALSE,col.names=NA)
}
