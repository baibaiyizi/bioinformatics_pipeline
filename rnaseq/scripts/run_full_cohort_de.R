#!/usr/bin/env Rscript
# 显式样本表驱动的完整 cohort DE；本入口不读取 PCA 选择名单。
args <- commandArgs(TRUE)
opts <- list(p = "0.05", lfc = "0.58")
if (length(args) %% 2L) stop("参数必须成对: --root MODULE --counts FILE --samples FILE --annotation FILE --out DIR")
for (i in seq.int(1L, length(args), by = 2L)) opts[[sub("^--", "", args[[i]])]] <- args[[i + 1L]]
for (key in c("root", "counts", "samples", "annotation", "out", "reference", "comparison")) if (is.null(opts[[key]])) stop("缺少 --", key)
for (pkg in c("DESeq2", "digest", "jsonlite")) if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少包: ", pkg)
root <- normalizePath(opts$root, mustWork = TRUE)
resolve <- function(p) normalizePath(if (grepl("^/", p)) p else file.path(root, p), mustWork = FALSE)
paths <- lapply(opts[c("counts", "samples", "annotation", "out")], resolve)
orgdb_mode <- identical(opts$annotation, "org.Mm.eg.db")
if (orgdb_mode) {
  for (pkg in c("org.Mm.eg.db", "AnnotationDbi")) if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少包: ", pkg)
  paths$annotation <- AnnotationDbi::dbfile(org.Mm.eg.db::org.Mm.eg.db)
}
out <- paths$out
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  if (!file.exists(file.path(out, "execution.tsv"))) stop("目标非空且不属于完整cohort producer，拒绝覆盖")
  archived <- file.path(root, "history", paste0("full_cohort_", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
  dir.create(dirname(archived), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(out, archived)) stop("不能原位保存前次producer产物")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
script <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg)), mustWork = TRUE)
producer_snapshot <- file.path(out, "producer_source.R")
if (!file.copy(script, producer_snapshot)) stop("无法冻结实际执行源码")
status_file <- file.path(out, "execution.tsv")
record <- function(status, message) write.table(data.frame(time = format(Sys.time(), tz = "UTC", usetz = TRUE), status = status, message = message), status_file, sep = "\t", quote = FALSE, row.names = FALSE)
record("STARTED", "完整cohort模型；不是FASTQ重分析")
tryCatch({
  if (opts$reference == opts$comparison) stop("reference 与 comparison 必须不同")
  for (p in paths[c("counts", "samples", "annotation")]) if (!file.exists(p)) stop("缺失输入: ", p)
  counts <- read.delim(paths$counts, row.names = 1L, check.names = FALSE)
  metadata <- read.delim(paths$samples, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample", "group") %in% names(metadata)) || anyDuplicated(metadata$sample)) stop("样本表需唯一 sample、显式 group")
  if (anyDuplicated(rownames(counts)) || anyDuplicated(names(counts))) stop("count 基因或样本身份重复")
  if (!setequal(metadata$sample, names(counts))) stop("完整count和样本表集合不一致；不得静默丢样本")
  metadata <- metadata[match(names(counts), metadata$sample), , drop = FALSE]
  if (anyNA(metadata$group) || !setequal(unique(metadata$group), c(opts$reference, opts$comparison))) stop("组别不完整或含其他组")
  if (any(table(metadata$group) < 2L)) stop("每组至少需要两个文库；生物学单位另行核验")
  m <- as.matrix(counts)
  if (any(!is.finite(m)) || any(m < 0)) stop("counts必须非负且有限")
  n_fractional <- sum(m != round(m))
  m <- round(m)
  condition <- factor(metadata$group, levels = c(opts$reference, opts$comparison))
  coldata <- data.frame(condition = condition, row.names = metadata$sample)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = coldata, design = ~ condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  result <- DESeq2::results(dds, contrast = c("condition", opts$comparison, opts$reference), independentFiltering = FALSE)
  norm <- DESeq2::counts(dds, normalized = TRUE)
  de <- data.frame(gene_id = rownames(result), sampleA = opts$comparison, sampleB = opts$reference,
                   baseMeanA = rowMeans(norm[, condition == opts$comparison, drop = FALSE]),
                   baseMeanB = rowMeans(norm[, condition == opts$reference, drop = FALSE]),
                   as.data.frame(result), row.names = NULL, check.names = FALSE)
  de$ensembl_id <- sub("\\.[0-9]+$", "", de$gene_id)
  if (orgdb_mode) {
    valid_keys <- intersect(de$ensembl_id, AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "ENSEMBL"))
    annotation <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = valid_keys, keytype = "ENSEMBL", columns = c("SYMBOL", "GENENAME", "ENTREZID"))
    names(annotation)[names(annotation) == "ENSEMBL"] <- "ensembl_id"
  } else annotation <- read.csv(paths$annotation, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"ensembl_id" %in% names(annotation)) stop("注释来源缺少ensembl_id")
  for (key in c("SYMBOL", "GENENAME", "ENTREZID")) {
    if (!key %in% names(annotation)) stop("注释缺少: ", key)
    values <- tapply(as.character(annotation[[key]]), annotation$ensembl_id, function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";"))
    de[[key]] <- unname(values[de$ensembl_id])
  }
  de$annotation_ambiguous <- grepl(";", de$SYMBOL) | grepl(";", de$ENTREZID)
  de$gene_label <- ifelse(is.na(de$SYMBOL) | !nzchar(de$SYMBOL), de$ensembl_id, de$SYMBOL)
  de$direction <- ifelse(!is.na(de$pvalue) & de$pvalue < as.numeric(opts$p) & abs(de$log2FoldChange) >= as.numeric(opts$lfc), ifelse(de$log2FoldChange > 0, "up", "down"), "ns")
  de$significance_value <- de$pvalue
  write.csv(de, file.path(out, "DE_results.csv"), row.names = FALSE, na = "")
  write.table(data.frame(gene_id = rownames(m), m, check.names = FALSE), file.path(out, "counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = rownames(norm), norm, check.names = FALSE), file.path(out, "normalized_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (key in c("biological_unit", "pair_id", "sire_id", "litter_id")) if (!key %in% names(metadata)) metadata[[key]] <- "unknown"
  write.table(metadata, file.path(out, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sha <- function(p) digest::digest(file = p, algo = "sha256")
  manifest <- data.frame(role = c("counts", "samples", "annotation", "producer"),
                         path = c(opts$counts, opts$samples, if (orgdb_mode) paths$annotation else opts$annotation, producer_snapshot),
                         sha256 = vapply(c(paths$counts, paths$samples, paths$annotation, producer_snapshot), sha, character(1)))
  write.table(manifest, file.path(out, "input_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  outputs <- c("DE_results.csv", "normalized_expression.tsv", "counts.tsv", "samples.tsv")
  provenance <- list(status = "completed", contrast = paste(opts$comparison, "-", opts$reference),
    samples = metadata$sample, group_counts = as.list(table(metadata$group)),
    design = "~condition", independentFiltering = FALSE, rounded_fractional_values = n_fractional,
    biological_unit_status = "unknown; library-level exploratory model, not independent sire/litter inference",
    annotation_source_status = if (orgdb_mode) paste0("local org.Mm.eg.db ", utils::packageVersion("org.Mm.eg.db")) else "historical annotation only; old statistical columns not consumed",
    primary_metric = "nominal pvalue, exploratory", p_cutoff = as.numeric(opts$p), lfc_cutoff = as.numeric(opts$lfc),
    source = manifest, outputs = as.list(setNames(vapply(file.path(out, outputs), sha, character(1)), outputs)),
    DESeq2_version = as.character(utils::packageVersion("DESeq2")), R_version = R.version.string)
  jsonlite::write_json(provenance, file.path(out, "provenance.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  record("COMPLETED", paste(ncol(m), "文库；", nrow(de), "基因；文库层探索性，父本/窝独立性未知"))
  cat("COMPLETED", root, ncol(m), nrow(de), "\n")
}, error = function(e) { record("FAILED", conditionMessage(e)); stop(e) })
