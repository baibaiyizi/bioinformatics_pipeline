#!/usr/bin/env Rscript

# Mouse bulk RNA-seq TF activity template using a locally installed DoRothEA
# A/B/C prior and decoupleR ULM.  It consumes a completed full-cohort matrix.

parse_args <- function(x) {
  defaults <- list(minsize = "5")
  out <- defaults
  i <- 1L
  valid <- c("--de", "--expression", "--samples", "--output", "--reference-group",
             "--comparison-group", "--species", "--assembly", "--minsize")
  while (i <= length(x)) {
    key <- x[[i]]
    if (!key %in% valid || i == length(x)) stop("未知或不完整参数: ", key)
    out[[sub("^--", "", key)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  required <- c("de", "expression", "samples", "output", "reference-group",
                "comparison-group", "species", "assembly")
  missing <- required[!required %in% names(out)]
  if (length(missing)) stop("缺少参数: ", paste(missing, collapse = ", "))
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!identical(tolower(args$species), "mouse")) {
  stop("当前模板仅资格验证 mouse DoRothEA prior；其他物种必须更换并冻结对应prior")
}
reference_group <- args[["reference-group"]]
comparison_group <- args[["comparison-group"]]
if (identical(reference_group, comparison_group)) stop("reference与comparison组必须不同")
minsize <- as.integer(args$minsize)
if (!is.finite(minsize) || minsize < 2L) stop("--minsize必须为>=2的整数")

required_packages <- c("decoupleR", "dorothea", "limma", "jsonlite", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("缺少已声明依赖: ", paste(missing_packages, collapse = ", "))

de_path <- normalizePath(args$de, mustWork = TRUE)
expression_path <- normalizePath(args$expression, mustWork = TRUE)
samples_path <- normalizePath(args$samples, mustWork = TRUE)
out_dir <- normalizePath(args$output, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sha256 <- function(path) as.character(digest::digest(file = path, algo = "sha256", serialize = FALSE))
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
clean_id <- function(x) sub("\\.[0-9]+$", "", as.character(x))

de <- read.csv(de_path, stringsAsFactors = FALSE, check.names = FALSE)
expression <- read.delim(expression_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.delim(samples_path, stringsAsFactors = FALSE, check.names = FALSE)
required_de <- c("gene_id", "SYMBOL", "baseMean", "stat", "log2FoldChange", "pvalue", "padj")
required_samples <- c("sample", "group")
if (length(setdiff(required_de, names(de)))) stop("DE表缺少字段: ", paste(setdiff(required_de, names(de)), collapse = ", "))
if (length(setdiff(required_samples, names(samples)))) stop("样本表缺少字段: ", paste(setdiff(required_samples, names(samples)), collapse = ", "))
if (anyDuplicated(names(expression))) stop("表达矩阵列名必须唯一")
if (anyNA(samples$sample) || any(samples$sample == "") || anyDuplicated(samples$sample)) {
  stop("samples$sample必须非空且唯一")
}
if (anyNA(samples$group) || any(samples$group == "")) stop("samples$group不得缺失")
if (!all(samples$sample %in% names(expression))) stop("表达矩阵缺少样本表中的列")
if (!setequal(unique(samples$group), c(reference_group, comparison_group))) {
  stop("样本表必须且只能包含指定reference/comparison两组")
}
group_counts <- table(samples$group)
if (any(group_counts < 2L)) stop("每组至少需要2个文库；独立生物学单位仍需由项目元数据确认")
for (column in c("biological_unit", "sire_id", "litter_id")) {
  if (!column %in% names(samples)) samples[[column]] <- "unknown"
}
for (column in c("baseMean", "stat", "log2FoldChange", "pvalue", "padj")) {
  de[[column]] <- suppressWarnings(as.numeric(de[[column]]))
}
if (all(!is.finite(de$baseMean)) || all(!is.finite(de$stat))) stop("DE表的baseMean/stat没有有限数值")

sample_order <- as.character(samples$sample)
numeric_expression <- as.matrix(expression[, sample_order, drop = FALSE])
storage.mode(numeric_expression) <- "double"
if (any(!is.finite(numeric_expression)) || any(numeric_expression < 0)) {
  stop("expression必须是非负、有限的归一化计数；本模板内部使用log2(x+1)")
}
de$gene_id_clean <- clean_id(de$gene_id)
expression$gene_id_clean <- clean_id(expression[[1L]])
de$SYMBOL <- trimws(as.character(de$SYMBOL))
map <- de[!is.na(de$SYMBOL) & de$SYMBOL != "", c("gene_id_clean", "SYMBOL", "baseMean")]
map <- map[order(map$gene_id_clean, -map$baseMean, map$SYMBOL), ]
map <- map[!duplicated(map$gene_id_clean), ]
joined <- merge(expression, map, by = "gene_id_clean", all = FALSE, sort = FALSE)
joined$row_mean <- rowMeans(joined[, sample_order, drop = FALSE], na.rm = TRUE)
joined <- joined[order(joined$SYMBOL, -joined$row_mean, joined$gene_id_clean), ]
joined <- joined[!duplicated(joined$SYMBOL), ]
mat <- as.matrix(joined[, sample_order, drop = FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- joined$SYMBOL
mat <- log2(mat + 1)
if (nrow(mat) < 100L) stop("映射后的symbol表达矩阵少于100行")

reg_env <- new.env(parent = emptyenv())
utils::data(list = "dorothea_mm", package = "dorothea", envir = reg_env)
if (!exists("dorothea_mm", envir = reg_env, inherits = FALSE)) stop("本地dorothea_mm对象不可用")
network_full <- as.data.frame(get("dorothea_mm", envir = reg_env), stringsAsFactors = FALSE)
network_full <- network_full[network_full$confidence %in% c("A", "B", "C"), c("tf", "confidence", "target", "mor")]
network_full$tf <- as.character(network_full$tf)
network_full$target <- as.character(network_full$target)
network_full$mor <- as.numeric(network_full$mor)
network_full <- unique(network_full[!is.na(network_full$tf) & network_full$tf != "" &
                                      !is.na(network_full$target) & network_full$target != "" &
                                      is.finite(network_full$mor), ])
network_tested <- network_full[network_full$target %in% rownames(mat), ]
if (!nrow(network_tested)) stop("DoRothEA A/B/C网络与表达矩阵没有交集")

activity <- as.data.frame(decoupleR::run_ulm(
  mat = mat, network = network_tested, .source = "tf", .target = "target",
  .mor = "mor", minsize = minsize
), stringsAsFactors = FALSE)
pick_column <- function(candidates, available, role) {
  found <- intersect(candidates, available)
  if (!length(found)) stop("decoupleR输出缺少", role, "字段")
  found[[1L]]
}
source_col <- pick_column(c("source", ".source", "tf"), names(activity), "TF")
condition_col <- pick_column(c("condition", "sample", ".condition"), names(activity), "sample")
score_col <- pick_column(c("score", "statistic"), names(activity), "score")
p_col <- intersect(c("p_value", "pvalue"), names(activity))
activity_long <- data.frame(
  tf = as.character(activity[[source_col]]), sample = as.character(activity[[condition_col]]),
  activity_score = as.numeric(activity[[score_col]]),
  ulm_target_association_pvalue = if (length(p_col)) as.numeric(activity[[p_col[[1L]]]]) else NA_real_,
  stringsAsFactors = FALSE
)
activity_long <- merge(activity_long, samples[, c("sample", "group", "biological_unit", "sire_id", "litter_id")],
                       by = "sample", all.x = TRUE, sort = FALSE)
activity_long$method <- paste0("decoupleR ULM; DoRothEA mouse A/B/C; minsize=", minsize)
activity_long$inference_scope <- "RNA_library_level; biological_unit_and_family_from_input_metadata"

activity_matrix <- reshape(activity_long[, c("tf", "sample", "activity_score")],
                           idvar = "tf", timevar = "sample", direction = "wide")
rownames(activity_matrix) <- activity_matrix$tf
activity_matrix$tf <- NULL
names(activity_matrix) <- sub("^activity_score\\.", "", names(activity_matrix))
activity_matrix <- as.matrix(activity_matrix[, sample_order, drop = FALSE])
storage.mode(activity_matrix) <- "double"
group <- factor(samples$group, levels = c(reference_group, comparison_group))
design <- stats::model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast_name <- paste0(comparison_group, " - ", reference_group)
contrast_vector <- rep(0, ncol(design)); names(contrast_vector) <- colnames(design)
contrast_vector[[comparison_group]] <- 1; contrast_vector[[reference_group]] <- -1
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(activity_matrix, design), contrast_vector))
differential <- limma::topTable(fit, number = Inf, sort.by = "none", adjust.method = "BH")
differential$tf <- rownames(differential); rownames(differential) <- NULL
names(differential)[names(differential) == "P.Value"] <- "pvalue"
names(differential)[names(differential) == "adj.P.Val"] <- "FDR"
full_counts <- aggregate(target ~ tf, data = network_full, FUN = function(x) length(unique(x)))
names(full_counts)[[2L]] <- "prior_target_n"
tested_counts <- aggregate(target ~ tf, data = network_tested, FUN = function(x) length(unique(x)))
names(tested_counts)[[2L]] <- "tested_target_n"
differential <- merge(merge(differential, full_counts, by = "tf", all.x = TRUE, sort = FALSE),
                      tested_counts, by = "tf", all.x = TRUE, sort = FALSE)
differential$contrast <- contrast_name
differential$raw_p_supported <- is.finite(differential$pvalue) & differential$pvalue < 0.05
differential$bh_supported <- is.finite(differential$FDR) & differential$FDR < 0.05
differential$experimental_unit <- "RNA_library; use input metadata to determine independent biological unit"
differential <- differential[order(differential$pvalue, -abs(differential$logFC), differential$tf), ]

de_target <- de[!is.na(de$SYMBOL) & de$SYMBOL != "", required_de]
de_target <- de_target[order(de_target$SYMBOL, -abs(de_target$stat), de_target$gene_id), ]
de_target <- de_target[!duplicated(de_target$SYMBOL), ]
names(de_target)[names(de_target) == "SYMBOL"] <- "target"
network_effects <- merge(network_full, de_target, by = "target", all.x = TRUE, sort = FALSE)
network_effects$target_tested_in_rna <- network_effects$target %in% rownames(mat)
network_effects$prior_version <- paste0("dorothea_", as.character(utils::packageVersion("dorothea")), "_mouse_ABC")
network_effects$evidence_role <- "signed prior edge plus target RNA effect; not TF binding or causal regulation"

output_files <- c(
  regulon_sample_activity = file.path(out_dir, "regulon_sample_activity.tsv"),
  regulon_group_effects = file.path(out_dir, "regulon_group_effects.tsv"),
  regulon_target_effects = file.path(out_dir, "regulon_target_effects.tsv"),
  dorothea_mouse_ABC_network = file.path(out_dir, "dorothea_mouse_ABC_network.tsv")
)
write_tsv(activity_long, output_files[["regulon_sample_activity"]])
write_tsv(differential, output_files[["regulon_group_effects"]])
write_tsv(network_effects, output_files[["regulon_target_effects"]])
write_tsv(network_full, output_files[["dorothea_mouse_ABC_network"]])

prior_files <- file.path(find.package("dorothea"), "data", c("Rdata.rdb", "Rdata.rdx", "Rdata.rds"))
prior_files <- prior_files[file.exists(prior_files)]
input_files <- c(DE = de_path, normalized_expression = expression_path, samples = samples_path,
                 setNames(prior_files, paste0("dorothea_serialized_data_", seq_along(prior_files))))
input_manifest <- data.frame(
  role = names(input_files), path = unname(input_files),
  sha256 = vapply(unname(input_files), sha256, character(1)), stringsAsFactors = FALSE
)
input_manifest_path <- file.path(out_dir, "input_manifest.tsv")
write_tsv(input_manifest, input_manifest_path)
output_files <- c(output_files, input_manifest = input_manifest_path)
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][[1L]]), mustWork = TRUE)
receipt <- list(
  status = "completed", generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  species = "mouse", assembly = args$assembly, contrast = contrast_name,
  group_counts = as.list(group_counts), biological_unit = unique(as.character(samples$biological_unit)),
  method = paste0("decoupleR::run_ulm, minsize=", minsize, ", log2(normalized_count+1)"),
  prior = "local installed dorothea_mm A/B/C; no network update",
  input_contract = list(
    de_columns = required_de, expression = "nonnegative normalized counts, first column gene_id, remaining columns sample IDs",
    samples_columns = c("sample", "group", "optional biological_unit/sire_id/litter_id")
  ),
  package_versions = as.list(vapply(required_packages, function(p) as.character(utils::packageVersion(p)), character(1))),
  input_sha256 = as.list(setNames(input_manifest$sha256, input_manifest$path)),
  code_sha256 = sha256(script_path),
  output_sha256 = as.list(vapply(output_files, sha256, character(1))),
  outputs = list(sample_activity_rows = nrow(activity_long), tf_effect_rows = nrow(differential),
                 prior_edges = nrow(network_full), tested_edges = nrow(network_tested)),
  claim_boundary = "TF activity is inferred from bulk RNA and a signed prior; it is not direct TF binding or causal regulation."
)
jsonlite::write_json(receipt, file.path(out_dir, "receipt.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
