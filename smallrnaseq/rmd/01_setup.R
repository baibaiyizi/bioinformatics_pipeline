# 这个文件是通用 small RNA 报告体系的公共底座：
# 1. 定义路径、颜色、分类规则和绘图辅助函数
# 2. 统一加载 SPORTS 结果并整理成后续模块可复用的数据对象
# 3. 为 02-08 这些对象中心模块提供共享环境，避免重复写同样的预处理代码
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 11, fig.height = 7, dpi = 300)
options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(knitr)
})

# 轻量级依赖检测：后面某些图依赖可选包，缺失时会优雅降级而不是直接报错。
pkg_has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# 常用的空值兜底操作符，便于在构造标题/注释时写出更紧凑的代码。
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

hellinger_transform <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  if (any(totals <= 0)) stop("Hellinger transform requires at least one non-zero feature per sample.")
  sqrt(x / totals)
}

# 通过当前 Rmd 文件的位置反推出项目根目录。
# 这样无论是直接 knit 还是通过 run_rmd.sh 调用，路径解析都尽量保持一致。
required_env <- function(key) {
  value <- Sys.getenv(key, "")
  if (!nzchar(value)) stop("Missing config variable: ", key, "; use run_rmd.sh to load project.env")
  value
}
ROOT_DIR <- normalizePath(required_env("SMALLRNA_ROOT_DIR"), mustWork = TRUE)
RAW_DIR <- file.path(ROOT_DIR, "inputs")
SPORTS_RESULT_DIR <- required_env("SPORTS_RESULT_DIR")
READS_INFO_FILE <- required_env("SAMPLES_TSV")
SPORTS_METADATA_DIR <- file.path(SPORTS_RESULT_DIR, "02_metadata")
SPORTS_QC_DIR <- file.path(SPORTS_RESULT_DIR, "03_qc")
SPORTS_ANNOTATION_DIR <- file.path(SPORTS_RESULT_DIR, "04_sports_annotation")
SPORTS_SAMPLE_ROOT <- file.path(SPORTS_ANNOTATION_DIR, "samples")
SPORTS_ENDOSIRNA_DIR <- file.path(SPORTS_RESULT_DIR, "05_endosirna")
SPORTS_ENDOSIRNA_SUMMARY_DIR <- file.path(SPORTS_ENDOSIRNA_DIR, "summary")
SPORTS_LOG_DIR <- file.path(SPORTS_RESULT_DIR, "06_logs")
METADATA_FILE <- READS_INFO_FILE
RESULT_DIR <- required_env("SMALLRNA_RESULT_DIR")
INTERNAL_CACHE_DIR <- file.path(RESULT_DIR, ".cache")
PROJECT_SPECIES <- required_env("PROJECT_SPECIES")
PROJECT_HOST_LABEL <- required_env("PROJECT_HOST_LABEL")
PROJECT_STAGE <- required_env("PROJECT_STAGE")
PROJECT_CONTRAST_BASELINE <- required_env("PROJECT_CONTRAST_BASELINE")
PROJECT_CONTRAST_CASE <- required_env("PROJECT_CONTRAST_CASE")
PROJECT_GROUP_LABELS <- setNames(c(PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE), c(PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE))
MIRNA_SPECIES <- PROJECT_SPECIES
MIRNA_PREFIX <- required_env("MIRNA_PREFIX")
MULTIMIR_ORG <- required_env("MULTIMIR_ORG")
KEGG_ORGANISM <- required_env("KEGG_ORGANISM")
MIRNA_TARGET_ORGDB_PACKAGE <- required_env("MIRNA_TARGET_ORGDB_PACKAGE")
MIRNA_TARGET_ORGDB_OBJECT <- required_env("MIRNA_TARGET_ORGDB_OBJECT")
set.seed(42)

SETUP_DIR         <- file.path(RESULT_DIR, "1_setup")
EXPRESSION_MATRIX_DIR <- file.path(SETUP_DIR, "1_expression_matrix")

SINGLE_DIR        <- file.path(RESULT_DIR, "2_single_sample_qc")
SINGLE_FIG        <- file.path(SINGLE_DIR, "plots")
SINGLE_TAB        <- file.path(SINGLE_DIR, "tables")

PCA_DIR           <- file.path(RESULT_DIR, "3_cohort_qc_pca")
PCA_FIG           <- file.path(PCA_DIR, "plots")
PCA_TAB           <- file.path(PCA_DIR, "tables")

SNC_COMPARE_DIR   <- file.path(RESULT_DIR, "4_global_sncRNA")
SNC_COMPARE_FIG   <- file.path(SNC_COMPARE_DIR, "plots")
SNC_COMPARE_TAB   <- file.path(SNC_COMPARE_DIR, "tables")

TSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "5_tsrna")
TSRNA_COMPARE_FIG <- file.path(TSRNA_COMPARE_DIR, "plots")
TSRNA_COMPARE_TAB <- file.path(TSRNA_COMPARE_DIR, "tables")

RSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "6_rsrna")
RSRNA_COMPARE_FIG <- file.path(RSRNA_COMPARE_DIR, "plots")
RSRNA_COMPARE_TAB <- file.path(RSRNA_COMPARE_DIR, "tables")

MIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "7_mirna")
MIRNA_COMPARE_FIG <- file.path(MIRNA_COMPARE_DIR, "plots")
MIRNA_COMPARE_TAB <- file.path(MIRNA_COMPARE_DIR, "tables")

PIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "8_pirna")
PIRNA_COMPARE_FIG <- file.path(PIRNA_COMPARE_DIR, "plots")
PIRNA_COMPARE_TAB <- file.path(PIRNA_COMPARE_DIR, "tables")

ENDOSIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "9_endosirna")
ENDOSIRNA_COMPARE_FIG <- file.path(ENDOSIRNA_COMPARE_DIR, "plots")
ENDOSIRNA_COMPARE_TAB <- file.path(ENDOSIRNA_COMPARE_DIR, "tables")

RUN_FASTQC <- tolower(Sys.getenv("RUN_FASTQC", "auto"))
RUN_FASTQC <- if (RUN_FASTQC %in% c("auto", "1", "0", "true", "false", "yes", "no")) RUN_FASTQC else "auto"
RUN_EXOGENOUS <- tolower(Sys.getenv("RUN_EXOGENOUS", "0"))
RUN_EXOGENOUS <- if (RUN_EXOGENOUS %in% c("1", "true", "yes")) "1" else "0"
EXOGENOUS_DB_ROOT <- Sys.getenv("EXOGENOUS_DB_ROOT", "")
MIRNA_TARGET_SOURCE <- tolower(Sys.getenv("MIRNA_TARGET_SOURCE", "validated"))
MIRNA_TARGET_SOURCE <- if (MIRNA_TARGET_SOURCE %in% c("validated", "predicted", "msigdb_fallback")) MIRNA_TARGET_SOURCE else "validated"
DE_SIGNIFICANCE <- toupper(Sys.getenv("DE_SIGNIFICANCE", "PVALUE"))
DE_SIGNIFICANCE <- if (DE_SIGNIFICANCE %in% c("FDR", "PVALUE", "P_VALUE", "P")) DE_SIGNIFICANCE else stop("DE_SIGNIFICANCE 必须为 PVALUE 或 FDR")
RAW_AUDIT_MAX_READS <- suppressWarnings(as.integer(Sys.getenv("RAW_AUDIT_MAX_READS", "100000")))
if (!is.finite(RAW_AUDIT_MAX_READS) || RAW_AUDIT_MAX_READS <= 0) RAW_AUDIT_MAX_READS <- 100000L
RAW_AUDIT_KMER <- suppressWarnings(as.integer(Sys.getenv("RAW_AUDIT_KMER", "12")))
if (!is.finite(RAW_AUDIT_KMER) || RAW_AUDIT_KMER < 4) RAW_AUDIT_KMER <- 12L
RAW_AUDIT_REFRESH <- tolower(Sys.getenv("RAW_AUDIT_REFRESH", "0"))
RAW_AUDIT_REFRESH <- if (RAW_AUDIT_REFRESH %in% c("1", "true", "yes")) "1" else "0"

# 当前报告默认只生成 sense 视图。
# 如果后面想扩成 merged / antisense，多数辅助函数已经支持，只需要放开 VIEW_LEVELS。
REPORT_VIEW <- "sense"
ALL_VIEW_LEVELS <- c("merged", "sense", "antisense")
VIEW_LEVELS <- c(REPORT_VIEW)

# 这两个变量会在具体模块里被覆盖，用于 save_plot_* 系列函数自动决定输出目录。
OUTPUT_FIG_DIR <- NULL
OUTPUT_TAB_DIR <- NULL

# 通用的 TSV 读取函数：
# - 优先用 data.table::fread 提速
# - 兼容只取部分列
# - 输出统一转成 tibble，便于后续 tidyverse 管道处理
read_tsv_flex <- function(file, select = NULL) {
  if (!file.exists(file)) return(tibble())
  if (pkg_has("data.table")) {
    data.table::fread(file, sep="\t", select=select, data.table=FALSE, showProgress=FALSE)
  } else {
    out <- read.delim(file, sep="\t", header=TRUE, check.names=FALSE, quote="", comment.char="", stringsAsFactors=FALSE)
    if (!is.null(select)) out <- out[, select, drop=FALSE]
    out
  } |> as_tibble()
}

# SPORTS 的某些中间比对文件没有表头，而且列数可能不稳定。
# 这里按给定的列名模板读取，并在列不够时自动补 NA，保证后续代码能稳定引用字段。
read_alignment_tsv_flex <- function(file, col_names) {
  if (!file.exists(file)) return(tibble())
  if (pkg_has("data.table")) {
    out <- data.table::fread(file, sep="\t", header=FALSE, data.table=FALSE, showProgress=FALSE, quote="", fill=TRUE)
  } else {
    out <- read.delim(file, sep="\t", header=FALSE, check.names=FALSE, quote="", comment.char="", fill=TRUE, stringsAsFactors=FALSE)
  }
  if (nrow(out) == 0) return(tibble())
  if (ncol(out) < length(col_names)) for (i in seq.int(ncol(out)+1L, length(col_names))) out[[i]] <- NA_character_
  out <- out[, seq_along(col_names), drop=FALSE]; names(out) <- col_names; as_tibble(out)
}

# 所有表格输出统一走这里，避免每个模块重复写目录创建逻辑。
write_csv_safe <- function(df, file_name, tab_dir = OUTPUT_TAB_DIR) {
  if (is.null(tab_dir) || !nzchar(tab_dir)) stop("tab_dir is not set for write_csv_safe().")
  target_file <- file.path(tab_dir, file_name)
  dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, target_file, row.names = FALSE)
}

write_tsv_safe <- function(df, file_path) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# 规范化视图参数，避免 save_plot_* 等函数收到不支持的 view 名称。
normalize_views <- function(view) {
  view <- unique(as.character(view)); invalid <- setdiff(view, VIEW_LEVELS)
  if (length(invalid) > 0) stop("Unsupported view: ", paste(invalid, collapse=", "))
  view
}

# 文件名清洗：把分组名/样本名里的空格和特殊字符转成可安全落盘的形式。
safe_dir_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

plot_png_sidecar_path <- function(file_path) {
  if (!grepl("[.]pdf$", file_path, ignore.case = TRUE)) return(NULL)
  sub("[.]pdf$", ".png", file_path, ignore.case = TRUE)
}

save_ggplot_file <- function(plot_obj, file_path, width = 11, height = 7, dpi = 300) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  png_path <- plot_png_sidecar_path(file_path)
  if (!is.null(png_path)) {
    ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  }
  invisible(file_path)
}

save_plot_device <- function(file_path, draw_fn, width = 11, height = 7, dpi = 300) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file_path, width = width, height = height)
  tryCatch(draw_fn(), finally = grDevices::dev.off())
  png_path <- plot_png_sidecar_path(file_path)
  if (!is.null(png_path)) {
    grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
    tryCatch(draw_fn(), finally = grDevices::dev.off())
  }
  invisible(file_path)
}

save_pheatmap_file <- function(..., filename, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  pheatmap::pheatmap(..., filename = filename, width = width, height = height)
  png_path <- plot_png_sidecar_path(filename)
  if (!is.null(png_path)) {
    pheatmap::pheatmap(..., filename = png_path, width = width, height = height)
  }
  invisible(filename)
}

# 初始化模块输出目录，并清理本次模块会重写的 section。
# 这样每次重跑单个模块时，旧图不会和新图混在一起。
prepare_module_outputs <- function(fig_dir, tab_dir = NULL, sections = character(), stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  if (is.null(fig_dir) || !nzchar(fig_dir)) stop("fig_dir is not set for prepare_module_outputs().")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(tab_dir) && nzchar(tab_dir)) dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(sections) > 0) cleanup_plot_sections(fig_dir, unique(sections), stale_views = unique(c(ALL_VIEW_LEVELS, stale_views)))
}

# 按 “模块目录 / section / 样本子目录 / view 子目录” 规则构造输出路径。
# 这是整个结果树目录结构保持整齐的关键。
build_nested_dir <- function(base_dir, section, sample = NULL, view = "merged") {
  parts <- c(base_dir)
  if (length(VIEW_LEVELS) > 1) parts <- c(parts, normalize_views(view)[[1]])
  if (!is.null(section)) parts <- c(parts, section)
  if (!is.null(sample)) parts <- c(parts, safe_dir_name(sample))
  out <- do.call(file.path, as.list(parts)); dir.create(out, recursive=TRUE, showWarnings=FALSE); out
}

# 保存简单图：当只有一个 view 时直接落盘；若有多个 view，则按 view 分目录保存。
save_plot <- function(plot_obj, file_name, width=11, height=7, view="merged") {
  if (length(VIEW_LEVELS) == 1) {
    dir.create(OUTPUT_FIG_DIR, recursive=TRUE, showWarnings=FALSE)
    save_ggplot_file(plot_obj, file.path(OUTPUT_FIG_DIR, file_name), width=width, height=height)
  } else {
    for (v in normalize_views(view)) {
      d <- file.path(OUTPUT_FIG_DIR, v); dir.create(d, recursive=TRUE, showWarnings=FALSE)
      save_ggplot_file(plot_obj, file.path(d, file_name), width=width, height=height)
    }
  }
}

# 保存需要 section/sample 层级目录的图，是本项目最常用的保存入口。
save_plot_nested <- function(plot_obj, section, file_name, sample=NULL, width=11, height=7, view="merged") {
  for (v in normalize_views(view)) {
    d <- build_nested_dir(OUTPUT_FIG_DIR, section, sample=sample, view=v)
    save_ggplot_file(plot_obj, file.path(d, file_name), width=width, height=height)
  }
}

# 保存 PDF 到指定绝对路径
save_pdf <- function(plot_obj, filepath, width=11, height=7) {
  save_ggplot_file(plot_obj, filepath, width=width, height=height)
}

filter_direction_view <- function(df, view, direction_col="Direction") {
  if (view == "merged" || !direction_col %in% names(df)) return(df)
  df %>% filter(as.character(.data[[direction_col]]) == view)
}

view_title_suffix <- function(view) {
  switch(as.character(view), merged="", sense="", antisense=" (antisense only)", "")
}

plot_note <- function(title, note) {
  ggplot() + annotate("text", x=0, y=0, label=paste(title, note, sep="\n"), size=5) +
    xlim(-1,1) + ylim(-1,1) + theme_void()
}

# 常用小工具函数：
# - clean_numeric：尽量把字符安全转成数字
# - group_percent：组内百分比统一算法
# - get_col：宽表取列时不存在则补 0，避免处处写 if 判断
clean_numeric <- function(x) suppressWarnings(as.numeric(x))

group_percent <- function(x) { t <- sum(x, na.rm=TRUE); if(t>0) 100*x/t else rep(0, length(x)) }

calculate_gini <- function(x) {
  x <- sort(as.numeric(x[is.finite(x) & x >= 0]))
  n <- length(x)
  total <- sum(x, na.rm = TRUE)
  if (n == 0 || total <= 0) return(NA_real_)
  (2 * sum(seq_len(n) * x, na.rm = TRUE) / (n * total)) - ((n + 1) / n)
}

sequence_entropy <- function(x) {
  x <- as.numeric(x)
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  p <- x[x > 0 & is.finite(x)] / total
  -sum(p * log(p), na.rm = TRUE)
}

top_read_sum <- function(x, n = 1L) {
  x <- sort(as.numeric(x), decreasing = TRUE)
  sum(head(x[is.finite(x) & x > 0], n), na.rm = TRUE)
}

get_col <- function(df, col_name) {
  if (col_name %in% colnames(df)) df[[col_name]] else rep(0, nrow(df))
}

summarise_group_mean_rpm <- function(sample_level_df, key_cols, sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$RPM <- numeric()
    return(out)
  }
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!"Reads" %in% names(sample_level_df)) sample_level_df$Reads <- 0
  if (!"RPM" %in% names(sample_level_df)) stop("sample_level_df must contain an RPM column.")

  sample_summary <- sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      RPM = sum(RPM, na.rm = TRUE),
      .groups = "drop"
    )
  key_tbl <- sample_level_df %>%
    select(all_of(key_cols)) %>%
    distinct()
  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")
  group_summary <- sample_summary %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      RPM = sum(RPM, na.rm = TRUE),
      .groups = "drop"
    )

  tidyr::expand_grid(group_grid, key_tbl) %>%
    left_join(group_summary, by = c("Group", key_cols)) %>%
    mutate(
      Reads = replace_na(Reads, 0),
      RPM = if_else(sample_n > 0, replace_na(RPM, 0) / sample_n, 0)
    ) %>%
    select(-sample_n)
}

summarise_sparse_group_mean_rpm <- function(sample_level_df, key_cols, sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$RPM <- numeric()
    return(out)
  }
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!"Reads" %in% names(sample_level_df)) sample_level_df$Reads <- 0
  if (!"RPM" %in% names(sample_level_df)) stop("sample_level_df must contain an RPM column.")

  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")

  sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    left_join(group_grid, by = "Group") %>%
    mutate(RPM = if_else(sample_n > 0, RPM / sample_n, 0)) %>%
    select(-sample_n)
}

summarise_group_mean_percent <- function(sample_level_df, key_cols, percent_col = "Percent", reads_col = "Reads", sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$Percent <- numeric()
    return(out)
  }
  if (!percent_col %in% names(sample_level_df)) stop("sample_level_df must contain a percent column.")
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!reads_col %in% names(sample_level_df)) sample_level_df[[reads_col]] <- 0

  sample_summary <- sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(.data[[reads_col]], na.rm = TRUE),
      Percent = sum(.data[[percent_col]], na.rm = TRUE),
      .groups = "drop"
    )
  key_tbl <- sample_level_df %>%
    select(all_of(key_cols)) %>%
    distinct()
  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")
  group_summary <- sample_summary %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      Percent = sum(Percent, na.rm = TRUE),
      .groups = "drop"
    )

  tidyr::expand_grid(group_grid, key_tbl) %>%
    left_join(group_summary, by = c("Group", key_cols)) %>%
    mutate(
      Reads = replace_na(Reads, 0),
      Percent = if_else(sample_n > 0, replace_na(Percent, 0) / sample_n, 0)
    ) %>%
    select(-sample_n)
}

resolve_raw_fastq_path <- function(path_value) {
  value <- trimws(as.character(path_value))
  if (is.na(value) || !nzchar(value)) return(NA_character_)
  if (grepl("^/", value)) value else file.path(ROOT_DIR, value)
}

capture_command_version <- function(command, args = c("--version")) {
  command_path <- Sys.which(command)
  if (!nzchar(command_path)) {
    return(tibble(tool = command, available = FALSE, path = "", version = "not available"))
  }

  version_text <- tryCatch(
    {
      out <- suppressWarnings(system2(command_path, args = args, stdout = TRUE, stderr = TRUE))
      paste(head(out, 3), collapse = " | ")
    },
    error = function(e) paste("version check failed:", conditionMessage(e))
  )

  tibble(tool = command, available = TRUE, path = unname(command_path), version = version_text)
}

build_runtime_config_table <- function() {
  tibble(
    key = c(
      "ROOT_DIR", "RAW_DIR", "SPORTS_RESULT_DIR", "SPORTS_METADATA_DIR", "SPORTS_QC_DIR",
      "SPORTS_SAMPLE_ROOT", "SPORTS_ENDOSIRNA_DIR", "SPORTS_LOG_DIR", "RESULT_DIR",
      "PROJECT_SPECIES", "PROJECT_STAGE", "PROJECT_CONTRAST_BASELINE", "PROJECT_CONTRAST_CASE",
      "RUN_FASTQC", "RUN_EXOGENOUS", "EXOGENOUS_DB_ROOT",
      "MIRNA_TARGET_SOURCE", "MIRNA_SPECIES", "MIRNA_PREFIX", "MULTIMIR_ORG",
      "MIRNA_TARGET_ORGDB_PACKAGE", "KEGG_ORGANISM", "DE_SIGNIFICANCE",
      "RAW_AUDIT_MAX_READS", "RAW_AUDIT_KMER", "RAW_AUDIT_REFRESH",
      "REPORT_VIEW"
    ),
    value = c(
      ROOT_DIR, RAW_DIR, SPORTS_RESULT_DIR, SPORTS_METADATA_DIR, SPORTS_QC_DIR,
      SPORTS_SAMPLE_ROOT, SPORTS_ENDOSIRNA_DIR, SPORTS_LOG_DIR, RESULT_DIR,
      PROJECT_SPECIES, PROJECT_STAGE, PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE,
      RUN_FASTQC, RUN_EXOGENOUS, EXOGENOUS_DB_ROOT,
      MIRNA_TARGET_SOURCE, MIRNA_SPECIES, MIRNA_PREFIX, MULTIMIR_ORG,
      MIRNA_TARGET_ORGDB_PACKAGE, KEGG_ORGANISM, DE_SIGNIFICANCE,
      RAW_AUDIT_MAX_READS, RAW_AUDIT_KMER, RAW_AUDIT_REFRESH,
      REPORT_VIEW
    )
  )
}

is_mito_label <- function(x) {
  str_detect(replace_na(as.character(x),""), regex("(^|[-_])(mt|mito)|mitochond|chrM|12S-rRNA|16S-rRNA", ignore_case=TRUE))
}

remove_path_if_exists <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  } else if (file.exists(path)) {
    unlink(path, force = TRUE)
  }
}

cleanup_stale_view_dirs <- function(base_dir, stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  walk(stale_views, ~remove_path_if_exists(file.path(base_dir, .x)))
}

cleanup_plot_sections <- function(base_dir, sections, views = VIEW_LEVELS, stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  cleanup_stale_view_dirs(base_dir, stale_views = stale_views)
  for (section_name in sections) {
    remove_path_if_exists(file.path(base_dir, section_name))
  }
  for (view_name in unique(c(views, ALL_VIEW_LEVELS))) {
    for (section_name in sections) {
      remove_path_if_exists(file.path(base_dir, view_name, section_name))
    }
  }
}

# 多图拼接帮助函数。优先用 patchwork；若环境里没有 patchwork，则回退成说明图。
build_plot_collage <- function(plot_list, title_text = NULL, cols = 4) {
  valid_plots <- compact(plot_list)
  if (length(valid_plots) == 0) {
    return(plot_note(title_text %||% "Plot collage", "No plots available."))
  }
  if (!pkg_has("patchwork")) {
    return(plot_note(title_text %||% "Plot collage", "patchwork package is not available."))
  }

  patchwork::wrap_plots(valid_plots, ncol = cols) +
    patchwork::plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
      )
    )
}

# 根据拼图数量粗略估算输出宽高，避免手工反复调图幅。
collage_dimensions <- function(n_plots, cols = 4, base_width = 4, base_height = 3.1) {
  rows <- max(1, ceiling(n_plots / cols))
  list(
    width = max(10, min(18, cols * base_width)),
    height = max(6, rows * base_height + 0.9)
  )
}

# 下面这组函数负责把 SPORTS 原始的 Class / SubClass 字符串，
# 映射成项目内部统一使用的类别标签。后续所有统计和作图都依赖这里的规则。

# 根据 class / subclass 文本判断方向。当前项目把带 antisense 的记录单独归类。
class_direction <- function(class_name, sub_class="") {
  if_else(str_detect(paste(class_name, sub_class), "antisense"), "antisense", "sense")
}

# SPORTS 输出同时包含 Match_Genome / Unmatch_Genome 两类结果，这里统一提取状态字段。
class_genome_status <- function(class_name) {
  case_when(str_detect(class_name,"_Match_Genome$")~"Match_Genome", str_detect(class_name,"_Unmatch_Genome$")~"Unmatch_Genome", TRUE~"Unknown")
}

# 识别条目来自哪个数据库，方便后面做来源分层统计。
class_source_db <- function(class_name) {
  case_when(
    str_starts(class_name,"miRBase-")~"miRBase", str_starts(class_name,"rRNAdb-")~"rRNAdb",
    str_starts(class_name,"YRNAdb-")~"YRNAdb", str_starts(class_name,"GtRNAdb-")~"GtRNAdb",
    str_starts(class_name,"mitotRNAdb-")~"mitotRNAdb", str_starts(class_name,"piRNAdb-")~"piRNA",
    str_starts(class_name,"ensembl-")~"Ensembl", str_starts(class_name,"Rfam-")~"Rfam",
    str_starts(class_name,"Unannotated_")~"Unannotated", TRUE~"Other")
}

# 将复杂的 SPORTS 注释压缩成我们关注的大类：miRNA / piRNA / tsRNA / rsRNA 等。
class_main_category <- function(class_name, sub_class="") {
  token <- paste(class_name, sub_class)
  case_when(
    str_detect(token,"^Unannotated_Match_Genome")~"unanno_MG", str_detect(token,"^Unannotated_Unmatch_Genome")~"unanno_UMG",
    str_detect(token,"YRNA")~"ysRNA", str_detect(token,"miRNA")~"miRNA", str_detect(token,"piRNA")~"piRNA",
    str_detect(token,"scaRNA")~"scaRNA", str_detect(token,"snoRNA")~"snoRNA", str_detect(token,"snRNA")~"snRNA",
    str_detect(token,"lncRNA|lincRNA")~"lncRNA", str_detect(token,"vault_RNA")~"vault_RNA",
    str_detect(token,"misc_RNA")~"misc_RNA", str_detect(token,"ribozyme")~"ribozyme",
    str_detect(token,"tRNA")~"tsRNA", str_detect(token,"rRNA")~"rsRNA", TRUE~"other")
}

# output.txt 里的 Annotation 格式和 summary.txt 不一样，这里单独给它准备一套归类规则。
annotation_main_category <- function(annotation) {
  token <- if_else(is.na(annotation)|annotation=="", "NO_Annotation", annotation)
  case_when(
    str_detect(token,"YRNA")~"ysRNA", str_detect(token,regex("mir|let-",ignore_case=TRUE))~"miRNA",
    str_detect(token,regex("piRNA",ignore_case=TRUE))~"piRNA", str_detect(token,"scaRNA")~"scaRNA",
    str_detect(token,"snoRNA")~"snoRNA", str_detect(token,"snRNA")~"snRNA",
    str_detect(token,"lncRNA|lincRNA")~"lncRNA", str_detect(token,"vault_RNA")~"vault_RNA",
    str_detect(token,"misc_RNA")~"misc_RNA", str_detect(token,"ribozyme")~"ribozyme",
    str_detect(token,"tRNA")~"tsRNA", str_detect(token,"rRNA")~"rsRNA",
    str_detect(token,"NO_Annotation")~"other", TRUE~"other")
}

# 长度分布图只需要比较粗的主类，因此按 class_name 直接映射到长度类别。
length_category_from_class <- function(class_name) {
  case_when(
    str_detect(class_name,"miRNA")~"miRNA", str_detect(class_name,"piRNA")~"piRNA",
    str_detect(class_name,"YRNA")~"ysRNA", str_detect(class_name,"^mitotRNAdb")~"mt_tsRNA",
    str_detect(class_name,"Mt_rRNA")~"mt_rsRNA", str_detect(class_name,"tRNA")~"tsRNA",
    str_detect(class_name,"rRNA")~"rsRNA", str_detect(class_name,"lncRNA")~"lncRNA",
    str_detect(class_name,"snRNA")~"snRNA", str_detect(class_name,"snoRNA")~"snoRNA",
    str_detect(class_name,"scaRNA")~"scaRNA", str_detect(class_name,"misc_RNA")~"misc_RNA",
    str_detect(class_name,"vault_RNA")~"vault_RNA", str_detect(class_name,"ribozyme")~"ribozyme",
    str_detect(class_name,"scRNA")~"scRNA", str_detect(class_name,"(^|-)sRNA_")~"sRNA",
    str_detect(class_name,"mt_DNA")~"mt_DNA", str_detect(class_name,"DNA")~"DNA", TRUE~NA_character_)
}

# rsRNA 细分类提取器：
# - 优先从 SubClass 取 subtype
# - 若没有 subtype，再回退到 Class
# - 统一去掉 antisense 后缀，便于和 sense 结果共用同一套标签
parse_rsrna_class <- function(class_name, sub_class) {
  token <- if_else(sub_class=="-"|is.na(sub_class), class_name, sub_class) |> str_remove("-antisense$")
  case_when(
    str_detect(token,"45S-rRNA")~"45S-rRNA", str_detect(token,"28S-rRNA")~"28S-rRNA",
    str_detect(token,"26S-rRNA")~"26S-rRNA", str_detect(token,"25S-rRNA")~"25S-rRNA",
    str_detect(token,"23S-rRNA")~"23S-rRNA", str_detect(token,"18S-rRNA")~"18S-rRNA",
    str_detect(token,"17S-rRNA")~"17S-rRNA", str_detect(token,"16S-rRNA")~"16S-rRNA",
    str_detect(token,"12S-rRNA")~"12S-rRNA", str_detect(token,"5\\.8S-rRNA")~"5.8S-rRNA",
    str_detect(token,"5\\.3S-rRNA")~"5.3S-rRNA", str_detect(token,"5S-rRNA")~"5S-rRNA",
    str_detect(token,"4\\.5S-rRNA")~"4.5S-rRNA", str_detect(token,"2S-rRNA")~"2S-rRNA",
    str_detect(token,"mt_rRNA|Mt_rRNA")~"mt-rRNA", str_detect(token,"other-rRNA")~"other-rRNA",
    str_starts(class_name,"ensembl-rRNA")~"Ensembl-rRNA", str_starts(class_name,"Rfam-")~"Rfam-rRNA-like",
    TRUE~"other-rRNA")
}

# tsRNA 相关的拆解函数：
# - origin：区分基因组 tRNA 和线粒体 tRNA
# - maturity：前体还是成熟体
# - cleavage：5' / 3' / CCA / internal
# - amino acid / anticodon：从 subtype 名称里解析
parse_tsrna_origin <- function(class_name, sub_class) if_else(str_detect(paste(class_name,sub_class),"mitotRNAdb|mt_tRNA"), "Mitochondrial", "Genomic")
parse_tsrna_maturity <- function(class_name, sub_class) if_else(str_detect(paste(class_name,sub_class),"pre"), "Pre-tRNA", "Mature-tRNA")
parse_tsrna_cleavage <- function(class_name, sub_class) {
  token <- paste(class_name, sub_class)
  case_when(str_detect(token,"5_end")~"5'end", str_detect(token,"3_end")~"3'end", str_detect(token,"CCA_end")~"CCA-end", TRUE~"Internal/full")
}
parse_tsrna_amino_acid <- function(sub_class) {
  out <- str_match(sub_class, "(?:mt_)?tRNA-([A-Za-z]+[0-9]*)-[A-Z]{3}")[, 2]
  if_else(is.na(out) | out == "", "Unknown", out)
}
parse_tsrna_anticodon <- function(sub_class) {
  out <- str_match(sub_class, "(?:mt_)?tRNA-[A-Za-z]+[0-9]*-([A-Z]{3})")[, 2]
  if_else(is.na(out) | out == "", "Unknown", out)
}

# SPORTS 的 tRNA 末端类和 rRNAdb subtype 是父级总量的子集，不能和父级相加。
is_tsrna_nested_class <- function(class_name) {
  str_detect(as.character(class_name), "(?:tRNA|mt_tRNA)_(?:5_end|3_end|CCA_end)_")
}

is_rsrna_nested_length_class <- function(class_name) {
  token <- as.character(class_name)
  str_starts(token, "rRNAdb-") & !str_starts(token, "rRNAdb-rRNA_")
}

strip_tsrna_end_from_class <- function(class_name) {
  str_replace(as.character(class_name), "_(?:5_end|3_end|CCA_end)(?=_)", "")
}

strip_tsrna_end_from_subclass <- function(sub_class) {
  str_replace(as.character(sub_class), "_(?:5_end|3_end|CCA_end)$", "")
}

# 一些低频 rsRNA subtype 在可视化层面会合并到 other，以控制图面复杂度。
collapse_rsrna_visual_class <- function(x) case_when(x %in% c("17S-rRNA","23S-rRNA","25S-rRNA","26S-rRNA","5.3S-rRNA")~"other-rRNA", TRUE~x)

# 调色板相关函数集中放在这里，保证整套报告的视觉语言统一。

# 通过降低 RGB 亮度生成更深的颜色，常用于同色系文字标签。
darken_color <- function(cols, amount = 0.22) {
  if (length(cols) == 0) return(character())
  rgb_mat <- grDevices::col2rgb(cols) / 255
  darker <- rgb_mat * (1 - amount)
  darker[darker < 0] <- 0
  grDevices::rgb(darker[1, ], darker[2, ], darker[3, ])
}

# 将实际出现的 levels 尽量按照参考顺序排序，未出现在参考顺序里的放到末尾。
order_by_reference <- function(levels, reference_levels) {
  levels <- unique(as.character(levels))
  c(reference_levels[reference_levels %in% levels], setdiff(levels, reference_levels))
}

# 雷达图外环的文字一般使用扇区颜色的加深版，便于保持关联性同时保证可读性。
make_sector_label_palette <- function(sector_palette, amount = 0.28) {
  sector_palette <- sector_palette[!is.na(sector_palette)]
  setNames(darken_color(unname(sector_palette), amount = amount), names(sector_palette))
}

# 每组一个颜色，采用清新美观的莫兰迪/马卡龙色系配色，协调且区分度好
make_group_palette <- function(groups) {
  base <- c("#8ECFC9", "#FFBE7A", "#FF9F9B", "#A1C9F4", "#BEB8DC", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4")
  groups <- as.character(groups)
  if (length(groups) > length(base)) {
    extra <- grDevices::hcl.colors(length(groups) - length(base), palette = "Pastel 1")
    base <- c(base, extra)
  }
  setNames(base[seq_along(groups)], groups)
}

# 各类主图使用的主类别顺序与颜色。
category_levels <- c("rsRNA","tsRNA","miRNA","piRNA","ysRNA","lncRNA","snRNA","snoRNA","scaRNA","misc_RNA","vault_RNA","ribozyme","other","unanno_MG","unanno_UMG")

# 各类主图使用的主类别顺序与颜色，采用清新美观的莫兰迪/马卡龙色系配色，协调且区分度好
# 保证雷达图顺时针颜色与 rsRNA 雷达图一致，按扇区出现顺序分配颜色
sports_main_palette <- c(
  unanno_UMG="#8ECFC9", miRNA="#FFBE7A", piRNA="#FF9F9B", tsRNA="#A1C9F4", 
  rsRNA="#BEB8DC", other="#B3E2CD", unanno_MG="#FDCDAC", antisense="#CBD5E8", 
  ysRNA="#F4CAE4"
)
main_palette <- c(
  sports_main_palette, 
  lncRNA="#E6F5C9", snRNA="#FFF2AE", snoRNA="#F1E2CC", scaRNA="#CCCCCC", 
  misc_RNA="#B2DF8A", ribozyme="#FB8072", mt_tsRNA="#FDB462", mt_rsRNA="#FCCDE5", 
  DNA="#D9D9D9", mt_DNA="#BC80BD", scRNA="#CCEBC5", sRNA="#FFED6F", vault_RNA="#8DD3C7", 
  mt="#FF9F9B", cyto="#A1C9F4", Genomic="#A1C9F4", Mitochondrial="#FF9F9B"
)
category_palette <- main_palette[category_levels]

# sncRNA 百分比雷达图专用扇区顺序与标签。
radar_levels  <- c("unanno_UMG","miRNA","piRNA","tsRNA","rsRNA","other","unanno_MG")
radar_labels  <- c(unanno_UMG="unanno UMG", miRNA="miRNA", piRNA="piRNA", tsRNA="tsRNA", rsRNA="rsRNA", other="other", unanno_MG="unanno MG")
radar_sector_palette <- main_palette[radar_levels]
radar_label_palette  <- make_sector_label_palette(radar_sector_palette, amount = 0.28)

# 长度分布图关注的类别顺序。
sports_length_categories <- c("miRNA","piRNA","ysRNA","tsRNA","mt_tsRNA","rsRNA","mt_rsRNA","lncRNA","snRNA","snoRNA","scaRNA","misc_RNA","vault_RNA","ribozyme","scRNA","sRNA","DNA","mt_DNA","other","unanno_MG","unanno_UMG")
sports_length_palette <- main_palette[sports_length_categories]

# rRNA / rsRNA 的 subtype 顺序和颜色，与主类别色系保持一致。
rrna_display_order <- c("5S","5.8S","12S","16S","18S","28S","45S","4.5S","other")
rrna_coverage_subtypes <- setdiff(rrna_display_order, "other")
rrna_palette <- c("12S"="#8ECFC9","16S"="#FFBE7A","5S"="#FF9F9B","5.8S"="#A1C9F4","18S"="#BEB8DC","28S"="#B3E2CD","45S"="#FDCDAC","4.5S"="#F4CAE4","other"="#CBD5E8")

# rsRNA 细分类可视化用的一整套 levels / color / label。
# 此处前8个颜色保持与 sncRNA雷达图顺时针颜色（unanno_UMG->...->unanno_MG）高度一致，实现清新统一。
rsrna_visual_levels  <- c("5S-rRNA","5.8S-rRNA","12S-rRNA","16S-rRNA","18S-rRNA","28S-rRNA","45S-rRNA","4.5S-rRNA","2S-rRNA","mt-rRNA","Ensembl-rRNA","Rfam-rRNA-like","other-rRNA")
rsrna_visual_palette <- c(
  "12S-rRNA" = "#8ECFC9",
  "16S-rRNA" = "#FFBE7A",
  "5S-rRNA" = "#FF9F9B",
  "5.8S-rRNA" = "#A1C9F4",
  "18S-rRNA" = "#BEB8DC",
  "28S-rRNA" = "#B3E2CD",
  "45S-rRNA" = "#FDCDAC",
  "other-rRNA" = "#CBD5E8",
  "4.5S-rRNA" = "#F4CAE4",
  "2S-rRNA" = "#E6F5C9",
  "mt-rRNA" = "#FFF2AE",
  "Ensembl-rRNA" = "#F1E2CC",
  "Rfam-rRNA-like" = "#CCCCCC"
)
rsrna_visual_labels  <- c("5S-rRNA"="5S","5.8S-rRNA"="5.8S","12S-rRNA"="12S","16S-rRNA"="16S","18S-rRNA"="18S","28S-rRNA"="28S","45S-rRNA"="45S","4.5S-rRNA"="4.5S","2S-rRNA"="2S","mt-rRNA"="mt-rRNA","Ensembl-rRNA"="Ensembl","Rfam-rRNA-like"="Rfam","other-rRNA"="other")
rsrna_radar_levels  <- c("12S-rRNA","16S-rRNA","5S-rRNA","5.8S-rRNA","18S-rRNA","28S-rRNA","45S-rRNA","other-rRNA")
rsrna_radar_palette <- rsrna_visual_palette[rsrna_radar_levels]
rsrna_radar_labels  <- rsrna_visual_labels[rsrna_radar_levels]

# rsRNA 雷达图按 12S / 16S / 5S / 5.8S / 18S / 28S / 45S / other 的顺序显示。
collapse_rsrna_radar_class <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% setdiff(rsrna_radar_levels, "other-rRNA") ~ x,
    TRUE ~ "other-rRNA"
  )
}

# tsRNA 相关图的专用配色，与主类别色系保持一致。
ts_origin_palette    <- c("Genomic tRNA"="#A1C9F4","Mitochondrial tRNA"="#FF9F9B")
ts_genomic_palette   <- c("Genomic tRNA 5'end"="#8ECFC9","Genomic tRNA 3'end"="#FFBE7A","Genomic tRNA CCA-end"="#FF9F9B","Genomic tRNA internal"="#A1C9F4")
ts_mito_palette      <- c("Mitochondrial tRNA 5'end"="#8ECFC9","Mitochondrial tRNA 3'end"="#FFBE7A","Mitochondrial tRNA CCA-end"="#FF9F9B","Mitochondrial tRNA internal"="#A1C9F4")
ts_cleavage_palette  <- c("5'end"="#8ECFC9","3'end"="#FFBE7A","CCA-end"="#FF9F9B","Internal/full"="#A1C9F4")
ts_maturity_palette  <- c("Genomic precursor tRNA"="#B3E2CD","Genomic mature tRNA"="#A1C9F4","Mito mature tRNA"="#FF9F9B")
mt_cyto_palette      <- c(mt="#FF9F9B", cyto="#A1C9F4")
rsrna_source_palette <- c("rRNAdb known subtypes"="#A1C9F4","rRNAdb other"="#CBD5E8","Ensembl rRNA"="#FFBE7A","Ensembl Mt rRNA"="#FF9F9B")

# 氨基酸相关图沿用固定的 canonical 顺序和相对稳定的配色，便于跨图比较。
canonical_aa_order <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","SeC","iMet")
amino_acid_palette <- c(
  Ala = "#8ECFC9", Arg = "#FFBE7A", Asn = "#FF9F9B", Asp = "#A1C9F4", Cys = "#BEB8DC",
  Gln = "#B3E2CD", Glu = "#FDCDAC", Gly = "#CBD5E8", His = "#F4CAE4", Ile = "#E6F5C9",
  Leu = "#FFF2AE", Lys = "#F1E2CC", Met = "#CCCCCC", Phe = "#B2DF8A", Pro = "#FB8072",
  Ser = "#FDB462", Thr = "#FCCDE5", Trp = "#D9D9D9", Tyr = "#BC80BD", Val = "#CCEBC5",
  SeC = "#FFED6F", iMet = "#8DD3C7", Unknown = "#E5D8BD"
)

# 统一管理整套报告里反复出现的基础色值与主题设置。
# 这样各模块只引用这里的对象，不再各自写死颜色或 theme 参数。
mii_blue <- "#A1C9F4"
mii_orange <- "#FFBE7A"
mii_green <- "#8ECFC9"
mii_red <- "#FF9F9B"
mii_purple <- "#BEB8DC"

# 中性色按用途收敛，避免维护多个视觉差异很小的近似灰色。
mii_neutral_line <- "#C9CED6"
mii_neutral_text <- "#9AA1AC"
mii_axis_label_color <- "#6D7784"
mii_text_primary <- "#1F2329"
mii_text_muted <- "#47596E"
mii_transition_yellow <- "#FDFFB6"

mii_qc_palette <- c(CleanReads = mii_blue, MatchGenome = mii_green)
mii_diverging_palette <- c(low = mii_blue, mid = mii_transition_yellow, high = mii_red)
mii_direction_palette <- c(up = mii_red, down = mii_blue, ns = mii_neutral_line)
mii_direction_binary_palette <- c(up = mii_red, down = mii_blue)
mii_target_group_palette <- c(all = mii_purple, up = mii_red, down = mii_blue)
mii_enrichment_gradient <- c(mii_blue, mii_transition_yellow, mii_red)
mii_heatmap_gradient <- c(mii_blue, mii_transition_yellow, mii_red)

mii_theme_bw <- function(...) {
  theme_bw() + theme(panel.grid.minor = element_blank(), ...)
}

mii_theme_minimal <- function(...) {
  theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    ...
  )
}

# 读取并标准化样本元数据。
# 这里同时会把各类 SPORTS 输出文件路径挂到 metadata 上，方便后续逐样本批量读取。
sports_sample_dir <- function(sample_id) {
  file.path(SPORTS_SAMPLE_ROOT, as.character(sample_id))
}

sports_sample_file <- function(sample_id, ...) {
  file.path(sports_sample_dir(sample_id), ...)
}

sports_sample_processed_dir <- function(sample_id) {
  sports_sample_file(sample_id, paste0(as.character(sample_id), "_processed"))
}

sports_endosirna_summary_file <- function() {
  file.path(SPORTS_ENDOSIRNA_SUMMARY_DIR, "endo_siRNA_summary.tsv")
}

load_metadata <- function() {
  meta <- read.delim(READS_INFO_FILE, sep="\t", check.names=FALSE, quote="", comment.char="", colClasses="character", na.strings=character())
  required <- c("sample", "group", "biological_sample_id", "unit_type", "batch", "raw_fq", "condition", "age_group", "species", "stage", "replicate")
  if (!all(required %in% names(meta)) || nrow(meta)==0) stop("samples.tsv 为空或缺必要字段")
  if (anyDuplicated(meta$sample)) stop("sample 必须唯一")
  groups <- c(PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE)
  if (!setequal(unique(meta$group), groups) || length(unique(groups)) != 2) stop("显式两组与样本表不一致")
  out <- tibble(Group=factor(meta$group, levels=groups), Sample=factor(meta$sample, levels=meta$sample),
    BiologicalSample=meta$biological_sample_id, UnitType=meta$unit_type, Batch=meta$batch,
    Condition=meta$condition, AgeGroup=meta$age_group, Species=meta$species, Stage=meta$stage,
    Replicate=meta$replicate, RawFastqOriginal=meta$raw_fq) %>%
    mutate(GroupOrder=match(as.character(Group), groups), RawFastq=vapply(RawFastqOriginal, resolve_raw_fastq_path, character(1)),
      RawFastqExists=file.exists(RawFastq), SampleSportsDir=sports_sample_dir(as.character(Sample)),
      SummaryFile=sports_sample_file(as.character(Sample),paste0(Sample,"_summary.txt")),
      LengthFile=sports_sample_file(as.character(Sample),paste0(Sample,"_length_distribution.txt")),
      OutputFile=sports_sample_file(as.character(Sample),paste0(Sample,"_output.txt")),
      PDF_sncRNA=sports_sample_file(as.character(Sample),paste0(Sample,"_sncRNA_distribution.pdf")),
      PDF_rsRNA=sports_sample_file(as.character(Sample),paste0(Sample,"_rsRNA_distribution.pdf")),
      PDF_rRNA=sports_sample_file(as.character(Sample),paste0(Sample,"_rRNA_mapping.pdf")),
      PDF_tRNA=sports_sample_file(as.character(Sample),paste0(Sample,"_tRNA_mature_mapping.pdf")),
      RRNAMatchFile=file.path(sports_sample_processed_dir(as.character(Sample)),paste0(Sample,"_output_rRNA_match_genome")),
      RRNAUnmatchFile=file.path(sports_sample_processed_dir(as.character(Sample)),paste0(Sample,"_output_rRNA_unmatch_genome")))
  if (any(!file.exists(out$SummaryFile)|!file.exists(out$LengthFile)|!file.exists(out$OutputFile))) stop("MISSING SPORTS input; run sports.sh or provide explicit tables")
  out
}

smallrna_inference_gate <- function(meta=load_metadata()) {
  known <- function(x) !is.na(x) & !tolower(trimws(x)) %in% c("", "unknown", "na", "none")
  reason <- character()
  if (!all(known(meta$BiologicalSample))) reason <- c(reason,"unknown biological_sample_id")
  if (anyDuplicated(meta$BiologicalSample)) reason <- c(reason,"technical/repeated libraries require a separately specified aggregation; not supported by this edgeR entry")
  if (!all(meta$UnitType %in% c("individual","pool"))) reason <- c(reason,"unknown unit type")
  if (!all(known(meta$Batch)) || length(unique(meta$Batch)) != 1) reason <- c(reason,"this independent two-group entry requires one known batch")
  if (length(unique(meta$Group)) != 2 || any(table(meta$Group)<3)) reason <- c(reason,"need at least 3 independent units per group")
  list(eligible=length(reason)==0, reason=paste(reason,collapse="; "))
}

# 批量读取每个样本的 summary / length / output 文件，并补上样本与分组字段。
# 这三个入口分别对应：总览统计、长度分布、逐条注释结果。
load_summary_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$SummaryFile[i]); if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]), Replicate=as.character(si$Replicate[i]),
           Class=as.character(df$Class), SubClass=as.character(df$Sub_Class), Reads=clean_numeric(df$Reads))
  })
}
load_length_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$LengthFile[i]); if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]),
           Class=as.character(df$Class), Length=clean_numeric(df$Length), Reads=clean_numeric(df$Reads))
  })
}
load_output_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$OutputFile[i], select=c("ID","Sequence","Length","Reads","Match_Genome","Annotation"))
    if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]), ID=as.character(df$ID),
           Sequence=as.character(df$Sequence), Length=clean_numeric(df$Length), Reads=clean_numeric(df$Reads),
           MatchGenome=as.character(df$Match_Genome), Annotation=as.character(df$Annotation))
  })
}

normalize_filter_sequence <- function(x) {
  toupper(gsub("U", "T", trimws(as.character(x))))
}

build_sequence_complexity_tables <- function(output_all_tbl, sample_info_tbl, clean_reads_tbl, top_rank_limit = 50L) {
  sample_base <- sample_info_tbl %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group))

  empty_result <- list(
    sequence_complexity_by_sample = sample_base %>%
      left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
      mutate(
        CleanReads = replace_na(CleanReads, 0),
        SequenceReads = 0,
        UniqueSequences = 0,
        Top1Reads = 0,
        Top5Reads = 0,
        Top10Reads = 0,
        Top1CleanFraction = 0,
        Top5CleanFraction = 0,
        Top10CleanFraction = 0,
        Top1OutputFraction = 0,
        Top5OutputFraction = 0,
        Top10OutputFraction = 0,
        ShannonEntropy = NA_real_,
        EffectiveSequences = NA_real_,
        GiniIndex = NA_real_,
        DominantSequence = NA_character_,
        DominantAnnotation = NA_character_,
        DominantCategory = NA_character_,
        DominantLength = NA_real_,
        DominantReads = 0,
        DominantRPM = 0
      ),
    sequence_complexity_by_group = tibble(),
    sequence_rank_by_sample = tibble(),
    sequence_length_complexity = tibble()
  )

  if (is.null(output_all_tbl) || nrow(output_all_tbl) == 0) return(empty_result)

  if (pkg_has("data.table")) {
    dt <- data.table::as.data.table(output_all_tbl)
    dt <- data.table::copy(dt)
    if (!"PrimaryAnnotation" %in% names(dt)) {
      dt[, Annotation := ifelse(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation))]
      dt[, PrimaryAnnotation := sub(";.*$", "", Annotation)]
    }
    if (!"FeatureCategory" %in% names(dt)) {
      dt[, FeatureCategory := annotation_main_category(PrimaryAnnotation)]
    }

    dt[, Sample := as.character(Sample)]
    dt[, Group := as.character(Group)]
    dt[, Sequence := normalize_filter_sequence(Sequence)]
    dt[, Length := clean_numeric(Length)]
    dt[, Reads := clean_numeric(Reads)]
    dt[, Annotation := ifelse(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation))]
    dt[, PrimaryAnnotation := as.character(PrimaryAnnotation)]
    dt[, FeatureCategory := as.character(FeatureCategory)]
    dt <- dt[!is.na(Sequence) & Sequence != "" & is.finite(Reads) & Reads > 0]
    if (nrow(dt) == 0) return(empty_result)
    dt[is.na(Length) | Length <= 0, Length := nchar(Sequence)]

    sequence_counts_dt <- dt[, .(Reads = sum(Reads, na.rm = TRUE)), by = .(Sample, Group, Sequence, Length, PrimaryAnnotation, FeatureCategory)]
    data.table::setorder(sequence_counts_dt, Sample, Group, Sequence, Length, -Reads, PrimaryAnnotation)
    sequence_counts_dt <- sequence_counts_dt[, .(
      Reads = sum(Reads, na.rm = TRUE),
      PrimaryAnnotation = PrimaryAnnotation[1L],
      FeatureCategory = FeatureCategory[1L]
    ), by = .(Sample, Group, Sequence, Length)]

    data.table::setorder(sequence_counts_dt, Sample, -Reads, Sequence)
    sequence_counts_dt[, Rank := seq_len(.N), by = Sample]
    sequence_counts_dt[, SequenceReads := sum(Reads, na.rm = TRUE), by = Sample]
    sequence_counts_dt[, CumulativeReads := cumsum(Reads), by = Sample]
    sequence_counts_dt[, FractionOfOutput := ifelse(SequenceReads > 0, Reads / SequenceReads, 0)]
    sequence_counts_dt[, CumulativeFractionOfOutput := ifelse(SequenceReads > 0, CumulativeReads / SequenceReads, 0)]

    clean_dt <- data.table::as.data.table(clean_reads_tbl %>% mutate(Sample = as.character(Sample)))
    sequence_rank_dt <- merge(
      sequence_counts_dt[Rank <= top_rank_limit],
      clean_dt,
      by = "Sample",
      all.x = TRUE,
      sort = FALSE
    )
    sequence_rank_dt[is.na(CleanReads), CleanReads := 0]
    sequence_rank_dt[, FractionOfClean := ifelse(CleanReads > 0, Reads / CleanReads, 0)]
    sequence_rank_dt[, CumulativeFractionOfClean := ifelse(CleanReads > 0, CumulativeReads / CleanReads, 0)]
    sequence_rank_dt[, RPM := ifelse(CleanReads > 0, 1e6 * Reads / CleanReads, 0)]
    data.table::setorder(sequence_rank_dt, Sample, Rank)

    dominant_dt <- sequence_rank_dt[Rank == 1, .(
      Sample,
      DominantSequence = Sequence,
      DominantAnnotation = PrimaryAnnotation,
      DominantCategory = FeatureCategory,
      DominantLength = Length,
      DominantReads = Reads,
      DominantRPM = RPM
    )]

    sequence_complexity_dt <- sequence_counts_dt[, .(
      SequenceReads = sum(Reads, na.rm = TRUE),
      UniqueSequences = .N,
      Top1Reads = top_read_sum(Reads, 1),
      Top5Reads = top_read_sum(Reads, 5),
      Top10Reads = top_read_sum(Reads, 10),
      ShannonEntropy = sequence_entropy(Reads),
      EffectiveSequences = exp(sequence_entropy(Reads)),
      GiniIndex = calculate_gini(Reads)
    ), by = .(Sample, Group)]

    sample_base_dt <- data.table::as.data.table(sample_base)
    sequence_complexity_dt <- merge(sample_base_dt, sequence_complexity_dt, by = c("Sample", "Group"), all.x = TRUE, sort = FALSE)
    sequence_complexity_dt <- merge(sequence_complexity_dt, clean_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    sequence_complexity_dt <- merge(sequence_complexity_dt, dominant_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    numeric_zero_cols <- c("CleanReads", "SequenceReads", "UniqueSequences", "Top1Reads", "Top5Reads", "Top10Reads", "DominantReads", "DominantRPM")
    for (col_name in numeric_zero_cols) {
      if (col_name %in% names(sequence_complexity_dt)) {
        data.table::set(sequence_complexity_dt, which(is.na(sequence_complexity_dt[[col_name]])), col_name, 0)
      }
    }
    sequence_complexity_dt[, Top1CleanFraction := ifelse(CleanReads > 0, Top1Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top5CleanFraction := ifelse(CleanReads > 0, Top5Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top10CleanFraction := ifelse(CleanReads > 0, Top10Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top1OutputFraction := ifelse(SequenceReads > 0, Top1Reads / SequenceReads, 0)]
    sequence_complexity_dt[, Top5OutputFraction := ifelse(SequenceReads > 0, Top5Reads / SequenceReads, 0)]
    sequence_complexity_dt[, Top10OutputFraction := ifelse(SequenceReads > 0, Top10Reads / SequenceReads, 0)]
    data.table::setorder(sequence_complexity_dt, Group, Sample)

    complexity_metric_cols <- c(
      "CleanReads", "SequenceReads", "UniqueSequences",
      "Top1CleanFraction", "Top5CleanFraction", "Top10CleanFraction",
      "Top1OutputFraction", "Top5OutputFraction", "Top10OutputFraction",
      "ShannonEntropy", "EffectiveSequences", "GiniIndex", "DominantRPM"
    )
    for (col_name in complexity_metric_cols) {
      if (col_name %in% names(sequence_complexity_dt)) {
        sequence_complexity_dt[, (col_name) := as.numeric(.SD[[1L]]), .SDcols = col_name]
      }
    }
    sequence_group_dt <- sequence_complexity_dt[, c(
      list(SampleN = .N),
      unlist(lapply(.SD, function(x) {
        list(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
      }), recursive = FALSE)
    ), by = Group, .SDcols = complexity_metric_cols]
    expected_group_names <- c(
      "Group", "SampleN",
      unlist(lapply(complexity_metric_cols, function(x) paste0(x, c("_mean", "_sd", "_median"))))
    )
    if (length(expected_group_names) == ncol(sequence_group_dt)) names(sequence_group_dt) <- expected_group_names

    length_complexity_dt <- sequence_counts_dt[, .(Reads = sum(Reads, na.rm = TRUE), UniqueSequences = .N), by = .(Sample, Group, Length)]
    length_complexity_dt <- merge(length_complexity_dt, clean_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    length_complexity_dt[is.na(CleanReads), CleanReads := 0]
    length_complexity_dt[, RPM := ifelse(CleanReads > 0, 1e6 * Reads / CleanReads, 0)]
    length_complexity_dt[, FractionOfClean := ifelse(CleanReads > 0, Reads / CleanReads, 0)]
    data.table::setorder(length_complexity_dt, Sample, Length)

    return(list(
      sequence_complexity_by_sample = as_tibble(sequence_complexity_dt),
      sequence_complexity_by_group = as_tibble(sequence_group_dt),
      sequence_rank_by_sample = as_tibble(sequence_rank_dt),
      sequence_length_complexity = as_tibble(length_complexity_dt)
    ))
  }

  input_tbl <- output_all_tbl
  if (!"PrimaryAnnotation" %in% names(input_tbl)) {
    input_tbl <- input_tbl %>%
      mutate(
        Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation)),
        PrimaryAnnotation = str_split_fixed(Annotation, ";", 2)[, 1]
      )
  }
  if (!"FeatureCategory" %in% names(input_tbl)) {
    input_tbl <- input_tbl %>%
      mutate(FeatureCategory = annotation_main_category(PrimaryAnnotation))
  }

  output_clean <- input_tbl %>%
    mutate(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Sequence = normalize_filter_sequence(Sequence),
      Length = clean_numeric(Length),
      Reads = clean_numeric(Reads),
      Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation)),
      PrimaryAnnotation = as.character(PrimaryAnnotation),
      FeatureCategory = as.character(FeatureCategory)
    ) %>%
    filter(!is.na(Sequence), Sequence != "", Reads > 0) %>%
    mutate(Length = if_else(is.na(Length) | Length <= 0, as.numeric(nchar(Sequence)), Length))

  if (nrow(output_clean) == 0) return(empty_result)

  sequence_counts <- output_clean %>%
    group_by(Sample, Group, Sequence, Length, PrimaryAnnotation, FeatureCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    arrange(Sample, Group, Sequence, Length, desc(Reads), PrimaryAnnotation) %>%
    group_by(Sample, Group, Sequence, Length) %>%
    summarise(
      PrimaryAnnotation = first(PrimaryAnnotation),
      FeatureCategory = first(FeatureCategory),
      Reads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      PrimaryAnnotation = replace_na(PrimaryAnnotation, "NO_Annotation"),
      FeatureCategory = replace_na(FeatureCategory, "unanno_UMG")
    )

  sequence_rank_by_sample <- sequence_counts %>%
    group_by(Sample) %>%
    arrange(desc(Reads), Sequence, .by_group = TRUE) %>%
    mutate(
      Rank = row_number(),
      SequenceReads = sum(Reads, na.rm = TRUE),
      CumulativeReads = cumsum(Reads),
      FractionOfOutput = if_else(SequenceReads > 0, Reads / SequenceReads, 0),
      CumulativeFractionOfOutput = if_else(SequenceReads > 0, CumulativeReads / SequenceReads, 0)
    ) %>%
    ungroup() %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      FractionOfClean = if_else(CleanReads > 0, Reads / CleanReads, 0),
      CumulativeFractionOfClean = if_else(CleanReads > 0, CumulativeReads / CleanReads, 0),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)
    ) %>%
    filter(Rank <= top_rank_limit) %>%
    arrange(Sample, Rank)

  dominant_by_sample <- sequence_rank_by_sample %>%
    filter(Rank == 1) %>%
    transmute(
      Sample,
      DominantSequence = Sequence,
      DominantAnnotation = PrimaryAnnotation,
      DominantCategory = FeatureCategory,
      DominantLength = Length,
      DominantReads = Reads,
      DominantRPM = RPM
    )

  sequence_complexity_by_sample <- sequence_counts %>%
    group_by(Sample, Group) %>%
    summarise(
      SequenceReads = sum(Reads, na.rm = TRUE),
      UniqueSequences = n(),
      Top1Reads = top_read_sum(Reads, 1),
      Top5Reads = top_read_sum(Reads, 5),
      Top10Reads = top_read_sum(Reads, 10),
      ShannonEntropy = sequence_entropy(Reads),
      EffectiveSequences = exp(ShannonEntropy),
      GiniIndex = calculate_gini(Reads),
      .groups = "drop"
    ) %>%
    right_join(sample_base, by = c("Sample", "Group")) %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    left_join(dominant_by_sample, by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      SequenceReads = replace_na(SequenceReads, 0),
      UniqueSequences = replace_na(UniqueSequences, 0),
      Top1Reads = replace_na(Top1Reads, 0),
      Top5Reads = replace_na(Top5Reads, 0),
      Top10Reads = replace_na(Top10Reads, 0),
      Top1CleanFraction = if_else(CleanReads > 0, Top1Reads / CleanReads, 0),
      Top5CleanFraction = if_else(CleanReads > 0, Top5Reads / CleanReads, 0),
      Top10CleanFraction = if_else(CleanReads > 0, Top10Reads / CleanReads, 0),
      Top1OutputFraction = if_else(SequenceReads > 0, Top1Reads / SequenceReads, 0),
      Top5OutputFraction = if_else(SequenceReads > 0, Top5Reads / SequenceReads, 0),
      Top10OutputFraction = if_else(SequenceReads > 0, Top10Reads / SequenceReads, 0),
      DominantReads = replace_na(DominantReads, 0),
      DominantRPM = replace_na(DominantRPM, 0)
    ) %>%
    arrange(Group, Sample)

  complexity_metric_cols <- c(
    "CleanReads", "SequenceReads", "UniqueSequences",
    "Top1CleanFraction", "Top5CleanFraction", "Top10CleanFraction",
    "Top1OutputFraction", "Top5OutputFraction", "Top10OutputFraction",
    "ShannonEntropy", "EffectiveSequences", "GiniIndex", "DominantRPM"
  )

  sequence_complexity_by_group <- sequence_complexity_by_sample %>%
    group_by(Group) %>%
    summarise(
      SampleN = n(),
      across(
        all_of(complexity_metric_cols),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~stats::sd(.x, na.rm = TRUE),
          median = ~stats::median(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  sequence_length_complexity <- sequence_counts %>%
    group_by(Sample, Group, Length) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), UniqueSequences = n(), .groups = "drop") %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0),
      FractionOfClean = if_else(CleanReads > 0, Reads / CleanReads, 0)
    ) %>%
    arrange(Sample, Length)

  list(
    sequence_complexity_by_sample = sequence_complexity_by_sample,
    sequence_complexity_by_group = sequence_complexity_by_group,
    sequence_rank_by_sample = sequence_rank_by_sample,
    sequence_length_complexity = sequence_length_complexity
  )
}

empty_raw_fastq_audit_tables <- function(sample_info_tbl = tibble()) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) {
    sample_base <- tibble(Sample = character(), Group = character(), RawFastq = character(), RawFastqExists = logical())
  } else {
    sample_info_tbl <- as_tibble(sample_info_tbl)
    if (!"RawFastq" %in% names(sample_info_tbl)) sample_info_tbl$RawFastq <- ""
    if (!"Group" %in% names(sample_info_tbl)) sample_info_tbl$Group <- ""
    if (!"Sample" %in% names(sample_info_tbl)) sample_info_tbl$Sample <- ""
  sample_base <- sample_info_tbl %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
        RawFastq = as.character(RawFastq),
      RawFastqExists = file.exists(RawFastq)
    )
  }

  list(
    raw_fastq_audit_by_sample = sample_base %>%
      mutate(
        Status = if_else(RawFastqExists, "not_run", "missing_fastq"),
        MaxReadsRequested = RAW_AUDIT_MAX_READS,
        ReadsExamined = 0L,
        KmerK = RAW_AUDIT_KMER,
        MinLength = NA_real_,
        Q1Length = NA_real_,
        MedianLength = NA_real_,
        MeanLength = NA_real_,
        Q3Length = NA_real_,
        MaxLength = NA_real_,
        FractionLt18 = NA_real_,
        Fraction18to32 = NA_real_,
        FractionGt32 = NA_real_,
        AdapterSeedDetectedFraction = NA_real_,
        PolyA3pFraction = NA_real_,
        PolyT3pFraction = NA_real_,
        PolyG3pFraction = NA_real_,
        PolyC3pFraction = NA_real_,
        MostCommon5pKmer = NA_character_,
        MostCommon5pKmerFraction = NA_real_,
        MostCommon3pKmer = NA_character_,
        MostCommon3pKmerFraction = NA_real_
      ),
    raw_fastq_audit_by_group = tibble(),
    raw_fastq_length_distribution = tibble(),
    raw_fastq_end_base = tibble(),
    raw_fastq_kmer = tibble()
  )
}

read_fastq_sequences_limited <- function(file_path, max_reads = RAW_AUDIT_MAX_READS, chunk_reads = 10000L) {
  if (!file.exists(file_path) || max_reads <= 0) return(character())
  con <- if (str_ends(file_path, fixed(".gz"))) gzfile(file_path, open = "rt") else file(file_path, open = "rt")
  on.exit(close(con), add = TRUE)

  out <- character()
  while (length(out) < max_reads) {
    read_target <- min(chunk_reads, max_reads - length(out))
    lines <- readLines(con, n = read_target * 4L, warn = FALSE)
    if (length(lines) == 0) break
    seq_idx <- seq(2L, length(lines), by = 4L)
    seq_lines <- lines[seq_idx]
    if (length(seq_lines) > 0) out <- c(out, seq_lines)
    if (length(lines) < read_target * 4L) break
  }

  normalize_filter_sequence(out)
}

top_kmer_table <- function(seqs, sample_id, group_id, side = c("5prime", "3prime"), k = RAW_AUDIT_KMER, top_n = 25L) {
  side <- match.arg(side)
  seq_len_vec <- nchar(seqs)
  seqs <- seqs[seq_len_vec >= k]
  if (length(seqs) == 0) {
    return(tibble(Sample = character(), Group = character(), KmerSide = character(), Kmer = character(), ReadN = integer(), Fraction = numeric()))
  }
  kmers <- if (side == "5prime") substr(seqs, 1L, k) else substr(seqs, nchar(seqs) - k + 1L, nchar(seqs))
  kmers <- kmers[str_detect(kmers, "^[ACGTN]+$")]
  if (length(kmers) == 0) {
    return(tibble(Sample = character(), Group = character(), KmerSide = character(), Kmer = character(), ReadN = integer(), Fraction = numeric()))
  }

  kmer_tbl <- tibble(Kmer = kmers) %>%
    count(Kmer, name = "ReadN") %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      KmerSide = side,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    arrange(desc(ReadN), Kmer) %>%
    slice_head(n = top_n) %>%
    select(Sample, Group, KmerSide, Kmer, ReadN, Fraction)

  kmer_tbl
}

audit_one_raw_fastq <- function(sample_id, group_id, raw_fastq, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  raw_fastq <- as.character(raw_fastq %||% "")
  exists_flag <- nzchar(raw_fastq) && file.exists(raw_fastq)
  adapter_seed <- substr(Sys.getenv("ADAPTER3", ""), 1L, max(8L, min(k, 12L)))

  if (!exists_flag) {
    return(list(
      summary = empty_raw_fastq_audit_tables(tibble(Sample = sample_id, Group = group_id, RawFastq = raw_fastq))$raw_fastq_audit_by_sample %>%
        mutate(Status = "missing_fastq"),
      length_distribution = tibble(),
      end_base = tibble(),
      kmer = tibble()
    ))
  }

  seqs <- tryCatch(
    read_fastq_sequences_limited(raw_fastq, max_reads = max_reads),
    error = function(e) character()
  )
  seqs <- seqs[nzchar(seqs)]
  read_n <- length(seqs)
  if (read_n == 0) {
    return(list(
      summary = tibble(
        Sample = sample_id,
        Group = group_id,
        RawFastq = raw_fastq,
        RawFastqExists = TRUE,
        Status = "empty_or_unreadable",
        MaxReadsRequested = max_reads,
        ReadsExamined = 0L,
        KmerK = k,
        MinLength = NA_real_,
        Q1Length = NA_real_,
        MedianLength = NA_real_,
        MeanLength = NA_real_,
        Q3Length = NA_real_,
        MaxLength = NA_real_,
        FractionLt18 = NA_real_,
        Fraction18to32 = NA_real_,
        FractionGt32 = NA_real_,
        AdapterSeedDetectedFraction = NA_real_,
        PolyA3pFraction = NA_real_,
        PolyT3pFraction = NA_real_,
        PolyG3pFraction = NA_real_,
        PolyC3pFraction = NA_real_,
        MostCommon5pKmer = NA_character_,
        MostCommon5pKmerFraction = NA_real_,
        MostCommon3pKmer = NA_character_,
        MostCommon3pKmerFraction = NA_real_
      ),
      length_distribution = tibble(),
      end_base = tibble(),
      kmer = tibble()
    ))
  }

  len <- nchar(seqs)
  first_base <- substr(seqs, 1L, 1L)
  last_base <- substr(seqs, len, len)
  kmer_tbl <- bind_rows(
    top_kmer_table(seqs, sample_id, group_id, side = "5prime", k = k),
    top_kmer_table(seqs, sample_id, group_id, side = "3prime", k = k)
  )
  top5 <- kmer_tbl %>% filter(KmerSide == "5prime") %>% arrange(desc(ReadN), Kmer) %>% slice(1)
  top3 <- kmer_tbl %>% filter(KmerSide == "3prime") %>% arrange(desc(ReadN), Kmer) %>% slice(1)

  summary_tbl <- tibble(
    Sample = sample_id,
    Group = group_id,
    RawFastq = raw_fastq,
    RawFastqExists = TRUE,
    Status = "completed",
    MaxReadsRequested = max_reads,
    ReadsExamined = read_n,
    KmerK = k,
    MinLength = min(len, na.rm = TRUE),
    Q1Length = as.numeric(stats::quantile(len, 0.25, na.rm = TRUE, names = FALSE)),
    MedianLength = stats::median(len, na.rm = TRUE),
    MeanLength = mean(len, na.rm = TRUE),
    Q3Length = as.numeric(stats::quantile(len, 0.75, na.rm = TRUE, names = FALSE)),
    MaxLength = max(len, na.rm = TRUE),
    FractionLt18 = mean(len < 18, na.rm = TRUE),
    Fraction18to32 = mean(len >= 18 & len <= 32, na.rm = TRUE),
    FractionGt32 = mean(len > 32, na.rm = TRUE),
    AdapterSeedDetectedFraction = if (nzchar(adapter_seed)) mean(str_detect(seqs, fixed(adapter_seed)), na.rm = TRUE) else NA_real_,
    PolyA3pFraction = mean(str_detect(seqs, "A{5,}$"), na.rm = TRUE),
    PolyT3pFraction = mean(str_detect(seqs, "T{5,}$"), na.rm = TRUE),
    PolyG3pFraction = mean(str_detect(seqs, "G{5,}$"), na.rm = TRUE),
    PolyC3pFraction = mean(str_detect(seqs, "C{5,}$"), na.rm = TRUE),
    MostCommon5pKmer = if (nrow(top5) > 0) top5$Kmer[[1]] else NA_character_,
    MostCommon5pKmerFraction = if (nrow(top5) > 0) top5$Fraction[[1]] else NA_real_,
    MostCommon3pKmer = if (nrow(top3) > 0) top3$Kmer[[1]] else NA_character_,
    MostCommon3pKmerFraction = if (nrow(top3) > 0) top3$Fraction[[1]] else NA_real_
  )

  length_distribution <- tibble(Length = len) %>%
    count(Length, name = "ReadN") %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    select(Sample, Group, Length, ReadN, Fraction) %>%
    arrange(Sample, Length)

  end_base <- bind_rows(
    tibble(End = "5prime", Base = first_base),
    tibble(End = "3prime", Base = last_base)
  ) %>%
    filter(str_detect(Base, "^[ACGTN]$")) %>%
    count(End, Base, name = "ReadN") %>%
    group_by(End) %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Sample, Group, End, Base, ReadN, Fraction)

  list(summary = summary_tbl, length_distribution = length_distribution, end_base = end_base, kmer = kmer_tbl)
}

build_raw_fastq_audit_tables <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) return(empty_raw_fastq_audit_tables())
  audits <- purrr::map(seq_len(nrow(sample_info_tbl)), function(i) {
    audit_one_raw_fastq(
      sample_id = as.character(sample_info_tbl$Sample[[i]]),
      group_id = as.character(sample_info_tbl$Group[[i]]),
      raw_fastq = as.character(sample_info_tbl$RawFastq[[i]]),
      max_reads = max_reads,
      k = k
    )
  })

  raw_fastq_audit_by_sample <- bind_rows(map(audits, "summary"))
  raw_fastq_length_distribution <- bind_rows(map(audits, "length_distribution"))
  raw_fastq_end_base <- bind_rows(map(audits, "end_base"))
  raw_fastq_kmer <- bind_rows(map(audits, "kmer"))

  raw_fastq_audit_by_group <- raw_fastq_audit_by_sample %>%
    filter(Status == "completed") %>%
    group_by(Group) %>%
    summarise(
      SampleN = n(),
      ReadsExamined_mean = mean(ReadsExamined, na.rm = TRUE),
      MeanLength_mean = mean(MeanLength, na.rm = TRUE),
      MeanLength_sd = stats::sd(MeanLength, na.rm = TRUE),
      MedianLength_mean = mean(MedianLength, na.rm = TRUE),
      FractionLt18_mean = mean(FractionLt18, na.rm = TRUE),
      Fraction18to32_mean = mean(Fraction18to32, na.rm = TRUE),
      FractionGt32_mean = mean(FractionGt32, na.rm = TRUE),
      AdapterSeedDetectedFraction_mean = mean(AdapterSeedDetectedFraction, na.rm = TRUE),
      PolyA3pFraction_mean = mean(PolyA3pFraction, na.rm = TRUE),
      PolyT3pFraction_mean = mean(PolyT3pFraction, na.rm = TRUE),
      PolyG3pFraction_mean = mean(PolyG3pFraction, na.rm = TRUE),
      PolyC3pFraction_mean = mean(PolyC3pFraction, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    raw_fastq_audit_by_sample = raw_fastq_audit_by_sample,
    raw_fastq_audit_by_group = raw_fastq_audit_by_group,
    raw_fastq_length_distribution = raw_fastq_length_distribution,
    raw_fastq_end_base = raw_fastq_end_base,
    raw_fastq_kmer = raw_fastq_kmer
  )
}

build_raw_fastq_audit_manifest <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) {
    return(tibble(Sample = character(), Group = character(), RawFastq = character(), RawFastqExists = logical()))
  }

  sample_info_tbl %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      RawFastq = as.character(RawFastq),
      RawFastqExists = file.exists(RawFastq),
      RawFastqSize = map_dbl(RawFastq, ~if (file.exists(.x)) file.info(.x)$size else NA_real_),
      RawFastqMTime = map_chr(RawFastq, ~if (file.exists(.x)) as.character(file.info(.x)$mtime) else ""),
      MaxReadsRequested = max_reads,
      KmerK = k
    ) %>%
    arrange(Group, Sample)
}

raw_fastq_audit_manifest_matches <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  a <- as.data.frame(a, stringsAsFactors = FALSE)
  b <- as.data.frame(b, stringsAsFactors = FALSE)
  rownames(a) <- NULL
  rownames(b) <- NULL
  identical(a, b)
}

load_or_build_raw_fastq_audit_tables <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  manifest <- build_raw_fastq_audit_manifest(sample_info_tbl, max_reads = max_reads, k = k)
  cache_dir <- file.path(INTERNAL_CACHE_DIR, "raw_fastq_audit")
  cache_file <- file.path(cache_dir, "raw_fastq_audit_cache.rds")

  if (RAW_AUDIT_REFRESH != "1" && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(cached) && raw_fastq_audit_manifest_matches(cached$manifest, manifest) && is.list(cached$tables)) {
      return(cached$tables)
    }
  }

  tables <- build_raw_fastq_audit_tables(sample_info_tbl, max_reads = max_reads, k = k)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    saveRDS(list(manifest = manifest, tables = tables), cache_file),
    error = function(e) NULL
  )
  tables
}

# 这是整个项目最重要的共享数据构建函数。
# 它会把原始 SPORTS 结果整理成：
# - 样本级和组级统计表
# - 各类图直接可用的宽/长表
# - 后续模块需要复用的颜色、分组和 legacy 索引信息
load_mii_context <- function(report_view = REPORT_VIEW) {
  sample_info <- load_metadata()
  group_levels <- levels(sample_info$Group)
  group_palette <- make_group_palette(group_levels)

  # 原始 summary 和 length 表先统一加上方向、数据库来源、主类别等派生字段。
  summary_all_raw <- load_summary_tables(sample_info) %>%
    mutate(
      Direction = class_direction(Class, SubClass),
      GenomeStatus = class_genome_status(Class),
      SourceDB = class_source_db(Class),
      MainCategory = class_main_category(Class, SubClass)
    )
  length_all_raw <- load_length_tables(sample_info) %>%
    mutate(
      Direction = class_direction(Class),
      GenomeStatus = class_genome_status(Class),
      SourceDB = class_source_db(Class),
      MainCategory = class_main_category(Class),
      LengthCategory = length_category_from_class(Class)
    )
  output_all <- load_output_tables(sample_info)

  # Clean_Reads / Match_Genome 是很多百分比与 RPM 计算的分母，因此单独提前抽出来。
  clean_reads_tbl <- summary_all_raw %>%
    filter(Class == "Clean_Reads") %>%
    transmute(Sample, CleanReads = Reads)
  match_reads_tbl <- summary_all_raw %>%
    filter(Class == "Match_Genome") %>%
    transmute(Sample, MatchGenome = Reads)

  # 根据 report_view 过滤方向。当前默认是 sense，但函数本身保留了 merged/antisense 扩展能力。
  summary_all <- filter_direction_view(summary_all_raw, report_view)
  summary_all_view <- summary_all
  length_all <- filter_direction_view(length_all_raw, report_view)
  length_all_view <- length_all

  # summary_top 仅保留互斥主类别；tRNA 末端类是基础 tRNA 总量的子集，不能重复相加。
  summary_top <- summary_all %>%
    filter(
      SubClass == "-",
      !Class %in% c("Clean_Reads", "Match_Genome"),
      !is_tsrna_nested_class(Class)
    ) %>%
    mutate(MainCategory = factor(MainCategory, levels = category_levels))
  summary_top_view <- summary_top

  # 组级 clean/match 读数在大量组间比较图中复用，因此单独聚合一次。
  group_clean_reads_tbl <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group)) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    group_by(Group) %>%
    summarise(CleanReads = sum(CleanReads, na.rm = TRUE), .groups = "drop")
  group_match_tbl <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group)) %>%
    left_join(match_reads_tbl, by = "Sample") %>%
    group_by(Group) %>%
    summarise(MatchGenome = sum(MatchGenome, na.rm = TRUE), .groups = "drop")

  # 样本级百分比先计算；组间百分比用样本等权均值，避免重复数或测序深度改变比较尺度。
  sample_comp_total <- summary_top %>%
    group_by(Sample, Group, MainCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(Percent = if_else(CleanReads > 0, 100 * Reads / CleanReads, 0))
  sample_comp_view <- sample_comp_total
  composition_integrity_check <- sample_comp_total %>%
    group_by(Sample) %>%
    summarise(
      ExclusiveClassReads = sum(Reads, na.rm = TRUE),
      CleanReads = first(CleanReads),
      ReadDifference = ExclusiveClassReads - CleanReads,
      PercentSum = sum(Percent, na.rm = TRUE),
      PercentDifference = PercentSum - 100,
      Pass = abs(PercentDifference) <= 0.01,
      .groups = "drop"
    )
  if (any(!composition_integrity_check$Pass)) {
    bad_samples <- composition_integrity_check %>% filter(!Pass) %>% pull(Sample)
    stop("Mutually exclusive RNA classes did not recover 100% for: ", paste(bad_samples, collapse = ", "))
  }
  sample_comp_direction <- summary_top %>%
    group_by(Sample, Group, Direction, MainCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(Percent = if_else(CleanReads > 0, 100 * Reads / CleanReads, 0))
  group_comp_total <- sample_comp_total %>%
    summarise_group_mean_percent("MainCategory") %>%
    left_join(group_clean_reads_tbl, by = "Group")
  group_comp_view <- group_comp_total
  group_comp_direction <- sample_comp_direction %>%
    summarise_group_mean_percent(c("Direction", "MainCategory")) %>%
    left_join(group_clean_reads_tbl, by = "Group")
  group_comp_direction_view <- group_comp_direction %>% select(-Direction)

  sample_qc_table <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), Replicate = Replicate) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    left_join(match_reads_tbl, by = "Sample") %>%
    left_join(
      sample_comp_view %>%
        arrange(Sample, desc(Percent)) %>%
        group_by(Sample) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        transmute(Sample, DominantCategory = as.character(MainCategory), DominantPercent = round(Percent, 2)),
      by = "Sample"
    ) %>%
    mutate(GenomeRate = round(if_else(CleanReads > 0, 100 * MatchGenome / CleanReads, 0), 2))

  group_qc_table <- group_clean_reads_tbl %>%
    left_join(group_match_tbl, by = "Group") %>%
    left_join(
      sample_qc_table %>%
        group_by(Group = as.character(Group)) %>%
        summarise(SampleN = n_distinct(Sample), GenomeRate = round(mean(GenomeRate, na.rm = TRUE), 2), .groups = "drop"),
      by = "Group"
    )

  # 长度表同时包含父级总量与 tRNA/rRNA 子集；全局长度构成只保留父级总量。
  length_detail <- length_all %>%
    filter(
      !Class %in% c("Clean_Reads", "Match_Genome"),
      !is_tsrna_nested_class(Class),
      !is_rsrna_nested_length_class(Class)
    ) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0))
  length_detail_view <- length_detail
  length_focus_categories <- length_detail_view %>%
    group_by(MainCategory) %>%
    summarise(TotalRPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalRPM)) %>%
    pull(MainCategory) %>%
    as.character()
  length_focus_categories <- unique(c(c("rsRNA", "tsRNA", "miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA"), length_focus_categories)) %>%
    intersect(category_levels)
  sample_length_display <- length_detail_view %>%
    filter(MainCategory %in% length_focus_categories) %>%
    group_by(Sample, Group, Direction, MainCategory, Length) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop")
  sample_length_view <- sample_length_display %>% select(-Direction)
  length_composition_integrity_check <- sample_length_display %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    group_by(Sample) %>%
    summarise(
      ExclusiveLengthReads = sum(Reads, na.rm = TRUE),
      CleanReads = first(CleanReads),
      ReadDifference = ExclusiveLengthReads - CleanReads,
      PercentSum = if_else(CleanReads > 0, 100 * ExclusiveLengthReads / CleanReads, NA_real_),
      Pass = abs(ReadDifference) <= 2,
      .groups = "drop"
    )
  if (any(!length_composition_integrity_check$Pass)) {
    bad_samples <- length_composition_integrity_check %>% filter(!Pass) %>% pull(Sample)
    stop("Mutually exclusive length classes did not recover clean reads for: ", paste(bad_samples, collapse = ", "))
  }
  group_length_display <- sample_length_display %>%
    summarise_group_mean_rpm(c("Direction", "MainCategory", "Length"))
  group_length_view <- group_length_display %>% select(-Direction)

  # rsRNA 先做 subtype 解析，再做样本级/组级占比统计。
  # visual 版本会把低频 subtype 合并，方便雷达图和 alluvial 展示。
  rsrna_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "rsRNA") %>%
    mutate(RsRNAClass = parse_rsrna_class(Class, SubClass))
  sample_rsrna <- rsrna_detail %>%
    group_by(Sample, Group, Direction, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_rsrna <- sample_rsrna %>%
    summarise_group_mean_percent(c("Direction", "RsRNAClass"))
  sample_rsrna_visual <- sample_rsrna %>%
    mutate(RsRNAClass = collapse_rsrna_visual_class(as.character(RsRNAClass))) %>%
    group_by(Sample, Group, Direction, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup() %>%
    mutate(RsRNAClass = factor(RsRNAClass, levels = rsrna_visual_levels))
  sample_rsrna_visual_by_view <- setNames(
    list(sample_rsrna_visual %>% select(-Direction)),
    report_view
  )
  sample_rsrna_visual_view <- sample_rsrna_visual_by_view[[report_view]]
  group_rsrna_visual <- sample_rsrna_visual %>%
    summarise_group_mean_percent(c("Direction", "RsRNAClass")) %>%
    mutate(RsRNAClass = factor(RsRNAClass, levels = rsrna_visual_levels))
  group_rsrna_visual_view <- group_rsrna_visual %>% select(-Direction)
  group_rsrna_visual_merged <- group_rsrna_visual_view %>%
    group_by(Group, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), Percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
    ungroup()
  group_rsrna_view <- group_rsrna %>% select(-Direction)

  # SPORTS 的 tRNA 末端行是基础总量的子集。基础表用于总量/origin/AA/anticodon；
  # cleavage 表以“基础总量 - 5' - 3' - CCA”重建 Internal/full，保证四类互斥。
  tsrna_base_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "tsRNA", !is_tsrna_nested_class(Class)) %>%
    mutate(
      BaseClass = as.character(Class),
      BaseSubClass = as.character(SubClass),
      Origin = parse_tsrna_origin(Class, SubClass),
      Maturity = parse_tsrna_maturity(Class, SubClass),
      AminoAcid = parse_tsrna_amino_acid(SubClass),
      Anticodon = parse_tsrna_anticodon(SubClass)
    )

  tsrna_end_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "tsRNA", is_tsrna_nested_class(Class)) %>%
    mutate(
      BaseClass = strip_tsrna_end_from_class(Class),
      BaseSubClass = strip_tsrna_end_from_subclass(SubClass),
      Origin = parse_tsrna_origin(BaseClass, BaseSubClass),
      Maturity = parse_tsrna_maturity(BaseClass, BaseSubClass),
      Cleavage = parse_tsrna_cleavage(Class, SubClass),
      AminoAcid = parse_tsrna_amino_acid(BaseSubClass),
      Anticodon = parse_tsrna_anticodon(BaseSubClass)
    )

  tsrna_reconstruction_keys <- c(
    "Sample", "Group", "Direction", "BaseClass", "BaseSubClass",
    "Origin", "Maturity", "AminoAcid", "Anticodon"
  )
  orphan_tsrna_ends <- tsrna_end_detail %>%
    anti_join(tsrna_base_detail, by = tsrna_reconstruction_keys)
  if (nrow(orphan_tsrna_ends) > 0) {
    stop("SPORTS tRNA end subsets lack matching base-total rows.")
  }

  tsrna_end_totals <- tsrna_end_detail %>%
    group_by(across(all_of(tsrna_reconstruction_keys))) %>%
    summarise(EndReads = sum(Reads, na.rm = TRUE), .groups = "drop")

  tsrna_cleavage_reconstruction_check <- tsrna_base_detail %>%
    transmute(across(all_of(tsrna_reconstruction_keys)), BaseReads = Reads) %>%
    left_join(tsrna_end_totals, by = tsrna_reconstruction_keys) %>%
    mutate(
      EndReads = coalesce(EndReads, 0),
      InternalReadsRaw = BaseReads - EndReads,
      ExcessEndReads = pmax(EndReads - BaseReads, 0),
      Pass = ExcessEndReads <= 1,
      InternalReads = pmax(InternalReadsRaw, 0),
      ReconstructedReads = EndReads + InternalReads,
      ResidualReads = ReconstructedReads - BaseReads
    )
  if (any(!tsrna_cleavage_reconstruction_check$Pass)) {
    bad_samples <- tsrna_cleavage_reconstruction_check %>% filter(!Pass) %>% pull(Sample) %>% unique()
    stop("tRNA end subsets exceed base totals by more than 1 read for: ", paste(bad_samples, collapse = ", "))
  }

  tsrna_internal_detail <- tsrna_cleavage_reconstruction_check %>%
    transmute(
      Sample, Group, Direction,
      Class = BaseClass, SubClass = BaseSubClass, Reads = InternalReads,
      MainCategory = "tsRNA", BaseClass, BaseSubClass, Origin, Maturity,
      Cleavage = "Internal/full", AminoAcid, Anticodon
    )
  tsrna_cleavage_detail <- bind_rows(
    tsrna_end_detail %>%
      transmute(
        Sample, Group, Direction, Class, SubClass, Reads,
        MainCategory = "tsRNA", BaseClass, BaseSubClass, Origin, Maturity,
        Cleavage, AminoAcid, Anticodon
      ),
    tsrna_internal_detail
  ) %>%
    filter(Reads > 0)

  sample_tsrna_origin <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, Origin) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_origin <- sample_tsrna_origin %>%
    summarise_group_mean_percent(c("Direction", "Origin"))
  group_tsrna_origin_view <- group_tsrna_origin %>% select(-Direction)
  sample_tsrna_cleavage <- tsrna_cleavage_detail %>%
    group_by(Sample, Group, Direction, Origin, Cleavage) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction, Origin) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_cleavage <- sample_tsrna_cleavage %>%
    summarise_group_mean_percent(c("Direction", "Origin", "Cleavage"))
  group_tsrna_cleavage_view <- group_tsrna_cleavage %>% select(-Direction)
  sample_tsrna_aa <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, AminoAcid) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_aa <- sample_tsrna_aa %>%
    summarise_group_mean_percent(c("Direction", "AminoAcid"))
  group_tsrna_aa_view <- group_tsrna_aa %>% select(-Direction)
  sample_tsrna_anticodon <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, Anticodon) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_anticodon <- sample_tsrna_anticodon %>%
    summarise_group_mean_percent(c("Direction", "Anticodon"))
  group_tsrna_anticodon_view <- group_tsrna_anticodon %>% select(-Direction)

  # output_all 保留逐条序列粒度。PrimaryAnnotation 只用于展示和逐序列分析，
  # tRNA 总量及 isoacceptor 统计使用上面的 SPORTS 分数分配结果，不由首个注释决定。
  output_all <- output_all %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(
      Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", Annotation),
      PrimaryAnnotation = str_split_fixed(Annotation, ";", 2)[, 1],
      FeatureCategory = annotation_main_category(PrimaryAnnotation),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)
    )

  sequence_complexity_tables <- build_sequence_complexity_tables(output_all, sample_info, clean_reads_tbl)
  sequence_complexity_by_sample <- sequence_complexity_tables$sequence_complexity_by_sample
  sequence_complexity_by_group <- sequence_complexity_tables$sequence_complexity_by_group
  sequence_rank_by_sample <- sequence_complexity_tables$sequence_rank_by_sample
  sequence_length_complexity <- sequence_complexity_tables$sequence_length_complexity
  raw_fastq_audit_tables <- load_or_build_raw_fastq_audit_tables(sample_info)
  raw_fastq_audit_by_sample <- raw_fastq_audit_tables$raw_fastq_audit_by_sample
  raw_fastq_audit_by_group <- raw_fastq_audit_tables$raw_fastq_audit_by_group
  raw_fastq_length_distribution <- raw_fastq_audit_tables$raw_fastq_length_distribution
  raw_fastq_end_base <- raw_fastq_audit_tables$raw_fastq_end_base
  raw_fastq_kmer <- raw_fastq_audit_tables$raw_fastq_kmer

  group_top_features <- output_all %>%
    group_by(Sample, Group, PrimaryAnnotation, FeatureCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    summarise_group_mean_rpm(c("PrimaryAnnotation", "FeatureCategory")) %>%
    filter(RPM > 0 | Reads > 0) %>%
    arrange(Group, desc(RPM), desc(Reads), PrimaryAnnotation) %>%
    group_by(Group) %>%
    mutate(Rank = row_number()) %>%
    filter(Rank <= 20) %>%
    ungroup()

  # legacy_pdf_index 用来保留 SPORTS 原生 PDF 的索引，便于在结果表里回链查看。
  ratio_all_sample_df <- sample_comp_view %>%
    filter(Reads > 0) %>%
    mutate(MainCategory = factor(as.character(MainCategory), levels = category_levels))
  legacy_pdf_index <- sample_info %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      sncRNA_distribution = PDF_sncRNA,
      rsRNA_distribution = PDF_rsRNA,
      rRNA_mapping = PDF_rRNA,
      tRNA_mature_mapping = PDF_tRNA
    )

  # keep_names 决定哪些对象会暴露给后续模块环境。
  # 这样可以明确共享边界，避免把太多临时变量泄露出去。
  keep_names <- c(
    "sample_info", "group_levels", "group_palette",
    "summary_all_raw", "summary_all", "summary_all_view", "summary_top", "summary_top_view",
    "length_all_raw", "length_all", "length_all_view",
    "clean_reads_tbl", "match_reads_tbl", "group_clean_reads_tbl", "group_match_tbl",
    "group_comp_total", "group_comp_view", "group_comp_direction_view",
    "sample_comp_view", "composition_integrity_check", "sample_qc_table", "group_qc_table",
    "length_detail", "length_detail_view", "length_focus_categories", "sample_length_display", "sample_length_view", "group_length_display", "group_length_view", "length_composition_integrity_check",
    "rsrna_detail", "sample_rsrna_visual_by_view", "sample_rsrna_visual_view", "group_rsrna_visual_view", "group_rsrna_visual_merged", "group_rsrna_view",
    "tsrna_base_detail", "tsrna_cleavage_detail", "tsrna_cleavage_reconstruction_check",
    "group_tsrna_origin", "group_tsrna_origin_view", "group_tsrna_cleavage", "group_tsrna_cleavage_view", "group_tsrna_aa_view", "group_tsrna_anticodon_view",
    "output_all", "group_top_features", "ratio_all_sample_df", "legacy_pdf_index",
    "sequence_complexity_by_sample", "sequence_complexity_by_group", "sequence_rank_by_sample", "sequence_length_complexity",
    "raw_fastq_audit_by_sample", "raw_fastq_audit_by_group", "raw_fastq_length_distribution", "raw_fastq_end_base", "raw_fastq_kmer"
  )
  mget(keep_names, inherits = FALSE)
}

# 把 load_mii_context 生成的一整套对象注入到指定环境中。
activate_mii_context <- function(env = parent.frame(), report_view = REPORT_VIEW) {
  list2env(load_mii_context(report_view = report_view), envir = env)
  invisible(env)
}

# 下面开始是整套报告共用的绘图核心函数。
# 顺序大致是：雷达图 -> 趋势图 -> 饼图。

# ===== 雷达图核心函数 =====

# 把稀疏的长度表补成连续长度序列，便于画线图时不因为缺失长度而断线。
complete_length_series <- function(df, len, category, clean_total) {
  if (nrow(df)==0) out <- tibble(Length=len,Reads=0) else
    out <- df %>% group_by(Length) %>% summarise(Reads=sum(Reads,na.rm=TRUE),.groups="drop") %>%
      right_join(tibble(Length=len),by="Length") %>% mutate(Reads=replace_na(Reads,0)) %>% arrange(Length)
  out %>% mutate(Category=category, RPM=if(clean_total>0) 1e6*Reads/clean_total else rep(0,n())) %>% select(Category,Length,Reads,RPM)
}

# 把样本/组级主类别统计整理成 sncRNA 总览雷达图需要的 pie/radar 输入表。
compute_pie_tables <- function(summary_top_df) {
  summary_top_df %>% filter(!is.na(MainCategory)) %>%
    mutate(PieCategory=case_when(MainCategory%in%c("rsRNA","tsRNA","piRNA","miRNA","unanno_UMG","unanno_MG")~as.character(MainCategory),TRUE~"other")) %>%
    group_by(Sample,Group,PieCategory) %>% summarise(Reads=sum(Reads,na.rm=TRUE),.groups="drop") %>%
    group_by(Sample,Group) %>% complete(PieCategory=radar_levels,fill=list(Reads=0)) %>%
    mutate(Percent={t<-sum(Reads);if(t>0)100*Reads/t else rep(0,n())}, PieCategory=factor(PieCategory,levels=radar_levels)) %>% ungroup()
}

# 单组/单样本扇区雷达图的数据构造器。
# 这里会给每个类别分配一个扇区角度、计算显示标签、并把原始百分比换成归一化半径。
build_sector_radar_df <- function(input_df, radar_max, sector_levels=radar_levels, sector_labels=radar_labels, sector_palette=radar_sector_palette, bottom_anchor_category="rsRNA", angle_offset=0) {
  n <- length(sector_levels); sw <- 2*pi/n
  bai <- match(bottom_anchor_category,sector_levels); if(is.na(bai)) bai <- ceiling(n/2)
  label_palette <- make_sector_label_palette(sector_palette[sector_levels], amount = 0.28)
  input_df %>% mutate(PieCategory=factor(as.character(PieCategory),levels=sector_levels)) %>%
    complete(PieCategory=factor(sector_levels,levels=sector_levels),fill=list(Reads=0,Percent=0)) %>% arrange(PieCategory) %>%
    mutate(SectorId=row_number(), CenterAngle=-pi/2+(bai-SectorId)*sw+angle_offset, StartAngle=CenterAngle+sw/2, EndAngle=CenterAngle-sw/2,
           Radius=Percent/radar_max, DisplayLabel=unname(sector_labels[as.character(PieCategory)]),
           DisplayValue=sprintf("%.2f",Percent), ValueColor=mii_text_primary, LabelColor=unname(label_palette[as.character(PieCategory)]),
           UpperHalf=sin(CenterAngle)>=0,
           TextAngle=((if_else(UpperHalf,CenterAngle*180/pi-90,CenterAngle*180/pi+90)+180)%%360)-180)
}

# 将扇区中心角转换成实际文本坐标；rv 控制文本落在离圆心多远的位置。
build_sector_text_df <- function(pdf,tcol,rv,ccol) tibble(Label=pdf[[tcol]],TextX=rv*cos(pdf$CenterAngle),TextY=rv*sin(pdf$CenterAngle),TextAngle=pdf$TextAngle,TextColor=pdf[[ccol]])

# 外圈类别名称沿着弧线逐字符摆放，保证长标签在圆周附近仍然可读。
build_sector_arc_text_df <- function(pdf,tcol,rv,sw,ccol) {
  pmap_dfr(list(pdf$CenterAngle,pdf[[tcol]],pdf[[ccol]],pdf$UpperHalf), function(ca,tv,tc,uh) {
    ch <- strsplit(as.character(tv),"",fixed=TRUE)[[1]]; ch <- ch[nzchar(ch)]; if(length(ch)==0) return(tibble())
    ts <- min(sw*0.72, sw*0.068*max(length(ch)-1,1))
    if(uh) { ca2<-seq(ca+ts/2,ca-ts/2,length.out=length(ch)); ta<-ca2*180/pi-90 }
    else   { ca2<-seq(ca-ts/2,ca+ts/2,length.out=length(ch)); ta<-ca2*180/pi+90 }
    tibble(Label=ch,TextX=rv*cos(ca2),TextY=rv*sin(ca2),TextAngle=((ta+180)%%360)-180,TextColor=tc)
  })
}

# 负责生成扇区多边形坐标：
# - use_data=TRUE 时画中心向外延伸的数据扇形
# - 否则画最外层的彩色标签环
build_sector_polygon_df <- function(pdf,ir,or,use_data=FALSE,pc=160) {
  pmap_dfr(list(pdf$SectorId,pdf$PieCategory,pdf$StartAngle,pdf$EndAngle,pdf$Radius), function(sid,pcat,sa,ea,rv) {
    cor <- if(use_data) rv else or; if(is.na(cor)||cor<=ir+1e-8) return(tibble())
    to <- seq(sa,ea,length.out=pc); od <- tibble(x=cor*cos(to),y=cor*sin(to))
    pd <- if(ir<=0) bind_rows(tibble(x=0,y=0),od,tibble(x=0,y=0)) else {
      ti <- seq(ea,sa,length.out=pc); id <- tibble(x=ir*cos(ti),y=ir*sin(ti)); bind_rows(od,id) }
    pd %>% mutate(Group=paste0("sector_",sid,"_",if(use_data)"wedge"else"ring"), PieCategory=pcat)
  })
}

# 生成同心圆的坐标路径，作为雷达图的背景网格。
build_circle_path_df <- function(rv,pc=240) { th<-seq(0,2*pi,length.out=pc); map_dfr(rv,function(r) tibble(r=r,x=r*cos(th),y=r*sin(th))) }

# 单个雷达图的主绘制函数。
# 这里同时画背景同心圆、中心扇形、外层标签环、扇区分割线、类别名和百分比数值。
plot_sector_radar_core <- function(plot_df, title_text, radar_max, radar_breaks, sector_levels=radar_levels, sector_palette=radar_sector_palette,
  axis_label="Relative Abundance (%)", label_radius=1.24, value_radius=NULL, coord_xlim=c(-1.92,1.72), coord_ylim=c(-1.54,1.54), label_size=3.35, value_size=3.55) {
  rir<-1.02; ror<-1.16; if(is.null(value_radius)) value_radius<-0.93
  sw <- 2*pi/length(sector_levels)
  cdf <- build_circle_path_df(radar_breaks[radar_breaks>0]/radar_max)
  wdf <- build_sector_polygon_df(plot_df,ir=0,or=NA_real_,use_data=TRUE)
  rdf <- build_sector_polygon_df(plot_df,ir=rir,or=ror)
  sdf <- tibble(x=0,y=0,xend=ror*cos(plot_df$StartAngle),yend=ror*sin(plot_df$StartAngle))
  ldf <- build_sector_arc_text_df(plot_df,"DisplayLabel",label_radius,sw,"LabelColor")
  vdf <- build_sector_text_df(plot_df,"DisplayValue",value_radius,"ValueColor")
  ggplot() +
    geom_path(data=cdf,aes(x=x,y=y,group=r),color=mii_neutral_line,linewidth=0.55) +
    geom_polygon(data=wdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.4,linejoin="round",alpha=0.86) +
    geom_polygon(data=rdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.45,linejoin="round",alpha=1) +
    geom_segment(data=sdf,aes(x=x,y=y,xend=xend,yend=yend),color=mii_neutral_line,linewidth=0.75,lineend="round") +
    geom_text(data=ldf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=label_size,fontface="bold",inherit.aes=FALSE) +
    geom_text(data=vdf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=value_size,fontface="bold",inherit.aes=FALSE) +
    annotate("text",x=-1.60,y=0,label=axis_label,angle=90,size=4.0,color=mii_axis_label_color,fontface="bold") +
    scale_fill_manual(values=sector_palette,drop=FALSE) + scale_color_identity() +
    coord_fixed(xlim=coord_xlim,ylim=coord_ylim,clip="off") + theme_void() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5,face="bold",size=15,color=mii_text_muted),
          plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
          plot.margin=margin(10,20,10,28)) + labs(title=title_text)
}

# 最基础的 sncRNA 样本雷达图封装，默认走公共的 radar_levels。
plot_sector_radar_sample <- function(sample_id, pie_df, rmax, rbrk, view="merged") {
  pdf <- pie_df %>% filter(as.character(Sample)==sample_id) %>% build_sector_radar_df(radar_max=rmax)
  plot_sector_radar_core(pdf, paste0("sncRNA Proportion of ",sample_id,view_title_suffix(view)), rmax, rbrk)
}

# 组间共享雷达图的数据构造器。
# 和单图版本的区别是：每个扇区会继续按组数切成更细的小扇片。
prepare_group_sector_radar <- function(input_df, radar_max, sector_levels = radar_levels, sector_labels = radar_labels,
  bottom_anchor_category = "rsRNA", sector_palette = radar_sector_palette, angle_offset = 0) {
  n_sector <- length(sector_levels)
  sector_width <- 2 * pi / n_sector
  bai <- match(bottom_anchor_category, sector_levels)
  if (is.na(bai)) bai <- ceiling(n_sector / 2)
  label_palette <- make_sector_label_palette(sector_palette[sector_levels], amount = 0.28)

  sector_df <- tibble(PieCategory = factor(sector_levels, levels = sector_levels)) %>%
    mutate(
      SectorId = row_number(),
      SectorCenterAngle = -pi / 2 + (bai - SectorId) * sector_width + angle_offset,
      SectorStartAngle = SectorCenterAngle + sector_width / 2,
      SectorEndAngle = SectorCenterAngle - sector_width / 2,
      CenterAngle = SectorCenterAngle,
      StartAngle = SectorStartAngle,
      EndAngle = SectorEndAngle,
      DisplayLabel = unname(sector_labels[as.character(PieCategory)]),
      LabelColor = unname(label_palette[as.character(PieCategory)]),
      UpperHalf = sin(CenterAngle) >= 0,
      TextAngle = ((if_else(UpperHalf, CenterAngle * 180 / pi - 90, CenterAngle * 180 / pi + 90) + 180) %% 360) - 180
    )

  radar_df <- expand_grid(
    Group = factor(group_levels, levels = group_levels),
    PieCategory = factor(sector_levels, levels = sector_levels)
  ) %>%
    left_join(
      input_df %>%
        mutate(
          Group = factor(as.character(Group), levels = group_levels),
          PieCategory = factor(as.character(PieCategory), levels = sector_levels)
        ) %>%
        select(Group, PieCategory, Percent),
      by = c("Group", "PieCategory")
    ) %>%
    mutate(Percent = replace_na(Percent, 0)) %>%
    left_join(sector_df, by = "PieCategory") %>%
    mutate(
      GroupIndex = match(as.character(Group), group_levels),
      GroupSliceWidth = sector_width / length(group_levels),
      StartAngle = SectorStartAngle - (GroupIndex - 1) * GroupSliceWidth,
      EndAngle = SectorStartAngle - GroupIndex * GroupSliceWidth,
      CenterAngle = (StartAngle + EndAngle) / 2,
      Radius = Percent / radar_max,
      DisplayValue = if_else(Percent > 0, sprintf("%.2f", Percent), ""),
      ValueColor = mii_text_primary,
      UpperHalf = sin(CenterAngle) >= 0,
      TextAngle = ((if_else(UpperHalf, CenterAngle * 180 / pi - 90, CenterAngle * 180 / pi + 90) + 180) %% 360) - 180
    )

  list(sectors = sector_df, slices = radar_df)
}

# 组间共享雷达图的主绘制函数。
# 外层彩色环表示类别，内层不同颜色的扇形表示不同组在该类别上的占比。
plot_group_sector_radar <- function(radar_data, title_text, radar_max, radar_breaks, axis_label = "Relative Abundance (%)",
  label_radius = 1.24, value_radius = 0.965, coord_xlim = c(-1.92, 1.72), coord_ylim = c(-1.54, 1.54),
  label_size = 3.35, value_size = 2.7, sector_palette = radar_sector_palette) {
  sector_df <- radar_data$sectors
  slice_df <- radar_data$slices
  ring_inner <- 1.02
  ring_outer <- 1.16

  cdf <- build_circle_path_df(radar_breaks[radar_breaks > 0] / radar_max)
  wedge_df <- build_sector_polygon_df(
    slice_df %>% mutate(SectorId = paste0(PieCategory, "_", Group), PieCategory = Group),
    ir = 0,
    or = NA_real_,
    use_data = TRUE
  ) %>%
    mutate(FillColor = unname(group_palette[as.character(PieCategory)]))
  ring_df <- build_sector_polygon_df(
    sector_df %>% mutate(Radius = ring_outer),
    ir = ring_inner,
    or = ring_outer
  )
  ring_df <- ring_df %>% mutate(FillColor = unname(sector_palette[as.character(PieCategory)]))
  sdf <- tibble(
    x = 0,
    y = 0,
    xend = ring_outer * cos(sector_df$SectorStartAngle),
    yend = ring_outer * sin(sector_df$SectorStartAngle)
  )
  ldf <- build_sector_arc_text_df(sector_df, "DisplayLabel", label_radius, 2 * pi / nrow(sector_df), "LabelColor")
  vdf <- build_sector_text_df(slice_df, "DisplayValue", value_radius, "ValueColor")

  ggplot() +
    geom_path(data = cdf, aes(x = x, y = y, group = r), color = mii_neutral_line, linewidth = 0.55) +
    geom_polygon(data = wedge_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.38, linejoin = "round", alpha = 0.86) +
    geom_polygon(data = ring_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.45, linejoin = "round", alpha = 1) +
    geom_segment(data = sdf, aes(x = x, y = y, xend = xend, yend = yend), color = mii_neutral_line, linewidth = 0.75, lineend = "round") +
    geom_text(data = ldf, aes(x = TextX, y = TextY, label = Label, angle = TextAngle, color = TextColor), size = label_size, fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = vdf %>% filter(Label != ""), aes(x = TextX, y = TextY, label = Label, angle = TextAngle), size = value_size, fontface = "bold", inherit.aes = FALSE) +
    annotate("text", x = -1.60, y = 0, label = axis_label, angle = 90, size = 4.0, color = mii_axis_label_color, fontface = "bold") +
    scale_fill_identity(
      guide = "legend",
      breaks = unname(group_palette[group_levels]),
      labels = group_levels,
      name = NULL
    ) +
    scale_color_identity() +
    coord_fixed(xlim = coord_xlim, ylim = coord_ylim, clip = "off") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, color = mii_text_muted),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 20, 10, 28)
    ) +
    labs(title = title_text, fill = NULL)
}

# ===== 趋势图函数（alluvial/bar） =====

# 趋势图统一先走这里，把类别和分组补齐，避免堆积柱图/alluvial 因缺失类别而顺序漂移。
prepare_trend_plot_df <- function(df, levels, group_levels) {
  if (nrow(df) == 0) {
    return(tibble(Group = factor(character(), levels = group_levels), TrendCategory = factor(character(), levels = levels), Reads = numeric(), Percent = numeric()))
  }
  df %>% mutate(Group = factor(as.character(Group), levels = group_levels), TrendCategory = factor(as.character(TrendCategory), levels = levels)) %>%
    filter(!is.na(TrendCategory)) %>%
    complete(Group = factor(group_levels, levels = group_levels), TrendCategory = factor(levels, levels = levels), fill = list(Reads = 0, Percent = 0)) %>%
    arrange(Group, TrendCategory)
}

# 堆积柱图：适合看各组内部组成。
build_trend_bar_plot <- function(df, levels, value_col, title, y_label, group_levels, fill_palette = main_palette) {
  plot_df <- prepare_trend_plot_df(df, levels, group_levels)
  display_df <- plot_df %>% filter(.data[[value_col]] > 0)
  if (nrow(display_df) == 0) return(plot_note(title, "No trend data available."))
  p <- ggplot(display_df, aes(x = Group, y = .data[[value_col]], fill = TrendCategory)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = fill_palette[levels], drop = TRUE) +
    labs(x = "Group", y = y_label, title = title) +
    mii_theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE))
  if (identical(value_col, "Reads")) p <- p + scale_y_continuous(labels = scales::comma)
  p
}

# 将 p value 转成图上常用的显著性星号。
pvalue_to_stars <- function(pvalue) {
  case_when(
    is.na(pvalue) ~ NA_character_,
    pvalue < 0.0001 ~ "****",
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

format_pvalue_label <- function(pvalue, show_pvalue = FALSE) {
  stars <- pvalue_to_stars(pvalue)
  if (is.na(stars)) return(NA_character_)
  if (!show_pvalue || is.na(pvalue)) return(stars)
  paste0(stars, "\n", format.pval(pvalue, digits = 2, eps = 0.001))
}

calculate_two_group_pvalue <- function(df, value_col = "Value", group_col = "Group", group_levels = NULL, min_replicates = 2) {
  if (Sys.getenv("RUN_INFERENCE","0") != "1" || !smallrna_inference_gate()$eligible) return(NA_real_)
  if (!"Sample" %in% names(df) || anyDuplicated(df$Sample)) return(NA_real_)
  min_replicates <- max(3,min_replicates)
  if (nrow(df) == 0 || !all(c(value_col, group_col) %in% names(df))) return(NA_real_)
  stat_df <- df %>%
    mutate(
      .group = as.character(.data[[group_col]]),
      .value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    filter(!is.na(.group), is.finite(.value))

  groups <- if (is.null(group_levels)) unique(stat_df$.group) else as.character(group_levels)
  groups <- groups[groups %in% unique(stat_df$.group)]
  if (length(groups) != 2) return(NA_real_)

  x <- stat_df$.value[stat_df$.group == groups[[1]]]
  y <- stat_df$.value[stat_df$.group == groups[[2]]]
  if (length(x) < min_replicates || length(y) < min_replicates) return(NA_real_)
  if (length(unique(c(x, y))) < 2) return(1)

  tryCatch(
    suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

complete_sample_metric <- function(df, value_col, sample_info, group_levels) {
  base_df <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group))

  if (nrow(df) == 0 || !value_col %in% names(df)) {
    return(base_df %>% mutate(Value = 0, Group = factor(Group, levels = group_levels)))
  }

  metric_df <- df %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    group_by(Sample, Group) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

  base_df %>%
    left_join(metric_df, by = c("Sample", "Group")) %>%
    mutate(
      Value = replace_na(Value, 0),
      Group = factor(Group, levels = group_levels)
    )
}

build_two_group_summary_bar <- function(sample_df, value_col, title, y_label, group_levels, fill_palette = group_palette,
                                        sample_info = sample_info, comma_labels = FALSE, label_digits = 1,
                                        show_pvalue = FALSE) {
  plot_df <- complete_sample_metric(sample_df, value_col, sample_info = sample_info, group_levels = group_levels)
  if (nrow(plot_df) == 0) return(plot_note(title, "No sample-level data available."))

  summary_df <- plot_df %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = stats::sd(Value, na.rm = TRUE),
      N = sum(is.finite(Value)),
      SEM = if_else(N > 1, SD / sqrt(N), 0),
      .groups = "drop"
    ) %>%
    mutate(
      SD = replace_na(SD, 0),
      SEM = replace_na(SEM, 0),
      Group = factor(as.character(Group), levels = group_levels)
    )

  max_y <- max(c(summary_df$Mean + summary_df$SEM, plot_df$Value), na.rm = TRUE)
  if (!is.finite(max_y)) max_y <- 0
  y_pad <- if (max_y > 0) max_y * 0.12 else 1
  bracket_y <- max_y + y_pad
  tip_y <- max_y + y_pad * 0.45
  label_y <- bracket_y + y_pad * 0.35
  upper_y <- label_y + y_pad * 0.45
  pvalue <- calculate_two_group_pvalue(plot_df, value_col = "Value", group_col = "Group", group_levels = group_levels)
  p_label <- format_pvalue_label(pvalue, show_pvalue = show_pvalue)

  bar_plot <- ggplot(summary_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.64, color = "white", linewidth = 0.25) +
    geom_errorbar(aes(ymin = pmax(Mean - SEM, 0), ymax = Mean + SEM), width = 0.16, linewidth = 0.45, color = mii_text_primary) +
    geom_point(
      data = plot_df,
      aes(x = Group, y = Value),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.065, height = 0, seed = 11),
      size = 2.2,
      shape = 21,
      stroke = 0.25,
      color = mii_text_primary,
      fill = "white",
      alpha = 0.9
    ) +
    geom_text(
      aes(label = sprintf(paste0("%.", label_digits, "f"), Mean), y = Mean + SEM + y_pad * 0.12),
      size = 3.2,
      color = mii_text_primary,
      vjust = 0
    ) +
    scale_fill_manual(values = fill_palette, drop = FALSE) +
    scale_y_continuous(
      labels = if (comma_labels) scales::comma else waiver(),
      limits = c(0, upper_y),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = y_label, title = title) +
    mii_theme_bw() +
    theme(legend.position = "none")

  if (length(group_levels) == 2 && !is.na(p_label)) {
    bracket_df <- tibble(
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      label = p_label,
      label_x = 1.5,
      label_y = label_y
    )
    tick_df <- tibble(x = c(1, 2), xend = c(1, 2), y = c(tip_y, tip_y), yend = c(bracket_y, bracket_y))
    bar_plot <- bar_plot +
      geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.45) +
      geom_segment(data = tick_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.45) +
      geom_text(data = bracket_df, aes(x = label_x, y = label_y, label = label), inherit.aes = FALSE, size = 4.0, fontface = "bold")
  }

  bar_plot
}

# alluvial：适合看类别在组与组之间的相对流动关系。
build_trend_alluvial_plot <- function(df, levels, value_col, title, y_label, group_levels, fill_palette = main_palette) {
  plot_df <- prepare_trend_plot_df(df, levels, group_levels)
  display_df <- plot_df %>% filter(.data[[value_col]] > 0)
  if (!pkg_has("ggalluvial")) return(plot_note(title, "ggalluvial package is not available."))
  if (nrow(display_df) == 0) return(plot_note(title, "No trend data available."))
  # 这里显式压缩组轴间距，让“柱宽 = 柱间距 = 0.1”。
  # flow 宽度与 stratum 宽度保持一致，这样流带会从柱子的边界出发，到下一根柱子的边界结束。
  stratum_width <- 0.1
  flow_width <- stratum_width
  group_axis_spacing <- stratum_width * 2
  group_axis_positions <- seq(1, by = group_axis_spacing, length.out = length(group_levels))
  names(group_axis_positions) <- group_levels
  display_df <- display_df %>%
    mutate(GroupX = unname(group_axis_positions[as.character(Group)]))

  left_group <- group_levels[[1]]
  if (identical(value_col, "Reads")) {
    label_layer <- geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(scales::comma(count))),
      size = 3,
      color = mii_text_primary,
      check_overlap = TRUE
    )
    axis_label_layer <- NULL
    y_scale <- scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0)))
  } else {
    label_layer <- geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(ifelse(count > 3.5, sprintf("%.1f%%", count), ""))),
      size = 3,
      color = mii_text_primary,
      check_overlap = FALSE
    )
    # 左侧柱子上的刻度标注 0/25/50/75/100
    axis_label_df <- data.frame(
      x = group_axis_positions[[1]],
      y = c(0, 25, 50, 75, 100),
      label = c("0", "25", "50", "75", "100")
    )
    axis_label_layer <- geom_text(
      data = axis_label_df,
      aes(x = x - stratum_width / 2 - 0.02, y = y, label = label),
      size = 3, color = mii_text_primary, hjust = 1, inherit.aes = FALSE
    )
    y_scale <- scale_y_continuous(
      breaks = seq(0, 100, 25),
      labels = NULL,
      expand = expansion(mult = c(0, 0))
    )
  }

  p <- ggplot(display_df, aes(x = GroupX, y = .data[[value_col]], stratum = TrendCategory, alluvium = TrendCategory, fill = TrendCategory)) +
    ggalluvial::geom_flow(width = flow_width, alpha = 0.32) +
    ggalluvial::geom_stratum(width = stratum_width, color = "white", linewidth = 0.25) +
    label_layer +
    axis_label_layer +
    scale_fill_manual(values = fill_palette[levels], drop = TRUE) +
    scale_x_continuous(
      breaks = unname(group_axis_positions),
      labels = names(group_axis_positions),
      expand = expansion(add = stratum_width / 2)
    ) +
    y_scale +
    labs(x = "Group", y = y_label, title = title) +
    mii_theme_bw() + theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE))
  p
}

# ===== 百分比饼图函数 =====

# 通用百分比饼图。
# 很多单样本 origin / cleavage 小图都直接复用这一个函数。
plot_percent_pie <- function(df, category_col, palette, title_text, empty_note) {
  df <- df %>% filter(.data[[category_col]] != "Unknown", Reads > 0 | Percent > 0)
  if (nrow(df) == 0) return(plot_note(title_text, empty_note))
  df <- df %>% mutate(!!category_col := factor(.data[[category_col]], levels = names(palette)))
  ggplot(df, aes(x = "", y = Percent, fill = .data[[category_col]])) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = palette, drop = FALSE) +
    geom_text(aes(label = ifelse(Percent >= 1, sprintf("%.1f%%", Percent), "")), position = position_stack(vjust = 0.5), size = 3.2) +
    labs(x = NULL, y = NULL, title = title_text, fill = NULL) +
    theme_void() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
}

# 对行做稳健标准化，常用于热图。
# 当某些行方差为 0 或只有一个值时，直接回填 0，避免 scale() 产生 NA。
safe_scale_rows <- function(x) {
  x <- as.matrix(x)
  if (length(x) == 0 || nrow(x) == 0 || ncol(x) == 0) return(x)
  scaled <- tryCatch(t(scale(t(x))), error = function(e) matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

# 把长表里的特征计数整理成宽矩阵，供 DE、相关性和热图复用。
# feature_col / sample_col / count_col 都用列名字符串传入。
build_feature_count_matrix <- function(df, feature_col, sample_col = "Sample", count_col = "Reads", sample_order = NULL, feature_order = NULL, min_total = 0) {
  if (is.null(df) || nrow(df) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = character(), feature_order = character()))
  }

  required_cols <- c(feature_col, sample_col, count_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns for build_feature_count_matrix(): ", paste(missing_cols, collapse = ", "))
  }

  count_tbl <- df %>%
    transmute(
      Sample = as.character(.data[[sample_col]]),
      Feature = as.character(.data[[feature_col]]),
      Reads = clean_numeric(.data[[count_col]])
    ) %>%
    filter(!is.na(Sample), Sample != "", !is.na(Feature), Feature != "", !is.na(Reads)) %>%
    group_by(Sample, Feature) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop")

  if (nrow(count_tbl) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = character(), feature_order = character()))
  }

  if (is.null(sample_order)) {
    sample_order <- count_tbl %>% distinct(Sample) %>% pull(Sample)
  } else {
    sample_order <- unique(as.character(sample_order))
  }

  if (is.null(feature_order)) {
    feature_order <- count_tbl %>%
      group_by(Feature) %>%
      summarise(TotalReads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
      filter(TotalReads > min_total) %>%
      arrange(desc(TotalReads), Feature) %>%
      pull(Feature)
  } else {
    feature_order <- unique(as.character(feature_order))
  }

  if (length(sample_order) == 0 || length(feature_order) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = sample_order, feature_order = feature_order))
  }

  full_tbl <- tidyr::expand_grid(Sample = sample_order, Feature = feature_order) %>%
    left_join(count_tbl, by = c("Sample", "Feature")) %>%
    mutate(Reads = replace_na(Reads, 0), Sample = factor(Sample, levels = sample_order), Feature = factor(Feature, levels = feature_order))

  count_matrix <- xtabs(Reads ~ Feature + Sample, data = full_tbl)
  count_matrix <- as.matrix(count_matrix)
  count_matrix <- count_matrix[feature_order, sample_order, drop = FALSE]

  list(
    matrix = count_matrix,
    table = full_tbl %>% arrange(Feature, Sample),
    sample_order = sample_order,
    feature_order = feature_order
  )
}

# 将 pairwise contrast 规格统一整理成可循环的数据表。
# 默认情况下会自动生成所有两两比较；也支持手工传入两列数据框。
resolve_pairwise_contrasts <- function(group_levels, contrast_specs = NULL) {
  group_levels <- unique(as.character(group_levels))
  if (length(group_levels) < 2) return(tibble(group1 = character(), group2 = character(), contrast = character()))

  if (is.null(contrast_specs)) {
    pairs <- utils::combn(group_levels, 2, simplify = FALSE)
    return(map_dfr(pairs, function(pair) tibble(group1 = pair[[1]], group2 = pair[[2]], contrast = paste0(pair[[2]], "_vs_", pair[[1]]))))
  }

  if (is.data.frame(contrast_specs)) {
    if (!all(c("group1", "group2") %in% names(contrast_specs))) {
      stop("contrast_specs data frame must contain group1 and group2 columns.")
    }
    return(as_tibble(contrast_specs) %>% transmute(group1 = as.character(group1), group2 = as.character(group2), contrast = if ("contrast" %in% names(.)) as.character(contrast) else paste0(group2, "_vs_", group1)))
  }

  if (is.list(contrast_specs)) {
    return(map_dfr(contrast_specs, function(pair) {
      pair <- as.character(pair)
      if (length(pair) < 2) return(tibble())
      tibble(group1 = pair[[1]], group2 = pair[[2]], contrast = if (length(pair) >= 3 && nzchar(pair[[3]])) pair[[3]] else paste0(pair[[2]], "_vs_", pair[[1]]))
    }))
  }

  if (is.character(contrast_specs)) {
    if (length(contrast_specs) == 2) {
      return(tibble(group1 = contrast_specs[[1]], group2 = contrast_specs[[2]], contrast = paste0(contrast_specs[[2]], "_vs_", contrast_specs[[1]])))
    }
  }

  stop("Unsupported contrast_specs format.")
}

# 统一定义 DE 显著性方向，后续火山图、摘要表和 top feature heatmap 都会用到。
classify_de_direction <- function(pvalue, log_fc, pvalue_cutoff = 0.05, lfc_cutoff = log2(1.5)) {
  case_when(
    is.na(pvalue) | is.na(log_fc) ~ "ns",
    pvalue < pvalue_cutoff & log_fc >= lfc_cutoff ~ "up",
    pvalue < pvalue_cutoff & log_fc <= -lfc_cutoff ~ "down",
    TRUE ~ "ns"
  )
}

de_significance_label <- function() {
  if (DE_SIGNIFICANCE %in% c("PVALUE", "P_VALUE", "P")) "PValue" else "FDR"
}

resolve_de_significance_values <- function(de_tbl) {
  if ("SignificanceMetric" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$SignificanceMetric))
  }
  preferred_col <- de_significance_label()
  if (preferred_col %in% names(de_tbl)) {
    return(as.numeric(de_tbl[[preferred_col]]))
  }
  if ("FDR" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$FDR))
  }
  if ("PValue" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$PValue))
  }
  rep(NA_real_, nrow(de_tbl))
}

# 计算一组 pairwise contrast 的 edgeR 差异结果。
# 若缺少足够重复或 edgeR 不可用，返回 skipped = TRUE 的结构，不强制报错。
run_edgeR_pairwise_de <- function(count_matrix, sample_groups, group1, group2, feature_ids = NULL, min_replicates = 3) {
  skipped_result <- function(message) list(skipped=TRUE, message=message, result=tibble(), sample_groups=sample_groups, feature_ids=feature_ids)
  if (Sys.getenv("RUN_INFERENCE", "0") != "1") return(skipped_result("disabled: descriptive default; no condition inference"))
  gate <- smallrna_inference_gate()
  if (!gate$eligible) return(skipped_result(paste("NO-GO:",gate$reason)))
  meta <- load_metadata()
  if (is.null(colnames(count_matrix)) || !all(colnames(count_matrix) %in% as.character(meta$Sample))) stop("Count matrix sample IDs must match explicit sample table")
  expected_groups <- as.character(meta$Group[match(colnames(count_matrix),as.character(meta$Sample))])
  if (!identical(as.character(sample_groups),expected_groups)) stop("count_matrix/sample_groups direction mismatch")
  if (any(!is.finite(count_matrix)) || any(count_matrix<0) || any(count_matrix!=round(count_matrix))) return(skipped_result("fractional/non-integer counts are descriptive only; no silent rounding"))
  min_replicates <- max(3,min_replicates)
  if (!pkg_has("edgeR")) {
    return(list(
      skipped = TRUE,
      message = "edgeR package is not available.",
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  sample_groups <- as.character(sample_groups)
  keep_samples <- sample_groups %in% c(group1, group2)
  if (sum(keep_samples) < 4) {
    return(list(
      skipped = TRUE,
      message = paste0("Not enough samples for contrast ", group2, " vs ", group1, "."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  count_matrix <- as.matrix(count_matrix[, keep_samples, drop = FALSE])
  sample_groups <- sample_groups[keep_samples]
  if (min(table(sample_groups)) < min_replicates) {
    return(list(
      skipped = TRUE,
      message = paste0("Contrast ", group2, " vs ", group1, " needs at least ", min_replicates, " replicates per group."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  group_factor <- factor(sample_groups, levels = c(group1, group2))
  base_mean_values <- rowMeans(count_matrix, na.rm = TRUE)
  base_mean_names <- rownames(count_matrix)
  if (is.null(base_mean_names) || anyNA(base_mean_names) || !all(nzchar(base_mean_names))) {
    if (!is.null(feature_ids) && length(feature_ids) == nrow(count_matrix)) {
      base_mean_names <- as.character(feature_ids)
    } else {
      base_mean_names <- as.character(seq_len(nrow(count_matrix)))
    }
  }
  names(base_mean_values) <- base_mean_names
  dge <- edgeR::DGEList(counts = count_matrix, group = group_factor)
  keep <- as.vector(edgeR::filterByExpr(dge, group = group_factor))
  dge <- edgeR::DGEList(counts = count_matrix[keep, , drop = FALSE], group = group_factor)

  if (nrow(dge) == 0) {
    return(list(
      skipped = TRUE,
      message = paste0("No features passed expression filtering for contrast ", group2, " vs ", group1, "."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ 0 + group_factor)
  safe_group_names <- setNames(make.names(levels(group_factor)), levels(group_factor))
  colnames(design) <- unname(safe_group_names)
  de_method <- "edgeR_robust_glmQLFTest"
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  contrast <- limma::makeContrasts(
    contrasts = paste0(safe_group_names[[as.character(group2)]], "-", safe_group_names[[as.character(group1)]]),
    levels = design
  )
  test <- edgeR::glmQLFTest(fit, contrast = contrast)
  de_table <- edgeR::topTags(test, n = Inf)$table
  result <- de_table
  result <- tibble::rownames_to_column(as.data.frame(result), var = "Feature")

  if (!is.null(feature_ids) && length(feature_ids) == nrow(count_matrix)) {
    result <- result %>% left_join(tibble(Feature = as.character(feature_ids), FeatureID = as.character(feature_ids)), by = "Feature")
  }

  result <- result %>%
    mutate(baseMean = unname(base_mean_values[Feature])) %>%
    mutate(
      SignificanceMetric = if (de_significance_label() == "FDR") FDR else PValue,
      SignificanceMethod = de_significance_label(),
      DEMethod = de_method,
      direction = classify_de_direction(SignificanceMetric, logFC),
      rank_score = sign(coalesce(logFC, 0)) * -log10(coalesce(SignificanceMetric, 1) + 1e-300),
      contrast = paste0(group2, "_vs_", group1)
    ) %>%
    arrange(SignificanceMetric, PValue, desc(abs(logFC)))

  list(
    skipped = FALSE,
    message = NULL,
    result = result,
    sample_groups = sample_groups,
    feature_ids = feature_ids
  )
}

# 生成 DE 火山图，必要时自动给最显著的条目加标签。
build_de_volcano_plot <- function(de_tbl, title_text, lfc_cutoff = log2(1.5), pvalue_cutoff = 0.05, label_n = 15) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No DE results available."))

  significance_values <- resolve_de_significance_values(de_tbl)
  significance_methods <- if ("SignificanceMethod" %in% names(de_tbl)) as.character(de_tbl$SignificanceMethod) else character(0)
  significance_methods <- significance_methods[!is.na(significance_methods) & nzchar(significance_methods)]
  significance_label <- if (length(significance_methods) > 0) {
    unique(significance_methods)[[1]]
  } else {
    de_significance_label()
  }

  plot_df <- de_tbl %>%
    mutate(
      SignificanceMetric = significance_values,
      direction = classify_de_direction(SignificanceMetric, logFC, pvalue_cutoff = pvalue_cutoff, lfc_cutoff = lfc_cutoff),
      label_score = abs(coalesce(logFC, 0)) * -log10(coalesce(SignificanceMetric, 1) + 1e-300),
      label = if_else(direction != "ns" & rank(-label_score, ties.method = "first") <= label_n, as.character(Feature), NA_character_)
    )

  p <- ggplot(plot_df, aes(x = logFC, y = -log10(SignificanceMetric + 1e-300), color = direction)) +
    geom_point(alpha = 0.75, size = 1.6) +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = mii_neutral_text) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = mii_neutral_text) +
    scale_color_manual(values = mii_direction_palette, drop = FALSE) +
    labs(x = "log2 Fold Change", y = paste0("-log10(", significance_label, ")"), title = title_text) +
    mii_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  if (pkg_has("ggrepel")) {
    p <- p + ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, na.rm = TRUE, show.legend = FALSE)
  } else {
    p <- p + geom_text(aes(label = label), size = 3, na.rm = TRUE, show.legend = FALSE, vjust = -0.4)
  }

  p
}

# 生成 DE MA 图。
build_de_ma_plot <- function(de_tbl, title_text) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No DE results available."))

  ggplot(de_tbl, aes(x = log10(baseMean + 1), y = logFC, color = direction)) +
    geom_point(alpha = 0.72, size = 1.5) +
    scale_color_manual(values = mii_direction_palette, drop = FALSE) +
    labs(x = "log10(baseMean + 1)", y = "log2 Fold Change", title = title_text) +
    mii_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")
}

build_de_robustness_table <- function(result_tbl, count_matrix, sample_order, sample_groups, group1, group2) {
  if (nrow(result_tbl) == 0 || nrow(count_matrix) == 0 || length(sample_order) == 0) return(tibble())
  sample_order <- as.character(sample_order)
  sample_groups <- as.character(sample_groups)
  if (!all(sample_order %in% colnames(count_matrix))) return(tibble())

  top_features <- result_tbl %>%
    arrange(SignificanceMetric, PValue, desc(abs(logFC))) %>%
    pull(Feature) %>%
    as.character() %>%
    unique()
  top_features <- top_features[top_features %in% rownames(count_matrix)]
  if (length(top_features) == 0) return(tibble())

  baseline_samples <- sample_order[sample_groups == group1]
  case_samples <- sample_order[sample_groups == group2]
  if (length(baseline_samples) == 0 || length(case_samples) == 0) return(tibble())

  purrr::map_dfr(top_features, function(feature_name) {
    values <- as.numeric(count_matrix[feature_name, sample_order, drop = TRUE])
    names(values) <- sample_order
    baseline_values <- values[baseline_samples]
    case_values <- values[case_samples]
    baseline_total <- sum(baseline_values, na.rm = TRUE)
    case_total <- sum(case_values, na.rm = TRUE)
    all_total <- baseline_total + case_total
    max_sample <- names(values)[which.max(values)]
    max_reads <- max(values, na.rm = TRUE)
    max_group <- if (max_sample %in% baseline_samples) group1 else group2
    max_baseline_fraction <- if (baseline_total > 0) max(baseline_values, na.rm = TRUE) / baseline_total else NA_real_
    max_case_fraction <- if (case_total > 0) max(case_values, na.rm = TRUE) / case_total else NA_real_
    max_group_fraction <- max(c(max_baseline_fraction, max_case_fraction), na.rm = TRUE)
    if (!is.finite(max_group_fraction)) max_group_fraction <- NA_real_

    de_row <- result_tbl %>% filter(as.character(Feature) == feature_name) %>% slice(1)
    direction <- as.character(de_row$direction[[1]])
    baseline_median <- stats::median(baseline_values, na.rm = TRUE)
    case_median <- stats::median(case_values, na.rm = TRUE)
    case_support_fraction <- case_when(
      direction == "up" ~ mean(case_values > baseline_median, na.rm = TRUE),
      direction == "down" ~ mean(case_values < baseline_median, na.rm = TRUE),
      TRUE ~ NA_real_
    )
    baseline_support_fraction <- case_when(
      direction == "up" ~ mean(baseline_values < case_median, na.rm = TRUE),
      direction == "down" ~ mean(baseline_values > case_median, na.rm = TRUE),
      TRUE ~ NA_real_
    )
    single_sample_dominated <- is.finite(max_group_fraction) & max_group_fraction >= 0.80
    weak_replicate_support <- direction %in% c("up", "down") & is.finite(case_support_fraction) & case_support_fraction < 0.50
    robustness_flag <- case_when(
      direction == "ns" ~ "not_significant",
      single_sample_dominated ~ "single_sample_dominated",
      weak_replicate_support ~ "weak_replicate_support",
      TRUE ~ "replicate_supported"
    )

    tibble(
      Feature = feature_name,
      direction = direction,
      logFC = de_row$logFC[[1]],
      SignificanceMetric = de_row$SignificanceMetric[[1]],
      PValue = de_row$PValue[[1]],
      FDR = if ("FDR" %in% names(de_row)) de_row$FDR[[1]] else NA_real_,
      BaselineGroup = group1,
      CaseGroup = group2,
      BaselineTotalReads = baseline_total,
      CaseTotalReads = case_total,
      TotalReads = all_total,
      BaselineDetectedSamples = sum(baseline_values > 0, na.rm = TRUE),
      CaseDetectedSamples = sum(case_values > 0, na.rm = TRUE),
      BaselineSampleN = length(baseline_samples),
      CaseSampleN = length(case_samples),
      BaselineMedianReads = baseline_median,
      CaseMedianReads = case_median,
      CaseDirectionSupportFraction = case_support_fraction,
      BaselineDirectionSupportFraction = baseline_support_fraction,
      MaxSample = max_sample,
      MaxSampleGroup = max_group,
      MaxSampleReads = max_reads,
      MaxSampleFractionTotal = if_else(all_total > 0, max_reads / all_total, NA_real_),
      MaxBaselineSampleFraction = max_baseline_fraction,
      MaxCaseSampleFraction = max_case_fraction,
      MaxGroupSampleFraction = max_group_fraction,
      SingleSampleDominated = single_sample_dominated,
      RobustnessFlag = robustness_flag,
      BaselineSampleReads = paste(paste0(baseline_samples, ":", baseline_values), collapse = ";"),
      CaseSampleReads = paste(paste0(case_samples, ":", case_values), collapse = ";")
    )
  })
}

build_de_robustness_plot <- function(robustness_tbl, title_text, label_n = 30) {
  if (nrow(robustness_tbl) == 0) return(plot_note(title_text, "No robustness table available."))
  plot_df <- robustness_tbl %>%
    arrange(desc(MaxGroupSampleFraction), SignificanceMetric) %>%
    slice_head(n = label_n) %>%
    mutate(
      Feature = factor(Feature, levels = rev(unique(Feature))),
      RobustnessFlag = factor(
        RobustnessFlag,
        levels = c("replicate_supported", "weak_replicate_support", "single_sample_dominated", "not_significant")
      )
    )

  ggplot(plot_df, aes(x = Feature, y = MaxGroupSampleFraction, fill = RobustnessFlag)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.15) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(
      values = c(
        replicate_supported = mii_green,
        weak_replicate_support = mii_orange,
        single_sample_dominated = mii_red,
        not_significant = mii_neutral_line
      ),
      drop = FALSE
    ) +
    labs(title = title_text, x = NULL, y = "Largest single-sample fraction within a group", fill = "Flag") +
    mii_theme_bw() +
    theme(legend.position = "top")
}

# ===== 对象中心 DE 共享入口 =====
resolve_mii_contrast_specs <- function(sample_info_tbl = sample_info, group_levels_vec = group_levels) {
  resolve_pairwise_contrasts(group_levels_vec, tibble(group1=PROJECT_CONTRAST_BASELINE, group2=PROJECT_CONTRAST_CASE,
    contrast=paste0(PROJECT_CONTRAST_CASE,"_vs_",PROJECT_CONTRAST_BASELINE)))
}

write_de_run_config <- function(contrast_specs, section = "00_de_config") {
  write_csv_safe(contrast_specs, file.path(section, "contrast_manifest.csv"))
  write_csv_safe(
    tibble(
      significance_method = de_significance_label(),
      significance_cutoff = 0.05,
      log2_fold_change_cutoff = log2(1.5),
      rule = paste0("direction != ns requires ", de_significance_label(), " < 0.05 and |logFC| >= log2(1.5)")
    ),
    file.path(section, "de_threshold_manifest.csv")
  )
  invisible(contrast_specs)
}

build_global_de_inputs <- function() {
  general_annotation_categories <- c("miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "ribozyme", "other")
  general_annotation_categories <- general_annotation_categories[general_annotation_categories %in% category_levels]
  if (length(general_annotation_categories) == 0) {
    general_annotation_categories <- c("miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "ribozyme", "other")
  }

  category_input <- summary_top %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), Feature = as.character(MainCategory), Reads = Reads)

  subtype_input <- summary_all %>%
    filter(SubClass != "-", !MainCategory %in% c("tsRNA", "rsRNA"), !Class %in% c("Clean_Reads", "Match_Genome")) %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), Feature = paste0(as.character(MainCategory), "::", as.character(SubClass)), Reads = Reads)

  annotation_input <- output_all %>%
    filter(
      FeatureCategory %in% general_annotation_categories,
      !is.na(PrimaryAnnotation),
      PrimaryAnnotation != "",
      PrimaryAnnotation != "NO_Annotation"
    ) %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = as.character(PrimaryAnnotation),
      FeatureCategory = as.character(FeatureCategory),
      Reads = Reads
    )

  list(
    category_input = category_input,
    subtype_input = subtype_input,
    annotation_input = annotation_input,
    category_meta = tibble(Feature = as.character(category_levels), FeatureCategory = as.character(category_levels)),
    subtype_meta = subtype_input %>% distinct(Feature) %>% mutate(FeatureCategory = sub("::.*$", "", Feature)),
    annotation_meta = annotation_input %>% distinct(Feature, FeatureCategory)
  )
}

build_annotation_de_input <- function(category_name) {
  output_all %>%
    filter(
      FeatureCategory == category_name,
      !is.na(PrimaryAnnotation),
      PrimaryAnnotation != "",
      PrimaryAnnotation != "NO_Annotation"
    ) %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = as.character(PrimaryAnnotation),
      FeatureCategory = as.character(FeatureCategory),
      Reads = Reads
    )
}

build_tsrna_annotation_de_input <- function() {
  tsrna_cleavage_detail %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = make_company_tsrna_feature(SubClass, Cleavage),
      FeatureCategory = if_else(Origin == "Mitochondrial", "mt_tsRNA", "tsRNA"),
      Reads = Reads
    )
}

build_rsrna_subtype_de_input <- function() {
  rsrna_detail %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = as.character(RsRNAClass),
      FeatureCategory = "rsRNA",
      Reads = Reads
    )
}

build_sequence_de_base <- function() {
  output_all %>%
    filter(!is.na(Sequence), Sequence != "", !is.na(Reads), Reads > 0) %>%
    mutate(
      Feature = as.character(Sequence),
      AnnotationLabel = if_else(
        is.na(PrimaryAnnotation) | PrimaryAnnotation == "" | PrimaryAnnotation == "NO_Annotation",
        "unannotated",
        as.character(PrimaryAnnotation)
      ),
      SequenceCategory = as.character(FeatureCategory),
      IsMito = is_mito_label(paste(Annotation, PrimaryAnnotation, sep = " "))
    ) %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature,
      Reads,
      AnnotationLabel,
      SequenceCategory,
      IsMito
    )
}

sequence_de_input_for <- function(sequence_de_base, category_name) {
  switch(
    category_name,
    miRNA = sequence_de_base %>% filter(SequenceCategory == "miRNA"),
    tsRNA = sequence_de_base %>% filter(SequenceCategory == "tsRNA", !IsMito),
    mt_tsRNA = sequence_de_base %>% filter(SequenceCategory == "tsRNA", IsMito),
    piRNA = sequence_de_base %>% filter(SequenceCategory == "piRNA"),
    rsRNA = sequence_de_base %>% filter(SequenceCategory == "rsRNA"),
    tibble()
  )
}

run_sequence_de_categories <- function(categories, layer_prefix = "sequence_de", sequence_de_base = build_sequence_de_base()) {
  categories <- unique(as.character(categories))
  category_labels <- c(miRNA = "miRNA", tsRNA = "tsRNA", mt_tsRNA = "mt-tsRNA", piRNA = "piRNA", rsRNA = "rsRNA")
  runs <- purrr::map(categories, function(category) {
    seq_input_raw <- sequence_de_input_for(sequence_de_base, category)
    seq_input <- seq_input_raw %>%
      group_by(Sample, Group, Feature) %>%
      summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop")

    seq_meta <- seq_input_raw %>%
      group_by(Feature) %>%
      summarise(
        FeatureCategory = category,
        TopAnnotation = {
          ann_counts <- sort(table(AnnotationLabel), decreasing = TRUE)
          if (length(ann_counts) == 0) "unannotated" else names(ann_counts)[[1]]
        },
        TotalReads = sum(Reads, na.rm = TRUE),
        .groups = "drop"
      )

    run_layer_de(
      layer_name = file.path(layer_prefix, paste0(category, "_seq")),
      layer_title = paste0(category_labels[[category]] %||% category, " sequence-level DE"),
      layer_input = seq_input,
      feature_order = NULL,
      feature_meta = seq_meta,
      min_total = 5
    )
  })
  names(runs) <- categories

  list(
    runs = runs,
    results = bind_rows(map(runs, "results")),
    summary = bind_rows(map(runs, "summary")),
    manifest = bind_rows(map(runs, "manifest")),
    robustness_summary = bind_rows(map(runs, "robustness_summary"))
  )
}

summarise_de_robustness_flags <- function(robustness_tbl) {
  if (is.null(robustness_tbl) || nrow(robustness_tbl) == 0 || !"RobustnessFlag" %in% names(robustness_tbl)) {
    return(tibble())
  }

  robustness_tbl %>%
    mutate(
      RobustnessFlag = as.character(RobustnessFlag),
      direction = as.character(direction),
      is_significant = direction != "ns"
    ) %>%
    group_by(layer, layer_title, contrast, RobustnessFlag) %>%
    summarise(
      FeatureN = n(),
      SignificantFeatureN = sum(is_significant, na.rm = TRUE),
      UpN = sum(direction == "up", na.rm = TRUE),
      DownN = sum(direction == "down", na.rm = TRUE),
      MedianMaxGroupSampleFraction = median(MaxGroupSampleFraction, na.rm = TRUE),
      MaxMaxGroupSampleFraction = max(MaxGroupSampleFraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(layer, contrast, desc(SignificantFeatureN), RobustnessFlag)
}

run_layer_de <- function(layer_name, layer_title, layer_input, feature_order = NULL, feature_meta = NULL, min_total = 0) {
  layer_results <- list()
  layer_summary <- list()
  layer_manifest <- list()
  layer_robustness <- list()
  layer_robustness_summary <- list()

  for (i in seq_len(nrow(contrast_specs))) {
    group1 <- as.character(contrast_specs$group1[[i]])
    group2 <- as.character(contrast_specs$group2[[i]])
    contrast_label <- as.character(contrast_specs$contrast[[i]])
    pair_samples <- sample_info %>%
      filter(as.character(Group) %in% c(group1, group2)) %>%
      mutate(Group = factor(as.character(Group), levels = c(group1, group2))) %>%
      arrange(Group, Sample)

    pair_sample_order <- as.character(pair_samples$Sample)
    pair_sample_groups <- as.character(pair_samples$Group)
    pair_section <- file.path(layer_name, safe_dir_name(contrast_label))

    layer_subset <- layer_input %>%
      filter(Sample %in% pair_sample_order) %>%
      select(Sample, Group, Feature, Reads)

    matrix_info <- build_feature_count_matrix(
      layer_subset,
      feature_col = "Feature",
      sample_order = pair_sample_order,
      feature_order = feature_order,
      min_total = min_total
    )

    if (nrow(matrix_info$matrix) == 0 || ncol(matrix_info$matrix) == 0) {
      skip_message <- paste0("No features available for ", layer_title, " in contrast ", contrast_label, ".")
      layer_manifest[[length(layer_manifest) + 1]] <- tibble(
        layer = layer_name,
        contrast = contrast_label,
        group1 = group1,
        group2 = group2,
        status = "skipped",
        message = skip_message,
        sample_n = length(pair_sample_order),
        feature_n = 0,
        significant_n = 0,
        up_n = 0,
        down_n = 0
      )
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), skip_message),
        pair_section,
        paste0("skipped_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
      next
    }

    feature_lookup <- tibble(
      FeatureInternal = sprintf("feature_%06d", seq_len(nrow(matrix_info$matrix))),
      Feature = rownames(matrix_info$matrix)
    )
    de_matrix <- matrix_info$matrix
    rownames(de_matrix) <- feature_lookup$FeatureInternal

    de_run <- run_edgeR_pairwise_de(
      count_matrix = de_matrix,
      sample_groups = pair_sample_groups,
      group1 = group1,
      group2 = group2,
      feature_ids = feature_lookup$FeatureInternal,
      min_replicates = 2
    )

    if (isTRUE(de_run$skipped)) {
      layer_manifest[[length(layer_manifest) + 1]] <- tibble(
        layer = layer_name,
        contrast = contrast_label,
        group1 = group1,
        group2 = group2,
        status = "skipped",
        message = de_run$message %||% "Skipped.",
        sample_n = length(pair_sample_order),
        feature_n = nrow(matrix_info$matrix),
        significant_n = 0,
        up_n = 0,
        down_n = 0
      )
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), de_run$message %||% "Skipped."),
        pair_section,
        paste0("skipped_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
      next
    }

    result_tbl <- de_run$result %>%
      mutate(FeatureInternal = as.character(Feature)) %>%
      select(-any_of(c("FeatureID", "Feature"))) %>%
      left_join(feature_lookup, by = "FeatureInternal") %>%
      select(Feature, everything(), -FeatureInternal) %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        group1 = group1,
        group2 = group2,
        contrast = contrast_label
      )

    if (!is.null(feature_meta) && nrow(feature_meta) > 0) {
      result_tbl <- result_tbl %>% left_join(feature_meta, by = c("Feature" = "Feature"))
    }

    summary_tbl <- result_tbl %>%
      count(direction, name = "feature_count") %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        group1 = group1,
        group2 = group2,
        contrast = contrast_label
      )

    sig_tbl <- result_tbl %>% filter(direction != "ns")
    up_tbl <- result_tbl %>% filter(direction == "up") %>% arrange(SignificanceMetric, PValue, desc(logFC))
    down_tbl <- result_tbl %>% filter(direction == "down") %>% arrange(SignificanceMetric, PValue, logFC)
    robustness_tbl <- build_de_robustness_table(
      result_tbl = result_tbl,
      count_matrix = matrix_info$matrix,
      sample_order = pair_sample_order,
      sample_groups = pair_sample_groups,
      group1 = group1,
      group2 = group2
    ) %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        contrast = contrast_label
      )
    robustness_summary_tbl <- summarise_de_robustness_flags(robustness_tbl)

    layer_results[[length(layer_results) + 1]] <- result_tbl
    layer_summary[[length(layer_summary) + 1]] <- summary_tbl
    layer_robustness[[length(layer_robustness) + 1]] <- robustness_tbl
    layer_robustness_summary[[length(layer_robustness_summary) + 1]] <- robustness_summary_tbl
    layer_manifest[[length(layer_manifest) + 1]] <- tibble(
      layer = layer_name,
      contrast = contrast_label,
      group1 = group1,
      group2 = group2,
      status = "completed",
      message = "",
      sample_n = length(pair_sample_order),
      feature_n = nrow(matrix_info$matrix),
      significant_n = nrow(sig_tbl),
      up_n = sum(result_tbl$direction == "up", na.rm = TRUE),
      down_n = sum(result_tbl$direction == "down", na.rm = TRUE)
    )

    write_csv_safe(result_tbl, file.path(pair_section, paste0("de_results_", safe_dir_name(contrast_label), ".csv")))
    write_csv_safe(summary_tbl, file.path(pair_section, paste0("de_summary_", safe_dir_name(contrast_label), ".csv")))
    if (nrow(up_tbl) > 0) {
      write_csv_safe(up_tbl, file.path(pair_section, paste0("de_up_", safe_dir_name(contrast_label), ".csv")))
    }
    if (nrow(down_tbl) > 0) {
      write_csv_safe(down_tbl, file.path(pair_section, paste0("de_down_", safe_dir_name(contrast_label), ".csv")))
    }
    if (nrow(robustness_tbl) > 0) {
      write_csv_safe(robustness_tbl, file.path(pair_section, paste0("de_robustness_", safe_dir_name(contrast_label), ".csv")))
      write_csv_safe(robustness_summary_tbl, file.path(pair_section, paste0("de_robustness_summary_", safe_dir_name(contrast_label), ".csv")))
      save_plot_nested(
        build_de_robustness_plot(
          robustness_tbl,
          paste0(layer_title, " robustness | ", contrast_label)
        ),
        pair_section,
        paste0("robustness_", safe_dir_name(contrast_label), ".pdf"),
        width = 8.8,
        height = 6.4,
        view = REPORT_VIEW
      )
    }

    volcano_plot <- build_de_volcano_plot(
      result_tbl,
      paste0(layer_title, " | ", contrast_label),
      label_n = 15
    )
    save_plot_nested(
      volcano_plot,
      pair_section,
      paste0("volcano_", safe_dir_name(contrast_label), ".pdf"),
      width = 8.6,
      height = 6.2,
      view = REPORT_VIEW
    )

    ma_plot <- build_de_ma_plot(
      result_tbl,
      paste0(layer_title, " MA | ", contrast_label)
    )
    save_plot_nested(
      ma_plot,
      pair_section,
      paste0("ma_", safe_dir_name(contrast_label), ".pdf"),
      width = 8.6,
      height = 6.2,
      view = REPORT_VIEW
    )

    top_heat_features <- result_tbl %>%
      arrange(SignificanceMetric, PValue, desc(abs(logFC))) %>%
      slice_head(n = 30) %>%
      pull(Feature) %>%
      unique()
    top_heat_features <- top_heat_features[top_heat_features %in% rownames(matrix_info$matrix)]

    if (length(top_heat_features) >= 2 && pkg_has("pheatmap")) {
      heat_dge <- edgeR::DGEList(counts = round(matrix_info$matrix[, pair_sample_order, drop = FALSE]))
      heat_dge <- edgeR::calcNormFactors(heat_dge)
      heat_mat <- edgeR::cpm(heat_dge, log = TRUE, prior.count = 1)[top_heat_features, , drop = FALSE]
      heat_mat <- safe_scale_rows(heat_mat)
      heat_mat <- pmax(pmin(heat_mat, 2.5), -2.5)
      heat_anno <- data.frame(Group = factor(pair_sample_groups, levels = c(group1, group2)))
      rownames(heat_anno) <- pair_sample_order

      heat_plot_path <- file.path(OUTPUT_FIG_DIR, pair_section, paste0("top_heatmap_", safe_dir_name(contrast_label), ".pdf"))
      save_pheatmap_file(
        heat_mat,
        annotation_col = heat_anno,
        annotation_colors = list(Group = group_palette[c(group1, group2)]),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = grDevices::colorRampPalette(mii_heatmap_gradient)(101),
        breaks = seq(-2.5, 2.5, length.out = 102),
        legend_breaks = seq(-2.5, 2.5, by = 1.25),
        fontsize_row = 7,
        fontsize_col = 9,
        main = paste0(layer_title, " top features | ", contrast_label),
        filename = heat_plot_path,
        width = 9,
        height = max(5.5, 0.25 * nrow(heat_mat) + 2.5)
      )
    } else {
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), "No heatmap available."),
        pair_section,
        paste0("heatmap_note_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
    }
  }

  combined_results <- bind_rows(layer_results)
  combined_summary <- bind_rows(layer_summary)
  combined_manifest <- bind_rows(layer_manifest)
  combined_robustness <- bind_rows(layer_robustness)
  combined_robustness_summary <- bind_rows(layer_robustness_summary)

  write_csv_safe(combined_results, paste0(layer_name, "_de_results.csv"))
  write_csv_safe(combined_summary, paste0(layer_name, "_de_summary.csv"))
  write_csv_safe(combined_manifest, paste0(layer_name, "_de_manifest.csv"))
  if (nrow(combined_robustness) > 0) {
    write_csv_safe(combined_robustness, paste0(layer_name, "_de_robustness.csv"))
    write_csv_safe(combined_robustness_summary, paste0(layer_name, "_de_robustness_summary.csv"))
  }

  list(
    results = combined_results,
    summary = combined_summary,
    manifest = combined_manifest,
    robustness = combined_robustness,
    robustness_summary = combined_robustness_summary
  )
}

normalize_mirna_name <- function(x) {
  x <- as.character(x)
  x <- tolower(trimws(x))
  x <- gsub("^[a-z]{2,3}[-_]+", "", x)
  x <- gsub("^mirna", "mir", x)
  x <- gsub("[-_. ]+(3p|5p)$", "", x)
  x <- gsub("[-_. ]+", "", x)
  x <- sub("^mirlet", "let", x)
  toupper(x)
}

extract_mirna_tokens_from_gs_name <- function(gs_name) {
  tokens <- unlist(strsplit(as.character(gs_name), "_", fixed = TRUE))
  tokens <- tokens[grepl("^(MIR|LET)", toupper(tokens))]
  unique(normalize_mirna_name(tokens))
}

empty_mirna_target_sets <- function(query_tokens = character(), source = "not_available", status = "skipped", message = "") {
  list(
    mapping = tibble(
      mirna_token = character(),
      target_gene = character(),
      source_database = character(),
      evidence_source = character(),
      gs_name = character()
    ),
    summary = tibble(
      mirna_token = character(),
      target_gene_count = integer(),
      target_set_count = integer(),
      total_links = integer(),
      target_genes = character(),
      target_sets = character()
    ),
    target_genes = tibble(target_gene = character()),
    query_tokens = query_tokens,
    source_summary = tibble(
      target_source = source,
      status = status,
      message = message,
      query_miRNA_n = length(query_tokens),
      target_gene_n = 0L,
      target_link_n = 0L
    )
  )
}

summarise_mirna_target_mapping <- function(mapping, query_tokens, source, status = "completed", message = "") {
  if (nrow(mapping) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = source, status = "empty", message = message))
  }

  mapping <- mapping %>%
    mutate(
      mirna_token = normalize_mirna_name(mirna_token),
      target_gene = trimws(as.character(target_gene)),
      gs_name = if ("gs_name" %in% names(.)) as.character(gs_name) else paste(source_database, mirna_token, sep = ":"),
      source_database = if ("source_database" %in% names(.)) as.character(source_database) else source,
      evidence_source = if ("evidence_source" %in% names(.)) as.character(evidence_source) else source
    ) %>%
    filter(mirna_token %in% query_tokens, !is.na(target_gene), nzchar(target_gene)) %>%
    group_by(mirna_token, target_gene, source_database, evidence_source) %>%
    summarise(
      target_set_n = n_distinct(gs_name),
      gs_name = paste(head(sort(unique(gs_name)), 20), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(mirna_token, target_gene, source_database)

  if (nrow(mapping) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = source, status = "empty", message = message))
  }

  summary_tbl <- mapping %>%
    group_by(mirna_token) %>%
    summarise(
      target_gene_count = n_distinct(target_gene),
      target_set_count = sum(target_set_n, na.rm = TRUE),
      total_links = n(),
      target_genes = paste(sort(unique(target_gene)), collapse = ";"),
      target_sets = paste(sort(unique(gs_name)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(target_gene_count), mirna_token)

  target_genes <- mapping %>%
    distinct(target_gene) %>%
    arrange(target_gene)

  list(
    mapping = mapping,
    summary = summary_tbl,
    target_genes = target_genes,
    query_tokens = query_tokens,
    source_summary = tibble(
      target_source = source,
      status = status,
      message = message,
      query_miRNA_n = length(query_tokens),
      target_gene_n = nrow(target_genes),
      target_link_n = nrow(mapping)
    )
  )
}

format_multimir_query <- function(mirna_names) {
  x <- trimws(as.character(mirna_names))
  x <- x[!is.na(x) & nzchar(x)]
  prefix <- paste0(MIRNA_PREFIX, "-")
  x <- str_replace(x, regex("^[a-z]{2,3}-mir", ignore_case = TRUE), paste0(prefix, "miR"))
  x <- str_replace(x, regex("^mir", ignore_case = TRUE), paste0(prefix, "miR"))
  x <- str_replace(x, regex("^[a-z]{2,3}-let", ignore_case = TRUE), paste0(prefix, "let"))
  x <- str_replace(x, regex("^let", ignore_case = TRUE), paste0(prefix, "let"))
  unique(x)
}

collect_multimir_target_sets <- function(mirna_names, query_tokens, target_source) {
  if (!pkg_has("multiMiR")) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "skipped", message = "multiMiR package is not available."))
  }

  target_tables <- if (identical(target_source, "predicted")) c("validated", "predicted") else "validated"
  query_mirnas <- format_multimir_query(mirna_names)
  if (length(query_mirnas) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "empty", message = "No miRNA names were available for multiMiR."))
  }

  mapping <- purrr::map_dfr(target_tables, function(target_table) {
    mm_obj <- tryCatch(
      multiMiR::get_multimir(org = MULTIMIR_ORG, mirna = query_mirnas, table = target_table, legacy.out = FALSE),
      error = function(e) e
    )
    if (inherits(mm_obj, "error")) {
      return(tibble(.message = conditionMessage(mm_obj), .table = target_table))
    }

    mm_df <- tryCatch(
      as_tibble(mm_obj@data),
      error = function(e) tryCatch(as_tibble(as.data.frame(mm_obj)), error = function(e2) tibble())
    )
    if (nrow(mm_df) == 0) return(tibble())

    mirna_col <- intersect(c("mature_mirna_id", "mirna", "mirna_id", "mature_mirna_acc"), names(mm_df))[1]
    target_col <- intersect(c("target_symbol", "target_gene", "symbol", "target"), names(mm_df))[1]
    db_col <- intersect(c("database", "source", "support_type"), names(mm_df))[1]
    if (is.na(mirna_col) || is.na(target_col)) return(tibble())

    tibble(
      mirna_token = normalize_mirna_name(mm_df[[mirna_col]]),
      target_gene = as.character(mm_df[[target_col]]),
      source_database = if (!is.na(db_col)) as.character(mm_df[[db_col]]) else "multiMiR",
      evidence_source = target_table,
      gs_name = paste0("multiMiR:", target_table, ":", if (!is.na(db_col)) as.character(mm_df[[db_col]]) else "database")
    )
  })

  if (".message" %in% names(mapping) && all(!all(c("mirna_token", "target_gene") %in% names(mapping)))) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "skipped", message = paste(unique(mapping$.message), collapse = "; ")))
  }
  mapping <- mapping %>% select(any_of(c("mirna_token", "target_gene", "source_database", "evidence_source", "gs_name")))

  summarise_mirna_target_mapping(mapping, query_tokens, source = paste0("multiMiR_", target_source), status = "completed")
}

collect_msigdb_mir_target_sets <- function(mirna_names, query_tokens, species = MIRNA_SPECIES, subcollections = c("MIR:MIRDB", "MIR:MIR_LEGACY"), fallback_message = "") {
  if (!pkg_has("msigdbr")) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "skipped", message = "msigdbr package is not available."))
  }

  target_db <- tryCatch(
    purrr::map_dfr(subcollections, function(subcollection) {
      msigdbr::msigdbr(species = species, collection = "C3", subcollection = subcollection) %>%
        transmute(
          gs_collection,
          gs_subcollection,
          gs_name,
          gs_description,
          target_gene = gene_symbol,
          mirna_token = list(extract_mirna_tokens_from_gs_name(gs_name)),
          source_database = "MSigDB_C3_MIR",
          evidence_source = "msigdb_fallback"
        ) %>%
        filter(purrr::map_lgl(mirna_token, ~any(.x %in% query_tokens)))
    }),
    error = function(e) tibble(.message = conditionMessage(e))
  )

  if (".message" %in% names(target_db)) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "skipped", message = unique(target_db$.message)[[1]]))
  }
  if (nrow(target_db) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "empty", message = "No MSigDB C3 MIR records were returned."))
  }

  mapping <- target_db %>%
    tidyr::unnest_longer(mirna_token, values_to = "mirna_token") %>%
    filter(!is.na(mirna_token), mirna_token != "", mirna_token %in% query_tokens) %>%
    transmute(
      mirna_token,
      target_gene,
      source_database,
      evidence_source,
      gs_name
    )

  status <- if (nzchar(fallback_message)) "fallback" else "completed"
  summarise_mirna_target_mapping(mapping, query_tokens, source = "msigdb_fallback", status = status, message = fallback_message)
}

collect_mirna_target_sets <- function(mirna_names, species = MIRNA_SPECIES, subcollections = c("MIR:MIRDB", "MIR:MIR_LEGACY"), target_source = MIRNA_TARGET_SOURCE) {
  query_tokens <- unique(normalize_mirna_name(mirna_names))
  query_tokens <- query_tokens[nzchar(query_tokens)]
  if (Sys.getenv("RUN_TARGET_ENRICHMENT","0") != "1") return(empty_mirna_target_sets(query_tokens, source=target_source, status="disabled", message="RUN_TARGET_ENRICHMENT=0; no network query"))
  if (length(query_tokens) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "empty", message = "No query miRNA tokens were available."))
  }

  if (!identical(target_source, "msigdb_fallback")) {
    multimir_sets <- collect_multimir_target_sets(mirna_names, query_tokens, target_source = target_source)
    if (nrow(multimir_sets$mapping) > 0) {
      return(multimir_sets)
    }
    multimir_sets$source_summary <- multimir_sets$source_summary %>%
      mutate(
        status = "skipped",
        message = paste0(
          "Primary source unavailable or empty: ",
          message %||% "multiMiR returned no targets.",
          " Exploratory MSigDB C3 MIR fallback is disabled unless MIRNA_TARGET_SOURCE=msigdb_fallback."
        )
      )
    return(multimir_sets)
  }

  collect_msigdb_mir_target_sets(mirna_names, query_tokens, species = species, subcollections = subcollections)
}

write_setup_provenance <- function(ctx, output_dir = SETUP_DIR) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  tool_manifest <- bind_rows(
    tibble(tool = "R", available = TRUE, path = R.home(), version = paste(R.version.string, collapse = " ")),
    capture_command_version("Rscript", "--version"),
    capture_command_version("perl", "-v"),
    capture_command_version("python3", "--version"),
    capture_command_version("cutadapt", "--version"),
    capture_command_version("fastqc", "--version"),
    capture_command_version("multiqc", "--version"),
    capture_command_version("bowtie", "--version"),
    capture_command_version("bowtie2", "--version"),
    capture_command_version("samtools", "--version")
  )

  input_audit <- ctx$sample_info %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      RawFastqOriginal = RawFastqOriginal,
      RawFastq = RawFastq,
      RawFastqExists = RawFastqExists,
      SummaryFile = SummaryFile,
      SummaryFileExists = file.exists(SummaryFile),
      LengthFile = LengthFile,
      LengthFileExists = file.exists(LengthFile),
      OutputFile = OutputFile,
      OutputFileExists = file.exists(OutputFile)
    )

  sports_script <- file.path(ROOT_DIR, "sports.sh")
  sports_config <- if (file.exists(sports_script)) {
    script_lines <- readLines(sports_script, warn = FALSE)
    config_keys <- c(
      "RESULT_DIR", "METADATA_DIR", "QC_DIR", "SPORTS_ANNOTATION_DIR", "SAMPLE_ROOT",
      "ENDOSIRNA_DIR", "ENDOSIRNA_SAMPLE_ROOT", "ENDOSIRNA_SUMMARY_DIR",
      "LOG_DIR", "SAMPLE_LOG_DIR", "CUTADAPT_LOG_DIR", "ENDOSIRNA_LOG_DIR", "WORK_DIR", "TMP_DIR",
      "DB_ROOT", "GENOME_BASE", "MIRNA_BASE", "RRNA_BASE", "TRNA_BASE",
      "TRNA_MITO_BASE", "PIRNA_BASE", "ENSEMBL_BASE", "RFAM_BASE",
      "ADAPTER3", "PRETRIM_FASTQ", "TRIM_ADAPTER", "CUTADAPT_NEXTSEQ_TRIM", "CUTADAPT_QUALITY_CUTOFF", "MISMATCH", "RUN_FASTQC",
      "RUN_ENDOSIRNA", "ENDOSIRNA_REFS", "ENDOSIRNA_MISMATCH", "ENDOSIRNA_MAX_ALIGNMENTS",
      "THREADS", "PARALLEL_JOBS", "MIN_LEN", "MAX_LEN"
    )
    purrr::map_dfr(config_keys, function(key) {
      hit <- script_lines[str_detect(script_lines, paste0("^", key, "="))][1]
      tibble(key = key, script_value = hit %||% "", environment_value = Sys.getenv(key, unset = ""))
    })
  } else {
    tibble(key = character(), script_value = character(), environment_value = character())
  }

  write_tsv_safe(build_runtime_config_table(), file.path(output_dir, "analysis_runtime_config.tsv"))
  write_tsv_safe(tool_manifest, file.path(output_dir, "tool_version_manifest.tsv"))
  write_tsv_safe(input_audit, file.path(output_dir, "input_path_audit.tsv"))
  write_tsv_safe(sports_config, file.path(output_dir, "sports_pipeline_config.tsv"))
  utils::write.csv(
    ctx$sample_comp_view %>%
      transmute(
        Sample = as.character(Sample), Group = as.character(Group),
        MainCategory = as.character(MainCategory), Reads, CleanReads, Percent
      ),
    file.path(output_dir, "exclusive_class_composition.csv"), row.names = FALSE
  )
  utils::write.csv(ctx$composition_integrity_check, file.path(output_dir, "composition_integrity_check.csv"), row.names = FALSE)
  utils::write.csv(ctx$length_composition_integrity_check, file.path(output_dir, "length_composition_integrity_check.csv"), row.names = FALSE)
  utils::write.csv(ctx$tsrna_cleavage_reconstruction_check, file.path(output_dir, "tsrna_cleavage_reconstruction_check.csv"), row.names = FALSE)
  write_tsv_safe(ctx$raw_fastq_audit_by_sample %||% tibble(), file.path(output_dir, "raw_fastq_audit_by_sample.tsv"))
  write_tsv_safe(ctx$raw_fastq_audit_by_group %||% tibble(), file.path(output_dir, "raw_fastq_audit_by_group.tsv"))
  invisible(output_dir)
}

# ===== 模块 1 表达矩阵 =====
collapse_unique_text <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else paste(sort(x), collapse = ";")
}

make_company_tsrna_feature <- function(sub_class, cleavage) {
  base_name <- as.character(sub_class) %>%
    str_remove("^mature-") %>%
    str_remove("^pre-") %>%
    str_replace("_(5_end|3_end|CCA_end)$", "")
  suffix <- case_when(
    cleavage == "5'end" ~ "_5_end",
    cleavage == "3'end" ~ "_3_end",
    cleavage == "CCA-end" ~ "_3_CCA",
    TRUE ~ "_inner"
  )
  paste0(base_name, suffix)
}

empty_expression_matrix <- function(output_col, sample_order) {
  out <- tibble(.feature = character())
  names(out)[1] <- output_col
  for (sample_name in sample_order) out[[sample_name]] <- numeric()
  out
}

build_expression_wide_matrix <- function(df, feature_col, value_col, output_col, sample_order) {
  if (nrow(df) == 0) return(empty_expression_matrix(output_col, sample_order))

  long_df <- df %>%
    transmute(
      .feature = as.character(.data[[feature_col]]),
      Sample = as.character(Sample),
      .value = clean_numeric(.data[[value_col]])
    ) %>%
    filter(!is.na(.feature), .feature != "", Sample %in% sample_order) %>%
    group_by(.feature, Sample) %>%
    summarise(.value = sum(.value, na.rm = TRUE), .groups = "drop")

  if (nrow(long_df) == 0) return(empty_expression_matrix(output_col, sample_order))

  wide_df <- long_df %>%
    mutate(Sample = factor(Sample, levels = sample_order)) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = .value, values_fill = 0)

  for (sample_name in sample_order) {
    if (!sample_name %in% names(wide_df)) wide_df[[sample_name]] <- 0
  }

  wide_df <- wide_df %>%
    mutate(.total = rowSums(across(all_of(sample_order)), na.rm = TRUE)) %>%
    arrange(desc(.total), .feature) %>%
    select(.feature, all_of(sample_order))
  names(wide_df)[1] <- output_col
  wide_df
}

prepare_annotation_expression_values <- function(annotation_df, clean_reads_tbl) {
  annotation_df %>%
    group_by(Feature, Sample) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0))
}

write_annotation_expression_set <- function(annotation_df, set_name, output_dir, sample_order, clean_reads_tbl) {
  values_df <- prepare_annotation_expression_values(annotation_df, clean_reads_tbl)
  count_mat <- build_expression_wide_matrix(values_df, "Feature", "Reads", "Sub_Class", sample_order)
  rpm_mat <- build_expression_wide_matrix(values_df, "Feature", "RPM", "Sub_Class", sample_order)
  write_tsv_safe(count_mat, file.path(output_dir, paste0("count_", set_name, ".txt")))
  write_tsv_safe(rpm_mat, file.path(output_dir, paste0("rpm_", set_name, ".txt")))
}

write_sequence_expression_set <- function(sequence_df, set_name, output_dir, seq_dir, sample_order, duplicate_root = FALSE) {
  if (pkg_has("data.table")) {
    sequence_dt <- data.table::as.data.table(sequence_df)
    if (nrow(sequence_dt) == 0) {
      count_mat <- empty_expression_matrix("Sequence", sample_order)
      rpm_mat <- empty_expression_matrix("Sequence", sample_order)
      type_tbl <- tibble(Sequence = character(), Annotation = character(), FeatureCategory = character())
    } else {
      sequence_dt <- sequence_dt[!is.na(Sequence) & Sequence != ""]
      sequence_dt[, Sample := as.character(Sample)]
      sequence_dt <- sequence_dt[Sample %in% sample_order]
      if (nrow(sequence_dt) == 0) {
        count_mat <- empty_expression_matrix("Sequence", sample_order)
        rpm_mat <- empty_expression_matrix("Sequence", sample_order)
        type_tbl <- tibble(Sequence = character(), Annotation = character(), FeatureCategory = character())
      } else {
        sequence_dt[, Reads := clean_numeric(Reads)]
        sequence_dt[, RPM := clean_numeric(RPM)]

        values_dt <- sequence_dt[, .(
          Reads = sum(Reads, na.rm = TRUE),
          RPM = sum(RPM, na.rm = TRUE)
        ), by = .(Sequence, Sample)]

        count_dt <- data.table::dcast(values_dt, Sequence ~ Sample, value.var = "Reads", fill = 0, fun.aggregate = sum)
        rpm_dt <- data.table::dcast(values_dt, Sequence ~ Sample, value.var = "RPM", fill = 0, fun.aggregate = sum)
        for (sample_name in sample_order) {
          if (!sample_name %in% names(count_dt)) count_dt[, (sample_name) := 0]
          if (!sample_name %in% names(rpm_dt)) rpm_dt[, (sample_name) := 0]
        }

        if (nrow(count_dt) > 0) {
          count_dt[, .total := rowSums(.SD, na.rm = TRUE), .SDcols = sample_order]
          data.table::setorderv(count_dt, c(".total", "Sequence"), c(-1L, 1L))
          sequence_order <- count_dt$Sequence
          rpm_dt <- rpm_dt[match(sequence_order, Sequence)]
          count_mat <- as_tibble(count_dt[, c("Sequence", sample_order), with = FALSE])
          rpm_mat <- as_tibble(rpm_dt[, c("Sequence", sample_order), with = FALSE])
        } else {
          sequence_order <- character()
          count_mat <- empty_expression_matrix("Sequence", sample_order)
          rpm_mat <- empty_expression_matrix("Sequence", sample_order)
        }

        collapse_dt_text <- function(x) {
          x <- unique(as.character(x))
          x <- x[!is.na(x) & nzchar(x)]
          if (length(x) == 0) "" else paste(sort(x), collapse = ";")
        }
        type_dt <- sequence_dt[, .(
          Annotation = collapse_dt_text(PrimaryAnnotation),
          FeatureCategory = collapse_dt_text(FeatureCategory)
        ), by = Sequence]
        type_dt[, .order := match(Sequence, sequence_order)]
        data.table::setorderv(type_dt, c(".order", "Sequence"), c(1L, 1L), na.last = TRUE)
        type_dt[, .order := NULL]
        type_tbl <- as_tibble(type_dt)
      }
    }
  } else {
    values_df <- sequence_df %>%
      group_by(Sequence, Sample) %>%
      summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop")

    count_mat <- build_expression_wide_matrix(values_df, "Sequence", "Reads", "Sequence", sample_order)
    rpm_mat <- build_expression_wide_matrix(values_df, "Sequence", "RPM", "Sequence", sample_order)
    sequence_order <- count_mat$Sequence

    type_tbl <- sequence_df %>%
      group_by(Sequence) %>%
      summarise(
        Annotation = collapse_unique_text(PrimaryAnnotation),
        FeatureCategory = collapse_unique_text(FeatureCategory),
        .groups = "drop"
      ) %>%
      mutate(.order = match(Sequence, sequence_order)) %>%
      arrange(.order, Sequence) %>%
      select(-.order)
  }

  write_tsv_safe(count_mat, file.path(seq_dir, paste0("seq_count_", set_name, ".txt")))
  write_tsv_safe(rpm_mat, file.path(seq_dir, paste0("seq_rpm_", set_name, ".txt")))
  write_tsv_safe(type_tbl, file.path(seq_dir, paste0("seq_type_", set_name, ".txt")))

  if (isTRUE(duplicate_root)) {
    write_tsv_safe(count_mat, file.path(output_dir, paste0("seq_count_", set_name, ".txt")))
    write_tsv_safe(rpm_mat, file.path(output_dir, paste0("seq_rpm_", set_name, ".txt")))
    write_tsv_safe(type_tbl, file.path(output_dir, paste0("seq_type_", set_name, ".txt")))
  }
}

write_mii_setup_expression_matrices <- function(report_view = REPORT_VIEW) {
  ctx <- load_mii_context(report_view = report_view)
  sample_order <- as.character(ctx$sample_info$Sample)
  output_dir <- EXPRESSION_MATRIX_DIR
  seq_dir <- file.path(output_dir, "seq_count")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)
  write_setup_provenance(ctx, output_dir = SETUP_DIR)

  non_tsrna_annotation <- ctx$summary_all %>%
    filter(SubClass != "-", MainCategory != "tsRNA", !Class %in% c("Clean_Reads", "Match_Genome")) %>%
    transmute(
      Sample = as.character(Sample),
      Feature = if_else(is.na(SubClass) | SubClass == "", Class, SubClass),
      FeatureCategory = as.character(MainCategory),
      Reads = Reads
    )

  tsrna_annotation <- ctx$tsrna_cleavage_detail %>%
    transmute(
      Sample = as.character(Sample),
      Feature = make_company_tsrna_feature(SubClass, Cleavage),
      FeatureCategory = if_else(Origin == "Mitochondrial", "mt_tsRNA", "tsRNA"),
      Reads = Reads
    )

  annotation_input <- bind_rows(non_tsrna_annotation, tsrna_annotation)
  annotation_sets <- list(
    all = annotation_input,
    miRNA = annotation_input %>% filter(FeatureCategory == "miRNA"),
    tsRNA = annotation_input %>% filter(FeatureCategory == "tsRNA"),
    mt_tsRNA = annotation_input %>% filter(FeatureCategory == "mt_tsRNA")
  )
  walk(names(annotation_sets), ~write_annotation_expression_set(annotation_sets[[.x]], .x, output_dir, sample_order, ctx$clean_reads_tbl))

  if (pkg_has("data.table")) {
    sequence_input <- data.table::as.data.table(ctx$output_all)
    sequence_input <- sequence_input[!is.na(Sequence) & Sequence != ""]
    sequence_input[, FeatureCategory := as.character(FeatureCategory)]
    sequence_input[, PrimaryAnnotation := as.character(PrimaryAnnotation)]
    sequence_input[is.na(PrimaryAnnotation) | PrimaryAnnotation == "", PrimaryAnnotation := "NO_Annotation"]
    sequence_input[, SequenceSet := FeatureCategory]
    sequence_input[
      FeatureCategory == "tsRNA" & grepl("mt_tRNA|mitochond", PrimaryAnnotation, ignore.case = TRUE),
      SequenceSet := "mt_tsRNA"
    ]
    sequence_input <- sequence_input[, .(Sample, Sequence, Reads, RPM, PrimaryAnnotation, FeatureCategory, SequenceSet)]
    sequence_sets <- list(
      all = sequence_input,
      miRNA = sequence_input[SequenceSet == "miRNA"],
      tsRNA = sequence_input[SequenceSet == "tsRNA"],
      mt_tsRNA = sequence_input[SequenceSet == "mt_tsRNA"],
      piRNA = sequence_input[SequenceSet == "piRNA"],
      rsRNA = sequence_input[SequenceSet == "rsRNA"]
    )
  } else {
    sequence_input <- ctx$output_all %>%
      filter(!is.na(Sequence), Sequence != "") %>%
      mutate(
        FeatureCategory = as.character(FeatureCategory),
        PrimaryAnnotation = if_else(is.na(PrimaryAnnotation) | PrimaryAnnotation == "", "NO_Annotation", PrimaryAnnotation),
        IsMtTsRNA = FeatureCategory == "tsRNA" & str_detect(PrimaryAnnotation, regex("mt_tRNA|mitochond", ignore_case = TRUE)),
        SequenceSet = case_when(
          IsMtTsRNA ~ "mt_tsRNA",
          TRUE ~ FeatureCategory
        )
      ) %>%
      select(Sample, Sequence, Reads, RPM, PrimaryAnnotation, FeatureCategory, SequenceSet)

    sequence_sets <- list(
      all = sequence_input,
      miRNA = sequence_input %>% filter(SequenceSet == "miRNA"),
      tsRNA = sequence_input %>% filter(SequenceSet == "tsRNA"),
      mt_tsRNA = sequence_input %>% filter(SequenceSet == "mt_tsRNA"),
      piRNA = sequence_input %>% filter(SequenceSet == "piRNA"),
      rsRNA = sequence_input %>% filter(SequenceSet == "rsRNA")
    )
  }
  walk(names(sequence_sets), ~write_sequence_expression_set(sequence_sets[[.x]], .x, output_dir, seq_dir, sample_order, duplicate_root = .x == "all"))

  invisible(output_dir)
}

should_write_setup_expression_matrices <- function() {
  exists("setup_file", inherits = TRUE) && !exists("input_file", inherits = TRUE)
}

if (should_write_setup_expression_matrices()) {
  write_mii_setup_expression_matrices(report_view = REPORT_VIEW)
}
