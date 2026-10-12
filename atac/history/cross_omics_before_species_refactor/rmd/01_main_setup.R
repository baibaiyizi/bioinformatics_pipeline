# Liver ATAC-seq shared setup
# Auto-generated from liveratac.Rmd setup and adjusted for modular execution.
first_nonempty_value <- function(...) {
  for(value in list(...)) {
    if(is.null(value) || length(value) == 0) next
    value_chr <- trimws(as.character(value[[1]]))
    if(!is.na(value_chr) && nzchar(value_chr)) return(value_chr)
  }
  ""
}

env_bool <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if(!nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "yes", "y", "on")
}

env_bool_checked <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if(!nzchar(value)) return(default)
  normalized <- tolower(trimws(value))
  if(!normalized %in% c("0", "1", "true", "false", "yes", "no", "y", "n", "on", "off")) {
    stop(name, " 必须是 0/1 或合法布尔值，当前为: ", value)
  }
  normalized %in% c("1", "true", "yes", "y", "on")
}

if(!exists("params", inherits = FALSE) || !is.list(params)) {
  params <- list(
    work_dir = first_nonempty_value(Sys.getenv("ATAC_ROOT_DIR", unset = ""), getwd()),
    show_code = env_bool("ATAC_SHOW_CODE", FALSE),
    cache_results = env_bool("ATAC_CACHE_RESULTS", FALSE),
    emit_light_plots = env_bool("ATAC_EMIT_LIGHT_PLOTS", TRUE),
    save_png_previews = env_bool("ATAC_SAVE_PNG_PREVIEWS", TRUE),
    png_dpi = first_nonempty_value(
      Sys.getenv("ATAC_PNG_DPI", unset = ""),
      Sys.getenv("FIG_PNG_DPI", unset = ""),
      "300"
    ),
    run_profile = first_nonempty_value(Sys.getenv("ATAC_RUN_PROFILE", unset = ""), "full"),
    reuse_step_tables = env_bool("ATAC_REUSE_STEP_TABLES", FALSE),
    strict_mode = env_bool("ATAC_STRICT_MODE", FALSE)
  )
}

light_setup <- env_bool("ATAC_LIGHT_SETUP", FALSE)

knitr::opts_chunk$set(
  echo = isTRUE(params$show_code),
  warning = FALSE,
  message = FALSE,
  error = FALSE,
  cache = isTRUE(params$cache_results),
  fig.align = "center"
)

# --- 基础 helper: 依赖、路径与状态打印 ---
`%||%` <- function(x, y) {
  if(is.null(x) || length(x) == 0) y else x
}

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs) > 0) {
    stop(
      "缺少以下 R 包，请先安装后再渲染: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

report_status <- function(level, text) {
  cat(sprintf("[%s] %s\n", level, text))
}

resolve_dir <- function(path_value, fallback) {
  chosen <- if(!is.null(path_value) && nzchar(path_value)) path_value else fallback
  if(!nzchar(chosen)) return("")
  normalizePath(chosen, mustWork = FALSE)
}

# --- 加载依赖包 ---
core_pkgs <- if(isTRUE(light_setup)) {
  c("knitr", "ggplot2", "dplyr", "tidyr", "jsonlite", "scales")
} else {
  c(
    "ggplot2", "dplyr", "tidyr", "jsonlite", "pheatmap", "ggrepel", "UpSetR", "patchwork", "scales",
    "GenomicRanges", "GenomeInfoDb", "Biostrings", "GenomicAlignments", "Rsamtools",
    "ChIPseeker", "DiffBind", "enrichplot", "clusterProfiler",
    "TxDb.Mmusculus.UCSC.mm39.refGene", "org.Mm.eg.db", "AnnotationDbi",
    "BSgenome.Mmusculus.UCSC.mm39", "GenomicScores",
    "motifmatchr", "JASPAR2020", "TFBSTools", "ggseqlogo",
    "igraph", "ggraph", "Gviz"
  )
}
load_required_packages(core_pkgs)

# Shared publication export contract. Width/height supplied by modules are
# retained only as aspect-ratio hints; all final figures use a 183 mm page.
atac_figure_width_mm <- 183
atac_figure_min_height_mm <- 70
atac_figure_max_height_mm <- 240
atac_figure_dpi <- 300
atac_figure_font_family <- "Liberation Sans"
atac_figure_base_size <- 8
atac_figure_tick_size <- 7

atac_figure_dimensions <- function(width = 8, height = 6) {
  width <- suppressWarnings(as.numeric(width)[1])
  height <- suppressWarnings(as.numeric(height)[1])
  ratio <- height / width
  if(!is.finite(ratio) || ratio <= 0) ratio <- 0.75
  height_mm <- min(atac_figure_max_height_mm, max(atac_figure_min_height_mm, atac_figure_width_mm * ratio))
  pdf_width_mm <- round(atac_figure_width_mm / 25.4 * 72) * 25.4 / 72
  pdf_height_mm <- round(height_mm / 25.4 * 72) * 25.4 / 72
  list(
    width_mm = atac_figure_width_mm,
    height_mm = height_mm,
    pdf_width_mm = pdf_width_mm,
    pdf_height_mm = pdf_height_mm,
    width_in = pdf_width_mm / 25.4,
    height_in = pdf_height_mm / 25.4
  )
}

atac_finalize_rgb_png <- function(path, dpi = atac_figure_dpi) {
  if(!file.exists(path) || !requireNamespace("png", quietly = TRUE)) return(invisible(path))
  raster <- png::readPNG(path, native = FALSE)
  if(length(dim(raster)) == 2L) {
    raster <- array(rep(raster, 3L), dim = c(dim(raster), 3L))
  } else if(dim(raster)[3L] == 4L) {
    alpha <- raster[, , 4L]
    rgb <- sweep(raster[, , 1:3, drop = FALSE], c(1L, 2L), alpha, `*`)
    raster <- rgb + array(rep(1 - alpha, 3L), dim = c(dim(alpha), 3L))
  }
  png::writePNG(raster[, , 1:3, drop = FALSE], target = path, dpi = dpi)
  invisible(path)
}

atac_publication_typography <- function(coordinate_system = c("axes", "none")) {
  coordinate_system <- match.arg(coordinate_system)
  out <- ggplot2::theme(
    text = ggplot2::element_text(family = atac_figure_font_family, size = atac_figure_base_size, color = "#333333"),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = atac_figure_base_size, color = "#333333"),
    axis.text = ggplot2::element_text(size = atac_figure_tick_size, color = "#333333"),
    strip.text = ggplot2::element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
    legend.title = ggplot2::element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
    legend.text = ggplot2::element_text(size = atac_figure_tick_size, color = "#333333"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
  )
  if(identical(coordinate_system, "none")) {
    out <- out + ggplot2::theme(
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank()
    )
  }
  out
}

atac_standardize_plot <- function(plot, coordinate_system = c("axes", "none")) {
  coordinate_system <- match.arg(coordinate_system)
  typography <- atac_publication_typography(coordinate_system)
  if(inherits(plot, "patchwork")) return(plot & typography)
  plot + typography
}

atac_open_pdf <- function(file, width, height, bg = "white") {
  dims <- atac_figure_dimensions(width, height)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if(file.exists(file)) unlink(file, force = TRUE)
  grDevices::cairo_pdf(
    file,
    width = dims$width_in,
    height = dims$height_in,
    pointsize = atac_figure_base_size,
    family = atac_figure_font_family,
    bg = "white"
  )
  invisible(grDevices::dev.cur())
}

atac_open_png <- function(file, width, height, bg = "white") {
  dims <- atac_figure_dimensions(width, height)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if(file.exists(file)) unlink(file, force = TRUE)
  grDevices::png(
    filename = file,
    width = dims$width_mm,
    height = dims$height_mm,
    units = "mm",
    res = atac_figure_dpi,
    pointsize = atac_figure_base_size,
    type = "cairo-png",
    bg = "white"
  )
  invisible(grDevices::dev.cur())
}

atac_close_device <- function(device) {
  device <- suppressWarnings(as.integer(device)[1])
  open_devices <- grDevices::dev.list()
  if(!is.finite(device) || is.null(open_devices) || !(device %in% unname(open_devices))) {
    return(invisible(FALSE))
  }
  grDevices::dev.off(which = device)
  invisible(TRUE)
}

# 绑定常用 dplyr 动词，避免被 AnnotationDbi / BiocGenerics 等同名函数覆盖。
select <- dplyr::select
count <- dplyr::count
filter <- dplyr::filter
mutate <- dplyr::mutate
transmute <- dplyr::transmute
arrange <- dplyr::arrange
distinct <- dplyr::distinct
left_join <- dplyr::left_join
inner_join <- dplyr::inner_join
full_join <- dplyr::full_join
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
summarise <- dplyr::summarise
summarize <- dplyr::summarize
slice_head <- dplyr::slice_head
slice_max <- dplyr::slice_max
pull <- dplyr::pull
rename <- dplyr::rename

# --- 路径配置 ---
work_dir <- resolve_dir(params$work_dir, getwd())
rna_seq_dir <- if(nzchar(Sys.getenv("ATAC_RNA_DE_FILE"))) dirname(Sys.getenv("ATAC_RNA_DE_FILE")) else ""
wgbs_dir <- if(nzchar(Sys.getenv("ATAC_WGBS_DMR_FILE"))) dirname(Sys.getenv("ATAC_WGBS_DMR_FILE")) else ""
grcm39_reference_dir <- Sys.getenv("REF_DIR")
if(!nzchar(grcm39_reference_dir)) stop("请通过 run_rmd.sh 加载 config/project.env")
reference_dir <- dirname(grcm39_reference_dir)
atac_mm10_to_mm39_chain_file <- first_nonempty_value(
  Sys.getenv("ATAC_MM10_TO_MM39_CHAIN", unset = ""),
  file.path(grcm39_reference_dir, "mm10ToMm39.over.chain")
)
atac_ensembl_gtf_file <- first_nonempty_value(
  Sys.getenv("ATAC_ENSEMBL_GTF", unset = ""),
  file.path(grcm39_reference_dir, "Mus_musculus.GRCm39.115.gtf")
)
atac_enable_liftover_related_gene_annotation <- env_bool(
  "ATAC_ENABLE_LIFTOVER_RELATED_GENE_ANNOTATION",
  FALSE
)
motif_reference_dir <- file.path(grcm39_reference_dir, "motif")
if(!dir.exists(motif_reference_dir)) {
  motif_reference_dir <- file.path(reference_dir, "motif")
}
default_jaspar_meme_db <- file.path(motif_reference_dir, "JASPAR2020_CORE_vertebrates.meme")
result_dir <- file.path(work_dir, "result")
rmd_dir <- resolve_dir(Sys.getenv("ATAC_RMD_DIR", unset = ""), file.path(work_dir, "rmd"))
source(file.path(rmd_dir, "significance_helpers.R"), local = TRUE)
cache_root_dir <- file.path(work_dir, ".cache", "rmd")
log_dir <- file.path(result_dir, "logs")
igv_dir <- file.path(result_dir, "IGV")
# Plot and table files are written directly into each module folder, e.g.
# result/07_peak_gene_enrichment/plots and result/07_peak_gene_enrichment/tables,
# not separate result/plots or result/tables trees.
plot_root_dir <- result_dir
table_root_dir <- result_dir
emit_light_plots <- isTRUE(params$emit_light_plots)
save_png_previews <- isTRUE(params$save_png_previews)
plot_png_dpi <- atac_figure_dpi
plot_pdf_preview_max_pages <- suppressWarnings(as.integer(first_nonempty_value(
  Sys.getenv("ATAC_PDF_PREVIEW_MAX_PAGES", unset = ""),
  "4"
)))
if(is.na(plot_pdf_preview_max_pages) || plot_pdf_preview_max_pages < 1) {
  plot_pdf_preview_max_pages <- 4L
}
plot_pdf_preview_max_pages <- min(plot_pdf_preview_max_pages, 10L)
run_profile <- tolower(as.character(params$run_profile %||% "full"))[1]
if(is.na(run_profile) || !run_profile %in% c("full", "quick", "report")) {
  run_profile <- "full"
}
reuse_step_tables <- if(is.null(params$reuse_step_tables)) TRUE else isTRUE(params$reuse_step_tables)
strict_mode <- isTRUE(params$strict_mode)
atac_significance_mode <- tolower(first_nonempty_value(Sys.getenv("ATAC_SIGNIFICANCE_MODE", unset = ""), "pvalue"))
if(!atac_significance_mode %in% c("pvalue", "padj", "fdr")) {
  stop("ATAC_SIGNIFICANCE_MODE 必须是 pvalue、padj 或 fdr，当前为: ", atac_significance_mode)
}
atac_significance_column <- if(atac_significance_mode == "pvalue") "pvalue" else "FDR"
atac_primary_require_fc <- env_bool_checked("ATAC_PRIMARY_REQUIRE_FC", FALSE)
atac_primary_p_cutoff <- suppressWarnings(as.numeric(first_nonempty_value(Sys.getenv("ATAC_PRIMARY_P_CUTOFF", unset = ""), "0.05")))
atac_diff_fdr_cutoff <- suppressWarnings(as.numeric(first_nonempty_value(Sys.getenv("ATAC_DIFF_FDR_CUTOFF", unset = ""), "0.05")))
atac_diff_log2fc_cutoff <- suppressWarnings(as.numeric(first_nonempty_value(Sys.getenv("ATAC_DIFF_LOG2FC_CUTOFF", unset = ""), "0.585")))
atac_diffbind_rpkm_filter <- suppressWarnings(as.numeric(first_nonempty_value(Sys.getenv("ATAC_DIFFBIND_RPKM_FILTER", unset = ""), "1")))
atac_diffbind_min_overlap <- suppressWarnings(as.integer(first_nonempty_value(Sys.getenv("ATAC_DIFFBIND_MIN_OVERLAP", unset = ""), "2")))
atac_diffbind_summits <- suppressWarnings(as.integer(first_nonempty_value(Sys.getenv("ATAC_DIFFBIND_SUMMITS", unset = ""), "200")))
atac_diffbind_min_count <- suppressWarnings(as.integer(first_nonempty_value(Sys.getenv("ATAC_DIFFBIND_MIN_COUNT", unset = ""), "0")))
atac_diffbind_filter_fun_name <- tolower(first_nonempty_value(Sys.getenv("ATAC_DIFFBIND_FILTER_FUN", unset = ""), "max"))
if(!atac_diffbind_filter_fun_name %in% c("max", "mean")) stop("ATAC_DIFFBIND_FILTER_FUN 必须是 max 或 mean")
atac_diffbind_filter_fun <- switch(atac_diffbind_filter_fun_name, max = base::max, mean = base::mean)
atac_nominal_p_cutoff <- suppressWarnings(as.numeric(first_nonempty_value(Sys.getenv("ATAC_NOMINAL_P_CUTOFF", unset = ""), as.character(atac_primary_p_cutoff))))
if(is.na(atac_primary_p_cutoff) || atac_primary_p_cutoff <= 0 || atac_primary_p_cutoff > 1) stop("ATAC_PRIMARY_P_CUTOFF 必须在 (0,1] 内")
if(is.na(atac_diff_fdr_cutoff) || atac_diff_fdr_cutoff <= 0 || atac_diff_fdr_cutoff > 1) stop("ATAC_DIFF_FDR_CUTOFF 必须在 (0,1] 内")
if(is.na(atac_diff_log2fc_cutoff) || atac_diff_log2fc_cutoff < 0) stop("ATAC_DIFF_LOG2FC_CUTOFF 必须是非负数")
if(is.na(atac_diffbind_rpkm_filter) || atac_diffbind_rpkm_filter < 0) stop("ATAC_DIFFBIND_RPKM_FILTER 必须是非负数")
if(is.na(atac_diffbind_min_overlap) || atac_diffbind_min_overlap < 1) stop("ATAC_DIFFBIND_MIN_OVERLAP 必须是正整数")
if(is.na(atac_diffbind_summits) || atac_diffbind_summits < 0) stop("ATAC_DIFFBIND_SUMMITS 必须是非负整数")
if(is.na(atac_diffbind_min_count) || atac_diffbind_min_count < 0) stop("ATAC_DIFFBIND_MIN_COUNT 必须是非负整数")
if(is.na(atac_nominal_p_cutoff) || atac_nominal_p_cutoff <= 0 || atac_nominal_p_cutoff > 1) stop("ATAC_NOMINAL_P_CUTOFF 必须在 (0,1] 内")

if(!dir.exists(work_dir)) {
  stop("ATAC 工作目录不存在: ", work_dir)
}
knitr::opts_knit$set(root.dir = work_dir)
options(stringsAsFactors = FALSE)
report_status("INFO", paste("ATAC work dir:", work_dir))
report_status(
  "INFO",
  paste(
    "Emit light plots in main report:",
    ifelse(emit_light_plots, "TRUE", "FALSE")
  )
)
report_status(
  "INFO",
  paste(
    "Save PNG previews for plots:",
    ifelse(save_png_previews, "TRUE", "FALSE")
  )
)
report_status("INFO", paste("PNG DPI:", plot_png_dpi))
report_status("INFO", paste("Run profile:", run_profile))
report_status("INFO", paste("Light setup:", light_setup))
report_status("INFO", paste("Reuse existing step tables:", reuse_step_tables))
report_status("INFO", paste("Strict mode:", strict_mode))
report_status(
  "INFO",
  paste0(
    "Primary differential peak threshold: ", atac_significance_column,
    if(atac_significance_column == "pvalue") " < " else " <= ",
    if(atac_significance_column == "pvalue") atac_primary_p_cutoff else atac_diff_fdr_cutoff,
    if(atac_primary_require_fc) paste0(" and |log2FC| >= ", atac_diff_log2fc_cutoff) else ""
  )
)
report_status(
  "INFO",
  paste(
    "Primary fold-change filter:",
    if(atac_primary_require_fc) paste0("enabled at |log2FC| >= ", round(atac_diff_log2fc_cutoff, 3)) else "disabled"
  )
)
report_status("INFO", paste("RNA-seq dir:", rna_seq_dir))
report_status("INFO", paste("WGBS dir:", wgbs_dir))

# 样本信息
resolve_atac_path <- function(path_value, root_dir = work_dir) {
  path_value <- trimws(as.character(path_value))
  if(length(path_value) == 0 || is.na(path_value) || !nzchar(path_value)) return(NA_character_)
  if(grepl("^/", path_value)) {
    normalizePath(path_value, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(root_dir, path_value), winslash = "/", mustWork = FALSE)
  }
}

read_atac_reads_info <- function(path, root_dir = work_dir, require_fastq_exists = TRUE) {
  if(!file.exists(path)) {
    stop("未找到 ATAC reads_info.txt: ", path)
  }
  first_line <- readLines(path, n = 1, warn = FALSE)
  if(length(first_line) == 0 || !nzchar(first_line)) {
    stop("ATAC reads_info.txt 为空: ", path)
  }
  info <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  if(!"sample" %in% colnames(info)) stop("ATAC reads_info.txt 缺少 sample 列")
  if(!"group" %in% colnames(info)) stop("ATAC reads_info.txt 缺少 group 列")
  if(!all(c("raw_fq1", "raw_fq2") %in% colnames(info))) {
    stop("ATAC reads_info.txt 缺少 raw_fq1/raw_fq2 列")
  }
  if(!"stage" %in% colnames(info) || anyNA(info$stage) || any(!nzchar(info$stage))) stop("必须显式提供 stage")
  if(!"group_stage" %in% colnames(info)) info$group_stage <- paste0(info$group, info$stage)

  trim_cols <- intersect(c("group_stage", "sample", "group", "stage", "raw_fq1", "raw_fq2"), colnames(info))
  for(col_name in trim_cols) {
    info[[col_name]] <- trimws(as.character(info[[col_name]]))
  }
  info <- info %>% filter(nzchar(sample) | nzchar(group) | nzchar(raw_fq1) | nzchar(raw_fq2))
  if(nrow(info) == 0) stop("ATAC reads_info.txt 中没有可用样本")
  if(any(!nzchar(info$sample))) stop("ATAC reads_info.txt 中存在空 sample")
  if(any(!nzchar(info$group))) stop("ATAC reads_info.txt 中存在空 group")
  if(any(!nzchar(info$raw_fq1) | !nzchar(info$raw_fq2))) {
    stop("ATAC reads_info.txt 中存在缺失 raw_fq1/raw_fq2 的样本")
  }
  if(anyDuplicated(info$sample)) {
    dup_samples <- unique(info$sample[duplicated(info$sample)])
    stop("ATAC reads_info.txt 中存在重复 sample: ", paste(dup_samples, collapse = ", "))
  }
  group_levels <- unique(info$group)
  if(length(group_levels) != 2) {
    stop("ATAC reads_info.txt 必须恰好包含两个分组，当前得到: ", paste(group_levels, collapse = ", "))
  }
  info$raw_fq1_resolved <- vapply(info$raw_fq1, resolve_atac_path, character(1), root_dir = root_dir)
  info$raw_fq2_resolved <- vapply(info$raw_fq2, resolve_atac_path, character(1), root_dir = root_dir)
  if(isTRUE(require_fastq_exists)) {
    missing_fastq <- c(info$raw_fq1_resolved[!file.exists(info$raw_fq1_resolved)], info$raw_fq2_resolved[!file.exists(info$raw_fq2_resolved)])
    if(length(missing_fastq) > 0) {
      stop("ATAC reads_info.txt 指向的 FASTQ 不存在: ", paste(missing_fastq, collapse = ", "))
    }
  }
  rownames(info) <- info$sample
  info
}

reads_info_file <- Sys.getenv("ATAC_SAMPLES")
# Rmd stages consume BAM/peak/count inputs; raw FASTQ presence is audited
# separately and remains required by the upstream FASTQ pipeline.
reads_info <- read_atac_reads_info(reads_info_file, work_dir, require_fastq_exists = FALSE)
if(!all(file.exists(c(reads_info$raw_fq1_resolved, reads_info$raw_fq2_resolved)))) {
  report_status("WARN", "Raw FASTQ inputs are missing: downstream BAM analysis only; upstream provenance remains STALE.")
}
manifest_file <- file.path(work_dir, "metadata", "sample_manifest.tsv")
sample_manifest <- reads_info %>%
  transmute(
    sample,
    group,
    stage,
    group_stage,
    raw_r1 = raw_fq1,
    raw_r2 = raw_fq2,
    raw_r1_resolved = raw_fq1_resolved,
    raw_r2_resolved = raw_fq2_resolved,
    clean_r1 = file.path("01rawdata", "cleaned_data", paste0(sample, "_1.fastq.gz")),
    clean_r2 = file.path("01rawdata", "cleaned_data", paste0(sample, "_2.fastq.gz")),
    bam = file.path("02mapping", "bam", paste0(sample, ".final.bam")),
    peak = file.path("03peaks", paste0(sample, "_peaks.clean.narrowPeak")),
    marked_bam = file.path("02mapping", "bam", paste0(sample, ".marked.bam")),
    bigwig = file.path("05bigwig", paste0(sample, "_RPGC.bw"))
  )

samples <- sample_manifest$sample
observed_groups <- unique(as.character(sample_manifest$group))
reference_group_config <- trimws(Sys.getenv("ATAC_REFERENCE_GROUP", unset = ""))
comparison_group_config <- trimws(Sys.getenv("ATAC_COMPARISON_GROUP", unset = ""))
if(xor(nzchar(reference_group_config), nzchar(comparison_group_config))) {
  stop("ATAC_REFERENCE_GROUP 与 ATAC_COMPARISON_GROUP 必须同时设置")
}
if(nzchar(reference_group_config)) {
  if(identical(reference_group_config, comparison_group_config)) stop("reference 与 comparison 分组不能相同")
  missing_groups <- setdiff(c(reference_group_config, comparison_group_config), observed_groups)
  if(length(missing_groups) > 0) stop("对比分组不在样本表中: ", paste(missing_groups, collapse = ", "))
  group_levels <- c(reference_group_config, comparison_group_config)
} else {
  stop("必须显式设置 comparison 和 reference，禁止按名称或行序推断")
}
baseline_group <- group_levels[1]
treatment_group <- group_levels[2]
reference_group <- baseline_group
comparison_group <- treatment_group
contrast_label <- paste(treatment_group, "vs", baseline_group)
conditions <- factor(sample_manifest$group, levels = group_levels)
names(conditions) <- samples
sample_group_map <- stats::setNames(as.character(sample_manifest$group), sample_manifest$sample)

# Peak 文件路径 - 使用 blacklist 过滤后的 clean peaks
peak_files <- file.path(work_dir, "03peaks", paste0(samples, "_peaks.clean.narrowPeak"))
names(peak_files) <- samples

# 同时读取原始 peaks 用于比较（可选）
raw_peak_files <- file.path(work_dir, "03peaks", paste0(samples, "_peaks.narrowPeak"))
names(raw_peak_files) <- samples

# Bam 文件路径 (用于 DiffBind)
bam_files <- file.path(work_dir, "02mapping/bam", paste0(samples, ".final.bam"))
names(bam_files) <- samples

# BigWig 路径仅接受明确标注的 RPGC 产物。
canonical_bigwig_files <- file.path(work_dir, "05bigwig", paste0(samples, "_RPGC.bw"))
bigwig_files <- canonical_bigwig_files
names(bigwig_files) <- samples
bigwig_method_status <- data.frame(
  SampleID = samples,
  BigWig = bigwig_files,
  Normalization = "RPGC",
  NamingStatus = ifelse(file.exists(canonical_bigwig_files), "CURRENT", "MISSING"),
  stringsAsFactors = FALSE
)

# 检查文件是否存在
if(!all(file.exists(peak_files))) {
  report_status("WARN", "部分 clean peak 文件未找到，请确认上游 Step6 是否完成。")
}
if(!all(file.exists(bam_files))) {
  report_status("WARN", "部分 BAM 文件未找到，请确认上游 mapping 步骤是否完成。")
}

file_freshness_record <- function(paths) {
  paths <- unique(normalizePath(as.character(paths), winslash = "/", mustWork = FALSE))
  info <- file.info(paths)
  data.frame(
    path = paths,
    exists = file.exists(paths),
    size = ifelse(file.exists(paths), info$size, NA_real_),
    mtime = ifelse(file.exists(paths), format(info$mtime, "%Y-%m-%dT%H:%M:%OS6%z"), NA_character_),
    stringsAsFactors = FALSE
  )
}

analysis_signature_inputs <- file_freshness_record(c(
  reads_info_file,
  bam_files,
  peak_files,
  file.path(work_dir, "04counts", c(
    "consensus_peak_counts.raw.tsv", "consensus_peak_counts.cpm.tsv",
    "consensus_peak_counts.library_size.tsv", "Consensus_Count_Matrix_Status.tsv",
    "Consensus_Count_Matrix_Method.tsv"
  )),
  file.path(work_dir, "run_upstream.sh"),
  file.path(work_dir, "run_rmd.sh"),
  list.files(rmd_dir, pattern = "\\.(R|Rmd)$", full.names = TRUE)
))
analysis_signature_config <- c(
  significance_mode = atac_significance_mode,
  significance_column = atac_significance_column,
  primary_p_cutoff = format(atac_primary_p_cutoff, digits = 16),
  fdr_cutoff = format(atac_diff_fdr_cutoff, digits = 16),
  log2fc_cutoff = format(atac_diff_log2fc_cutoff, digits = 16),
  primary_require_fc = as.character(atac_primary_require_fc),
  reference_group = baseline_group,
  comparison_group = treatment_group,
  diffbind_rpkm_filter = format(atac_diffbind_rpkm_filter, digits = 16),
  diffbind_min_overlap = as.character(atac_diffbind_min_overlap),
  diffbind_summits = as.character(atac_diffbind_summits),
  diffbind_min_count = as.character(atac_diffbind_min_count),
  diffbind_filter_fun = atac_diffbind_filter_fun_name,
  diffbind_version = if(requireNamespace("DiffBind", quietly = TRUE)) as.character(utils::packageVersion("DiffBind")) else "MISSING"
)
analysis_signature_text <- paste(
  paste(names(analysis_signature_config), analysis_signature_config, sep = "="),
  apply(analysis_signature_inputs, 1, function(x) paste(x, collapse = "|")),
  collapse = "\n"
)
current_analysis_signature <- if(requireNamespace("digest", quietly = TRUE)) {
  digest::digest(analysis_signature_text, algo = "sha256", serialize = FALSE)
} else {
  analysis_signature_text
}
analysis_freshness_tbl <- transform(
  analysis_signature_inputs,
  AnalysisSignature = current_analysis_signature,
  Status = ifelse(exists, "CURRENT_INPUT", "MISSING_INPUT")
)
report_status("INFO", paste("Analysis signature:", current_analysis_signature))

# 检查 blacklist 是否已应用
blacklist_file <- file.path(reference_dir, "GRCm39", "mm39.blacklist.bed")
if(file.exists(blacklist_file)) {
  report_status("INFO", paste("Blacklist 文件存在:", blacklist_file))
} else {
  report_status("WARN", "Blacklist 文件未找到；不能声称已完成 blacklist 过滤。")
}

step_dirs <- paste0("Step", 2:16)
step_output_dirs <- stats::setNames(
  c(
    "02_project_inventory",
    "03_atac_qc",
    "04_peak_catalog_annotation",
    "05_differential_accessibility",
    "06_genome_tracks",
    "07_peak_gene_enrichment",
    "08_motif_activity",
    "09_rna_multiomics",
    "10_footprinting_tf_evidence",
    "11_wgbs_integration",
    "12_tf_gene_network",
    "13_sequence_conservation_gc",
    "14_machine_learning",
    "15_keygenes_locus_reports",
    "16_export_summary_manifest"
  ),
  step_dirs
)

step_result_dir <- function(step) {
  step <- as.character(step)[1]
  step_dir_name <- unname(step_output_dirs[step])
  if(is.na(step_dir_name) || !nzchar(step_dir_name)) {
    stop("Unknown downstream result step: ", step)
  }
  dir_path <- file.path(result_dir, step_dir_name)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

# Result module directories are created just in time by step_table_dir() and
# step_plot_dir(); run_rmd.sh cleans only the active module before rendering.
dir.create(cache_root_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 建立 TxDb 对象 (mm39). Inventory/export-only refreshes can use light setup.
txdb <- if(isTRUE(light_setup)) NULL else TxDb.Mmusculus.UCSC.mm39.refGene

# --- feature levels; concrete colors are defined by the template palette layer below. ---
feature_levels_main <- c("Promoter", "UTR", "Exon", "Intron", "Downstream", "Intergenic", "Other")
feature_layer_levels <- c("Promoter", "Gene Body", "Distal", "Other")

save_csv_df <- function(x, file) {
  if(is.null(x)) return(invisible(NULL))
  write.csv(as.data.frame(x), file, row.names = FALSE)
  invisible(file)
}

plot_png_file <- function(file) {
  sub("\\.pdf$", ".png", file, ignore.case = TRUE)
}

atac_save_object_dual <- function(plot_obj, file, width = 10, height = 8,
                                  standardize = FALSE, bg = "white",
                                  coordinate_system = c("axes", "none")) {
  if(is.null(plot_obj) || is.null(file)) return(invisible(NULL))
  coordinate_system <- match.arg(coordinate_system)
  dims <- atac_figure_dimensions(width, height)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if(isTRUE(standardize) && inherits(plot_obj, c("ggplot", "patchwork"))) {
    plot_obj <- atac_standardize_plot(plot_obj, coordinate_system = coordinate_system)
  }
  ggplot2::ggsave(
    filename = file,
    plot = plot_obj,
    width = dims$pdf_width_mm,
    height = dims$pdf_height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    family = atac_figure_font_family,
    limitsize = FALSE,
    bg = "white"
  )
  if(isTRUE(save_png_previews)) {
    png_file <- plot_png_file(file)
    ggplot2::ggsave(
      filename = png_file,
      plot = plot_obj,
      width = dims$width_mm,
      height = dims$height_mm,
      units = "mm",
      dpi = atac_figure_dpi,
      device = "png",
      type = "cairo-png",
      limitsize = FALSE,
      bg = "white"
    )
    atac_finalize_rgb_png(png_file)
  }
  invisible(file)
}

save_gg_dual <- function(plot_obj, file = NULL, width = 10, height = 8,
                         dpi = plot_png_dpi, bg = "white",
                         coordinate_system = c("axes", "none")) {
  if(is.character(plot_obj) && length(plot_obj) == 1 && !is.null(file)) {
    tmp <- plot_obj
    plot_obj <- file
    file <- tmp
  }
  if(is.null(plot_obj) || is.null(file)) return(invisible(NULL))
  coordinate_system <- match.arg(coordinate_system)
  atac_save_object_dual(
    plot_obj = plot_obj,
    file = file,
    width = width,
    height = height,
    standardize = TRUE,
    bg = bg,
    coordinate_system = coordinate_system
  )
}

draw_base_plot_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi, bg = "white") {
  expr_sub <- substitute(expr)
  pdf_device <- atac_open_pdf(file, width = width, height = height)
  tryCatch({
    graphics::par(
      family = atac_figure_font_family,
      fg = "#333333", col.axis = "#333333", col.lab = "#333333", col.main = "#333333",
      cex.axis = atac_figure_tick_size / atac_figure_base_size,
      cex.lab = 1, cex.main = 1, bty = "l"
    )
    eval.parent(expr_sub)
  }, finally = atac_close_device(pdf_device))
  if(isTRUE(save_png_previews)) {
    png_file <- plot_png_file(file)
    png_device <- atac_open_png(png_file, width = width, height = height)
    tryCatch({
      graphics::par(
        family = atac_figure_font_family,
        fg = "#333333", col.axis = "#333333", col.lab = "#333333", col.main = "#333333",
        cex.axis = atac_figure_tick_size / atac_figure_base_size,
        cex.lab = 1, cex.main = 1, bty = "l"
      )
      eval.parent(expr_sub)
    }, finally = atac_close_device(png_device))
    atac_finalize_rgb_png(png_file)
  }
  invisible(file)
}

call_head_name <- function(expr) {
  if(!is.call(expr)) return("")
  head <- expr[[1]]
  if(is.symbol(head)) return(as.character(head))
  paste(as.character(head), collapse = "::")
}

add_call_arg_if_missing <- function(expr, target_names, arg_name, arg_value) {
  if(!is.call(expr)) return(expr)
  if(call_head_name(expr) %in% target_names && !arg_name %in% names(as.list(expr))[-1]) {
    expr[[arg_name]] <- arg_value
    return(expr)
  }
  for(i in seq_along(expr)[-1]) {
    expr[[i]] <- add_call_arg_if_missing(expr[[i]], target_names, arg_name, arg_value)
  }
  expr
}

force_pheatmap_silent <- function(expr) {
  expr <- add_call_arg_if_missing(expr, c("pheatmap"), "silent", TRUE)
  expr <- add_call_arg_if_missing(expr, c("pheatmap"), "fontfamily", atac_figure_font_family)
  add_call_arg_if_missing(expr, c("pheatmap"), "fontsize", atac_figure_base_size)
}

draw_plot_object <- function(plot_obj) {
  if(is.null(plot_obj)) return(invisible(FALSE))
  if(inherits(plot_obj, "ggplot")) {
    print(atac_standardize_plot(plot_obj))
    return(invisible(TRUE))
  }
  if(inherits(plot_obj, "upset")) {
    print(plot_obj)
    return(invisible(TRUE))
  }
  if(inherits(plot_obj, "pheatmap") && !is.null(plot_obj$gtable)) {
    grid::grid.newpage()
    grid::grid.draw(plot_obj$gtable)
    return(invisible(TRUE))
  }
  if(inherits(plot_obj, c("gtable", "gTree", "grob"))) {
    grid::grid.newpage()
    grid::grid.draw(plot_obj)
    return(invisible(TRUE))
  }
  if(inherits(plot_obj, "recordedplot")) {
    replayPlot(plot_obj)
    return(invisible(TRUE))
  }
  if(inherits(plot_obj, c("Heatmap", "HeatmapList")) && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    ComplexHeatmap::draw(plot_obj)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

png_rgb_sd <- function(png_file) {
  if(!file.exists(png_file) || !requireNamespace("png", quietly = TRUE)) return(NA_real_)
  img <- tryCatch(png::readPNG(png_file), error = function(e) NULL)
  if(is.null(img)) return(NA_real_)
  if(length(dim(img)) == 3 && dim(img)[3] >= 3) {
    img <- img[, , 1:3, drop = FALSE]
  }
  stats::sd(as.numeric(img), na.rm = TRUE)
}

write_pdf_first_page_png <- function(pdf_file, png_file, density = plot_png_dpi, page = 1L) {
  if(!file.exists(pdf_file)) return(FALSE)
  page <- suppressWarnings(as.integer(page)[1])
  if(is.na(page) || page < 1) page <- 1L
  ok <- FALSE
  if(nzchar(Sys.which("gs"))) {
    gs_args <- c(
      "-dSAFER", "-dBATCH", "-dNOPAUSE",
      "-sDEVICE=png16m",
      paste0("-r", density),
      paste0("-dFirstPage=", page), paste0("-dLastPage=", page),
      paste0("-sOutputFile=", png_file),
      pdf_file
    )
    ok <- tryCatch(
      identical(system2(Sys.which("gs"), gs_args, stdout = FALSE, stderr = FALSE), 0L),
      error = function(e) FALSE
    )
  }

  if(!ok && requireNamespace("magick", quietly = TRUE)) {
    ok <- tryCatch({
      img <- magick::image_read_pdf(pdf_file, density = density)
      if(length(img) >= page) {
        magick::image_write(img[page], path = png_file, format = "png")
        TRUE
      } else {
        FALSE
      }
    }, error = function(e) FALSE)
  }

  if(ok && file.exists(png_file) && requireNamespace("magick", quietly = TRUE)) {
    ok <- tryCatch({
      img <- magick::image_read(png_file)
      info <- magick::image_info(img)[1, , drop = FALSE]
      target_width <- as.integer(round(atac_figure_width_mm / 25.4 * density))
      target_height <- as.integer(round(info$height[[1]] * target_width / info$width[[1]]))
      if(info$width[[1]] != target_width) {
        img <- magick::image_resize(img, sprintf("%dx%d!", target_width, target_height))
      }
      img <- magick::image_background(img, color = "white", flatten = TRUE)
      magick::image_write(
        img,
        path = png_file,
        format = "png",
        depth = 8,
        density = sprintf("%dx%d", density, density)
      )
      TRUE
    }, error = function(e) FALSE)
  }
  isTRUE(ok) && file.exists(png_file) && file.info(png_file)$size > 0
}

pdf_render_pixel_sd <- function(pdf_file, density = 72, max_pages = plot_pdf_preview_max_pages,
                                min_png_sd = 0.001) {
  if(!file.exists(pdf_file)) return(NA_real_)
  tmp_files <- character()
  on.exit(unlink(tmp_files, force = TRUE), add = TRUE)
  sd_values <- numeric()
  for(page in seq_len(max_pages)) {
    tmp_png <- tempfile(fileext = ".png")
    tmp_files <- c(tmp_files, tmp_png)
    if(!write_pdf_first_page_png(pdf_file, tmp_png, density = density, page = page)) next
    page_sd <- png_rgb_sd(tmp_png)
    if(!is.na(page_sd)) sd_values <- c(sd_values, page_sd)
    if(!is.na(page_sd) && page_sd >= min_png_sd) return(page_sd)
  }
  if(length(sd_values) > 0 && any(!is.na(sd_values))) return(max(sd_values, na.rm = TRUE))
  NA_real_
}

render_pdf_preview_png <- function(pdf_file, png_file, density = plot_png_dpi, min_png_sd = 0.001,
                                   max_pages = plot_pdf_preview_max_pages) {
  tmp_files <- character()
  on.exit(unlink(tmp_files, force = TRUE), add = TRUE)
  for(page in seq_len(max_pages)) {
    tmp_png <- tempfile(fileext = ".png")
    tmp_files <- c(tmp_files, tmp_png)
    if(!write_pdf_first_page_png(pdf_file, tmp_png, density = density, page = page)) next
    rendered_sd <- png_rgb_sd(tmp_png)
    if(is.na(rendered_sd) || rendered_sd < min_png_sd) next
    return(isTRUE(file.copy(tmp_png, png_file, overwrite = TRUE)))
  }
  FALSE
}

atac_exportable_plot <- function(plot_obj, width, height) {
  if(is.null(plot_obj)) stop("Plot expression returned NULL.", call. = FALSE)
  if(inherits(plot_obj, "pheatmap") && !is.null(plot_obj$gtable)) return(plot_obj$gtable)
  if(inherits(plot_obj, c("ggplot", "patchwork", "gtable", "gTree", "grob"))) return(plot_obj)
  if(inherits(plot_obj, c("upset", "recordedplot", "Heatmap", "HeatmapList"))) {
    dims <- atac_figure_dimensions(width, height)
    captured <- grid::grid.grabExpr(
      draw_plot_object(plot_obj),
      width = dims$width_in,
      height = dims$height_in,
      wrap = TRUE,
      wrap.grobs = TRUE
    )
    if(!is.null(captured)) return(captured)
  }
  stop("Unsupported plot object class: ", paste(class(plot_obj), collapse = ", "), call. = FALSE)
}

draw_object_call_dual <- function(file, width = 10, height = 8, expr_sub, expr_env, dpi = plot_png_dpi,
                                  bg = "white") {
  plot_obj <- eval(expr_sub, envir = expr_env)
  export_obj <- atac_exportable_plot(plot_obj, width = width, height = height)
  atac_save_object_dual(
    plot_obj = export_obj,
    file = file,
    width = width,
    height = height,
    standardize = inherits(export_obj, c("ggplot", "patchwork")),
    bg = bg
  )
}

draw_object_expression_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi,
                                        bg = "white", expr_transform = identity) {
  draw_object_call_dual(
    file = file,
    width = width,
    height = height,
    expr_sub = expr_transform(substitute(expr)),
    expr_env = parent.frame(),
    dpi = dpi,
    bg = bg
  )
}

draw_pheatmap_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi, bg = "white") {
  expr_sub <- force_pheatmap_silent(substitute(expr))
  plot_obj <- eval(expr_sub, envir = parent.frame())
  export_obj <- atac_exportable_plot(plot_obj, width = width, height = height)
  atac_save_object_dual(
    plot_obj = export_obj,
    file = file,
    width = width,
    height = height,
    standardize = inherits(export_obj, c("ggplot", "patchwork")),
    bg = bg
  )
}

draw_grid_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi, bg = "white") {
  draw_object_call_dual(
    file = file,
    width = width,
    height = height,
    expr_sub = substitute(expr),
    expr_env = parent.frame(),
    dpi = dpi,
    bg = bg
  )
}

draw_upset_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi, bg = "white") {
  draw_object_call_dual(
    file = file,
    width = width,
    height = height,
    expr_sub = substitute(expr),
    expr_env = parent.frame(),
    dpi = dpi,
    bg = bg
  )
}

draw_gviz_dual <- function(file, width = 10, height = 8, expr, dpi = plot_png_dpi, bg = "white") {
  draw_base_plot_dual(file = file, width = width, height = height, expr = expr, dpi = dpi, bg = bg)
  if(isTRUE(save_png_previews) && !file.exists(plot_png_file(file))) {
    sync_pdf_preview_png(file, force = TRUE, density = dpi)
  }
  invisible(file)
}

pal_pick <- function(x) {
  if(length(x) == 0 || all(is.na(x))) return(NA_character_)
  unname(as.character(x)[1])
}

pal_named <- function(...) {
  vals <- list(...)
  stats::setNames(vapply(vals, pal_pick, character(1)), names(vals))
}

sync_pdf_preview_png <- function(pdf_file, force = FALSE, density = plot_png_dpi) {
  if(!isTRUE(save_png_previews) || !file.exists(pdf_file)) return(invisible(""))
  png_file <- plot_png_file(pdf_file)
  if(
    !force &&
    file.exists(png_file) &&
    !is.na(file.info(png_file)$mtime) &&
    !is.na(file.info(pdf_file)$mtime) &&
    file.info(png_file)$mtime >= file.info(pdf_file)$mtime
  ) {
    return(invisible(png_file))
  }

  old_png_sd <- png_rgb_sd(png_file)
  ok <- render_pdf_preview_png(pdf_file, png_file, density = density)
  if(!ok && !is.na(old_png_sd) && old_png_sd >= 0.001 && file.exists(png_file)) {
    return(invisible(png_file))
  }

  if(!ok && file.exists(png_file)) unlink(png_file, force = TRUE)
  invisible(if(ok) png_file else "")
}

sync_plot_previews <- function(plot_root_dir, force = FALSE, density = plot_png_dpi) {
  if(!isTRUE(save_png_previews) || !dir.exists(plot_root_dir)) return(invisible(character(0)))
  pdf_files <- list.files(plot_root_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
  if(length(pdf_files) == 0) return(invisible(character(0)))
  png_files <- vapply(pdf_files, sync_pdf_preview_png, character(1), force = force, density = density)
  invisible(unname(png_files[nzchar(png_files)]))
}

audit_plot_preview_pairs <- function(plot_root = plot_root_dir, step = "Step16") {
  empty_audit <- data.frame(
    PdfPath = character(),
    PngPath = character(),
    Step = character(),
    PdfFile = character(),
    PngFile = character(),
    PdfSizeBytes = numeric(),
    PngSizeBytes = numeric(),
    PngNewerOrSame = logical(),
    Status = character(),
    Reason = character(),
    stringsAsFactors = FALSE
  )
  if(!dir.exists(plot_root)) {
    save_step_table(empty_audit, step, "Plot_Preview_Pairing_Audit.csv")
    return(invisible(empty_audit))
  }

  pdf_files <- list.files(plot_root, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  pdf_files <- unique(pdf_files[file.exists(pdf_files)])
  if(length(pdf_files) == 0) {
    save_step_table(empty_audit, step, "Plot_Preview_Pairing_Audit.csv")
    return(invisible(empty_audit))
  }

  audit <- bind_rows(lapply(pdf_files, function(pdf_file) {
    png_file <- plot_png_file(pdf_file)
    pdf_info <- file.info(pdf_file)
    png_exists <- file.exists(png_file)
    png_info <- if(png_exists) file.info(png_file) else NULL
    png_size <- if(png_exists) as.numeric(png_info$size) else NA_real_
    png_newer <- if(png_exists) {
      !is.na(pdf_info$mtime) && !is.na(png_info$mtime) && png_info$mtime >= pdf_info$mtime
    } else {
      FALSE
    }
    status <- dplyr::case_when(
      !isTRUE(save_png_previews) ~ "SKIPPED_PREVIEWS_DISABLED",
      !png_exists ~ "MISSING_PNG",
      is.na(png_size) || png_size < 5000 ~ "FAIL_PNG_TOO_SMALL",
      !png_newer ~ "WARN_PNG_OLDER_THAN_PDF",
      TRUE ~ "PASS"
    )
    reason <- dplyr::case_when(
      status == "SKIPPED_PREVIEWS_DISABLED" ~ "ATAC_SAVE_PNG_PREVIEWS is FALSE",
      status == "MISSING_PNG" ~ "same-basename PNG preview not found",
      status == "FAIL_PNG_TOO_SMALL" ~ "PNG preview is below 5000 bytes",
      status == "WARN_PNG_OLDER_THAN_PDF" ~ "PNG exists but is older than the PDF",
      TRUE ~ "ok"
    )
    data.frame(
      PdfPath = normalizePath(pdf_file, winslash = "/", mustWork = FALSE),
      PngPath = normalizePath(png_file, winslash = "/", mustWork = FALSE),
      Step = basename(dirname(pdf_file)),
      PdfFile = basename(pdf_file),
      PngFile = basename(png_file),
      PdfSizeBytes = as.numeric(pdf_info$size),
      PngSizeBytes = png_size,
      PngNewerOrSame = png_newer,
      Status = status,
      Reason = reason,
      stringsAsFactors = FALSE
    )
  }))
  audit <- audit %>% arrange(Status, Step, PdfFile)
  save_step_table(audit, step, "Plot_Preview_Pairing_Audit.csv")
  invisible(audit)
}

save_plot_pdf <- function(plot_obj, file, width = 10, height = 8,
                          dpi = plot_png_dpi,
                          coordinate_system = c("axes", "none")) {
  coordinate_system <- match.arg(coordinate_system)
  save_gg_dual(
    plot_obj = plot_obj,
    file = file,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    coordinate_system = coordinate_system
  )
}

remove_output_with_preview <- function(file) {
  unlink(c(file, plot_png_file(file)), force = TRUE)
  invisible(file)
}

save_step_table <- function(x, step, file_name) {
  file <- step_table_file(step, file_name)
  if(is.null(x)) return(invisible(NULL))
  write.csv(as.data.frame(x), file, row.names = FALSE)
  invisible(file)
}

read_step_table <- function(step, file_name, ...) {
  bindetect_derived <- c(
    "Footprint_Regulator_Evidence.csv", "Integrated_Footprint_Regulator_Summary.csv",
    "Footprint_Supported_Motif_Target_Edges.csv", "Footprint_Supported_Motif_Summary.csv",
    "Footprint_Supported_Gene_Summary.csv", "TF_Gene_Concordance_Summary.csv",
    "TF_Evidence_Matrix.csv", "TF_Target_Module_Edges.csv"
  )
  if(identical(step, "Step10") && file_name %in% bindetect_derived &&
     !identical(validate_bindetect_provenance()$Status[1], "VERIFIED")) return(data.frame())
  file <- step_table_file(step, file_name)
  if(!file.exists(file)) return(data.frame())
  header <- readLines(file, n = 5L, warn = FALSE)
  if(length(header) == 0L || all(trimws(header) %in% c("", "\"\""))) return(data.frame())
  read.csv(file, stringsAsFactors = FALSE, ...)
}

read_tsv_if_exists <- function(file, ...) {
  if(!file.exists(file)) return(data.frame())
  utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

read_consensus_count_matrix <- function(type = c("raw", "cpm")) {
  type <- match.arg(type)
  file_name <- if(type == "raw") "consensus_peak_counts.raw.tsv" else "consensus_peak_counts.cpm.tsv"
  status <- read_tsv_if_exists(file.path(work_dir, "04counts", "Consensus_Count_Matrix_Status.tsv"))
  if("Status" %in% names(status) && any(grepl("FAILED|SKIPPED|STALE|RUNNING", status$Status, ignore.case = TRUE))) {
    report_status("WARN", "Consensus counts are not admitted: upstream count status is failed, skipped, stale or running")
    return(data.frame())
  }
  file <- file.path(work_dir, "04counts", file_name)
  mat <- read_tsv_if_exists(file)
  if(nrow(mat) == 0) return(mat)

  if(!"CHR" %in% colnames(mat) && "Chr" %in% colnames(mat)) mat$CHR <- mat$Chr
  if(!"START" %in% colnames(mat) && "Start" %in% colnames(mat)) mat$START <- mat$Start
  if(!"END" %in% colnames(mat) && "End" %in% colnames(mat)) mat$END <- mat$End
  if(!"peak_id" %in% colnames(mat) && all(c("CHR", "START", "END") %in% colnames(mat))) {
    mat$peak_id <- canonical_peak_id(mat$CHR, mat$START, mat$END)
  }
  if(!"peak_id_join" %in% colnames(mat) && all(c("CHR", "START", "END") %in% colnames(mat))) {
    mat$peak_id_join <- canonical_peak_id(mat$CHR, mat$START, mat$END)
  }
  if("CHR" %in% colnames(mat)) mat$seqnames <- ensure_ucsc_seqnames(mat$CHR)
  if("START" %in% colnames(mat)) {
    mat$START <- suppressWarnings(as.integer(mat$START))
    mat$start <- mat$START
    mat$START_1BASED <- mat$START + 1L
  }
  if("END" %in% colnames(mat)) {
    mat$END <- suppressWarnings(as.integer(mat$END))
    mat$end <- mat$END
  }

  sample_cols <- intersect(samples, colnames(mat))
  for(nm in sample_cols) mat[[nm]] <- suppressWarnings(as.numeric(mat[[nm]]))
  attr(mat, "source_file") <- file
  attr(mat, "matrix_type") <- type
  attr(mat, "coordinate_system") <- "BED_0_based_half_open"
  mat
}

read_consensus_count_status <- function() {
  read_tsv_if_exists(file.path(work_dir, "04counts", "Consensus_Count_Matrix_Status.tsv"))
}

read_consensus_count_method <- function() {
  read_tsv_if_exists(file.path(work_dir, "04counts", "Consensus_Count_Matrix_Method.tsv"))
}

read_consensus_count_library_size <- function() {
  read_tsv_if_exists(file.path(work_dir, "04counts", "consensus_peak_counts.library_size.tsv"))
}

build_stepwise_result_view <- function(result_dir = result_dir) {
  step_info <- data.frame(
    Module = sprintf("Module%02d", 2:16),
    Directory = unname(step_output_dirs[paste0("Step", 2:16)]),
    Title = c(
      "Project inventory and run context",
      "ATAC sample QC and integrated grading",
      "Peak catalog, annotation and coverage profiles",
      "Differential accessibility",
      "Genome tracks",
      "Peak-associated gene enrichment",
      "Motif enrichment and activity",
      "RNA multi-omics integration",
      "Footprinting and TF evidence",
      "WGBS integration",
      "TF-gene network evidence",
      "Sequence conservation and GC context",
      "Machine-learning feature ranking",
      "Key genes and locus reports",
      "Export and delivery summary"
    ),
    Note = c(
      "Input inventory, external dependency context, run configuration, and analysis route.",
      "Sample-level ATAC quality outputs: insert-size/NFR summaries, FRiP, per-sample peak counts, library complexity, TSS metadata, sequencing QC, and alignment/effective-read tables.",
      "Peak genomic annotation, per-sample length/depth/element distributions, group overlap, chromosome distribution, and coverage profiles.",
      "DiffBind count/PCA outputs with configurable primary differential peaks plus fixed p-value and FDR layers.",
      "Representative genome-browser style tracks for top differential accessibility loci.",
      "Gene annotation and GO/KEGG enrichment for differential peak-associated genes.",
      "Motif enrichment, motif activity, and direction-aware TF motif evidence.",
      "ATAC-RNA overlap, distal peak-gene links, WGCNA/splicing support, and evidence-weighted candidate scores.",
      "TOBIAS/internal footprint profiles, motif binding-site summaries, and regulator evidence.",
      "ATAC-WGBS overlap, gene-centric methylation support, Ensembl-first joins, and WGBS evidence scores.",
      "TF-target network outputs and footprint-supported regulator/gene summaries.",
      "Conservation, GC-bias, and high-GC peak context for differential regions.",
      "Supervised modelling eligibility and NO-GO status; no model fitted under the current sample design.",
      "Gene-centric evidence matrix, ranked key genes, representative loci, and locus report figures.",
      "Final delivery tables, input audits, track manifests, IGV files, and summary text."
    ),
    stringsAsFactors = FALSE
  )

  safe_symlink <- function(target_abs, link_path, target_rel) {
    if(!dir.exists(target_abs) && !file.exists(target_abs)) return(FALSE)
    existing_target <- suppressWarnings(Sys.readlink(link_path))
    if((file.exists(link_path) || dir.exists(link_path) || nzchar(existing_target)) && nzchar(existing_target)) {
      unlink(link_path, recursive = FALSE, force = TRUE)
    } else if(file.exists(link_path) || dir.exists(link_path)) {
      return(FALSE)
    }
    isTRUE(file.symlink(target_rel, link_path))
  }

  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  for(i in seq_len(nrow(step_info))) {
    info <- step_info[i, , drop = FALSE]
    step_dir <- file.path(result_dir, info$Directory)
    dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
    unlink(file.path(step_dir, "README.md"), force = TRUE)

    table_abs <- file.path(step_dir, "tables")
    plot_abs <- file.path(step_dir, "plots")
    has_tables <- dir.exists(table_abs)
    has_plots <- dir.exists(plot_abs)
    has_igv <- identical(info$Module, "Module16") && dir.exists(file.path(result_dir, "IGV"))

    if(!has_tables) dir.create(table_abs, recursive = TRUE, showWarnings = FALSE)
    if(!has_plots) dir.create(plot_abs, recursive = TRUE, showWarnings = FALSE)
    if(has_igv) {
      safe_symlink(file.path(result_dir, "IGV"), file.path(step_dir, "igv"), file.path("..", "IGV"))
    }
  }

  root_lines <- c(
    "# Liver ATAC result index",
    "",
    "Module folders are provided at the top level of this result directory.",
    "Plots are stored directly under each module folder's plots/ directory.",
    "Machine-readable tables are stored directly under each module folder's tables/ directory.",
    "",
    "Module folders:",
    paste0("- ", step_info$Directory, " (", step_info$Module, "): ", step_info$Title)
  )
  writeLines(root_lines, con = file.path(result_dir, "README.md"))

  step_info$Path <- file.path(result_dir, step_info$Directory)
  step_info
}

step_tables_ready <- function(step, files) {
  if(length(files) == 0) return(FALSE)
  all(file.exists(vapply(files, function(x) step_table_file(step, x), character(1))))
}

should_compute_step <- function(step, required_files = character()) {
  if(run_profile == "full") return(TRUE)
  if(run_profile == "report") {
    if(length(required_files) == 0) return(FALSE)
    if(!reuse_step_tables) return(TRUE)
    return(!step_tables_ready(step, required_files))
  }
  if(!reuse_step_tables) return(TRUE)
  !step_tables_ready(step, required_files)
}

should_render_plots <- function() {
  run_profile %in% c("full", "quick", "report")
}

write_skipped_state <- function(step, module_name, reason, clear_plot_files = character(), extra = list()) {
  if(length(clear_plot_files) > 0) {
    invisible(lapply(clear_plot_files, function(plot_name) remove_output_with_preview(step_plot_file(step, plot_name))))
  }
  skipped_df <- data.frame(
    module = module_name,
    status = "skipped",
    reason = reason,
    run_profile = run_profile,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  if(length(extra) > 0) {
    for(nm in names(extra)) skipped_df[[nm]] <- as.character(extra[[nm]])
  }
  save_step_table(skipped_df, step, paste0(module_name, "_Skipped_State.csv"))
  invisible(skipped_df)
}

write_module_status_plot <- function(step, title = NULL, file_name = "Module_Status.pdf", only_if_no_plots = TRUE) {
  if(!should_render_plots()) return(invisible(FALSE))

  plot_dir <- step_plot_dir(step)
  if(isTRUE(only_if_no_plots)) {
    existing_pdfs <- list.files(plot_dir, pattern = "\\.pdf$", recursive = FALSE, full.names = TRUE)
    existing_pdfs <- existing_pdfs[basename(existing_pdfs) != file_name]
    if(length(existing_pdfs) > 0) return(invisible(FALSE))
  }

  skip_files <- list.files(step_table_dir(step), pattern = "_Skipped_State\\.csv$", full.names = TRUE)
  if(length(skip_files) == 0) return(invisible(FALSE))

  skipped_rows <- bind_rows(lapply(skip_files, function(path) {
    tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  }))
  if(nrow(skipped_rows) == 0) return(invisible(FALSE))

  skipped_rows$module <- ifelse(is.na(skipped_rows$module) | skipped_rows$module == "", "Module component", skipped_rows$module)
  skipped_rows$reason <- ifelse(is.na(skipped_rows$reason) | skipped_rows$reason == "", "No reason recorded.", skipped_rows$reason)
  status_lines <- paste0(seq_len(nrow(skipped_rows)), ". ", skipped_rows$module, ": ", skipped_rows$reason)
  status_text <- paste(strwrap(paste(status_lines, collapse = "\n"), width = 100), collapse = "\n")
  if(is.null(title) || !nzchar(title)) title <- paste(step, "status")

  p_status <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 1,
      label = status_text,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      color = viz_pal$neutral["dark"],
      lineheight = 1.08
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 11) +
    labs(title = title, subtitle = "No biological plot was generated for this module because required inputs were unavailable.") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 13, color = viz_pal$neutral["dark"]),
      plot.subtitle = element_text(size = 10.5, color = viz_pal$neutral["dark"]),
      plot.margin = ggplot2::margin(18, 18, 18, 18)
    )

  plot_height <- max(4.5, min(10, 2.5 + 0.75 * nrow(skipped_rows)))
  save_gg_dual(p_status, step_plot_file(step, file_name), width = 10, height = plot_height)
  invisible(TRUE)
}

handle_nonfatal <- function(message_text, step = "GLOBAL") {
  if(strict_mode) {
    stop(sprintf("[%s] %s", step, message_text))
  }
  report_status("WARN", paste0("[", step, "] ", message_text))
}

step_plot_dir <- function(step) {
  dir_path <- file.path(step_result_dir(step), "plots")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

step_table_dir <- function(step) {
  dir_path <- file.path(step_result_dir(step), "tables")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  dir_path
}

step_plot_file <- function(step, file_name) {
  file.path(step_plot_dir(step), file_name)
}

step_table_file <- function(step, file_name) {
  file.path(step_table_dir(step), file_name)
}

plot_pixel_sd <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if(identical(ext, "png") && requireNamespace("png", quietly = TRUE)) {
    return(png_rgb_sd(file))
  }
  if(identical(ext, "pdf")) {
    preview_png <- plot_png_file(file)
    if(file.exists(preview_png)) {
      preview_sd <- png_rgb_sd(preview_png)
      if(!is.na(preview_sd)) return(preview_sd)
    }
    sd_value <- pdf_render_pixel_sd(file, density = 72)
    if(!is.na(sd_value)) return(sd_value)
    if(requireNamespace("magick", quietly = TRUE)) {
      sd_value <- tryCatch({
        img <- magick::image_read_pdf(file, density = 72)
        if(length(img) == 0) {
          NA_real_
        } else {
          page_sds <- vapply(seq_len(min(length(img), plot_pdf_preview_max_pages)), function(i) {
            raw_img <- magick::image_data(img[i], channels = "rgb")
            stats::sd(as.integer(raw_img), na.rm = TRUE) / 255
          }, numeric(1))
          nonblank <- page_sds[!is.na(page_sds) & page_sds >= 0.001]
          if(length(nonblank) > 0) nonblank[1] else max(page_sds, na.rm = TRUE)
        }
      }, error = function(e) NA_real_)
    }
    return(sd_value)
  }
  NA_real_
}

audit_plot_integrity <- function(plot_root = plot_root_dir, step = "Step16") {
  deep_plot_audit <- env_bool("ATAC_DEEP_PLOT_AUDIT", FALSE)
  if(!dir.exists(plot_root)) {
    audit <- data.frame(
      Path = character(),
      Step = character(),
      File = character(),
      Extension = character(),
      SizeBytes = numeric(),
      PixelSD = numeric(),
      MinExpectedBytes = numeric(),
      Status = character(),
      Reason = character(),
      stringsAsFactors = FALSE
    )
    save_step_table(audit, step, "Plot_Integrity_Audit.csv")
    return(invisible(audit))
  }
  plot_files <- list.files(plot_root, pattern = "\\.(pdf|png)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  plot_files <- unique(plot_files[file.exists(plot_files)])
  audit <- bind_rows(lapply(plot_files, function(path) {
    ext <- tolower(tools::file_ext(path))
    size_bytes <- as.numeric(file.info(path)$size)
    min_expected <- dplyr::case_when(
      ext == "pdf" ~ 4096,
      ext == "png" ~ 5000,
      TRUE ~ 1
    )
    pixel_sd <- if(deep_plot_audit && (size_bytes <= 250000 || ext == "png")) {
      plot_pixel_sd(path)
    } else {
      NA_real_
    }
    small_file <- is.na(size_bytes) || size_bytes < min_expected
    low_pixels <- !is.na(pixel_sd) && pixel_sd < 0.001
    status <- ifelse(small_file || low_pixels, "FAIL_POSSIBLE_BLANK", "PASS")
    reason <- paste(c(
      if(small_file) paste0("file_size_below_", min_expected, "_bytes") else character(),
      if(low_pixels) "near_blank_pixel_variance" else character()
    ), collapse = ";")
    if(!nzchar(reason)) reason <- "ok"
    data.frame(
      Path = normalizePath(path, winslash = "/", mustWork = FALSE),
      Step = basename(dirname(path)),
      File = basename(path),
      Extension = ext,
      SizeBytes = size_bytes,
      PixelSD = pixel_sd,
      MinExpectedBytes = min_expected,
      Status = status,
      Reason = reason,
      stringsAsFactors = FALSE
    )
  }))
  audit <- audit %>% arrange(Status, Step, File)
  save_step_table(audit, step, "Plot_Integrity_Audit.csv")
  invisible(audit)
}

ensure_columns <- function(df, spec) {
  if(is.null(df) || !is.data.frame(df)) {
    df <- data.frame(stringsAsFactors = FALSE)
  }
  n <- nrow(df)
  for(nm in names(spec)) {
    if(!nm %in% colnames(df)) {
      df[[nm]] <- rep(spec[[nm]], n)
    }
  }
  df
}

first_existing_file <- function(candidates) {
  existing <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if(length(existing) == 0) return(NA_character_)
  existing[1]
}

get_rna_input_specs <- function(...) {
  list(
    RNA_DE_Annotated = list(paths = Sys.getenv("ATAC_RNA_DE_FILE"), required = c("gene_id", "ensembl_id", "log2FoldChange", "pvalue", "padj", "ENTREZID", "sampleA", "sampleB")),
    RNA_Expression_Matrix = list(paths = Sys.getenv("ATAC_RNA_EXPRESSION_FILE"), required = c()),
    RNA_WGCNA_Module = list(paths = Sys.getenv("ATAC_RNA_WGCNA_MODULE_FILE"), required = c("ENTREZID", "module_label")),
    RNA_WGCNA_Trait = list(paths = Sys.getenv("ATAC_RNA_WGCNA_TRAIT_FILE"), required = c("Module", "correlation", "pvalue", "module_size")),
    RNA_Splicing_Candidates = list(paths = Sys.getenv("ATAC_RNA_SPLICING_FILE"), required = c("candidate_priority", "candidate_status")),
    RNA_GSEA_KEGG = list(paths = Sys.getenv("ATAC_RNA_GSEA_FILE"), required = c("Description", "NES", "pvalue", "core_enrichment"))
  )
}
get_wgbs_input_specs <- function(...) {
  list(
    WGBS_DMR_Annotated = list(paths = Sys.getenv("ATAC_WGBS_DMR_FILE"), required = c("seqnames", "start", "end", "diffMethy", "group1_label", "group2_label")),
    WGBS_Gene_Centric_Wide = list(paths = Sys.getenv("ATAC_WGBS_GENE_FILE"), required = c("gene_id_clean", "SYMBOL", "diffMethy_Promoter")),
    WGBS_RNA_KEGG_Comparison = list(paths = Sys.getenv("ATAC_WGBS_KEGG_FILE"), required = c("Description", "NES", "p.adjust", "core_enrichment")),
    WGBS_Representative_Loci = list(paths = Sys.getenv("ATAC_WGBS_LOCI_FILE"), required = c("SYMBOL"))
  )
}

external_input_path <- function(specs, input_name) {
  if(is.null(specs[[input_name]])) return(NA_character_)
  first_existing_file(specs[[input_name]]$paths)
}

get_rna_input_path <- function(input_name) {
  external_input_path(get_rna_input_specs(), input_name)
}

get_wgbs_input_path <- function(input_name) {
  external_input_path(get_wgbs_input_specs(), input_name)
}

read_rna_expression_matrix <- function(path) {
  if(is.na(path) || !file.exists(path) || dir.exists(path)) return(NULL)
  expr_mat <- read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- sub("\\..*$", "", rownames(expr_mat))
  expr_mat
}

summarise_rna_expression_by_group <- function(expr_mat, baseline = baseline_group, treatment = treatment_group) {
  if(is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    return(data.frame(
      ensembl_id = character(),
      RNA_MeanExpr_Reference = numeric(),
      RNA_MeanExpr_Comparison = numeric(),
      RNA_Expr_Delta = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  expr_groups <- sample_group(colnames(expr_mat))
  baseline_cols <- colnames(expr_mat)[expr_groups == baseline]
  treatment_cols <- colnames(expr_mat)[expr_groups == treatment]
  if(length(baseline_cols) == 0 || length(treatment_cols) == 0) {
    return(data.frame(
      ensembl_id = character(),
      RNA_MeanExpr_Reference = numeric(),
      RNA_MeanExpr_Comparison = numeric(),
      RNA_Expr_Delta = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    ensembl_id = rownames(expr_mat),
    RNA_MeanExpr_Reference = rowMeans(expr_mat[, baseline_cols, drop = FALSE], na.rm = TRUE),
    RNA_MeanExpr_Comparison = rowMeans(expr_mat[, treatment_cols, drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  ) %>%
    mutate(RNA_Expr_Delta = RNA_MeanExpr_Comparison - RNA_MeanExpr_Reference)
}

read_table_header <- function(path) {
  if(is.na(path) || !file.exists(path) || dir.exists(path)) return(character())
  header <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) character())
  if(length(header) == 0) return(character())
  delimiter <- if(grepl("\t", header, fixed = TRUE)) "\t" else ","
  out <- strsplit(header, delimiter, fixed = TRUE)[[1]]
  gsub('^"|"$', "", out)
}

count_file_rows <- function(path) {
  if(is.na(path) || !file.exists(path) || dir.exists(path)) return(NA_integer_)
  out <- tryCatch(system2("wc", c("-l", path), stdout = TRUE, stderr = FALSE), error = function(e) character())
  if(length(out) == 0) return(NA_integer_)
  suppressWarnings(as.integer(strsplit(trimws(out[1]), "\\s+")[[1]][1]))
}

audit_external_input_specs <- function(source_name, specs) {
  bind_rows(lapply(names(specs), function(input_name) {
    spec <- specs[[input_name]]
    selected <- first_existing_file(spec$paths)
    if(is.na(selected) && length(spec$paths) > 0) selected <- spec$paths[1]
    header <- read_table_header(selected)
    required <- unique(as.character(spec$required %||% character()))
    missing <- setdiff(required, header)
    data.frame(
      Source = source_name,
      Input = input_name,
      AnalysisUnit = as.character(spec$analysis_unit %||% "See assay-specific source metadata"),
      Path = selected,
      Exists = !is.na(selected) && file.exists(selected) && !dir.exists(selected),
      RowCount = count_file_rows(selected),
      RequiredColumns = paste(required, collapse = ";"),
      MissingColumns = paste(missing, collapse = ";"),
      KeyColumnsOK = length(missing) == 0,
      stringsAsFactors = FALSE
    )
  }))
}

build_external_input_inventory <- function() {
  bind_rows(
    audit_external_input_specs("RNA", get_rna_input_specs()),
    audit_external_input_specs("WGBS", get_wgbs_input_specs())
  )
}


filter_enrich_result <- function(enrich_obj, p_cutoff = 0.05) {
  if(is.null(enrich_obj)) return(NULL)
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) data.frame())
  p_col <- if("pvalue" %in% colnames(enrich_df)) "pvalue" else NA_character_
  if(nrow(enrich_df) == 0 || is.na(p_col)) {
    enrich_obj@result <- enrich_df[0, , drop = FALSE]
    return(enrich_obj)
  }
  p_vec <- suppressWarnings(as.numeric(enrich_df[[p_col]]))
  enrich_obj@result <- enrich_df[is.finite(p_vec) & p_vec >= 0 & p_vec < p_cutoff, , drop = FALSE]
  enrich_obj
}

filter_compare_cluster_result <- function(compare_obj, p_cutoff = 0.05) {
  if(is.null(compare_obj)) return(NULL)
  compare_df <- tryCatch(as.data.frame(compare_obj), error = function(e) data.frame())
  p_col <- if("pvalue" %in% colnames(compare_df)) "pvalue" else NA_character_
  if(nrow(compare_df) == 0 || is.na(p_col)) {
    compare_obj@compareClusterResult <- compare_df[0, , drop = FALSE]
    return(compare_obj)
  }
  p_vec <- suppressWarnings(as.numeric(compare_df[[p_col]]))
  compare_obj@compareClusterResult <- compare_df[is.finite(p_vec) & p_vec >= 0 & p_vec < p_cutoff, , drop = FALSE]
  compare_obj
}

run_go_enrichment <- function(genes, label, ont, universe_genes = NULL) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  universe_genes <- unique(as.character(universe_genes))
  universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
  if(length(genes) < 5) {
    cat("Too few genes for GO analysis:", label, ont, "\n")
    return(NULL)
  }
  tryCatch({
    enrichGO(
      gene = genes,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      universe = if(length(universe_genes) > 0) universe_genes else NULL,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )
  }, error = function(e) {
    cat("GO analysis failed for", label, ont, ":", e$message, "\n")
    NULL
  })
}

run_kegg_enrichment <- function(genes, label, universe_genes = NULL) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  universe_genes <- unique(as.character(universe_genes))
  universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]
  if(length(genes) < 5) {
    cat("Too few genes for KEGG analysis:", label, "\n")
    return(NULL)
  }
  tryCatch({
    enrichKEGG(
      gene = genes,
      organism = "mmu",
      universe = if(length(universe_genes) > 0) universe_genes else NULL,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    cat("KEGG analysis failed for", label, ":", e$message, "\n")
    NULL
  })
}

safe_save_enrichment_plot <- function(plot_obj, file, width = 10, height = 8,
                                      coordinate_system = c("axes", "none")) {
  coordinate_system <- match.arg(coordinate_system)
  tryCatch({
    plot_obj <- atac_standardize_plot(plot_obj, coordinate_system = coordinate_system)
    atac_save_object_dual(
      plot_obj = plot_obj,
      file = file,
      width = width,
      height = height,
      standardize = FALSE,
      bg = "white"
    )
    TRUE
  }, error = function(e) {
    report_status("WARN", paste("Skipping enrichment plot", basename(file), ":", e$message))
    remove_output_with_preview(file)
    FALSE
  })
}

format_enrichment_display_label <- function(x, width = 38) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  x <- trimws(gsub("\\s+", " ", x))
  stringr::str_trunc(x, width = width, side = "right", ellipsis = "...")
}

prepare_enrichment_plot_object <- function(enrich_obj, width = 38) {
  if(is.null(enrich_obj) || !methods::is(enrich_obj, "enrichResult")) return(enrich_obj)
  plot_obj <- enrich_obj
  plot_obj@result$Description <- make.unique(
    format_enrichment_display_label(plot_obj@result$Description, width = width),
    sep = " "
  )
  plot_obj
}

plot_enrichment_suite <- function(enrich_obj, prefix, title, plot_dir, show_n = 15) {
  enrich_obj <- filter_enrich_result(enrich_obj, p_cutoff = 0.05)
  enrich_df <- tryCatch(as.data.frame(enrich_obj), error = function(e) data.frame())
  if(is.null(enrich_obj) || nrow(enrich_df) == 0) return(invisible(NULL))

  plot_obj <- prepare_enrichment_plot_object(enrich_obj)
  show_n <- min(show_n, nrow(enrich_df))
  # pvalue is plotted directly: reverse the RNA-seq enrichment palette so
  # smaller (more significant) values are coral and larger values are blue.
  pvalue_colors <- rev(unname(viz_pal$enrich_p))
  safe_save_enrichment_plot(
    barplot(plot_obj, showCategory = show_n, color = "pvalue") +
      labs(title = NULL) +
      scale_fill_gradientn(colors = pvalue_colors) +
      theme_atac_classic(base_size = 8),
    file.path(plot_dir, paste0(prefix, "_barplot.pdf")),
    width = 10,
    height = 8
  )
  safe_save_enrichment_plot(
    enrichplot::dotplot(plot_obj, showCategory = show_n, color = "pvalue") +
      labs(title = NULL) +
      scale_color_gradientn(colors = pvalue_colors) +
      scale_size(range = c(1.5, 4)) +
      theme_atac_classic(base_size = 8),
    file.path(plot_dir, paste0(prefix, "_dotplot.pdf")),
    width = 10,
    height = 8
  )

  if(nrow(enrich_df) >= 3) {
    tree_obj <- tryCatch(enrichplot::pairwise_termsim(plot_obj), error = function(e) NULL)
    if(!is.null(tree_obj) && nrow(as.data.frame(tree_obj)) > 1) {
      safe_save_enrichment_plot(
        enrichplot::treeplot(
          tree_obj,
          color = "pvalue",
          showCategory = min(8L, show_n, nrow(as.data.frame(tree_obj))),
          nCluster = min(3L, 8L, show_n, nrow(as.data.frame(tree_obj))),
          label_format = 30,
          fontsize_tiplab = 2.4,
          fontsize_cladelab = 2.4,
          hilight = FALSE,
          align = "none",
          hexpand = 0.45,
          tiplab_offset = 0.1
        ) +
          scale_color_gradientn(colors = pvalue_colors) +
          labs(title = NULL, x = NULL, y = NULL) +
          theme_atac_void(base_size = 8),
        file.path(plot_dir, paste0(prefix, "_treeplot.pdf")),
        width = 12,
        height = 10,
        coordinate_system = "none"
      )
    }
  }

  if(nrow(enrich_df) >= 5) {
    cnet_obj <- tryCatch(
      enrichplot::cnetplot(plot_obj, showCategory = min(5, nrow(enrich_df)), node_label = "category"),
      error = function(e) NULL
    )
    if(!is.null(cnet_obj)) {
      for(i in seq_along(cnet_obj$layers)) {
        geom_class <- class(cnet_obj$layers[[i]]$geom)
        if(any(geom_class %in% c("GeomText", "GeomTextRepel", "GeomLabel", "GeomLabelRepel"))) {
          cnet_obj$layers[[i]]$aes_params$size <- 2.4
        }
      }
      safe_save_enrichment_plot(
        cnet_obj + labs(title = NULL, x = NULL, y = NULL) + theme_atac_void(base_size = 8),
        file.path(plot_dir, paste0(prefix, "_cnetplot.pdf")),
        width = 10,
        height = 8,
        coordinate_system = "none"
      )
    }
  }

  invisible(enrich_obj)
}

plot_compare_cluster_suite <- function(compare_obj, prefix, title, plot_dir, show_n = 15) {
  compare_obj <- filter_compare_cluster_result(compare_obj, p_cutoff = 0.05)
  compare_df <- tryCatch(as.data.frame(compare_obj), error = function(e) data.frame())
  if(is.null(compare_obj) || nrow(compare_df) == 0) return(invisible(NULL))
  ratio_to_numeric <- function(x) {
    vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
      nums <- suppressWarnings(as.numeric(parts))
      if(length(nums) >= 2 && !any(is.na(nums[1:2])) && nums[2] != 0) {
        nums[1] / nums[2]
      } else {
        suppressWarnings(as.numeric(parts[1]))
      }
    }, numeric(1))
  }

  compare_df <- compare_df %>%
    mutate(
      Description = as.character(Description),
      DescriptionShort = ifelse(nchar(Description) > 55, paste0(substr(Description, 1, 52), "..."), Description),
      GeneRatioNum = if("GeneRatio" %in% colnames(compare_df)) ratio_to_numeric(GeneRatio) else NA_real_,
      pvalue = suppressWarnings(as.numeric(pvalue)),
      p.adjust.num = if("p.adjust" %in% colnames(compare_df)) suppressWarnings(as.numeric(.data[["p.adjust"]])) else NA_real_,
      Count = suppressWarnings(as.numeric(Count))
    ) %>%
    filter(!is.na(DescriptionShort), DescriptionShort != "") %>%
    group_by(Cluster) %>%
    arrange(pvalue, p.adjust.num, desc(Count), .by_group = TRUE) %>%
    slice_head(n = show_n) %>%
    ungroup()

  if(nrow(compare_df) == 0) return(invisible(NULL))

  cluster_cols <- resolve_named_colors(unique(compare_df$Cluster), c(viz_pal$group, viz_pal$de), fallback_colors = viz_pal$feature_accent)
  # pvalue is plotted directly: reverse the RNA-seq enrichment palette so
  # smaller (more significant) values are coral and larger values are blue.
  pvalue_colors <- rev(unname(viz_pal$enrich_p))

  dotplot_obj <- ggplot(
    compare_df,
    aes(x = GeneRatioNum, y = reorder(DescriptionShort, GeneRatioNum), size = Count, color = pvalue)
  ) +
    geom_point(alpha = 0.85) +
    facet_wrap(~Cluster, scales = "free_y") +
    scale_color_gradientn(colors = pvalue_colors) +
    scale_size(range = c(1.5, 4.5)) +
    theme_atac_classic(base_size = 12) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      color = "pvalue"
    )
  safe_save_enrichment_plot(
    dotplot_obj,
    file.path(plot_dir, paste0(prefix, "_dotplot.pdf")),
    width = 11,
    height = 8
  )

  barplot_obj <- ggplot(
    compare_df,
    aes(x = GeneRatioNum, y = reorder(DescriptionShort, GeneRatioNum), fill = Cluster)
  ) +
    geom_col(width = 0.75) +
    facet_wrap(~Cluster, scales = "free_y") +
    scale_fill_manual(values = cluster_cols) +
    theme_atac_classic(base_size = 12) +
    labs(
      title = paste(title, "Barplot"),
      x = "Gene ratio",
      y = NULL,
      fill = "Cluster"
    )
  safe_save_enrichment_plot(
    barplot_obj,
    file.path(plot_dir, paste0(prefix, "_barplot.pdf")),
    width = 11,
    height = 8
  )

  invisible(compare_obj)
}
normalize_feature_label <- function(feature_vec) {
  feature_vec <- as.character(feature_vec)
  feature_vec <- sub(" \\(.*", "", feature_vec)
  dplyr::case_when(
    is.na(feature_vec) | feature_vec == "" ~ "Other",
    grepl("^Promoter|^Upstream$", feature_vec, ignore.case = TRUE) ~ "Promoter",
    grepl("UTR", feature_vec, ignore.case = TRUE) ~ "UTR",
    grepl("Exon", feature_vec, ignore.case = TRUE) ~ "Exon",
    grepl("Intron", feature_vec, ignore.case = TRUE) ~ "Intron",
    grepl("Downstream", feature_vec, ignore.case = TRUE) ~ "Downstream",
    grepl("Distal Intergenic|Intergenic|^Distal$", feature_vec, ignore.case = TRUE) ~ "Intergenic",
    TRUE ~ "Other"
  )
}

derive_feature_layer <- function(feature_vec) {
  feature_simple <- normalize_feature_label(feature_vec)
  dplyr::case_when(
    feature_simple == "Promoter" ~ "Promoter",
    feature_simple %in% c("UTR", "Exon", "Intron") ~ "Gene Body",
    feature_simple %in% c("Downstream", "Intergenic") ~ "Distal",
    TRUE ~ "Other"
  )
}

resolve_named_colors <- function(levels_vec, palette_map, fallback_colors = NULL, default_color = "#D4D9DE") {
  levels_vec <- unique(as.character(levels_vec))
  cols <- unname(palette_map[levels_vec])
  names(cols) <- levels_vec
  missing_idx <- is.na(cols)
  if(any(missing_idx)) {
    if(is.null(fallback_colors) || length(fallback_colors) == 0) {
      cols[missing_idx] <- default_color
    } else {
      cols[missing_idx] <- rep(fallback_colors, length.out = sum(missing_idx))
    }
  }
  cols
}

parse_alignment_rate <- function(log_file) {
  if(!file.exists(log_file)) return(NA_real_)
  hit <- grep("overall alignment rate", readLines(log_file, warn = FALSE), value = TRUE)
  if(length(hit) == 0) return(NA_real_)
  as.numeric(sub("%.*$", "", trimws(hit[1])))
}

extract_numeric_regex <- function(text, pattern, default = NA_real_) {
  if(length(text) == 0 || is.na(text[1])) return(default)
  hit <- regexec(pattern, text[1], perl = TRUE)
  match <- regmatches(text[1], hit)[[1]]
  if(length(match) < 2) return(default)
  suppressWarnings(as.numeric(match[2]))
}

parse_bowtie2_summary <- function(log_file) {
  out <- list(
    Bowtie2InputPairs = NA_real_,
    Bowtie2Concordant0Pairs = NA_real_,
    Bowtie2Concordant1Pairs = NA_real_,
    Bowtie2ConcordantMultiPairs = NA_real_,
    Bowtie2OverallAlignmentRatePct = NA_real_
  )
  if(!file.exists(log_file)) return(out)
  lines <- readLines(log_file, warn = FALSE)
  input_line <- grep("reads; of these:", lines, value = TRUE)
  concordant0_line <- grep("aligned concordantly 0 times", lines, value = TRUE)
  concordant1_line <- grep("aligned concordantly exactly 1 time", lines, value = TRUE)
  concordant_multi_line <- grep("aligned concordantly >1 times", lines, value = TRUE)
  rate_line <- grep("overall alignment rate", lines, value = TRUE)
  out$Bowtie2InputPairs <- extract_numeric_regex(input_line, "^\\s*([0-9]+)\\s+reads;")
  out$Bowtie2Concordant0Pairs <- extract_numeric_regex(concordant0_line, "^\\s*([0-9]+)\\s+\\(")
  out$Bowtie2Concordant1Pairs <- extract_numeric_regex(concordant1_line, "^\\s*([0-9]+)\\s+\\(")
  out$Bowtie2ConcordantMultiPairs <- extract_numeric_regex(concordant_multi_line, "^\\s*([0-9]+)\\s+\\(")
  out$Bowtie2OverallAlignmentRatePct <- if(length(rate_line) > 0) as.numeric(sub("%.*$", "", trimws(rate_line[1]))) else NA_real_
  out
}

parse_mito_stats <- function(stats_file) {
  out <- list(
    MappedQ30Reads = NA_real_,
    MitoReads = NA_real_,
    MitoPct = NA_real_,
    FilteredNonMitoProperPairReads = NA_real_,
    FinalDedupReads = NA_real_,
    DupRatePct = NA_real_
  )
  if(!file.exists(stats_file)) return(out)
  lines <- readLines(stats_file, warn = FALSE)
  total_line <- grep("Total=", lines, value = TRUE)
  filter_line <- grep("After filtering=", lines, value = TRUE)
  dedup_line <- grep("After dedup=", lines, value = TRUE)
  if(length(total_line) > 0) {
    out$MappedQ30Reads <- extract_numeric_regex(total_line, "Total=([0-9]+)")
    out$MitoReads <- extract_numeric_regex(total_line, "MT=([0-9]+)")
    out$MitoPct <- extract_numeric_regex(total_line, "\\(([0-9.]+)%\\)")
  }
  if(length(filter_line) > 0) {
    out$FilteredNonMitoProperPairReads <- extract_numeric_regex(filter_line, "After filtering=([0-9]+)")
  }
  if(length(dedup_line) > 0) {
    out$FinalDedupReads <- extract_numeric_regex(dedup_line, "After dedup=([0-9]+)")
    out$DupRatePct <- extract_numeric_regex(dedup_line, "DupRate=([0-9.]+)%")
  }
  out
}

parse_picard_dup_metrics <- function(metrics_file) {
  out <- list(
    PicardDupPct = NA_real_,
    EstimatedLibrarySize = NA_real_
  )
  if(!file.exists(metrics_file)) return(out)
  lines <- readLines(metrics_file, warn = FALSE)
  header_idx <- grep("^LIBRARY\tUNPAIRED_READS_EXAMINED", lines)
  if(length(header_idx) == 0 || header_idx[1] >= length(lines)) return(out)
  vals <- strsplit(lines[header_idx[1] + 1], "\t", fixed = TRUE)[[1]]
  if(length(vals) >= 10) {
    out$PicardDupPct <- as.numeric(vals[9]) * 100
    out$EstimatedLibrarySize <- as.numeric(vals[10])
  }
  out
}

parse_fastp_summary <- function(json_file) {
  out <- list(
    FastpRawReads = NA_real_,
    FastpRawBases = NA_real_,
    FastpCleanReads = NA_real_,
    FastpCleanBases = NA_real_,
    FastpPassedReads = NA_real_,
    FastpRawRead1MeanLength = NA_real_,
    FastpRawRead2MeanLength = NA_real_,
    FastpCleanRead1MeanLength = NA_real_,
    FastpCleanRead2MeanLength = NA_real_,
    FastpQ30RatePct = NA_real_,
    FastpGCContentPct = NA_real_,
    FastpSurvivalPct = NA_real_
  )
  if(!file.exists(json_file)) return(out)
  payload <- jsonlite::fromJSON(json_file)
  passed <- payload$filtering_result$passed_filter_reads
  before_reads <- payload$summary$before_filtering$total_reads
  after_reads <- payload$summary$after_filtering$total_reads
  out$FastpRawReads <- as.numeric(before_reads)
  out$FastpRawBases <- as.numeric(payload$summary$before_filtering$total_bases)
  out$FastpCleanReads <- as.numeric(after_reads)
  out$FastpCleanBases <- as.numeric(payload$summary$after_filtering$total_bases)
  out$FastpPassedReads <- as.numeric(passed)
  out$FastpRawRead1MeanLength <- as.numeric(payload$summary$before_filtering$read1_mean_length)
  out$FastpRawRead2MeanLength <- as.numeric(payload$summary$before_filtering$read2_mean_length)
  out$FastpCleanRead1MeanLength <- as.numeric(payload$summary$after_filtering$read1_mean_length)
  out$FastpCleanRead2MeanLength <- as.numeric(payload$summary$after_filtering$read2_mean_length)
  out$FastpQ30RatePct <- as.numeric(payload$summary$after_filtering$q30_rate) * 100
  out$FastpGCContentPct <- as.numeric(payload$summary$after_filtering$gc_content) * 100
  out$FastpSurvivalPct <- ifelse(before_reads > 0, 100 * passed / before_reads, NA_real_)
  out
}

primary_chr_ucsc <- paste0("chr", c(1:19, "X", "Y"))

granges_genome_keep_idx <- function(gr_obj, genome = BSgenome.Mmusculus.UCSC.mm39) {
  if(length(gr_obj) == 0) return(logical(0))
  genome_lengths <- tryCatch(GenomeInfoDb::seqlengths(genome), error = function(e) NULL)
  if(is.null(genome_lengths) || length(genome_lengths) == 0) {
    return(rep(TRUE, length(gr_obj)))
  }
  seq_chr <- as.character(seqnames(gr_obj))
  seq_len_vec <- genome_lengths[seq_chr]
  !is.na(seq_len_vec) &
    !is.na(start(gr_obj)) &
    !is.na(end(gr_obj)) &
    start(gr_obj) >= 1L &
    end(gr_obj) >= start(gr_obj) &
    end(gr_obj) <= seq_len_vec
}

drop_out_of_bound_granges <- function(gr_obj, genome = BSgenome.Mmusculus.UCSC.mm39) {
  if(length(gr_obj) == 0) return(gr_obj)
  keep_idx <- granges_genome_keep_idx(gr_obj, genome = genome)
  dropped_n <- sum(!keep_idx, na.rm = TRUE)
  gr_obj <- gr_obj[keep_idx]
  if(length(gr_obj) > 0) {
    seqlevels(gr_obj) <- seqlevelsInUse(gr_obj)
  }
  attr(gr_obj, "dropped_out_of_bounds") <- dropped_n
  gr_obj
}

keep_primary_chromosomes <- function(gr_obj, primary_chr = primary_chr_ucsc) {
  if(length(gr_obj) == 0) return(gr_obj)
  tryCatch({
    seqlevelsStyle(gr_obj) <- "UCSC"
  }, error = function(e) NULL)
  keep_idx <- as.character(seqnames(gr_obj)) %in% primary_chr
  gr_obj <- gr_obj[keep_idx]
  if(length(gr_obj) > 0) {
    gr_obj <- drop_out_of_bound_granges(gr_obj)
    seqlevels(gr_obj) <- seqlevelsInUse(gr_obj)
  }
  gr_obj
}

approx_tss_enrichment <- function(matrix_file, flank_bins = 100L, chunk_size = 4000L, max_rows = 10000L) {
  if(!file.exists(matrix_file)) return(NA_real_)
  con <- gzfile(matrix_file, open = "rt")
  on.exit(close(con), add = TRUE)
  sum_vec <- NULL
  n_rows <- 0L
  repeat {
    if(!is.null(max_rows) && !is.na(max_rows) && n_rows >= max_rows) break
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if(length(lines) == 0) break
    lines <- lines[!grepl("^[@#]", lines)]
    if(length(lines) == 0) next
    split_lines <- strsplit(lines, "\t", fixed = TRUE)
    for(parts in split_lines) {
      if(!is.null(max_rows) && !is.na(max_rows) && n_rows >= max_rows) break
      if(length(parts) <= 6) next
      signal <- suppressWarnings(as.numeric(parts[-(1:6)]))
      if(length(signal) == 0 || all(is.na(signal))) next
      if(is.null(sum_vec)) sum_vec <- numeric(length(signal))
      sum_vec <- sum_vec + signal
      n_rows <- n_rows + 1L
    }
  }
  if(is.null(sum_vec) || n_rows == 0) return(NA_real_)
  mean_signal <- sum_vec / n_rows
  center_idx <- floor(length(mean_signal) / 2) + 1L
  flank <- c(
    mean_signal[seq_len(min(flank_bins, length(mean_signal)))],
    tail(mean_signal, min(flank_bins, length(mean_signal)))
  )
  denom <- mean(flank, na.rm = TRUE)
  if(is.na(denom) || denom <= 0) return(NA_real_)
  mean_signal[center_idx] / denom
}

normalize_tss_profile <- function(signal, positions = NULL, flank_abs_min = 1500) {
  signal <- suppressWarnings(as.numeric(signal))
  if(length(signal) == 0 || all(is.na(signal))) {
    return(data.frame(
      signal = signal,
      baseline = NA_real_,
      signal_fold_over_flank = NA_real_,
      signal_centered = NA_real_,
      signal_zscore = NA_real_
    ))
  }

  if(is.null(positions) || length(positions) != length(signal)) {
    flank_idx <- c(
      seq_len(max(1L, floor(length(signal) * 0.1))),
      seq.int(max(1L, ceiling(length(signal) * 0.9)), length(signal))
    )
  } else {
    flank_idx <- which(abs(positions) >= flank_abs_min)
  }
  if(length(flank_idx) == 0) flank_idx <- seq_along(signal)

  flank_signal <- signal[flank_idx]
  baseline <- mean(flank_signal, na.rm = TRUE)
  baseline_sd <- sd(flank_signal, na.rm = TRUE)
  if(is.na(baseline) || baseline <= 0) baseline <- mean(signal, na.rm = TRUE)
  if(is.na(baseline) || baseline <= 0) baseline <- NA_real_
  if(is.na(baseline_sd) || baseline_sd <= 0) baseline_sd <- NA_real_

  data.frame(
    signal = signal,
    baseline = baseline,
    signal_fold_over_flank = if(is.na(baseline)) NA_real_ else signal / baseline,
    signal_centered = if(is.na(baseline)) NA_real_ else signal - baseline,
    signal_zscore = if(is.na(baseline) || is.na(baseline_sd)) NA_real_ else (signal - baseline) / baseline_sd
  )
}

qc_status <- function(value, cutoff, direction = c("ge", "le")) {
  direction <- match.arg(direction)
  status <- rep("NA", length(value))
  keep <- !is.na(value)
  if(direction == "ge") {
    status[keep] <- ifelse(value[keep] >= cutoff, "PASS", "WARN")
  } else {
    status[keep] <- ifelse(value[keep] <= cutoff, "PASS", "WARN")
  }
  status
}

score_by_abs_effect <- function(effect, full_effect = 1) {
  effect <- abs(suppressWarnings(as.numeric(effect)))
  pmin(1, replace_na(effect, 0) / full_effect)
}

score_by_pvalue <- function(pvalue, strong = 0.05, nominal = 0.1, missing_weight = 0.15) {
  pvalue <- suppressWarnings(as.numeric(pvalue))
  dplyr::case_when(
    is.na(pvalue) ~ missing_weight,
    pvalue <= strong ~ 1,
    pvalue <= nominal ~ 0.5,
    TRUE ~ 0.15
  )
}

build_sample_qc_summary <- function(sample_names) {
  frip_tbl <- if(file.exists(file.path(work_dir, "04counts", "FRiP_scores.txt"))) {
    read.delim(file.path(work_dir, "04counts", "FRiP_scores.txt"), sep = "\t", stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  bind_rows(lapply(sample_names, function(sample) {
    mito <- parse_mito_stats(file.path(work_dir, "02mapping", "qc", paste0(sample, ".mito_stats.txt")))
    picard <- parse_picard_dup_metrics(file.path(work_dir, "02mapping", "qc", paste0(sample, ".picard_dup_metrics.txt")))
    fastp <- parse_fastp_summary(file.path(work_dir, "01rawdata", "qc_results", paste0(sample, "_fastp.json")))
    bowtie2 <- parse_bowtie2_summary(file.path(work_dir, "02mapping", "qc", paste0(sample, ".bowtie2.log")))
    frip_row <- frip_tbl %>% filter(Sample == sample)
    sample_meta <- sample_manifest[sample_manifest$sample == sample, , drop = FALSE]
    if(nrow(sample_meta) > 1) sample_meta <- sample_meta[1, , drop = FALSE]
    raw_r1 <- if(nrow(sample_meta) > 0) sample_meta$raw_r1[1] else file.path("01rawdata", paste0(sample, "_1.fastq.gz"))
    raw_r2 <- if(nrow(sample_meta) > 0) sample_meta$raw_r2[1] else file.path("01rawdata", paste0(sample, "_2.fastq.gz"))
    clean_r1 <- file.path("01rawdata", "cleaned_data", paste0(sample, "_1.fastq.gz"))
    clean_r2 <- file.path("01rawdata", "cleaned_data", paste0(sample, "_2.fastq.gz"))
    reads_in_peaks <- ifelse(nrow(frip_row) > 0 && "Reads_in_Peaks" %in% colnames(frip_row), as.numeric(frip_row$Reads_in_Peaks[1]), NA_real_)
    frip_total_reads <- ifelse(nrow(frip_row) > 0 && "Total_Reads" %in% colnames(frip_row), as.numeric(frip_row$Total_Reads[1]), NA_real_)
    data.frame(
      Sample = sample,
      Group = sample_group(sample),
      Stage = if(nrow(sample_meta) > 0) sample_meta$stage[1] else Sys.getenv("ATAC_TISSUE"),
      RawFastq1 = raw_r1,
      RawFastq2 = raw_r2,
      CleanFastq1 = clean_r1,
      CleanFastq2 = clean_r2,
      FastpRawReads = fastp$FastpRawReads,
      FastpRawBases = fastp$FastpRawBases,
      FastpCleanReads = fastp$FastpCleanReads,
      FastpCleanBases = fastp$FastpCleanBases,
      FastpPassedReads = fastp$FastpPassedReads,
      FastpSurvivalPct = fastp$FastpSurvivalPct,
      FastpRawRead1MeanLength = fastp$FastpRawRead1MeanLength,
      FastpRawRead2MeanLength = fastp$FastpRawRead2MeanLength,
      FastpCleanRead1MeanLength = fastp$FastpCleanRead1MeanLength,
      FastpCleanRead2MeanLength = fastp$FastpCleanRead2MeanLength,
      Bowtie2InputPairs = bowtie2$Bowtie2InputPairs,
      Bowtie2Concordant0Pairs = bowtie2$Bowtie2Concordant0Pairs,
      Bowtie2Concordant1Pairs = bowtie2$Bowtie2Concordant1Pairs,
      Bowtie2ConcordantMultiPairs = bowtie2$Bowtie2ConcordantMultiPairs,
      AlignmentRatePct = bowtie2$Bowtie2OverallAlignmentRatePct,
      MappedQ30Reads = mito$MappedQ30Reads,
      MitoReads = mito$MitoReads,
      MitoPct = mito$MitoPct,
      FilteredNonMitoProperPairReads = mito$FilteredNonMitoProperPairReads,
      FinalDedupReads = mito$FinalDedupReads,
      DupRatePct = mito$DupRatePct,
      PicardDupPct = picard$PicardDupPct,
      EstimatedLibrarySize = picard$EstimatedLibrarySize,
      FastpQ30RatePct = fastp$FastpQ30RatePct,
      FastpGCContentPct = fastp$FastpGCContentPct,
      FRiPTotalReads = frip_total_reads,
      ReadsInPeaks = reads_in_peaks,
      FRiP = ifelse(nrow(frip_row) > 0, as.numeric(frip_row$FRiP[1]), NA_real_),
      CleanPeakCount = ifelse(file.exists(file.path(work_dir, "03peaks", paste0(sample, "_peaks.clean.narrowPeak"))),
                              length(readLines(file.path(work_dir, "03peaks", paste0(sample, "_peaks.clean.narrowPeak")))),
                              NA_real_),
      ApproxTSSEnrichment = approx_tss_enrichment(file.path(work_dir, "06qc_report", paste0(sample, "_TSS_matrix.gz"))),
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      QC_PassCount =
        (qc_status(AlignmentRatePct, 90, "ge") == "PASS") +
        (qc_status(MitoPct, 20, "le") == "PASS") +
        (qc_status(DupRatePct, 50, "le") == "PASS") +
        (qc_status(FRiP, 0.2, "ge") == "PASS") +
        (qc_status(ApproxTSSEnrichment, 5, "ge") == "PASS") +
        (qc_status(CleanPeakCount, 30000, "ge") == "PASS")
    )
}

sample_group <- function(sample_name) {
  mapped <- unname(sample_group_map[as.character(sample_name)])
  if(anyNA(mapped) || any(!nzchar(mapped))) stop("样本未在样本表声明，禁止按名称推断分组")
  mapped
}

sample_display_order <- function(sample_names) {
  sample_names <- as.character(sample_names)
  groups <- sample_group(sample_names)
  known_sample_order <- match(sample_names, samples)
  known_sample_order[is.na(known_sample_order)] <- seq_along(sample_names)[is.na(known_sample_order)] + length(samples)
  order(match(groups, group_levels), known_sample_order, sample_names, na.last = TRUE)
}

order_sample_columns <- function(sample_names) {
  sample_names <- as.character(sample_names)
  sample_names[sample_display_order(sample_names)]
}

canonical_peak_id <- function(seqnames, start, end) {
  seqnames <- sub("^chr", "", as.character(seqnames), ignore.case = TRUE)
  paste(seqnames, as.integer(start), as.integer(end), sep = ":")
}

# 安全地将 seqnames 转为 UCSC 风格 (chr 前缀)
# 注意: ifelse() 在 factor 输入下会返回整数索引而非字符值，
# 因此必须先用 as.character() 转换
ensure_ucsc_seqnames <- function(seqnames) {
  seqnames <- as.character(seqnames)
  ifelse(grepl("^chr", seqnames), seqnames, paste0("chr", seqnames))
}

comparison_open_label <- paste0(treatment_group, "-open")
reference_open_label <- paste0(baseline_group, "-open")

direction_from_fold <- function(fold_vec, pos_label = comparison_open_label, neg_label = reference_open_label, zero_label = "Neutral") {
  dplyr::case_when(
    is.na(fold_vec) ~ zero_label,
    fold_vec > 0 ~ pos_label,
    fold_vec < 0 ~ neg_label,
    TRUE ~ zero_label
  )
}

direction_factor <- function(fold_vec, pos_label = comparison_open_label, neg_label = reference_open_label, zero_label = "Neutral") {
  factor(
    direction_from_fold(fold_vec, pos_label = pos_label, neg_label = neg_label, zero_label = zero_label),
    levels = c(neg_label, pos_label, zero_label)
  )
}

format_threshold <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))[1]
  if(is.na(x)) return("NA")
  format(signif(x, digits), scientific = FALSE, trim = TRUE)
}

primary_diff_peak_label <- function() {
  metric <- if(atac_significance_column == "pvalue") "Exploratory nominal p-value" else "FDR/padj"
  cutoff <- if(atac_significance_column == "pvalue") atac_primary_p_cutoff else atac_diff_fdr_cutoff
  paste0(
    metric, if(atac_significance_column == "pvalue") " < " else " <= ", format_threshold(cutoff),
    if(atac_primary_require_fc) paste0(" and |log2FC| >= ", format_threshold(atac_diff_log2fc_cutoff)) else ""
  )
}

official_diff_peak_label <- function() {
  primary_diff_peak_label()
}

primary_pvalue_diff_peak_label <- function() {
  paste0("Exploratory nominal p-value < ", format_threshold(atac_primary_p_cutoff))
}

primary_diff_class_label <- function() {
  paste0("Primary: ", primary_diff_peak_label())
}

nominal_diff_peak_label <- function() {
  paste0(
    "Exploratory nominal p-value < ", format_threshold(atac_nominal_p_cutoff)
  )
}

standardize_diff_peak_df <- function(df) {
  if(is.null(df)) df <- data.frame(stringsAsFactors = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  colnames(df) <- make.unique(colnames(df), sep = "_dup")
  if(!"seqnames" %in% colnames(df) && "Chr" %in% colnames(df)) df$seqnames <- df$Chr
  if(!"seqnames" %in% colnames(df) && "CHR" %in% colnames(df)) df$seqnames <- df$CHR
  if(!"start" %in% colnames(df) && "Start" %in% colnames(df)) df$start <- df$Start
  if(!"start" %in% colnames(df) && "START" %in% colnames(df)) df$start <- df$START
  if(!"end" %in% colnames(df) && "End" %in% colnames(df)) df$end <- df$End
  if(!"end" %in% colnames(df) && "END" %in% colnames(df)) df$end <- df$END
  if(!"Fold" %in% colnames(df) && "log2Fold" %in% colnames(df)) df$Fold <- df$log2Fold
  if(!"pvalue" %in% colnames(df) && "p-value" %in% colnames(df)) df$pvalue <- df$`p-value`
  if(!"pvalue" %in% colnames(df) && "p.value" %in% colnames(df)) df$pvalue <- df$`p.value`
  if(!"FDR" %in% colnames(df) && "fdr" %in% colnames(df)) df$FDR <- df$fdr

  df <- ensure_columns(df, list(
    seqnames = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    Fold = NA_real_,
    pvalue = NA_real_,
    FDR = NA_real_
  ))

  seq_chr <- as.character(df$seqnames)
  seq_chr[is.na(seq_chr) | !nzchar(seq_chr)] <- NA_character_
  df$seqnames <- ifelse(is.na(seq_chr), NA_character_, ensure_ucsc_seqnames(seq_chr))
  df$start <- suppressWarnings(as.integer(df$start))
  df$end <- suppressWarnings(as.integer(df$end))
  df$Fold <- suppressWarnings(as.numeric(df$Fold))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$FDR <- suppressWarnings(as.numeric(df$FDR))
  # Drop the retired independent tier when reading legacy cache/table files.
  df$HighConfidenceDiffPeak <- NULL

  valid_coord <- !is.na(df$seqnames) & !is.na(df$start) & !is.na(df$end)
  df$peak_id <- NA_character_
  df$peak_id[valid_coord] <- paste(df$seqnames[valid_coord], df$start[valid_coord], df$end[valid_coord], sep = ":")
  df$peak_id_join <- NA_character_
  df$peak_id_join[valid_coord] <- canonical_peak_id(df$seqnames[valid_coord], df$start[valid_coord], df$end[valid_coord])

  flag_tbl <- compute_atac_diff_peak_flags(
    pvalue = df$pvalue,
    fdr = df$FDR,
    fold = df$Fold,
    significance_mode = atac_significance_mode,
    primary_p_cutoff = atac_primary_p_cutoff,
    fdr_cutoff = atac_diff_fdr_cutoff,
    log2fc_cutoff = atac_diff_log2fc_cutoff,
    primary_require_fc = atac_primary_require_fc
  )
  for(flag_name in colnames(flag_tbl)) df[[flag_name]] <- flag_tbl[[flag_name]]
  df$PrimaryPvalueDiffPeak <- df$PvalueDiffPeak
  df$NominalDiffPeak <- is.finite(df$pvalue) & df$pvalue >= 0 & df$pvalue < atac_nominal_p_cutoff
  primary_class <- primary_diff_class_label()
  df$DiffPeakClass <- dplyr::case_when(
    df$PrimaryDiffPeak ~ primary_class,
    TRUE ~ "Not significant"
  )
  df$DiffPeakClass <- factor(df$DiffPeakClass, levels = c("Not significant", primary_class))
  df$DiffDirection <- direction_factor(df$Fold)
  df
}

filter_diff_peaks <- function(df, mode = c("primary", "nominal")) {
  mode <- match.arg(mode)
  df <- standardize_diff_peak_df(df)
  if(mode == "primary") return(df %>% filter(PrimaryDiffPeak))
  df %>% filter(NominalDiffPeak)
}

diff_peak_granges <- function(df) {
  df <- standardize_diff_peak_df(df) %>%
    filter(!is.na(seqnames), !is.na(start), !is.na(end), start >= 1, end >= start)
  if(nrow(df) == 0) return(GRanges())
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end)
  )
  meta_cols <- setdiff(colnames(df), c("seqnames", "start", "end"))
  if(length(meta_cols) > 0) {
    mcols(gr) <- S4Vectors::DataFrame(df[, meta_cols, drop = FALSE])
  }
  gr
}

build_footprint_supported_tables <- function() {
  motif_candidates <- c(step_table_file("Step10", "Motif_Binding_Sites.csv"))
  anno_candidates <- c(
    step_table_file("Step7", "Annotated_DiffPeaks.csv")
  )
  motif_file <- motif_candidates[file.exists(motif_candidates)][1]
  anno_file <- anno_candidates[file.exists(anno_candidates)][1]
  if(!file.exists(motif_file) || !file.exists(anno_file)) return(NULL)

  motif_sites <- read.csv(motif_file, stringsAsFactors = FALSE)
  peak_gene_raw <- read.csv(anno_file, stringsAsFactors = FALSE)
  if("SYMBOL" %in% colnames(peak_gene_raw) && !"geneSymbol" %in% colnames(peak_gene_raw)) {
    colnames(peak_gene_raw)[colnames(peak_gene_raw) == "SYMBOL"] <- "geneSymbol"
  }
  if(!"geneSymbol" %in% colnames(peak_gene_raw)) {
    peak_gene_raw$geneSymbol <- NA_character_
  }
  if(!"peak_id" %in% colnames(peak_gene_raw) && all(c("seqnames", "start", "end") %in% colnames(peak_gene_raw))) {
    peak_gene_raw$peak_id <- paste(peak_gene_raw$seqnames, peak_gene_raw$start, peak_gene_raw$end, sep = ":")
  }
  if(!all(c("peak_id", "geneSymbol") %in% colnames(peak_gene_raw))) return(NULL)

  peak_gene <- peak_gene_raw %>%
    transmute(
      peak_id = peak_id,
      geneSymbol = geneSymbol,
      annotation = if("annotation" %in% colnames(peak_gene_raw)) annotation else NA_character_
    ) %>%
    filter(!is.na(geneSymbol) & geneSymbol != "")
  if(nrow(peak_gene) == 0) return(NULL)

  edges <- motif_sites %>%
    distinct(motif_name, peak_id) %>%
    inner_join(peak_gene, by = "peak_id") %>%
    dplyr::count(motif_name, geneSymbol, name = "SupportingPeakCount", sort = TRUE)
  if(nrow(edges) == 0) return(NULL)

  motif_summary <- edges %>%
    group_by(motif_name) %>%
    summarise(
      TargetGeneCount = n_distinct(geneSymbol),
      SupportingPeakCount = sum(SupportingPeakCount),
      .groups = "drop"
    ) %>%
    arrange(desc(TargetGeneCount), desc(SupportingPeakCount))

  gene_summary <- edges %>%
    group_by(geneSymbol) %>%
    summarise(
      MotifCount = n_distinct(motif_name),
      Motifs = paste(sort(unique(motif_name)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(MotifCount), geneSymbol)

  write.csv(edges, step_table_file("Step10", "Footprint_Supported_Motif_Target_Edges.csv"), row.names = FALSE)
  write.csv(motif_summary, step_table_file("Step10", "Footprint_Supported_Motif_Summary.csv"), row.names = FALSE)
  write.csv(gene_summary, step_table_file("Step10", "Footprint_Supported_Gene_Summary.csv"), row.names = FALSE)

  list(edges = edges, motif_summary = motif_summary, gene_summary = gene_summary)
}

detect_first_column <- function(col_names, pattern) {
  hits <- grep(pattern, col_names, value = TRUE)
  if(length(hits) == 0) NA_character_ else hits[1]
}

infer_bindetect_labels <- function(change_col) {
  if(identical(change_col, "comparison_reference_change")) return(c(positive="Comparison_gain", negative="Reference_gain"))
  if(identical(change_col, "reference_comparison_change")) return(c(positive="Reference_gain", negative="Comparison_gain"))
  stop("BINDetect change 列方向未经确认: ", change_col)
}

clean_motif_name <- function(x) {
  x <- ifelse(is.na(x), NA_character_, x)
  x <- gsub("^_+", "", x)
  x <- sub("_.*$", "", x)
  x
}

normalize_tf_symbol <- function(x) {
  x <- ifelse(is.na(x), NA_character_, x)
  x <- sub("::.*$", "", x)
  x <- gsub("[^A-Za-z0-9]", "", x)
  toupper(x)
}

validate_bindetect_provenance <- function(root_dir = file.path(work_dir, "07footprinting")) {
  result <- function(status, reason) data.frame(Status = status, Reason = reason, stringsAsFactors = FALSE)
  provenance_file <- file.path(root_dir, "BINDetect_Provenance.tsv")
  status_file <- file.path(root_dir, "BINDetect_Status.tsv")
  if(!file.exists(provenance_file) || !file.exists(status_file)) {
    return(result("UNVERIFIED", "Missing input/output SHA256 provenance or execution status; historical BINDetect tables are not current evidence"))
  }
  provenance <- tryCatch(read.delim(provenance_file, stringsAsFactors = FALSE), error = function(e) data.frame())
  execution <- tryCatch(read.delim(status_file, stringsAsFactors = FALSE), error = function(e) data.frame())
  if(!all(c("Status", "Reason") %in% names(execution)) || nrow(execution) != 1L) return(result("UNVERIFIED", "Malformed BINDetect execution status"))
  if(!identical(execution$Status[1], "COMPLETED")) return(result(execution$Status[1], execution$Reason[1]))
  roles <- c("signal_comparison", "signal_reference", "genome", "motifs", "peaks", "results")
  if(!all(c("Role", "Path", "SHA256") %in% names(provenance)) ||
     nrow(provenance) != length(roles) || !setequal(provenance$Role, roles)) return(result("UNVERIFIED", "Incomplete or duplicated BINDetect provenance roles"))
  provenance <- provenance[match(roles, provenance$Role), , drop = FALSE]
  # Expected paths must still be the paths declared by the footprint manifests.
  declared_inputs <- read_tsv_if_exists(file.path(root_dir, "Footprinting_Input_Manifest.tsv"))
  declared_tracks <- read_tsv_if_exists(file.path(root_dir, "Footprinting_Track_Manifest.tsv"))
  if(!all(c("type", "path") %in% names(declared_inputs)) ||
     !all(c("group", "track_type", "path") %in% names(declared_tracks))) return(result("UNVERIFIED", "Footprint input/track declarations are missing or malformed"))
  expected_paths <- c(
    declared_tracks$path[declared_tracks$group == comparison_group & declared_tracks$track_type == "footprints_bw"],
    declared_tracks$path[declared_tracks$group == reference_group & declared_tracks$track_type == "footprints_bw"],
    declared_inputs$path[declared_inputs$type == "genome_fasta"],
    declared_inputs$path[declared_inputs$type == "motif_db"],
    declared_inputs$path[declared_inputs$type == "consensus_peaks"],
    file.path(root_dir, "bindetect_results.txt")
  )
  if(length(expected_paths) != length(roles) || anyNA(expected_paths) || anyNA(provenance$Path) ||
     any(normalizePath(provenance$Path, mustWork = FALSE) != normalizePath(expected_paths, mustWork = FALSE))) return(result("UNVERIFIED", "BINDetect input/output paths disagree with current declarations"))
  if(any(!file.exists(provenance$Path))) return(result("UNVERIFIED", "A recorded BINDetect input or output is missing"))
  if(!requireNamespace("digest", quietly = TRUE)) return(result("UNVERIFIED", "digest package unavailable; SHA256 verification cannot be performed"))
  observed <- vapply(provenance$Path, function(path) digest::digest(file = path, algo = "sha256"), character(1))
  if(any(is.na(provenance$SHA256)) || any(observed != provenance$SHA256)) return(result("STALE", "Current input/output SHA256 differs from the completed BINDetect run"))
  result("VERIFIED", "Completed BINDetect inputs and results match all recorded SHA256 values; upstream FASTQ provenance is assessed separately")
}

read_bindetect_results <- function(root_dir = file.path(work_dir, "07footprinting")) {
  provenance_status <- validate_bindetect_provenance(root_dir)
  if(!identical(provenance_status$Status[1], "VERIFIED")) {
    report_status("WARN", paste("BINDetect evidence rejected:", provenance_status$Reason[1]))
    return(NULL)
  }
  candidates <- c(
    file.path(root_dir, "bindetect_results.txt"),
    file.path(root_dir, "BINDetect", "bindetect_results.txt"),
    file.path(root_dir, "bindetect", "bindetect_results.txt")
  )
  bindetect_file <- candidates[file.exists(candidates)][1]
  if(length(bindetect_file) == 0 || is.na(bindetect_file)) return(NULL)

  bindetect_df <- read.delim(bindetect_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if(nrow(bindetect_df) == 0) return(NULL)

  change_col <- detect_first_column(colnames(bindetect_df), "_change$")
  pvalue_col <- detect_first_column(colnames(bindetect_df), "_pvalue$")
  highlight_col <- detect_first_column(colnames(bindetect_df), "_highlighted$")

  if(is.na(change_col) || is.na(pvalue_col)) return(NULL)

  motif_id_source <- if("motif_id" %in% colnames(bindetect_df)) {
    bindetect_df$motif_id
  } else if("output_prefix" %in% colnames(bindetect_df)) {
    parsed_id <- sub("^.*_(MA[0-9]+\\.[0-9]+)$", "\\1", bindetect_df$output_prefix)
    ifelse(grepl("^MA[0-9]+\\.[0-9]+$", parsed_id), parsed_id, NA_character_)
  } else {
    rep(NA_character_, nrow(bindetect_df))
  }

  motif_name_source <- if("name" %in% colnames(bindetect_df)) {
    bindetect_df$name
  } else {
    motif_id_source
  }

  direction_labels <- infer_bindetect_labels(change_col)
  bindetect_df$motif_name <- vapply(motif_name_source, clean_motif_name, character(1))
  if(!"motif_id" %in% colnames(bindetect_df)) {
    bindetect_df$motif_id <- motif_id_source
  }
  bindetect_df$change_value <- suppressWarnings(as.numeric(bindetect_df[[change_col]]))
  bindetect_df$pvalue_value <- suppressWarnings(as.numeric(bindetect_df[[pvalue_col]]))
  bindetect_df$highlighted_value <- is.finite(bindetect_df$pvalue_value) &
    bindetect_df$pvalue_value >= 0 & bindetect_df$pvalue_value < 0.05
  bindetect_df$direction <- ifelse(
    !bindetect_df$highlighted_value | is.na(bindetect_df$highlighted_value) | is.na(bindetect_df$change_value),
    "NS",
    ifelse(bindetect_df$change_value >= 0, direction_labels["positive"], direction_labels["negative"])
  )
  bindetect_df$neg_log10_p <- -log10(pmax(bindetect_df$pvalue_value, 1e-300))
  bindetect_df$total_tfbs <- if("total_tfbs" %in% colnames(bindetect_df)) {
    suppressWarnings(as.numeric(bindetect_df$total_tfbs))
  } else {
    NA_real_
  }

  attr(bindetect_df, "bindetect_file") <- bindetect_file
  attr(bindetect_df, "bindetect_change_col") <- change_col
  bindetect_df
}

build_footprint_regulator_evidence <- function(bindetect_df, motif_df = NULL, de_gene_df = NULL, site_summary_df = NULL) {
  evidence_df <- bindetect_df %>%
    mutate(
      TF_SYMBOL_KEY = normalize_tf_symbol(motif_name)
    ) %>%
    dplyr::select(motif_name, motif_id, total_tfbs, change_value, pvalue_value,
                  highlighted_value, direction, neg_log10_p, TF_SYMBOL_KEY) %>%
    distinct()

  if(!is.null(motif_df) && nrow(motif_df) > 0) {
    evidence_df <- evidence_df %>%
      left_join(
        motif_df %>% dplyr::select(motif_name, count, freq),
        by = "motif_name"
      )
  }

  if(!is.null(site_summary_df) && nrow(site_summary_df) > 0) {
    evidence_df <- evidence_df %>%
      left_join(site_summary_df, by = "motif_name")
  }

  if(!"count" %in% colnames(evidence_df)) evidence_df$count <- NA_real_
  if(!"freq" %in% colnames(evidence_df)) evidence_df$freq <- NA_real_
  if(!"SupportingSiteCount" %in% colnames(evidence_df)) evidence_df$SupportingSiteCount <- NA_real_

  if(!is.null(de_gene_df) && nrow(de_gene_df) > 0) {
    tf_rna_df <- de_gene_df %>%
      mutate(TF_SYMBOL_KEY = toupper(GeneSymbol)) %>%
      dplyr::select(TF_SYMBOL_KEY, TF_RNA_log2FC = log2FoldChange, TF_RNA_pvalue = pvalue) %>%
      distinct()
    evidence_df <- evidence_df %>% left_join(tf_rna_df, by = "TF_SYMBOL_KEY")
  } else {
    evidence_df$TF_RNA_log2FC <- NA_real_
    evidence_df$TF_RNA_pvalue <- NA_real_
  }

  max_freq <- if("freq" %in% colnames(evidence_df) && any(!is.na(evidence_df$freq))) {
    max(evidence_df$freq, na.rm = TRUE)
  } else {
    NA_real_
  }

  evidence_df <- evidence_df %>%
    mutate(
      BindetectEvidence = ifelse(is.na(highlighted_value), 0, ifelse(highlighted_value, 2, 0)) +
        pmin(3, replace_na(neg_log10_p, 0) / 2) +
        pmin(3, abs(replace_na(change_value, 0))),
      MotifEvidence = if(!is.na(max_freq) && max_freq > 0) 2 * replace_na(freq, 0) / max_freq else 0,
      SiteEvidence = if("SupportingSiteCount" %in% colnames(evidence_df)) {
        pmin(2, log10(replace_na(SupportingSiteCount, 0) + 1))
      } else {
        0
      },
      RNAEvidence = ifelse(is.finite(TF_RNA_pvalue) & TF_RNA_pvalue >= 0 & TF_RNA_pvalue < 0.05,
                           1 + pmin(2, abs(TF_RNA_log2FC) / 2),
                           0),
      IntegratedEvidenceScore = round(BindetectEvidence + MotifEvidence + SiteEvidence + RNAEvidence, 2),
      EvidenceTier = case_when(
        IntegratedEvidenceScore >= 7 ~ "Tier1",
        IntegratedEvidenceScore >= 4.5 ~ "Tier2",
        TRUE ~ "Tier3"
      )
    ) %>%
    arrange(desc(IntegratedEvidenceScore), pvalue_value, desc(total_tfbs))

  evidence_df
}

# ------------------------------------------------------------------
# Modular state cache
# ------------------------------------------------------------------

atac_state_file <- function(file_name) {
  file.path(cache_root_dir, file_name)
}

cache_04_peak_context_path <- atac_state_file("04_peak_context.rds")
cache_05_differential_accessibility_state_path <- atac_state_file("05_differential_accessibility_state.rds")
cache_07_peak_gene_enrichment_state_path <- atac_state_file("07_peak_gene_enrichment_state.rds")
cache_08_motif_state_path <- atac_state_file("08_motif_state.rds")
cache_09_multiomics_state_path <- atac_state_file("09_multiomics_state.rds")
cache_10_footprint_state_path <- atac_state_file("10_footprint_state.rds")
cache_11_wgbs_state_path <- atac_state_file("11_wgbs_state.rds")

save_state <- function(file_name, object_names, env = parent.frame()) {
  object_names <- unique(as.character(object_names))
  object_names <- object_names[nzchar(object_names)]
  payload <- list()
  for(object_name in object_names) {
    if(exists(object_name, envir = env, inherits = TRUE)) {
      payload[[object_name]] <- get(object_name, envir = env, inherits = TRUE)
    }
  }
  if(length(payload) == 0) {
    report_status("WARN", paste("No objects available for state cache:", file_name))
    return(invisible(NULL))
  }
  payload$.__analysis_signature__ <- current_analysis_signature
  payload$.__analysis_signature_config__ <- analysis_signature_config
  dir.create(cache_root_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(payload, atac_state_file(file_name))
  report_status("INFO", paste("Saved state cache:", atac_state_file(file_name)))
  invisible(payload)
}

load_state <- function(file_name, env = parent.frame(), required = FALSE) {
  path <- atac_state_file(file_name)
  if(!file.exists(path)) {
    msg <- paste("Missing state cache:", path)
    if(isTRUE(required)) stop(msg)
    report_status("WARN", msg)
    return(invisible(NULL))
  }
  payload <- readRDS(path)
  if(!is.list(payload)) {
    stop("State cache is not a named list: ", path)
  }
  cached_signature <- payload$.__analysis_signature__ %||% ""
  if(!identical(cached_signature, current_analysis_signature)) {
    msg <- paste(
      "STALE state cache (analysis signature mismatch):", path,
      "cached=", ifelse(nzchar(cached_signature), cached_signature, "missing"),
      "current=", current_analysis_signature
    )
    if(isTRUE(required)) stop(msg)
    report_status("WARN", msg)
    return(invisible(NULL))
  }
  payload$.__analysis_signature__ <- NULL
  payload$.__analysis_signature_config__ <- NULL
  for(object_name in names(payload)) {
    assign(object_name, payload[[object_name]], envir = env)
  }
  report_status("INFO", paste("Loaded state cache:", path))
  invisible(payload)
}

load_atac_context <- function(module_id, env = parent.frame()) {
  module_id <- as.integer(module_id)
  if(is.na(module_id) || module_id <= 1) {
    return(invisible(NULL))
  }

  if(module_id %in% c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)) {
    load_state("04_peak_context.rds", env = env, required = module_id %in% c(5, 6, 7, 8))
  }
  if(module_id %in% c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)) {
    load_state("05_differential_accessibility_state.rds", env = env, required = TRUE)
  }
  if(module_id %in% c(9, 10, 11, 12, 15, 16)) {
    load_state("07_peak_gene_enrichment_state.rds", env = env, required = module_id %in% c(9, 11, 15))
  }
  if(module_id %in% c(10, 12, 15, 16)) {
    load_state("08_motif_state.rds", env = env, required = module_id %in% c(10, 12))
  }
  if(module_id %in% c(10, 11, 12, 15, 16)) {
    load_state("09_multiomics_state.rds", env = env, required = env_bool("ATAC_ENABLE_MULTIOMICS", FALSE) && module_id %in% c(11, 15))
  }
  if(module_id %in% c(11, 12, 15, 16)) {
    load_state("10_footprint_state.rds", env = env, required = FALSE)
  }
  if(module_id %in% c(15, 16)) {
    load_state("11_wgbs_state.rds", env = env, required = FALSE)
  }

  invisible(NULL)
}

# ------------------------------------------------------------------
# 绘图与组别配色
# ------------------------------------------------------------------

rna_group_palette <- stats::setNames(c("#A1C9F4", "#FF9F9B"), c(baseline_group, treatment_group))


rna_group_palette <- rna_group_palette[!duplicated(names(rna_group_palette))]
group1_color <- unname(rna_group_palette[[baseline_group]])
group2_color <- unname(rna_group_palette[[treatment_group]])
rna_de_palette <- c(
  "up" = "#FF9F9B",
  "down" = "#A1C9F4",
  "ns" = "#D4D9DE",
  "Neutral" = "#D4D9DE"
)
rna_de_palette[comparison_open_label] <- "#FF9F9B"
rna_de_palette[reference_open_label] <- "#A1C9F4"

rna_de_palette <- rna_de_palette[!duplicated(names(rna_de_palette))]
rna_sig_palette <- c("significant" = "#FF9F9B", "Significant" = "#FF9F9B", "ns" = "#D4D9DE", "Not Significant" = "#D4D9DE")
rna_neutral_palette <- c("light" = "#FFFDF8", "mid" = "#D4D9DE", "dark" = "#5C6B73")
rna_accent_palette <- c("warm" = "#FFB482", "green" = "#8DE5A1", "purple" = "#B39FDB", "yellow" = "#FDFFB6", "teal" = "#9ED9CC", "peach" = "#FFD3A8")
rna_red_yellow_blue_palette <- c("low" = "#A1C9F4", "mid" = "#FDFFB6", "high" = "#FF9F9B")
rna_red_yellow_blue_binary_palette <- c(
  "zero" = unname(rna_red_yellow_blue_palette["low"]),
  "mid" = unname(rna_red_yellow_blue_palette["mid"]),
  "one" = unname(rna_red_yellow_blue_palette["high"])
)

feature_core_pal <- c(
  "Promoter" = "#FF9F9B",
  "UTR" = "#FFD3A8",
  "Exon" = "#FFB482",
  "Intron" = "#A1C9F4",
  "Downstream" = "#B39FDB",
  "Intergenic" = "#8DE5A1",
  "Other" = "#D4D9DE"
)
feature_alias_pal <- c(
  "Upstream" = feature_core_pal["Promoter"],
  "5' UTR" = feature_core_pal["UTR"],
  "3' UTR" = feature_core_pal["UTR"],
  "Gene Body" = "#9ED9CC",
  "Distal" = feature_core_pal["Intergenic"]
)

viz_pal <- list(
  group = rna_group_palette,
  de = rna_de_palette,
  sig_status = rna_sig_palette,
  sig_gradient = c("low" = group1_color, "high" = group2_color),
  neutral = rna_neutral_palette,
  neutral_light = unname(rna_neutral_palette["light"]),
  neutral_mid = unname(rna_neutral_palette["mid"]),
  neutral_dark = unname(rna_neutral_palette["dark"]),
  accent = rna_accent_palette,
  feature_core = feature_core_pal,
  feature_alias = feature_alias_pal,
  feature_accent = unname(rna_accent_palette),
  feature = c(feature_core_pal, feature_alias_pal),
  dmr = c("Hyper" = unname(rna_accent_palette["warm"]), "Hypo" = group1_color, "Neutral" = unname(rna_neutral_palette["mid"])),
  quadrant = c(
    "Hyper-Up" = unname(rna_red_yellow_blue_palette["high"]),
    "Hyper-Down" = unname(rna_accent_palette["purple"]),
    "Hypo-Down" = unname(rna_accent_palette["green"]),
    "Hypo-Up" = unname(rna_red_yellow_blue_palette["low"]),
    "NS" = unname(rna_neutral_palette["mid"])
  ),
  cgi = c(
    "Island" = unname(rna_red_yellow_blue_palette["low"]),
    "N_Shore" = unname(rna_accent_palette["green"]),
    "S_Shore" = unname(rna_red_yellow_blue_palette["mid"]),
    "N_Shelf" = unname(rna_accent_palette["warm"]),
    "S_Shelf" = unname(rna_red_yellow_blue_palette["high"]),
    "Open_Sea" = unname(rna_neutral_palette["mid"]),
    "Others" = unname(rna_neutral_palette["mid"]),
    "NS" = unname(rna_neutral_palette["mid"])
  ),
  status = c("PASS" = "#8DE5A1", "WARN" = "#FFB482", "NA" = "#D4D9DE"),
  overlap = c(
    "Activated" = "#FF9F9B",
    "Repressed" = "#A1C9F4",
    "ATAC_up_RNA_down" = "#B39FDB",
    "ATAC_down_RNA_up" = "#8DE5A1",
    "Neutral" = "#D4D9DE"
  ),
  candidate_status = c(
    "DE_and_DAS_splicing" = "#FF9F9B",
    "DE_splicing" = "#A1C9F4",
    "DAS_splicing" = "#8DE5A1",
    "Background_splicing" = "#D4D9DE"
  ),
  tier = c("Tier1" = "#FF9F9B", "Tier2" = "#FFB482", "Tier3" = "#D4D9DE"),
  go_category = c("BP" = "#FF9F9B", "CC" = "#A1C9F4", "MF" = "#B39FDB"),
  enrich_p = rna_red_yellow_blue_palette,
  score = c("pos" = "#FF9F9B", "neg" = "#A1C9F4"),
  upset = c("sets" = group2_color, "main" = group1_color),
  upset_sets = group2_color,
  upset_main = group1_color,
  motif = unname(rna_accent_palette["green"]),
  burden_density = group1_color,
  burden_ratio = unname(rna_red_yellow_blue_palette["mid"]),
  chr = c(group1_color, unname(rna_accent_palette["warm"])),
  track_cgi = unname(rna_red_yellow_blue_palette["mid"]),
  diverging = rna_red_yellow_blue_palette,
  heatmap_diverging = rna_red_yellow_blue_palette,
  heatmap_sequential = rna_red_yellow_blue_palette,
  heatmap_binary = rna_red_yellow_blue_binary_palette,
  heat = colorRampPalette(c("#A1C9F4", "#FDFFB6", "#FF9F9B"))(100),
  heat_meth = colorRampPalette(c("#A1C9F4", "#FDFFB6", "#FF9F9B"))(100),
  heat_warm = colorRampPalette(c("#FFFDF8", "#FDFFB6", "#FF9F9B"))(50),
  heat_compare = colorRampPalette(c("#FFFDF8", "#A1C9F4", "#B39FDB"))(60),
  transition = colorRampPalette(c("#FFFDF8", "#FDFFB6", "#FFB482", "#FF9F9B"))(100),
  feature_fun = function(n) colorRampPalette(c(unname(rna_red_yellow_blue_palette), unname(rna_accent_palette)))(n),
  cgi_fun = function(n) colorRampPalette(c(group1_color, unname(rna_accent_palette["green"]), unname(rna_red_yellow_blue_palette["mid"]), unname(rna_accent_palette["warm"]), group2_color))(n),
  meth_bin_fun = function(n) colorRampPalette(unname(rna_red_yellow_blue_palette))(n)
)

group_colors <- viz_pal$group
extra_groups <- setdiff(unique(as.character(sample_manifest$group)), names(group_colors))
if(length(extra_groups) > 0) {
  extra_cols <- colorRampPalette(viz_pal$accent)(length(extra_groups))
  names(extra_cols) <- extra_groups
  group_colors <- c(group_colors, extra_cols)
}

build_color_palette_audit <- function() {
  semantic_palettes <- c(
    "group", "neutral", "accent", "feature_core", "feature_alias", "dmr",
    "quadrant", "cgi", "status", "overlap", "candidate_status",
    "go_category", "enrich_p", "score", "upset", "de", "sig_status",
    "heatmap_diverging", "heatmap_sequential", "heatmap_binary"
  )
  bind_rows(lapply(semantic_palettes, function(palette_name) {
    pal <- viz_pal[[palette_name]]
    if(is.null(pal) || !is.atomic(pal)) return(data.frame())
    pal <- as.character(pal)
    names(pal) <- names(viz_pal[[palette_name]]) %||% seq_along(pal)
    data.frame(
      Palette = palette_name,
      Key = names(pal),
      Hex = unname(pal),
      SourceStyle = "template figure contract",
      stringsAsFactors = FALSE
    ) %>%
      filter(grepl("^#[0-9A-Fa-f]{6}$", Hex))
  }))
}

write_color_palette_audit <- function(step = "Step16") {
  save_step_table(build_color_palette_audit(), step, "Color_Palette_Audit.csv")
}

viz_style <- list(
  base_size = atac_figure_base_size,
  compact_size = atac_figure_tick_size,
  line = 0.5,
  point = 2.8,
  alpha = 0.85,
  legend = "right"
)

theme_rnaseq_classic <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_classic(base_size = atac_figure_base_size, base_family = atac_figure_font_family) +
    theme(
      text = element_text(family = atac_figure_font_family, size = atac_figure_base_size, color = "#333333"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = atac_figure_base_size, color = "#333333"),
      axis.text = element_text(size = atac_figure_tick_size, color = "#333333"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "#333333", linewidth = 0.35),
      axis.ticks = element_line(color = "#333333", linewidth = 0.35),
      strip.background = element_rect(fill = viz_pal$neutral["light"], color = NA),
      strip.text = element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
      legend.position = legend,
      legend.title = element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
      legend.text = element_text(size = atac_figure_tick_size, color = "#333333"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
    )
}

theme_rnaseq_minimal <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_rnaseq_classic(base_size = base_size, legend = legend)
}

theme_rnaseq_bw <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_rnaseq_classic(base_size = base_size, legend = legend)
}

theme_atac_classic <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_rnaseq_classic(base_size = base_size, legend = legend)
}

theme_atac_void <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_void(base_size = atac_figure_base_size, base_family = atac_figure_font_family) +
    theme(
      text = element_text(family = atac_figure_font_family, size = atac_figure_base_size, color = "#333333"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      strip.background = element_rect(fill = viz_pal$neutral["light"], color = NA),
      strip.text = element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
      legend.position = legend,
      legend.title = element_text(size = atac_figure_base_size, face = "bold", color = "#333333"),
      legend.text = element_text(size = atac_figure_tick_size, color = "#333333"),
      legend.background = element_blank(), legend.key = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
    )
}

if(requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  ComplexHeatmap::ht_opt(
    heatmap_row_names_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_tick_size),
    heatmap_column_names_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_tick_size),
    heatmap_row_title_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_base_size),
    heatmap_column_title_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_base_size),
    legend_title_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_base_size, fontface = "bold"),
    legend_labels_gp = grid::gpar(fontfamily = atac_figure_font_family, fontsize = atac_figure_tick_size)
  )
}

viz_colors <- function(palette, keys) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & keys != ""]
  if(length(keys) == 0) return(setNames(character(0), character(0)))
  palette_values <- unname(palette)
  if(length(palette_values) == 0) palette_values <- unname(viz_pal$neutral["mid"])
  out <- setNames(rep(palette_values, length.out = length(keys)), keys)
  known_keys <- intersect(names(palette), keys)
  out[known_keys] <- unname(palette[known_keys])
  out
}

viz_color <- function(palette, key) {
  unname(viz_colors(palette, key)[[1]])
}

viz_rekey_colors <- function(palette, key_map) {
  out <- vapply(unname(key_map), function(one_key) viz_color(palette, one_key), character(1))
  names(out) <- names(key_map)
  out
}

heatmap_group_palette <- function(groups, palette = group_colors) {
  group_levels <- unique(as.character(groups))
  resolve_named_colors(group_levels, palette)
}

heatmap_group_annotation <- function(groups, palette = group_colors, annotation_name = "Group", which = c("column", "row"), ...) {
  which <- match.arg(which)
  group_colors_local <- heatmap_group_palette(groups, palette = palette)
  annotation_data <- setNames(list(as.character(groups)), annotation_name)
  annotation_colors <- setNames(list(group_colors_local), annotation_name)
  if(identical(which, "row")) {
    return(do.call(ComplexHeatmap::rowAnnotation, c(annotation_data, list(col = annotation_colors), list(...))))
  }
  do.call(ComplexHeatmap::HeatmapAnnotation, c(annotation_data, list(col = annotation_colors), list(...)))
}

heatmap_diverging_col_fun <- function(x = NULL, midpoint = 0, breaks = NULL, colors = unname(viz_pal$heatmap_diverging)) {
  if(!is.null(breaks)) return(circlize::colorRamp2(breaks, colors))
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if(length(values) == 0) return(circlize::colorRamp2(c(-1, midpoint, 1), colors))
  low_value <- min(values, na.rm = TRUE)
  high_value <- max(values, na.rm = TRUE)
  if(!is.finite(low_value) || !is.finite(high_value) || low_value == high_value) {
    low_value <- midpoint - 1
    high_value <- midpoint + 1
  } else if(midpoint <= low_value || midpoint >= high_value) {
    midpoint <- mean(c(low_value, high_value))
  }
  circlize::colorRamp2(c(low_value, midpoint, high_value), colors)
}

heatmap_sequential_col_fun <- function(x = NULL, breaks = NULL, colors = unname(viz_pal$heatmap_sequential)) {
  if(!is.null(breaks)) return(circlize::colorRamp2(breaks, colors))
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if(length(values) == 0) return(circlize::colorRamp2(c(0, 0.5, 1), colors))
  low_value <- min(values, na.rm = TRUE)
  high_value <- max(values, na.rm = TRUE)
  if(!is.finite(low_value) || !is.finite(high_value) || low_value == high_value) {
    low_value <- 0
    high_value <- max(1, high_value)
  }
  mid_value <- mean(c(low_value, high_value))
  circlize::colorRamp2(c(low_value, mid_value, high_value), colors)
}

heatmap_binary_col_fun <- function() {
  circlize::colorRamp2(c(0, 0.5, 1), unname(viz_pal$heatmap_binary))
}

scale_fill_heatmap_diverging <- function(midpoint = 0, limits = NULL, ...) {
  ggplot2::scale_fill_gradient2(
    low = unname(viz_pal$heatmap_diverging["low"]),
    mid = unname(viz_pal$heatmap_diverging["mid"]),
    high = unname(viz_pal$heatmap_diverging["high"]),
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

scale_fill_heatmap_sequential <- function(...) {
  ggplot2::scale_fill_gradientn(colors = unname(viz_pal$heatmap_sequential), ...)
}

scale_color_heatmap_sequential <- function(...) {
  ggplot2::scale_color_gradientn(colors = unname(viz_pal$heatmap_sequential), ...)
}

scale_color_heatmap_diverging <- function(midpoint = 0, limits = NULL, ...) {
  ggplot2::scale_color_gradient2(
    low = unname(viz_pal$heatmap_diverging["low"]),
    mid = unname(viz_pal$heatmap_diverging["mid"]),
    high = unname(viz_pal$heatmap_diverging["high"]),
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

scale_fill_heatmap_binary <- function(limits = c(0, 1), ...) {
  ggplot2::scale_fill_gradientn(
    colors = unname(viz_pal$heatmap_binary),
    values = c(0, 0.5, 1),
    limits = limits,
    ...
  )
}
