#!/usr/bin/env Rscript

metabo_run <- function(opts) {
  required <- c("yaml", "readr", "dplyr", "tidyr", "ggplot2", "pheatmap", "ggrepel", "knitr", "tibble", "scales", "limma")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)[1]) y else x
  root_dir <- normalizePath(opts$root_dir, winslash = "/", mustWork = TRUE)
  config_file <- normalizePath(opts$config_file, winslash = "/", mustWork = TRUE)
  cfg <- yaml::read_yaml(config_file)
  safe <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9_.+-]+", "_", x))
  num <- function(x, default = NA_real_) {
    y <- suppressWarnings(as.numeric(x))
    ifelse(is.na(y), default, y)
  }
  bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0 || is.na(x)[1]) return(default)
    tolower(as.character(x)[1]) %in% c("1", "true", "yes", "on")
  }
  cfg_get <- function(path, default = NULL) {
    out <- cfg
    for (part in strsplit(path, ".", fixed = TRUE)[[1]]) {
      if (is.null(out[[part]])) return(default)
      out <- out[[part]]
    }
    out
  }
  resolve_path <- function(path) {
    if (is.null(path) || !nzchar(as.character(path))) return("")
    path <- as.character(path)[1]
    if (grepl("^/", path)) normalizePath(path, winslash = "/", mustWork = FALSE)
    else normalizePath(file.path(root_dir, path), winslash = "/", mustWork = FALSE)
  }
  ensure_dir <- function(...) {
    for (d in unlist(list(...), use.names = FALSE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  write_csv2 <- function(x, path) {
    ensure_dir(dirname(path))
    readr::write_csv(as.data.frame(x), path, na = "NA")
  }
  write_tsv2 <- function(x, path) {
    ensure_dir(dirname(path))
    readr::write_tsv(as.data.frame(x), path, na = "NA")
  }
  write_matrix <- function(x, path, id = "feature_id") {
    out <- data.frame(setNames(list(rownames(x)), id), x, check.names = FALSE)
    write_csv2(out, path)
  }
  save_plot <- function(plot, base, width = 7, height = 5) {
    ensure_dir(dirname(base))
    ggplot2::ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf)
    ggplot2::ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300)
  }
  write_matrix_txt <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.table(x, path, sep = "\t", quote = FALSE, col.names = NA)
  }

  result_dir <- if (!is.null(opts$result_dir) && nzchar(opts$result_dir)) resolve_path(opts$result_dir) else resolve_path(cfg_get("result_dir", "result"))
  dirs <- list(
    audit = file.path(result_dir, "01_audit"),
    peak = file.path(result_dir, "02_peak_table"),
    qc = file.path(result_dir, "03_qc"),
    stat = file.path(result_dir, "04_statistics"),
    annotation = file.path(result_dir, "05_annotation"),
    pathway = file.path(result_dir, "06_pathway"),
    plot = file.path(result_dir, "plots")
  )
  ensure_dir(unlist(dirs, use.names = FALSE))

  raw_dir <- resolve_path(cfg_get("raw_dir", "01rawdata"))
  sample_info_file <- resolve_path(cfg_get("sample_info_file", "01rawdata/sample_info.tsv"))
  mzml_dir <- resolve_path(cfg_get("mzml_dir", "02mzML"))
  peak_dir <- resolve_path(cfg_get("peak_dir", "02peak_tables"))
  case_group <- as.character(cfg_get("case_group"))
  control_group <- as.character(cfg_get("control_group"))
  qc_group <- as.character(cfg_get("qc_group", "QC"))
  comparison_prefix <- paste0(safe(case_group), "_vs_", safe(control_group))
  direction_case <- paste0("up_in_", safe(case_group))
  direction_control <- paste0("up_in_", safe(control_group))
  palette <- c(setNames("#A1C9F4", control_group), setNames("#FF9F9B", case_group), setNames("#8DE5A1", qc_group), not_sig = "#D4D9DE")
  group_palette <- function(x) {
    z <- palette[as.character(x)]
    z[is.na(z)] <- "#D4D9DE"
    unname(z)
  }
  direction_palette <- c(setNames("#FF9F9B", direction_case), setNames("#A1C9F4", direction_control), not_sig = "#D4D9DE", exploratory = "#FFB482", robust = "#8DE5A1")

  get_opt <- function(name, path, default) {
    if (!is.null(opts[[name]]) && length(opts[[name]]) && !is.na(opts[[name]])) return(opts[[name]])
    cfg_get(path, default)
  }
  qc_role <- as.character(get_opt("qc_role", "pipeline.default_qc_role", "unverified"))
  qc_filter <- as.character(get_opt("qc_filter", "pipeline.default_qc_filter", "off"))
  drift_correction <- as.character(get_opt("drift_correction", "pipeline.default_drift_correction", "none"))
  normalization <- as.character(get_opt("normalization", "pipeline.default_normalization", "median"))
  run_order_adjustment <- as.character(get_opt("run_order_adjustment", "pipeline.default_run_order_adjustment", "off"))
  significance <- as.character(get_opt("significance", "pipeline.default_significance", "nominal"))
  company_comparison <- as.character(get_opt("company_comparison", "pipeline.default_company_comparison", "off"))
  detection_cutoff <- num(get_opt("detection_rate", "filter.min_group_detection_rate", 0.8), 0.8)
  qc_rsd_cutoff <- num(get_opt("qc_rsd_cutoff", "filter.qc_rsd_cutoff", 0.3), 0.3)
  fc_cutoff <- num(get_opt("fc_cutoff", "filter.fold_change_cutoff", 1.5), 1.5)
  p_cutoff <- num(get_opt("p_cutoff", "filter.p_cutoff", 0.05), 0.05)
  fdr_cutoff <- num(get_opt("fdr_cutoff", "filter.fdr_cutoff", 0.05), 0.05)
  if (!qc_role %in% c("unverified", "pooled")) stop("qc-role must be unverified or pooled", call. = FALSE)
  if (!qc_filter %in% c("off", "rsd")) stop("qc-filter must be off or rsd", call. = FALSE)
  if (!drift_correction %in% c("none", "linear")) stop("drift-correction must be none or linear", call. = FALSE)
  if (!normalization %in% c("median", "pqn", "none")) stop("normalization must be median, pqn, or none", call. = FALSE)
  if (!run_order_adjustment %in% c("off", "on")) stop("run-order-adjustment must be off or on", call. = FALSE)
  if (!significance %in% c("dual", "fdr", "nominal")) stop("significance must be dual, fdr, or nominal", call. = FALSE)
  if (!company_comparison %in% c("off", "on")) stop("company-comparison must be off or on", call. = FALSE)
  if (qc_filter == "rsd" && qc_role != "pooled") stop("qc-filter=rsd requires qc-role=pooled", call. = FALSE)
  if (drift_correction == "linear" && qc_role != "pooled") stop("drift-correction=linear requires qc-role=pooled", call. = FALSE)

  # 样本标识完全来自输入表，不能根据名称改变组别或删除字符。
  clean_sample <- function(x) as.character(x)
  read_sample_info <- function(samples = character()) {
    info <- readr::read_tsv(sample_info_file, show_col_types = FALSE, progress = FALSE)
    if ("include" %in% names(info)) info <- info[tolower(as.character(info$include)) %in% c("1","true","yes"), , drop=FALSE]
    if (length(samples)) {
      if (!setequal(samples, info$sample)) stop("峰表样本与样本表不一致")
      info <- info[match(samples, info$sample), , drop=FALSE]
    }
    info
  }
  read_mzml_manifest <- function() {
    m <- readr::read_tsv(resolve_path(cfg_get("input.mzml_manifest")), show_col_types=FALSE)
    m$file <- vapply(m$file, resolve_path, character(1))
    m
  }
  inventory <- function() {
    info <- read_sample_info()
    mz <- if (cfg_get("input.mode") == "mzml") read_mzml_manifest() else data.frame(file=character(), sample=character(), polarity=character())
    if (nrow(mz)) mz <- left_join(mz, info |> select(any_of(c("sample", "run_order"))), by="sample")
    if (!"run_order" %in% names(mz)) mz$run_order <- rep(NA_real_,nrow(mz))
    raw_df <- data.frame(raw_file=character(),sample=character(),polarity=character())
    write_csv2(mz, file.path(dirs$audit,"mzml_manifest.csv"))
    write_tsv2(info, file.path(dirs$audit,"sample_info_used.tsv"))
    write_csv2(data.frame(metric=c("mzml_files","samples",case_group,control_group,qc_group), value=c(nrow(mz),nrow(info),sum(info$group==case_group),sum(info$group==control_group),sum(info$group==qc_group))), file.path(dirs$audit,"inventory_summary.csv"))
    list(raw=raw_df,mzml=mz,samples=info)
  }
  tool_audit <- function() {
    message("[audit] checking packages")
    packages <- c("yaml", "readr", "dplyr", "ggplot2", "pheatmap", "ggrepel", "limma", "xcms", "BiocParallel")
    pkg_version <- vapply(packages, function(p) {
      d <- suppressWarnings(utils::packageDescription(p, fields = "Version"))
      if (is.null(d) || is.na(d)) "" else as.character(d)
    }, character(1))
    pkg <- data.frame(tool = packages, available = nzchar(pkg_version), version = pkg_version)
    bins <- data.frame(tool="Rscript",path=Sys.which("Rscript"))
    write_tsv2(pkg, file.path(dirs$audit, "software_versions.tsv"))
    write_tsv2(bins, file.path(dirs$audit, "external_tools.tsv"))
    writeLines(capture.output(sessionInfo()), file.path(dirs$audit, "R_sessionInfo.txt"))
    message("[audit] package check complete")
    invisible(pkg)
  }
  effective_config <- function() {
    message("[audit] writing effective config")
    vals <- data.frame(parameter = c("config_file", "result_dir", "qc_role", "qc_filter", "qc_rsd_cutoff", "drift_correction", "normalization", "run_order_adjustment", "significance", "company_comparison", "detection_cutoff", "fc_cutoff", "p_cutoff", "fdr_cutoff", "case_group", "control_group", "qc_group"), value = c(config_file, result_dir, qc_role, qc_filter, qc_rsd_cutoff, drift_correction, normalization, run_order_adjustment, significance, company_comparison, detection_cutoff, fc_cutoff, p_cutoff, fdr_cutoff, case_group, control_group, qc_group))
    write_tsv2(vals, file.path(dirs$audit, "effective_config.tsv"))
    yaml::write_yaml(as.list(setNames(as.character(vals$value), vals$parameter)), file.path(dirs$audit, "effective_config.yml"))
    message("[audit] effective config complete")
  }
  manifest <- function(paths) {
    message("[audit] building input manifest")
    paths <- unique(paths[file.exists(paths)])
    if (!length(paths)) return(invisible(NULL))
    sizes <- file.info(paths)$size
    sha <- vapply(seq_along(paths), function(i) {
      if (!grepl("\\.(yml|yaml|tsv|csv)$", paths[i], ignore.case = TRUE)) return("deferred_large_file")
      value <- tryCatch(system2("sha256sum", c("--", shQuote(paths[i])), stdout = TRUE, stderr = FALSE), error = function(e) character())
      if (length(value) != 1L || !nzchar(value)) "unavailable" else sub("  .*", "", value)
    }, character(1))
    write_tsv2(data.frame(path = paths, sha256 = sha, size_bytes = sizes, modified = file.info(paths)$mtime), file.path(dirs$audit, "input_manifest.tsv"))
    message("[audit] input manifest complete")
  }

  run_xcms_modern <- function() {
    if (!requireNamespace("xcms", quietly = TRUE) || !requireNamespace("BiocParallel", quietly = TRUE)) stop("XCMS and BiocParallel are required to regenerate the peak table from mzML.", call. = FALSE)
    if (!requireNamespace("MsExperiment", quietly=TRUE)) stop("mzML 模式需要 MsExperiment")
    fmap <- read_mzml_manifest()
    if (!nrow(fmap)) stop("mzML manifest 为空")
    write_csv2(fmap, file.path(dirs$peak, "mzml_file_map.csv"))
    xcfg <- cfg_get("xcms", list())
    getx <- function(n, d) xcfg[[n]] %||% d
    workers <- as.integer(opts$threads %||% getx("workers", 4))
    if (workers > 1 && .Platform$OS.type != "windows") BiocParallel::register(BiocParallel::MulticoreParam(workers = workers)) else BiocParallel::register(BiocParallel::SerialParam())
    tables <- list()
    histories <- list()
    for (pol in unique(fmap$polarity)) {
      sub <- fmap[fmap$polarity == pol, , drop = FALSE]
      message("[xcms] processing ", pol, " (", nrow(sub), " files)")
      if (anyDuplicated(sub$sample)) stop("Duplicated sample names in ", pol, call. = FALSE)
      pd <- data.frame(sample_name = sub$sample, sample_group = read_sample_info()$group[match(sub$sample, read_sample_info()$sample)], stringsAsFactors = FALSE, row.names = basename(sub$file))
      xe <- MsExperiment::readMsExperiment(spectraFiles = sub$file, sampleData = pd)
      cwp <- xcms::CentWaveParam(ppm = num(getx("ppm", 15), 15), peakwidth = as.numeric(unlist(getx("peakwidth", c(5, 40)))), snthresh = num(getx("snthresh", 10), 10), prefilter = as.numeric(unlist(getx("prefilter", c(3, 100)))), mzdiff = num(getx("mzdiff", -0.001), -0.001))
      xe <- xcms::findChromPeaks(xe, param = cwp, msLevel = 1L, BPPARAM = BiocParallel::bpparam())
      if (nrow(sub) > 1 && !identical(getx("retention_correction", "obiwarp"), "none")) xe <- xcms::adjustRtime(xe, param = xcms::ObiwarpParam())
      pdp <- xcms::PeakDensityParam(sampleGroups = pd$sample_group, minFraction = num(getx("group_min_fraction", 0.5), 0.5), bw = num(getx("group_bw", 5), 5))
      xe <- xcms::groupChromPeaks(xe, param = pdp)
      detected <- xcms::featureValues(xe, value = "into", filled = FALSE, method = "sum")
      defs <- as.data.frame(xcms::featureDefinitions(xe), stringsAsFactors = FALSE)
      xe <- xcms::fillChromPeaks(xe, param = xcms::ChromPeakAreaParam())
      filled <- xcms::featureValues(xe, value = "into", filled = TRUE, method = "sum")
      if (is.null(dim(detected)) || !nrow(defs)) stop("No XCMS features detected for ", pol, call. = FALSE)
      colnames(detected) <- sub$sample; colnames(filled) <- sub$sample
      mz <- num(defs$mzmed, NA_real_); rt <- num(defs$rtmed, NA_real_)
      ids <- sprintf("%s_M%.6f_T%.3f", pol, mz, rt)
      rownames(detected) <- ids; rownames(filled) <- ids
      defs_out <- data.frame(feature_id = ids, mz = mz, rt = rt, polarity = pol, stringsAsFactors = FALSE)
      tables[[pol]] <- list(detected = detected, filled = filled, defs = defs_out)
      histories[[pol]] <- data.frame(polarity = pol, step = c("readMsExperiment", "findChromPeaks", "adjustRtime", "groupChromPeaks", "featureValues_filled_false", "fillChromPeaks", "featureValues_filled_true"), package = as.character(utils::packageVersion("xcms")), stringsAsFactors = FALSE)
      saveRDS(xe, file.path(dirs$peak, paste0("xcms_", pol, ".rds")))
      message("[xcms] completed ", pol, " (", nrow(defs_out), " features)")
    }
    all_samples <- unique(unlist(lapply(tables, function(x) colnames(x$filled))))
    bind_assay <- function(name) {
      out <- do.call(rbind, lapply(tables, function(x) {
        m <- x[[name]]; miss <- setdiff(all_samples, colnames(m)); if (length(miss)) m <- cbind(m, matrix(NA_real_, nrow(m), length(miss), dimnames = list(rownames(m), miss))); m[, all_samples, drop = FALSE]
      }))
      out
    }
    detected <- bind_assay("detected"); filled <- bind_assay("filled"); defs <- bind_rows(lapply(tables, `[[`, "defs"))
    write_matrix(detected, file.path(dirs$peak, "abundance_matrix_detected.csv"))
    write_matrix(filled, file.path(dirs$peak, "abundance_matrix_filled.csv"))
    write_csv2(defs, file.path(dirs$peak, "feature_annotation.csv"))
    write_csv2(data.frame(feature_id = rownames(filled), mz = defs$mz[match(rownames(filled), defs$feature_id)], rt = defs$rt[match(rownames(filled), defs$feature_id)], polarity = defs$polarity[match(rownames(filled), defs$feature_id)], filled = TRUE), file.path(peak_dir, "xcms_peak_table_combined.csv"))
    write_tsv2(bind_rows(histories), file.path(dirs$peak, "process_history.tsv"))
    write_csv2(data.frame(source="modern_xcms", detection_available=TRUE), file.path(dirs$peak,"peak_source.csv"))
    invisible(TRUE)
  }
  import_peak_table <- function() {
    filled <- resolve_path(cfg_get("input.filled_table"))
    detected <- resolve_path(cfg_get("input.detected_table", ""))
    annotation <- resolve_path(cfg_get("input.feature_table"))
    file.copy(filled,file.path(dirs$peak,"abundance_matrix_filled.csv"),overwrite=TRUE)
    if (nzchar(detected)) file.copy(detected,file.path(dirs$peak,"abundance_matrix_detected.csv"),overwrite=TRUE)
    file.copy(annotation,file.path(dirs$peak,"feature_annotation.csv"),overwrite=TRUE)
    write_csv2(data.frame(source="explicit_peak_table", detection_available=nzchar(detected)), file.path(dirs$peak,"peak_source.csv"))
  }
  load_peak <- function() {
    read_matrix <- function(path) {
      d <- readr::read_csv(path, show_col_types=FALSE, progress=FALSE)
      if (names(d)[1] != "feature_id" || anyDuplicated(d$feature_id)) stop("峰表需要唯一 feature_id")
      m <- as.matrix(d[,-1,drop=FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- d$feature_id
      if (any(is.infinite(m)) || any(m < 0,na.rm=TRUE)) stop("峰面积必须非负，缺失使用 NA")
      m
    }
    filled <- read_matrix(file.path(dirs$peak,"abundance_matrix_filled.csv"))
    source <- readr::read_csv(file.path(dirs$peak,"peak_source.csv"),show_col_types=FALSE)
    detected <- if (isTRUE(source$detection_available[[1]])) read_matrix(file.path(dirs$peak,"abundance_matrix_detected.csv")) else NULL
    if (!is.null(detected) && (!identical(dimnames(detected),dimnames(filled)))) stop("填峰前/后矩阵行列必须一致")
    ann <- readr::read_csv(file.path(dirs$peak,"feature_annotation.csv"),show_col_types=FALSE)
    if (!all(c("feature_id","mz","rt","polarity") %in% names(ann)) || anyDuplicated(ann$feature_id) || !setequal(ann$feature_id,rownames(filled))) stop("feature_table 与峰矩阵不一致")
    list(detected=detected,filled=filled,annotation=ann[match(rownames(filled),ann$feature_id),],source=source)
  }
  detection <- function(m) is.finite(m) & !is.na(m) & m > 0
  half_min <- function(m) {
    out <- m
    for (i in seq_len(nrow(out))) {
      x <- out[i, ]; obs <- x[is.finite(x) & !is.na(x) & x > 0]; fill <- if (length(obs)) min(obs) / 2 else 1; x[!is.finite(x) | is.na(x) | x <= 0] <- fill; out[i, ] <- x
    }
    out
  }
  fit_linear_drift <- function(m, si) {
    q <- which(si$group == qc_group & is.finite(si$run_order))
    if (length(q) < 3) stop("Linear drift correction requires at least 3 pooled QC samples", call. = FALSE)
    out <- m
    for (i in seq_len(nrow(m))) {
      y <- log2(pmax(m[i, ], .Machine$double.xmin)); qrun <- si$run_order[q]; fit <- tryCatch(stats::lm(y[q] ~ qrun), error = function(e) NULL)
      if (!is.null(fit) && all(is.finite(stats::coef(fit)))) {
        target <- stats::median(qrun); correction <- stats::predict(fit, newdata = data.frame(qrun = si$run_order)) - stats::predict(fit, newdata = data.frame(qrun = target)); out[i, ] <- 2^(y - correction)
      }
    }
    out
  }
  preprocess <- function() {
    peak <- load_peak(); si <- read_sample_info(colnames(peak$filled))
    if (!"run_order" %in% names(si)) si$run_order <- NA_real_
    si$run_order <- num(si$run_order)
    if ((drift_correction != "none" || run_order_adjustment == "on") && any(!is.finite(si$run_order))) stop("校正需要真实 run_order，不允许用行序代替")
    si$group <- as.character(si$group)
    det <- if (!is.null(peak$detected)) detection(peak$detected) else matrix(NA, nrow(peak$filled), ncol(peak$filled), dimnames=dimnames(peak$filled)); fill <- peak$filled
    feature <- peak$annotation |> select(any_of(c("feature_id", "mz", "rt", "polarity"))) |> distinct(feature_id, .keep_all = TRUE)
    feature$detection_rate <- rowMeans(det[match(feature$feature_id, rownames(det)), , drop = FALSE], na.rm = TRUE)
    for (g in c(case_group, control_group, qc_group)) {
      idx <- which(si$group == g); feature[[paste0("detection_rate_", g)]] <- if (length(idx)) rowMeans(det[, idx, drop = FALSE], na.rm = TRUE) else NA_real_
    }
    bpc_proxy <- apply(fill, 2, max, na.rm = TRUE); bpc_proxy[!is.finite(bpc_proxy)] <- NA_real_
    sample_qc <- data.frame(sample = colnames(fill), group = si$group[match(colnames(fill), si$sample)], run_order = si$run_order[match(colnames(fill), si$sample)], detected_features = colSums(det), missing_rate = colMeans(!det), total_signal = colSums(fill, na.rm = TRUE), tic_proxy_from_peak_area = colSums(fill, na.rm = TRUE), bpc_proxy_from_peak_area = bpc_proxy, stringsAsFactors = FALSE)
    qidx <- which(si$group == qc_group); bidx <- which(si$group %in% c(case_group, control_group))
    qc_cv <- data.frame(feature_id = rownames(fill), qc_mean = NA_real_, qc_rsd = NA_real_, stringsAsFactors = FALSE)
    if (length(qidx) >= 2) { qmean <- rowMeans(fill[, qidx, drop = FALSE], na.rm = TRUE); qc_cv$qc_mean <- qmean; qc_cv$qc_rsd <- apply(fill[, qidx, drop = FALSE], 1, stats::sd, na.rm = TRUE) / qmean }
    feature <- left_join(feature, qc_cv, by = "feature_id")
    feature$pass_detection <- if (is.null(peak$detected)) rep(TRUE,nrow(feature)) else feature[[paste0("detection_rate_", case_group)]] >= detection_cutoff | feature[[paste0("detection_rate_", control_group)]] >= detection_cutoff
    feature$detection_status <- if (is.null(peak$detected)) "not_available_no_prefill_matrix" else "prefill_measured"
    feature$detection_rate[!is.finite(feature$detection_rate)] <- NA_real_
    feature$pass_qc_filter <- if (qc_filter == "rsd") is.finite(feature$qc_rsd) & feature$qc_rsd <= qc_rsd_cutoff else TRUE
    feature$keep <- feature$pass_detection & feature$pass_qc_filter
    feature$filter_reason <- ifelse(feature$keep, "keep", ifelse(!feature$pass_detection, "low_group_detection", "qc_rsd"))
    write_csv2(sample_qc, file.path(dirs$qc, "sample_qc_summary.csv")); write_csv2(feature, file.path(dirs$qc, "feature_filter_decision.csv")); write_csv2(qc_cv, file.path(dirs$qc, "qc_rsd_summary.csv"))
    keep_ids <- feature$feature_id[feature$keep]; fill <- fill[keep_ids, , drop = FALSE]
    if (!nrow(fill)) stop("过滤后无特征")
    fill <- fill[rowSums(is.finite(fill) & fill>0)>0,,drop=FALSE]
    if (!nrow(fill)) stop("无正峰面积特征")
    if (drift_correction == "linear") fill <- fit_linear_drift(fill, si)
    imputed <- half_min(fill)
    pol_by_feature <- peak$annotation$polarity[match(rownames(imputed), peak$annotation$feature_id)]
    pols <- unique(pol_by_feature[!is.na(pol_by_feature)])
    factors <- bind_rows(lapply(pols, function(pol) data.frame(sample = colnames(imputed), polarity = pol, normalization_factor = 1, method = normalization, stringsAsFactors = FALSE)))
    norm <- imputed
    for (pol in pols) {
      rows <- which(peak$annotation$polarity[match(rownames(imputed), peak$annotation$feature_id)] == pol); if (!length(rows)) next
      x <- imputed[rows, , drop = FALSE]
      fac <- if (normalization == "median") apply(x, 2, median, na.rm = TRUE) / stats::median(x, na.rm = TRUE) else if (normalization == "pqn") { ref <- apply(log2(x), 1, median, na.rm = TRUE); q <- sweep(x, 1, 2^ref, "/"); apply(q, 2, median, na.rm = TRUE) / stats::median(q, na.rm = TRUE) } else rep(1, ncol(x))
      fac[!is.finite(fac) | fac <= 0] <- 1; norm[rows, ] <- sweep(x, 2, fac, "/"); factors$normalization_factor[factors$polarity == pol] <- fac[match(factors$sample[factors$polarity == pol], colnames(x))]
    }
    log2m <- log2(pmax(norm, .Machine$double.xmin)); centered <- sweep(log2m, 1, rowMeans(log2m), "-"); sdv <- apply(centered, 1, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1; scaled <- sweep(centered, 1, sdv, "/")
    rla_summary <- data.frame(sample = colnames(log2m), median_rla = apply(centered, 2, stats::median, na.rm = TRUE), iqr_rla = apply(centered, 2, stats::IQR, na.rm = TRUE), mean_abs_rla = colMeans(abs(centered), na.rm = TRUE), stringsAsFactors = FALSE) |> left_join(sample_qc |> select(sample, group, run_order), by = "sample")
    signal_summary <- sample_qc |> select(sample, group, run_order, tic_proxy_from_peak_area, bpc_proxy_from_peak_area) |> mutate(source = "peak-area proxy; not raw chromatographic TIC/BPC")
    source_label <- if (nrow(peak$source)) as.character(peak$source$source[[1]]) else "modern_xcms"
    spectra_file <- file.path(dirs$audit, "mzml_spectra_summary.csv")
    spectra_summary <- if (file.exists(spectra_file)) {
      sp <- readr::read_csv(spectra_file, show_col_types = FALSE) |> group_by(sample) |> summarise(ms1_spectra = sum(ms1_spectra, na.rm = TRUE), ms2_spectra = sum(ms2_spectra, na.rm = TRUE), .groups = "drop")
      data.frame(sample = colnames(fill), stringsAsFactors = FALSE) |> left_join(sp, by = "sample") |> mutate(status = ifelse(is.finite(ms1_spectra) & is.finite(ms2_spectra), "mzml_metadata", "not_available"))
    } else data.frame(sample = colnames(fill), ms1_spectra = NA_integer_, ms2_spectra = NA_integer_, status = ifelse(grepl("explicit_peak_table", source_label, ignore.case = TRUE), "not_available_from_cached_peak_table", "not_collected_in_peak_export"), stringsAsFactors = FALSE)
    rt_method <- as.character(cfg_get("xcms.retention_correction", "obiwarp"))
    rt_summary <- data.frame(polarity = c("POS", "NEG"), retention_correction = rt_method, status = ifelse(grepl("explicit_peak_table", source_label, ignore.case = TRUE), "not_run_external_peak_table", "executed_in_xcms"), stringsAsFactors = FALSE)
    bio_sd <- if (length(bidx)) apply(log2m[, bidx, drop = FALSE], 1, stats::sd, na.rm = TRUE) else rep(NA_real_, nrow(log2m)); qc_sd <- if (length(qidx)) apply(log2m[, qidx, drop = FALSE], 1, stats::sd, na.rm = TRUE) else rep(NA_real_, nrow(log2m)); dratio <- qc_sd / bio_sd
    write_csv2(factors |> left_join(sample_qc, by = "sample"), file.path(dirs$qc, "normalization_factors.csv")); write_csv2(signal_summary, file.path(dirs$qc, "sample_signal_summary.csv")); write_csv2(rla_summary, file.path(dirs$qc, "rla_summary.csv")); write_csv2(spectra_summary, file.path(dirs$qc, "spectra_summary.csv")); write_csv2(rt_summary, file.path(dirs$qc, "rt_correction_summary.csv")); write_matrix(imputed, file.path(dirs$qc, "matrix_imputed.csv")); write_matrix(norm, file.path(dirs$qc, "matrix_normalized.csv")); write_matrix(log2m, file.path(dirs$qc, "matrix_log2.csv")); write_matrix(scaled, file.path(dirs$qc, "matrix_scaled.csv")); write_csv2(data.frame(feature_id = rownames(log2m), d_ratio = dratio), file.path(dirs$qc, "d_ratio_summary.csv")); write_tsv2(si, file.path(dirs$qc, "sample_info_used.tsv"))
    saveRDS(list(raw_detected = peak$detected, filled = peak$filled, filtered = fill, imputed = imputed, norm = norm, log2 = log2m, scaled = scaled, annotation = peak$annotation, feature_qc = feature, qc_cv = qc_cv, d_ratio = dratio, sample_qc = sample_qc, signal_summary = signal_summary, rla_summary = rla_summary, spectra_summary = spectra_summary, rt_summary = rt_summary, sample_info = si, source = peak$source), file.path(dirs$qc, "processed_data.rds"))
    invisible(TRUE)
  }

  plot_qc <- function() {
    p <- readRDS(file.path(dirs$qc, "processed_data.rds")); sq <- p$sample_qc
    p1 <- ggplot(sq |> arrange(run_order), aes(run_order, total_signal, color = group, label = sample)) + geom_line(aes(group = 1), color = "grey75", linewidth = .35) + geom_point(size = 2.5) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = "真实进样顺序", y = "总峰面积", color = NULL)
    save_plot(p1, file.path(dirs$plot, "sample_run_order_total_signal"), 7.2, 4.8)
    long <- as.data.frame(p$log2, check.names = FALSE) |> tibble::rownames_to_column("feature_id") |> pivot_longer(-feature_id, names_to = "sample", values_to = "log2_intensity") |> left_join(p$sample_info, by = "sample")
    p2 <- ggplot(long, aes(sample, log2_intensity, fill = group)) + geom_boxplot(outlier.size = .2, linewidth = .2) + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = palette) + labs(x = NULL, y = "log2 normalized intensity", fill = NULL)
    save_plot(p2, file.path(dirs$plot, "sample_log2_intensity_boxplot"), 8, 4.8)
    if (any(is.finite(p$qc_cv$qc_rsd))) {
    p3 <- ggplot(p$qc_cv, aes(qc_rsd)) + geom_histogram(bins = 40, fill = "#8DE5A1", color = "white") + geom_vline(xintercept = qc_rsd_cutoff, linetype = 2, color = "#5C6B73") + theme_classic(base_size = 10) + labs(x = "QC RSD", y = "Feature count")
    save_plot(p3, file.path(dirs$plot, "qc_rsd_distribution"), 6, 4)
    }
    write_csv2(sq, file.path(dirs$qc, "sample_qc_summary.csv"))
    rla <- p$rla_summary; if (!is.null(rla) && nrow(rla)) { prla <- ggplot(rla, aes(run_order, median_rla, color = group, label = sample)) + geom_hline(yintercept = 0, linetype = 2, color = "grey70") + geom_line(aes(group = 1), color = "grey75", linewidth = .35) + geom_point(size = 2.5) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = "真实进样顺序", y = "RLA 中位数", color = NULL); save_plot(prla, file.path(dirs$plot, "sample_rla_summary"), 7.2, 4.8) }
    for (pol in c("POS", "NEG")) {
      ids <- p$annotation$feature_id[p$annotation$polarity == pol]; ids <- intersect(ids, rownames(p$scaled)); if (length(ids) < 2) next
      x <- p$scaled[ids, , drop = FALSE]; keep <- apply(x, 1, stats::sd) > 0; x <- x[keep, , drop = FALSE]; if (nrow(x) < 2) next
      pc <- stats::prcomp(t(x), center = FALSE, scale. = FALSE); vr <- pc$sdev^2 / sum(pc$sdev^2); sc <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2]) |> left_join(p$sample_info, by = "sample"); write_csv2(sc, file.path(dirs$qc, paste0("pca_scores_", pol, ".csv")))
      pp <- ggplot(sc, aes(PC1, PC2, color = group, label = sample)) + geom_point(size = 2.7) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = sprintf("PC1 (%.1f%%)", vr[1] * 100), y = sprintf("PC2 (%.1f%%)", vr[2] * 100), color = NULL, title = paste0(pol, " polarity PCA"))
      save_plot(pp, file.path(dirs$plot, paste0("pca_", pol)), 6.2, 5)
      cor_mat <- cor(p$log2[ids, , drop = FALSE], method = "spearman", use = "pairwise.complete.obs"); ann <- data.frame(group = p$sample_info$group); rownames(ann) <- p$sample_info$sample
      grDevices::pdf(file.path(dirs$plot, paste0("sample_spearman_correlation_", pol, ".pdf")), width = 7, height = 6); pheatmap::pheatmap(cor_mat, annotation_col = ann, annotation_row = ann, annotation_colors = list(group = palette), border_color = NA); grDevices::dev.off()
    }
    invisible(TRUE)
  }

  load_annotation_table <- function() {
    path <- resolve_path(cfg_get("annotation.table", "")); if (!nzchar(path) || !file.exists(path)) return(data.frame())
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
  annotate <- function() {
    peak <- load_peak(); ft <- peak$annotation |> distinct(feature_id, .keep_all = TRUE); id <- load_annotation_table()
    if (!nrow(id)) { ft$name <- NA_character_; ft$confidence_level <- "unannotated"; write_csv2(ft, file.path(dirs$annotation, "annotation_full.csv")); writeLines("No MetDNA identification table was available.", file.path(dirs$annotation, "annotation_status.txt")); return(invisible(TRUE)) }
    if ("confidence_level" %in% names(id)) {
      level1 <- as.character(id$confidence_level) %in% c("1","level 1","level1")
      confirmed <- if ("standard_confirmed" %in% names(id)) tolower(as.character(id$standard_confirmed)) %in% c("true","1","yes") else rep(FALSE,nrow(id))
      if (any(level1 & !confirmed)) stop("level 1 注释必须有 standard_confirmed 证据")
    }
    if (!"name" %in% names(id)) id$name <- NA_character_
    if (!"confidence_level" %in% names(id)) id$confidence_level <- "unverified"
    id_name <- intersect(c("peak_name", "metabolite", "name"), names(id))[1]; if (is.na(id_name)) stop("MetDNA table has no peak identifier column", call. = FALSE)
    id$match_feature_id <- as.character(id[[id_name]]); id <- id[!duplicated(id$match_feature_id), , drop = FALSE]
    out <- ft; out$annotation_match_type <- "unmatched"; out$annotation_mz_error_ppm <- NA_real_; out$annotation_rt_error_seconds <- NA_real_; out$annotation_ambiguous <- FALSE
    fields <- setdiff(names(id), c("match_feature_id", "mz", "rt")); for (f in fields) out[[f]] <- NA
    ## Resolve exact IDs in one vectorized lookup, then use a sorted m/z window
    ## for the remaining features. This preserves the matching rule while
    ## avoiding a full 26k x 1.6k distance matrix for every feature.
    match_idx <- match(as.character(out$feature_id), id$match_feature_id)
    exact_ok <- !is.na(match_idx); out$annotation_match_type[exact_ok] <- "exact_feature_id"
    mz_tol <- num(cfg_get("annotation.mz_tolerance_ppm", 15), 15); rt_tol <- num(cfg_get("annotation.rt_tolerance_seconds", 2), 2); ambig_delta <- num(cfg_get("annotation.match_ambiguity_delta", .05), .05)
    id_mz <- num(id$mz, NA_real_); id_rt <- num(id$rt, NA_real_); id_pol <- if ("polarity" %in% names(id)) as.character(id$polarity) else rep(NA_character_, nrow(id)); out_mz <- num(out$mz, NA_real_); out_rt <- num(out$rt, NA_real_); out_pol <- as.character(out$polarity)
    for (pol in unique(out_pol[!exact_ok])) {
      ii <- which(!exact_ok & out_pol == pol & is.finite(out_mz) & is.finite(out_rt)); if (!length(ii)) next
      cand_pool <- which(is.finite(id_mz) & is.finite(id_rt) & (is.na(id_pol) | id_pol == pol)); if (!length(cand_pool)) next
      ord_mz <- cand_pool[order(id_mz[cand_pool])]; mz_sorted <- id_mz[ord_mz]
      for (i in ii) {
        mz_window <- out_mz[i] * mz_tol / 1e6; lo <- findInterval(out_mz[i] - mz_window, mz_sorted); hi <- findInterval(out_mz[i] + mz_window, mz_sorted); if (hi < lo + 1L) next
        cand <- ord_mz[seq.int(max(1L, lo), min(length(ord_mz), hi))]; mzerr <- abs(id_mz[cand] - out_mz[i]) / out_mz[i] * 1e6; rterr <- abs(id_rt[cand] - out_rt[i]); ok <- is.finite(mzerr) & is.finite(rterr) & mzerr <= mz_tol & rterr <= rt_tol; if (!any(ok)) next
        cand <- cand[ok]; score <- mzerr[ok] / mz_tol + rterr[ok] / rt_tol; oo <- order(score); j <- cand[oo[1]]; match_idx[i] <- j; out$annotation_match_type[i] <- "nearest_mz_rt"; out$annotation_mz_error_ppm[i] <- mzerr[ok][oo[1]]; out$annotation_rt_error_seconds[i] <- rterr[ok][oo[1]]; if (length(oo) > 1L && (score[oo[2]] - score[oo[1]]) <= ambig_delta) out$annotation_ambiguous[i] <- TRUE
      }
    }
    for (f in fields) out[[f]] <- id[[f]][match_idx]
    out$annotated <- !is.na(out$name) & nzchar(as.character(out$name)); write_csv2(out, file.path(dirs$annotation, "annotation_full.csv")); write_csv2(out |> mutate(annotation_class = case_when(!annotated ~ "unannotated", !is.na(confidence_level) ~ as.character(confidence_level), TRUE ~ "annotated_no_confidence")) |> count(annotation_class, name = "n_features"), file.path(dirs$annotation, "annotation_summary.csv")); writeLines("MetDNA annotation was matched by exact feature ID first, then polarity/mz/RT tolerance; no correlation-based direction flipping was used.", file.path(dirs$annotation, "annotation_status.txt")); invisible(TRUE)
  }

  differential <- function() {
    p <- readRDS(file.path(dirs$qc, "processed_data.rds")); ann <- readr::read_csv(file.path(dirs$annotation, "annotation_full.csv"), show_col_types = FALSE); si <- p$sample_info; bio <- si$group %in% c(case_group, control_group); si_bio <- si[bio, , drop = FALSE]; mat <- p$log2[, si_bio$sample, drop = FALSE]; grp <- factor(si_bio$group, levels = c(control_group, case_group)); design <- model.matrix(~0 + grp); colnames(design) <- make.names(levels(grp)); if (run_order_adjustment == "on") { design <- cbind(design, run_order_scaled = as.numeric(scale(si_bio$run_order))) }
    if (qr(design)$rank < ncol(design) || nrow(design) <= ncol(design)) stop("统计设计不满秩或无残差自由度")
    fit <- limma::lmFit(mat, design); contrast_name <- paste0(make.names(case_group), "-", make.names(control_group)); cm <- limma::makeContrasts(contrasts = contrast_name, levels = design); fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), robust = TRUE); tt <- limma::topTable(fit2, number = Inf, sort.by = "none", adjust.method = "BH")
    n1 <- sum(grp == control_group); n2 <- sum(grp == case_group); auc <- function(v) { a <- v[grp == case_group]; b <- v[grp == control_group]; if (any(!is.finite(c(a, b)))) return(NA_real_); r <- rank(c(a, b), ties.method = "average"); u <- sum(r[seq_len(n2)]) - n2 * (n2 + 1) / 2; u / (n1 * n2) }
    aucs <- apply(mat, 1, auc); logfc <- tt$logFC; fdr <- tt$adj.P.Val; pval <- tt$P.Value; feature <- p$annotation[match(rownames(mat), p$annotation$feature_id), , drop = FALSE]
    res <- data.frame(feature_id = rownames(mat), mz = feature$mz, rt = feature$rt, polarity = feature$polarity, mean_case = rowMeans(mat[, si_bio$group == case_group, drop = FALSE]), mean_control = rowMeans(mat[, si_bio$group == control_group, drop = FALSE]), log2FC = logfc, fold_change = 2^logfc, limma_p = pval, fdr = fdr, auc_case_vs_control = aucs, rank_biserial = 2 * aucs - 1, stringsAsFactors = FALSE)
    res <- left_join(res, p$qc_cv |> select(feature_id, qc_rsd), by = "feature_id")
    ## Keep the feature coordinates from the statistical matrix and avoid
    ## duplicate mz/rt/polarity columns from the annotation table.
    res <- left_join(res, ann |> select(-any_of(c("mz", "rt", "polarity"))), by = "feature_id")
    robust <- is.finite(res$fdr) & res$fdr < fdr_cutoff & abs(res$log2FC) >= log2(fc_cutoff); exploratory <- is.finite(res$limma_p) & res$limma_p < p_cutoff & abs(res$log2FC) >= log2(fc_cutoff)
    res$classification <- if (significance == "fdr") ifelse(robust, "robust", "not_sig") else if (significance == "nominal") ifelse(exploratory, "exploratory", "not_sig") else ifelse(robust, "robust", ifelse(exploratory, "exploratory", "not_sig"))
    res$direction <- ifelse(res$classification == "not_sig", "not_sig", ifelse(res$log2FC > 0, direction_case, direction_control)); res$qc_rsd_pass <- is.na(res$qc_rsd) | res$qc_rsd <= qc_rsd_cutoff; res$candidate_score <- -log10(pmax(res$limma_p, .Machine$double.xmin)) * abs(res$log2FC)
    if (company_comparison == "on") {
      cpath <- resolve_path(cfg_get("external_statistics_table", "")); if (file.exists(cpath)) { cdat <- readr::read_csv(cpath, show_col_types = FALSE); cid <- intersect(c("name", "metabolite", "peak_name"), names(cdat))[1]; if (!is.na(cid) && all(c("p_values", "fold_changes") %in% names(cdat))) { cs <- cdat |> transmute(feature_id = as.character(.data[[cid]]), company_p_value = num(p_values, NA_real_), company_fold_change_raw = num(fold_changes, NA_real_)); res <- left_join(res, cs, by = "feature_id"); res$company_direction_convention <- as.character(cfg_get("company_comparison_direction", "unverified")) } }
    }
    res <- arrange(res, limma_p); robust_out <- res |> filter(is.finite(fdr), fdr < fdr_cutoff, abs(log2FC) >= log2(fc_cutoff)); exploratory_out <- res |> filter(classification == "exploratory"); write_csv2(res, file.path(dirs$stat, paste0(comparison_prefix, "_differential_features.csv"))); write_csv2(robust_out, file.path(dirs$stat, paste0(comparison_prefix, "_robust_differential_features.csv"))); write_csv2(exploratory_out, file.path(dirs$stat, paste0(comparison_prefix, "_exploratory_candidates.csv"))); write_csv2(res |> count(classification, direction, name = "n_features"), file.path(dirs$stat, paste0(comparison_prefix, "_summary.csv")))
    company_cols <- intersect(c("feature_id", "name", "polarity", "log2FC", "fold_change", "limma_p", "fdr", "classification", "direction", "company_p_value", "company_fold_change_raw", "company_direction_convention"), names(res))
    company_match <- if (company_comparison == "on" && "company_p_value" %in% names(res)) res |> filter(!is.na(company_p_value)) |> select(all_of(company_cols)) else res[0, company_cols, drop = FALSE]
    write_csv2(company_match, file.path(dirs$stat, paste0(comparison_prefix, "_company_comparison.csv")))
    writeLines(if (company_comparison == "on") paste0("Company statistics were joined by the supplied identifier without direction reversal. Matched rows: ", nrow(company_match), ".") else "Company comparison disabled by configuration.", file.path(dirs$stat, "company_comparison_status.txt"))
    labels <- head(res |> filter(classification != "not_sig") |> arrange(limma_p), 12); volcano <- ggplot(res, aes(log2FC, -log10(pmax(limma_p, .Machine$double.xmin)), color = classification)) + geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey65") + geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), linetype = 2, color = "grey75") + geom_point(alpha = .75, size = 1.4) + ggrepel::geom_text_repel(data = labels, aes(label = ifelse(is.na(name), feature_id, name)), size = 2.7, max.overlaps = 20) + scale_color_manual(values = c(robust = "#5C6B73", exploratory = "#FFB482", not_sig = "#D4D9DE"), drop = FALSE) + theme_classic(base_size = 10) + labs(x = paste0("log2 fold change (", case_group, " - ", control_group, ")"), y = "-log10 limma P", color = NULL)
    save_plot(volcano, file.path(dirs$plot, paste0(comparison_prefix, "_volcano")), 7, 5.5)
    top <- head(res |> filter(classification != "not_sig") |> arrange(desc(candidate_score)) |> pull(feature_id), 20); if (length(top)) { long <- as.data.frame(p$log2[top, , drop = FALSE], check.names = FALSE) |> tibble::rownames_to_column("feature_id") |> pivot_longer(-feature_id, names_to = "sample", values_to = "log2_intensity") |> left_join(si, by = "sample") |> filter(group != qc_group); bx <- ggplot(long, aes(group, log2_intensity, fill = group)) + geom_boxplot(outlier.shape = NA, linewidth = .2) + geom_jitter(width = .1, size = 1, alpha = .85) + facet_wrap(~feature_id, scales = "free_y", ncol = 4) + theme_classic(base_size = 8) + scale_fill_manual(values = palette) + labs(x = NULL, y = "log2 normalized intensity", fill = NULL); save_plot(bx, file.path(dirs$plot, paste0(comparison_prefix, "_top_feature_boxplots")), 10, 8); hm <- p$scaled[top, si_bio$sample, drop = FALSE]; ann_col <- data.frame(group = si_bio$group); rownames(ann_col) <- si_bio$sample; grDevices::pdf(file.path(dirs$plot, paste0(comparison_prefix, "_top_feature_heatmap.pdf")), width = 8, height = 9); pheatmap::pheatmap(hm, annotation_col = ann_col, annotation_colors = list(group = palette), border_color = NA, show_colnames = TRUE); grDevices::dev.off() }
    invisible(TRUE)
  }

  pathway <- function() {
    qpath <- resolve_path(cfg_get("pathway.quantitative_table", "")); epath <- resolve_path(cfg_get("pathway_enrichment_file", "")); if (nzchar(epath) && file.exists(epath)) write_csv2(readr::read_csv(epath, show_col_types = FALSE), file.path(dirs$pathway, "external_metdna_pathway_enrichment.csv"))
    if (!nzchar(qpath) || !file.exists(qpath)) { writeLines("No external pathway quantitative table was available.", file.path(dirs$pathway, "pathway_status.txt")); return(invisible(TRUE)) }
    dat <- readr::read_csv(qpath,show_col_types=FALSE); pc <- intersect(c("pathway_name","pathway"),names(dat))[1]
    si <- read_sample_info(); cases <- si$sample[si$group==case_group]; controls <- si$sample[si$group==control_group]
    if (is.na(pc) || !all(c(cases,controls) %in% names(dat))) stop("通路定量表需要所有生物学样本列")
    vals <- lapply(seq_len(nrow(dat)),function(i) {
      x <- as.numeric(unlist(dat[i,cases])); y <- as.numeric(unlist(dat[i,controls])); z <- c(x,y)
      if (any(!is.finite(z)) || any(z < 0) || !any(z>0)) return(data.frame(pathway=dat[[pc]][i],log2FC=NA_real_,p=NA_real_))
      shift <- min(z[z>0])/2
      data.frame(pathway=dat[[pc]][i],log2FC=mean(log2(x+shift))-mean(log2(y+shift)),p=tryCatch(t.test(log2(x+shift),log2(y+shift))$p.value,error=function(e) NA_real_))
    }) |> bind_rows()
    vals$fdr <- p.adjust(vals$p,"BH"); vals$interpretation <- "exploratory_external_pathway_quantification"
    write_csv2(vals,file.path(dirs$pathway,"pathway_quantitative_statistics.csv"))
    membership <- resolve_path(cfg_get("pathway_membership_file", "")); if (!nzchar(membership) || !file.exists(membership)) writeLines("KEGG pathway membership file not supplied; external MetDNA enrichment is reported separately and no new ORA is claimed.", file.path(dirs$pathway, "pathway_ora_status.txt"))
    invisible(TRUE)
  }

  report <- function() {
    tables <- list.files(result_dir,"[.](csv|tsv)$",recursive=TRUE)
    figures <- list.files(result_dir,"[.]png$",recursive=TRUE)
    lines <- c(paste0("# ",cfg_get("project_name","LC-MS")," 分析报告"),"",
      paste0("比较方向：",case_group," - ",control_group,"；统计单位：独立生物学样本。"),
      paste0("默认筛选：",significance,"；原始 P 值候选为探索性，完整结果保留 BH/FDR。"),
      paste0("输入模式：",cfg_get("input.mode"),"；QC身份：",qc_role,"；漂移校正：",drift_correction,"。"),
      "未提供填峰前矩阵时，原始检出率及其过滤为 not_available，不用填峰值代替。",
      "注释证据不足时保留未鉴定特征；level 1 要求标准品确认。",
      "上游外部峰表并不证明本流程已执行 XCMS；真实 mzML 上游须独立验证。","","## 表格与执行状态","",paste0("- [",tables,"](",tables,")"),"","## 图形","",if(length(figures)) paste0("![",figures,"](",figures,")") else "本次无图，参见表格和模块状态。")
    writeLines(lines,file.path(result_dir,"report.md"))
  }

  module <- opts$module %||% "all"
  stages <- c("inventory","xcms","preprocess","plots","annotate","stats","pathway","report")
  selected <- switch(module, all=stages, upstream="xcms", downstream=stages[3:8], audit="inventory", module)
  status_file <- file.path(result_dir,"module_status.tsv")
  states <- if (file.exists(status_file)) read.delim(status_file,header=FALSE,col.names=c("module","status","time"),stringsAsFactors=FALSE) else data.frame(module=character(),status=character(),time=character())
  states <- states[!duplicated(states$module,fromLast=TRUE),,drop=FALSE]
  record <- function(stage,status) {
    cat(paste(stage,status,format(Sys.time(),tz="UTC"),sep="\t"),"\n",file=status_file,append=TRUE)
    states <<- rbind(states[states$module!=stage,,drop=FALSE],data.frame(module=stage,status=status,time=as.character(Sys.time())))
  }
  dependencies <- list(preprocess="xcms",plots="preprocess",annotate="preprocess",stats=c("preprocess","annotate"))
  for (stage in selected) {
    position <- match(stage,stages)
    if (is.na(position)) stop("未知模块: ",stage)
    if (stage != "report") {
      prior <- stages[seq_len(position-1L)]
      if (any(states$module %in% prior & states$status %in% c("failed","running"))) stop("上游失败或未完成，拒绝使用旧结果")
      required <- dependencies[[stage]]
      completed <- states$module[states$status=="completed"]
      if (length(required) && !all(required %in% completed)) stop("输入生产模块未完成: ",paste(setdiff(required,completed),collapse=","))
      later <- stages[seq_along(stages)>position]
      if (any(states$module %in% later & states$status=="completed")) stop("下游已完成；重跑上游需新的result_dir")
    }
    record(stage,"running")
    tryCatch({
      switch(stage,
        inventory={inv <- inventory();tool_audit();effective_config();manifest(c(inv$raw$raw_file,inv$mzml$file,sample_info_file,config_file))},
        xcms={if(cfg_get("input.mode")=="peak_table") import_peak_table() else run_xcms_modern()},
        preprocess=preprocess(), plots=plot_qc(), annotate=annotate(), stats=differential(), pathway=pathway(), report=report())
      record(stage,"completed")
    },error=function(e){record(stage,"failed");stop(e)})
  }
  message("Completed module: ", module)
  invisible(TRUE)
}
