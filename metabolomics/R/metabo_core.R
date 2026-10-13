#!/usr/bin/env Rscript

metabo_run <- function(opts) {
  required <- c("yaml", "readr", "dplyr", "tidyr", "ggplot2", "pheatmap", "ggrepel", "knitr", "tibble", "scales", "limma", "MsExperiment")
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
  plots_only <- identical(opts$module, "plots-only")
  source(file.path(root_dir, "R", "figure_style.R"), local = TRUE)
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
    if (plots_only && file.exists(path)) return(invisible(path))
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
    plot <- plot + theme_classic(base_size = 8, base_family = pfos_figure_style$family) +
      theme(plot.title = element_blank(), strip.text = element_text(size = 7))
    pfos_save(plot, paste0(base, ".pdf"), height_mm = max(height, 4) * 25.4)
  }
  save_heatmap <- function(mat, base, height, scale_kind = "row_z", ...) {
    extra <- list(...)
    if (!is.null(extra$annotation_colors)) for (field in names(extra$annotation_colors)) {
      observed <- unique(c(as.character(extra$annotation_col[[field]]), as.character(extra$annotation_row[[field]])))
      extra$annotation_colors[[field]] <- extra$annotation_colors[[field]][names(extra$annotation_colors[[field]]) %in% observed]
    }
    limit <- if (scale_kind == "correlation") 1 else pfos_profile$scales$row_zscore[[2]]
    if (!is.finite(limit) || limit == 0) limit <- 1
    breaks <- seq(-limit, limit, length.out = 101)
    draw <- function() {
      hm <- do.call(pheatmap::pheatmap, c(list(mat = mat, silent = TRUE, fontsize = 8,
        fontsize_row = 7, fontsize_col = 7, fontfamily = pfos_figure_style$family,
        color = grDevices::colorRampPalette(pfos_figure_style$diverging)(100), breaks = breaks,
        na_col = pfos_profile$scales$missing_color), extra))
      hm$gtable$widths[3] <- hm$gtable$widths[3] - grid::unit(12, "pt")
      grid::grid.newpage(); grid::grid.draw(hm$gtable)
    }
    pfos_save(draw, paste0(base, ".pdf"), height_mm = height * 25.4)
  }
  write_matrix_txt <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.table(x, path, sep = "\t", quote = FALSE, col.names = NA)
  }

  result_dir <- if (!is.null(opts$result_dir) && nzchar(opts$result_dir)) resolve_path(opts$result_dir) else resolve_path(cfg_get("result_dir", "result"))
  dirs <- list(
    audit = file.path(result_dir, "01_audit", "tables"),
    peak = file.path(result_dir, "02_peak_table", "tables"),
    qc = file.path(result_dir, "03_qc", "tables"),
    stat = file.path(result_dir, "04_statistics", "tables"),
    annotation = file.path(result_dir, "05_annotation", "tables"),
    pathway = file.path(result_dir, "06_pathway", "tables"),
    plot = file.path(result_dir, "03_qc", "plots"),
    audit_logs = file.path(result_dir, "01_audit", "logs"),
    annotation_logs = file.path(result_dir, "05_annotation", "logs"),
    stat_logs = file.path(result_dir, "04_statistics", "logs"),
    pathway_logs = file.path(result_dir, "06_pathway", "logs"),
    qc_objects = file.path(result_dir, "03_qc", "objects"),
    peak_objects = file.path(result_dir, "02_peak_table", "objects")
  )
  ensure_dir(unlist(dirs, use.names = FALSE))

  raw_dir <- resolve_path(cfg_get("raw_dir", "01rawdata"))
  sample_info_file <- resolve_path(cfg_get("sample_info_file", "01rawdata/sample_info.tsv"))
  mzml_dir <- resolve_path(cfg_get("mzml_dir", "02mzML"))
  peak_dir <- resolve_path(cfg_get("peak_dir", "02peak_tables"))
  case_group <- as.character(cfg_get("case_group", "PFOS"))
  control_group <- as.character(cfg_get("control_group", "WT"))
  qc_group <- as.character(cfg_get("qc_group", "QC"))
  comparison_prefix <- paste0(safe(case_group), "_vs_", safe(control_group))
  direction_case <- paste0("up_in_", safe(case_group))
  direction_control <- paste0("up_in_", safe(control_group))
  palette <- c(setNames(pfos_palette("group")[["WT"]], control_group), setNames(pfos_palette("group")[["PFOS"]], case_group), setNames(pfos_palette("accent")[["green"]], qc_group), Unknown = pfos_palette("de")[["ns"]], not_sig = pfos_palette("de")[["ns"]])
  group_palette <- function(x) {
    z <- palette[as.character(x)]
    z[is.na(z)] <- pfos_palette("de")[["ns"]]
    unname(z)
  }
  direction_palette <- c(setNames(pfos_palette("de")[["up"]], direction_case), setNames(pfos_palette("de")[["down"]], direction_control), not_sig = pfos_palette("de")[["ns"]], exploratory = pfos_palette("accent")[["warm"]], robust = pfos_palette("accent")[["green"]])

  get_opt <- function(name, path, default) {
    if (!is.null(opts[[name]]) && length(opts[[name]]) && !is.na(opts[[name]])) return(opts[[name]])
    cfg_get(path, default)
  }
  qc_role <- as.character(get_opt("qc_role", "pipeline.default_qc_role", "unverified"))
  qc_filter <- as.character(get_opt("qc_filter", "pipeline.default_qc_filter", "off"))
  drift_correction <- as.character(get_opt("drift_correction", "pipeline.default_drift_correction", "none"))
  normalization <- as.character(get_opt("normalization", "pipeline.default_normalization", "median"))
  run_order_adjustment <- as.character(get_opt("run_order_adjustment", "pipeline.default_run_order_adjustment", "off"))
  significance <- as.character(get_opt("significance", "pipeline.default_significance", "dual"))
  company_comparison <- as.character(get_opt("company_comparison", "pipeline.default_company_comparison", "on"))
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

  clean_stale <- bool(cfg_get("pipeline.clean_stale_outputs", TRUE), TRUE)
  if (isTRUE(opts$force) && clean_stale && identical(opts$module, "audit")) {
    old <- c(file.path(result_dir, "00_inventory"), file.path(result_dir, "01_peak_table"), file.path(result_dir, "02_qc"), file.path(result_dir, "03_statistics"), file.path(result_dir, "04_pathway"))
    old <- old[dir.exists(old)]
    if (length(old)) unlink(old, recursive = TRUE, force = TRUE)
    stale_plots <- if (dir.exists(dirs$plot)) list.files(dirs$plot, full.names = TRUE) else character()
    if (length(stale_plots)) unlink(stale_plots, recursive = TRUE, force = TRUE)
    ensure_dir(unlist(dirs, use.names = FALSE))
  }

  clean_sample <- function(x) {
    x <- tools::file_path_sans_ext(basename(x))
    x <- sub("\\.wiff$", "", x, ignore.case = TRUE)
    x <- sub("_(POS|NEG|Positive|Negative)$", "", x, ignore.case = TRUE)
    gsub("-", "", x)
  }
  infer_polarity <- function(path) {
    x <- toupper(path)
    if (grepl("(^|[/_.-])POS(ITIVE)?([/_.-]|$)", x)) return("POS")
    if (grepl("(^|[/_.-])NEG(ATIVE)?([/_.-]|$)", x)) return("NEG")
    "ALL"
  }
  detect_group <- function(x) {
    vapply(as.character(x), function(one) {
      one <- toupper(one)
      if (grepl("^PFOS", one)) case_group else if (grepl("^WT", one)) control_group else if (grepl("^QC", one)) qc_group else "Unknown"
    }, character(1))
  }
  read_sample_info <- function(samples = character()) {
    if (!file.exists(sample_info_file)) stop("Missing sample_info.tsv: ", sample_info_file, call. = FALSE)
    info <- readr::read_tsv(sample_info_file, show_col_types = FALSE, progress = FALSE)
    if (!all(c("sample", "group") %in% names(info))) stop("sample_info.tsv needs sample and group columns", call. = FALSE)
    info$sample <- clean_sample(info$sample)
    info$group <- as.character(info$group)
    info$include <- if ("include" %in% names(info)) !tolower(as.character(info$include)) %in% c("0", "false", "no") else TRUE
    info <- info[info$include, , drop = FALSE]
    if (length(samples)) {
      miss <- setdiff(samples, info$sample)
      if (length(miss)) stop("Samples missing from sample_info.tsv: ", paste(miss, collapse = ", "), call. = FALSE)
      info <- info[match(samples, info$sample), , drop = FALSE]
    }
    info
  }
  parse_timestamp <- function(path) {
    lines <- tryCatch(readLines(path, n = 400, warn = FALSE), error = function(e) character())
    hit <- grep("startTimeStamp=", lines, value = TRUE)
    if (!length(hit)) return(as.POSIXct(NA))
    value <- sub(".*startTimeStamp=\\\"([^\\\"]+)\\\".*", "\\1", hit[[1]])
    as.POSIXct(value, tz = "UTC", tryFormats = c("%Y-%m-%dT%H:%M:%OS", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"))
  }
  count_mzml_spectra <- function(path) {
    count_one <- function(accession) {
      temp <- tempfile(fileext = ".spectra")
      grep_bin <- Sys.which("grep"); if (!nzchar(grep_bin)) return(NA_integer_)
      tryCatch(system2(grep_bin, c("-o", accession, "--", path), stdout = temp, stderr = FALSE), error = function(e) NA_integer_)
      n <- if (file.exists(temp)) length(readLines(temp, warn = FALSE)) else NA_integer_
      if (file.exists(temp)) unlink(temp)
      n
    }
    data.frame(ms1_spectra = count_one("MS:1000511"), ms2_spectra = count_one("MS:1000515"), stringsAsFactors = FALSE)
  }
  inventory <- function() {
    message("[inventory] scanning WIFF")
    raw <- if (dir.exists(raw_dir)) list.files(raw_dir, "\\.wiff$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character()
    raw <- raw[!grepl("\\.scan$", raw, ignore.case = TRUE)]
    raw_df <- data.frame(raw_file = sort(raw), sample = clean_sample(raw), polarity = vapply(raw, infer_polarity, character(1)), stringsAsFactors = FALSE)
    sample_names <- unique(raw_df$sample)
    info <- read_sample_info(sample_names)
    message("[inventory] reading sample metadata")
    mzml <- if (dir.exists(mzml_dir)) list.files(mzml_dir, "\\.(mzML|mzXML)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character()
    mz <- data.frame(file = sort(mzml), sample = clean_sample(mzml), polarity = vapply(mzml, infer_polarity, character(1)), stringsAsFactors = FALSE)
    message("[inventory] parsing mzML timestamps: ", nrow(mz), " files")
    if (nrow(mz)) {
      mz$acquired_at <- as.character(vapply(mz$file, function(x) as.character(parse_timestamp(x)), character(1)))
      mz$timestamp_ok <- !is.na(as.POSIXct(mz$acquired_at, tz = "UTC"))
      mz$run_order <- NA_integer_
      if (any(mz$timestamp_ok)) mz$run_order <- rank(as.POSIXct(mz$acquired_at, tz = "UTC"), ties.method = "first", na.last = "keep")
      mz$file_size_bytes <- file.info(mz$file)$size
    }
    if (nrow(mz)) {
      message("[inventory] counting MS1/MS2 spectra")
      spectra <- bind_rows(lapply(mz$file, count_mzml_spectra))
      spectra_summary <- bind_cols(mz |> select(file, sample, polarity), spectra)
      write_csv2(spectra_summary, file.path(dirs$audit, "mzml_spectra_summary.csv"))
    } else spectra_summary <- data.frame(ms1_spectra = integer(), ms2_spectra = integer())
    sample_order <- if (nrow(mz) && any(!is.na(mz$run_order))) mz |> filter(!is.na(run_order)) |> arrange(run_order) |> distinct(sample, .keep_all = TRUE) |> select(sample, run_order, acquired_at) else data.frame(sample = sample_names, run_order = seq_along(sample_names), acquired_at = NA_character_)
    if (nrow(mz) && any(mz$timestamp_ok)) {
      by_pol <- split(mz, mz$polarity)
      if (length(by_pol) > 1 && any(vapply(by_pol, nrow, integer(1)) > 1)) {
        ords <- lapply(by_pol, function(x) paste(x |> arrange(run_order) |> distinct(sample) |> pull(sample), collapse = ","))
        if (length(unique(ords)) > 1) stop("POS/NEG mzML acquisition orders disagree", call. = FALSE)
      }
    }
    info <- info |> left_join(sample_order, by = "sample")
    message("[inventory] writing audit tables")
    write_csv2(raw_df |> left_join(info |> select(sample, group), by = "sample") |> mutate(has_scan = file.exists(paste0(raw_file, ".scan"))), file.path(dirs$audit, "wiff_inventory.csv"))
    write_csv2(mz, file.path(dirs$audit, "mzml_manifest.csv"))
    write_tsv2(info, file.path(dirs$audit, "sample_info_used.tsv"))
    write_csv2(data.frame(metric = c("wiff_files", "mzml_files", "samples", "PFOS", "WT", "QC", "ms1_spectra", "ms2_spectra", "polarity_groups"), value = c(nrow(raw_df), nrow(mz), nrow(info), sum(info$group == case_group), sum(info$group == control_group), sum(info$group == qc_group), sum(spectra_summary$ms1_spectra), sum(spectra_summary$ms2_spectra), paste(sort(unique(mz$polarity)), collapse = ";"))), file.path(dirs$audit, "inventory_summary.csv"))
    message("[inventory] complete")
    list(raw = raw_df, mzml = mz, samples = info)
  }
  tool_audit <- function() {
    message("[audit] checking packages")
    packages <- c("yaml", "readr", "dplyr", "ggplot2", "pheatmap", "ggrepel", "limma", "xcms", "MsExperiment", "BiocParallel")
    pkg_version <- vapply(packages, function(p) {
      d <- suppressWarnings(utils::packageDescription(p, fields = "Version"))
      if (is.null(d) || is.na(d)) "" else as.character(d)
    }, character(1))
    pkg <- data.frame(tool = packages, available = nzchar(pkg_version), version = pkg_version)
    bins <- data.frame(tool = c("Rscript", "singularity", "msconvert"), path = c(Sys.which("Rscript"), Sys.which("singularity"), Sys.which(as.character(cfg_get("msconvert_cmd", "msconvert")))))
    write_tsv2(pkg, file.path(dirs$audit, "software_versions.tsv"))
    write_tsv2(bins, file.path(dirs$audit, "external_tools.tsv"))
    writeLines(capture.output(sessionInfo()), file.path(dirs$audit_logs, "R_sessionInfo.txt"))
    message("[audit] package check complete")
    invisible(pkg)
  }
  effective_config <- function() {
    message("[audit] writing effective config")
    vals <- data.frame(parameter = c("config_file", "result_dir", "qc_role", "qc_filter", "qc_rsd_cutoff", "drift_correction", "normalization", "run_order_adjustment", "significance", "company_comparison", "detection_cutoff", "fc_cutoff", "p_cutoff", "fdr_cutoff", "case_group", "control_group", "qc_group"), value = c(config_file, result_dir, qc_role, qc_filter, qc_rsd_cutoff, drift_correction, normalization, run_order_adjustment, significance, company_comparison, detection_cutoff, fc_cutoff, p_cutoff, fdr_cutoff, case_group, control_group, qc_group))
    write_tsv2(vals, file.path(dirs$audit, "effective_config.tsv"))
    yaml::write_yaml(as.list(setNames(as.character(vals$value), vals$parameter)), file.path(dirs$audit_logs, "effective_config.yml"))
    message("[audit] effective config complete")
  }
  manifest <- function(paths) {
    message("[audit] building input manifest")
    paths <- unique(paths[file.exists(paths)])
    if (!length(paths)) return(invisible(NULL))
    sizes <- file.info(paths)$size
    sha <- vapply(seq_along(paths), function(i) {
      if (!grepl("\\.(yml|yaml|tsv|csv)$", paths[i], ignore.case = TRUE)) return("deferred_large_file")
      value <- tryCatch(system2("sha256sum", c("--", paths[i]), stdout = TRUE, stderr = FALSE), error = function(e) character())
      if (length(value) != 1L || !nzchar(value)) "unavailable" else sub("  .*", "", value)
    }, character(1))
    write_tsv2(data.frame(path = paths, sha256 = sha, size_bytes = sizes, modified = file.info(paths)$mtime), file.path(dirs$audit, "input_manifest.tsv"))
    message("[audit] input manifest complete")
  }

  convert_wiff <- function(inv) {
    if (nrow(inv$raw) == 0) return(invisible(FALSE))
    sif <- resolve_path(Sys.getenv("MSCONVERT_SIF", unset = cfg_get("msconvert_sif", "")))
    backend <- as.character(cfg_get("msconvert_backend", "singularity"))
    cmd <- as.character(cfg_get("msconvert_cmd", "msconvert"))
    available <- if (backend == "singularity") nzchar(Sys.which("singularity")) && file.exists(sif) else nzchar(Sys.which(cmd))
    if (!available) {
      writeLines("WIFF conversion not run: converter unavailable; existing mzML may still be used.", file.path(dirs$audit_logs, "wiff_conversion_status.txt"))
      return(invisible(FALSE))
    }
    logs <- lapply(inv$raw$raw_file, function(raw) {
      out <- file.path(mzml_dir, paste0(tools::file_path_sans_ext(basename(raw)), ".mzML"))
      if (file.exists(out) && file.info(out)$size > 0 && bool(cfg_get("msconvert_skip_existing", TRUE), TRUE) && !isTRUE(opts$force)) return(data.frame(raw_file = raw, mzml_file = out, status = "skipped_existing"))
      ensure_dir(dirname(out))
      if (backend == "singularity") {
        args <- c("exec", "--writable-tmpfs", sif, "wine", cfg_get("msconvert_container_exe", "/wineprefix64/drive_c/pwiz/msconvert.exe"), raw, "--mzML", "--64", "--zlib", "--filter", "peakPicking true 1-", "-o", dirname(out), "--outfile", basename(out))
        status <- system2("singularity", args)
      } else {
        status <- system2(cmd, c(raw, "--mzML", "--64", "--zlib", "--filter", "peakPicking true 1-", "-o", dirname(out), "--outfile", basename(out)))
      }
      data.frame(raw_file = raw, mzml_file = out, status = ifelse(identical(status, 0L), "converted", paste0("failed_", status)))
    }) |> bind_rows()
    write_csv2(logs, file.path(dirs$audit, "wiff_conversion_log.csv"))
    writeLines(paste0("WIFF conversion backend: ", backend), file.path(dirs$audit_logs, "wiff_conversion_status.txt"))
    invisible(TRUE)
  }

  run_xcms_modern <- function() {
    if (!requireNamespace("xcms", quietly = TRUE) || !requireNamespace("BiocParallel", quietly = TRUE)) stop("XCMS and BiocParallel are required to regenerate the peak table from mzML.", call. = FALSE)
    files <- list.files(mzml_dir, "\\.(mzML|mzXML)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    if (!length(files)) stop("No mzML/mzXML files found in ", mzml_dir, call. = FALSE)
    fmap <- data.frame(file = sort(files), sample = clean_sample(files), polarity = vapply(files, infer_polarity, character(1)), stringsAsFactors = FALSE)
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
      pd <- data.frame(sample_name = sub$sample, sample_group = detect_group(sub$sample), stringsAsFactors = FALSE, row.names = basename(sub$file))
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
      saveRDS(xe, file.path(dirs$peak_objects, paste0("xcms_", pol, ".rds")))
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
    invisible(TRUE)
  }
  hydrate_cached_peak_table <- function() {
    cached <- file.path(peak_dir, "xcms_peak_table_combined.csv")
    if (!file.exists(cached)) stop("No cached XCMS peak table is available; run the forced modern XCMS step.", call. = FALSE)
    dat <- readr::read_csv(cached, show_col_types = FALSE, progress = FALSE)
    id_col <- intersect(c("feature_id", "metabolite", "peak_name"), names(dat))[1]
    if (is.na(id_col)) stop("Cached XCMS peak table has no feature identifier column", call. = FALSE)
    meta_cols <- intersect(c(id_col, "mz", "rt", "polarity"), names(dat)); sample_cols <- setdiff(names(dat), meta_cols)
    ids <- make.unique(as.character(dat[[id_col]])); m <- as.matrix(dat[, sample_cols, drop = FALSE]); suppressWarnings(storage.mode(m) <- "numeric"); rownames(m) <- ids
    ann <- data.frame(feature_id = ids, mz = num(dat$mz, NA_real_), rt = num(dat$rt, NA_real_), polarity = if ("polarity" %in% names(dat)) as.character(dat$polarity) else "ALL", stringsAsFactors = FALSE)
    write_matrix(m, file.path(dirs$peak, "abundance_matrix_detected.csv")); write_matrix(m, file.path(dirs$peak, "abundance_matrix_filled.csv")); write_csv2(ann, file.path(dirs$peak, "feature_annotation.csv")); write_csv2(data.frame(source = "legacy_cached_peak_table_gap_filled_unknown", detected_matrix_source = "cached_table_not_raw_detection"), file.path(dirs$peak, "peak_source.csv")); write_tsv2(data.frame(step = "hydrate_cached_peak_table", note = "The existing cached table was used because modern XCMS was not forced; detected and filled assays are identical and must be replaced by --force for strict raw-detection QC."), file.path(dirs$peak, "process_history.tsv")); invisible(TRUE)
  }
  load_peak <- function() {
    detected_path <- file.path(dirs$peak, "abundance_matrix_detected.csv")
    filled_path <- file.path(dirs$peak, "abundance_matrix_filled.csv")
    annotation_path <- file.path(dirs$peak, "feature_annotation.csv")
    if (!file.exists(detected_path) || !file.exists(filled_path) || !file.exists(annotation_path)) stop("Missing XCMS detected/filled matrices. Run the XCMS step first.", call. = FALSE)
    read_matrix <- function(path) {
      d <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE); id <- names(d)[1]; m <- as.matrix(d[, -1, drop = FALSE]); suppressWarnings(storage.mode(m) <- "numeric"); rownames(m) <- d[[id]]; m
    }
    source_path <- file.path(dirs$peak, "peak_source.csv")
    source <- if (file.exists(source_path)) readr::read_csv(source_path, show_col_types = FALSE) else data.frame(source = "modern_xcms")
    list(detected = read_matrix(detected_path), filled = read_matrix(filled_path), annotation = readr::read_csv(annotation_path, show_col_types = FALSE), source = source)
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
    manifest_path <- file.path(dirs$audit, "mzml_manifest.csv")
    if (file.exists(manifest_path)) {
      mz <- readr::read_csv(manifest_path, show_col_types = FALSE) |> filter(!is.na(run_order)) |> arrange(run_order) |> distinct(sample, .keep_all = TRUE) |> select(sample, run_order, acquired_at)
      si <- left_join(si, mz, by = "sample")
    }
    if (!"run_order" %in% names(si)) si$run_order <- seq_len(nrow(si))
    si$run_order <- num(si$run_order, seq_len(nrow(si)))
    si$group <- as.character(si$group)
    det <- detection(peak$detected); fill <- peak$filled
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
    feature$pass_detection <- feature[[paste0("detection_rate_", case_group)]] >= detection_cutoff | feature[[paste0("detection_rate_", control_group)]] >= detection_cutoff
    feature$pass_qc_filter <- if (qc_filter == "rsd") is.finite(feature$qc_rsd) & feature$qc_rsd <= qc_rsd_cutoff else TRUE
    feature$keep <- feature$pass_detection & feature$pass_qc_filter
    feature$filter_reason <- ifelse(feature$keep, "keep", ifelse(!feature$pass_detection, "low_group_detection", "qc_rsd"))
    write_csv2(sample_qc, file.path(dirs$qc, "sample_qc_summary.csv")); write_csv2(feature, file.path(dirs$qc, "feature_filter_decision.csv")); write_csv2(qc_cv, file.path(dirs$qc, "qc_rsd_summary.csv"))
    keep_ids <- feature$feature_id[feature$keep]; fill <- fill[keep_ids, , drop = FALSE]
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
    } else data.frame(sample = colnames(fill), ms1_spectra = NA_integer_, ms2_spectra = NA_integer_, status = ifelse(grepl("legacy", source_label, ignore.case = TRUE), "not_available_from_cached_peak_table", "not_collected_in_peak_export"), stringsAsFactors = FALSE)
    rt_method <- as.character(cfg_get("xcms.retention_correction", "obiwarp"))
    rt_summary <- data.frame(polarity = c("POS", "NEG"), retention_correction = rt_method, status = ifelse(grepl("legacy", source_label, ignore.case = TRUE), "unknown_cached", "requested_in_xcms"), stringsAsFactors = FALSE)
    bio_sd <- if (length(bidx)) apply(log2m[, bidx, drop = FALSE], 1, stats::sd, na.rm = TRUE) else rep(NA_real_, nrow(log2m)); qc_sd <- if (length(qidx)) apply(log2m[, qidx, drop = FALSE], 1, stats::sd, na.rm = TRUE) else rep(NA_real_, nrow(log2m)); dratio <- qc_sd / bio_sd
    write_csv2(factors |> left_join(sample_qc, by = "sample"), file.path(dirs$qc, "normalization_factors.csv")); write_csv2(signal_summary, file.path(dirs$qc, "sample_signal_summary.csv")); write_csv2(rla_summary, file.path(dirs$qc, "rla_summary.csv")); write_csv2(spectra_summary, file.path(dirs$qc, "spectra_summary.csv")); write_csv2(rt_summary, file.path(dirs$qc, "rt_correction_summary.csv")); write_matrix(imputed, file.path(dirs$qc, "matrix_imputed.csv")); write_matrix(norm, file.path(dirs$qc, "matrix_normalized.csv")); write_matrix(log2m, file.path(dirs$qc, "matrix_log2.csv")); write_matrix(scaled, file.path(dirs$qc, "matrix_scaled.csv")); write_csv2(data.frame(feature_id = rownames(log2m), d_ratio = dratio), file.path(dirs$qc, "d_ratio_summary.csv")); write_tsv2(si, file.path(dirs$qc, "sample_info_used.tsv"))
    saveRDS(list(raw_detected = peak$detected, filled = peak$filled, filtered = fill, imputed = imputed, norm = norm, log2 = log2m, scaled = scaled, annotation = peak$annotation, feature_qc = feature, qc_cv = qc_cv, d_ratio = dratio, sample_qc = sample_qc, signal_summary = signal_summary, rla_summary = rla_summary, spectra_summary = spectra_summary, rt_summary = rt_summary, sample_info = si, source = peak$source), file.path(dirs$qc_objects, "processed_data.rds"))
    invisible(TRUE)
  }

  plot_qc <- function() {
    dirs$plot <- file.path(result_dir, "03_qc", "plots")
    p <- readRDS(file.path(dirs$qc_objects, "processed_data.rds")); sq <- p$sample_qc
    p1 <- ggplot(sq |> arrange(run_order), aes(run_order, total_signal, color = group, label = sample)) + geom_line(aes(group = 1), color = "grey75", linewidth = .35) + geom_point(size = 2.5) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30, show.legend = FALSE) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = "Acquisition order", y = "Total peak area", color = NULL)
    save_plot(p1, file.path(dirs$plot, "sample_run_order_total_signal"), 7.2, 4.8)
    long <- as.data.frame(p$log2, check.names = FALSE) |> tibble::rownames_to_column("feature_id") |> pivot_longer(-feature_id, names_to = "sample", values_to = "log2_intensity") |> left_join(p$sample_info, by = "sample")
    p2 <- ggplot(long, aes(sample, log2_intensity, fill = group)) + geom_boxplot(outlier.size = .2, linewidth = .2) + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = palette) + labs(x = NULL, y = "log2 normalized intensity", fill = NULL)
    save_plot(p2, file.path(dirs$plot, "sample_log2_intensity_boxplot"), 8, 4.8)
    p3 <- ggplot(p$qc_cv, aes(qc_rsd)) + geom_histogram(bins = 40, fill = "#8DE5A1", color = "white") + geom_vline(xintercept = qc_rsd_cutoff, linetype = 2, color = "#5C6B73") + theme_classic(base_size = 10) + labs(x = "QC RSD", y = "Feature count")
    save_plot(p3, file.path(dirs$plot, "qc_rsd_distribution"), 6, 4)
    write_csv2(sq, file.path(dirs$qc, "sample_qc_summary.csv"))
    rla <- p$rla_summary; if (!is.null(rla) && nrow(rla)) { prla <- ggplot(rla, aes(run_order, median_rla, color = group, label = sample)) + geom_hline(yintercept = 0, linetype = 2, color = "grey70") + geom_line(aes(group = 1), color = "grey75", linewidth = .35) + geom_point(size = 2.5) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30, show.legend = FALSE) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = "Acquisition order", y = "Median RLA", color = NULL); save_plot(prla, file.path(dirs$plot, "sample_rla_summary"), 7.2, 4.8) }
    for (pol in c("POS", "NEG")) {
      ids <- p$annotation$feature_id[p$annotation$polarity == pol]; ids <- intersect(ids, rownames(p$scaled)); if (length(ids) < 2) next
      x <- p$scaled[ids, , drop = FALSE]; keep <- apply(x, 1, stats::sd) > 0; x <- x[keep, , drop = FALSE]; if (nrow(x) < 2) next
      pc <- stats::prcomp(t(x), center = FALSE, scale. = FALSE); vr <- pc$sdev^2 / sum(pc$sdev^2); sc <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2]) |> left_join(p$sample_info, by = "sample"); write_csv2(sc, file.path(dirs$qc, paste0("pca_scores_", pol, ".csv")))
      pp <- ggplot(sc, aes(PC1, PC2, color = group, label = sample)) + geom_point(size = 2.7) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30, show.legend = FALSE) + theme_classic(base_size = 10) + scale_color_manual(values = palette) + labs(x = sprintf("PC1 (%.1f%%)", vr[1] * 100), y = sprintf("PC2 (%.1f%%)", vr[2] * 100), color = NULL, title = paste0(pol, " polarity PCA"))
      save_plot(pp, file.path(dirs$plot, paste0("pca_", pol)), 6.2, 5)
      cor_mat <- cor(p$log2[ids, , drop = FALSE], method = "spearman", use = "pairwise.complete.obs"); ann <- data.frame(group = p$sample_info$group); rownames(ann) <- p$sample_info$sample
      write_matrix(cor_mat, file.path(dirs$qc, paste0("sample_spearman_correlation_", pol, ".csv")), "sample")
      save_heatmap(cor_mat, file.path(dirs$plot, paste0("sample_spearman_correlation_", pol)), 6, scale_kind = "correlation",
        annotation_col = ann, annotation_row = ann, annotation_colors = list(group = palette), border_color = NA)
    }
    invisible(TRUE)
  }

  load_annotation_table <- function() {
    path <- resolve_path(cfg_get("fallback_identification_table", "")); if (!nzchar(path) || !file.exists(path)) return(data.frame())
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
  annotate <- function() {
    peak <- load_peak(); ft <- peak$annotation |> distinct(feature_id, .keep_all = TRUE); id <- load_annotation_table()
    if (!nrow(id)) { write_csv2(ft, file.path(dirs$annotation, "annotation_full.csv")); writeLines("No MetDNA identification table was available.", file.path(dirs$annotation_logs, "annotation_status.txt")); return(invisible(TRUE)) }
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
      cand_pool <- which(is.finite(id_mz) & is.finite(id_rt) & (is.na(id_pol) | id_pol == pol)); if (!length(cand_pool)) cand_pool <- which(is.finite(id_mz) & is.finite(id_rt)); if (!length(cand_pool)) next
      ord_mz <- cand_pool[order(id_mz[cand_pool])]; mz_sorted <- id_mz[ord_mz]
      for (i in ii) {
        mz_window <- out_mz[i] * mz_tol / 1e6; lo <- findInterval(out_mz[i] - mz_window, mz_sorted); hi <- findInterval(out_mz[i] + mz_window, mz_sorted); if (hi < lo + 1L) next
        cand <- ord_mz[seq.int(max(1L, lo), min(length(ord_mz), hi))]; mzerr <- abs(id_mz[cand] - out_mz[i]) / out_mz[i] * 1e6; rterr <- abs(id_rt[cand] - out_rt[i]); ok <- is.finite(mzerr) & is.finite(rterr) & mzerr <= mz_tol & rterr <= rt_tol; if (!any(ok)) next
        cand <- cand[ok]; score <- mzerr[ok] / mz_tol + rterr[ok] / rt_tol; oo <- order(score); j <- cand[oo[1]]; match_idx[i] <- j; out$annotation_match_type[i] <- "nearest_mz_rt"; out$annotation_mz_error_ppm[i] <- mzerr[ok][oo[1]]; out$annotation_rt_error_seconds[i] <- rterr[ok][oo[1]]; if (length(oo) > 1L && (score[oo[2]] - score[oo[1]]) <= ambig_delta) out$annotation_ambiguous[i] <- TRUE
      }
    }
    for (f in fields) out[[f]] <- id[[f]][match_idx]
    out$annotated <- !is.na(out$name) & nzchar(as.character(out$name)); write_csv2(out, file.path(dirs$annotation, "annotation_full.csv")); write_csv2(out |> mutate(annotation_class = case_when(!annotated ~ "unannotated", !is.na(confidence_level) ~ as.character(confidence_level), TRUE ~ "annotated_no_confidence")) |> count(annotation_class, name = "n_features"), file.path(dirs$annotation, "annotation_summary.csv")); writeLines("MetDNA annotation was matched by exact feature ID first, then polarity/mz/RT tolerance; no correlation-based direction flipping was used.", file.path(dirs$annotation_logs, "annotation_status.txt")); invisible(TRUE)
  }

  draw_differential <- function(res, p) {
    dirs$plot <- file.path(result_dir, "04_statistics", "plots")
    si <- p$sample_info
    si_bio <- si[si$group %in% c(case_group, control_group), , drop = FALSE]
    labels <- head(res |> filter(classification != "not_sig") |> arrange(limma_p), 12); volcano <- ggplot(res, aes(log2FC, -log10(pmax(limma_p, .Machine$double.xmin)), color = direction, shape = classification)) + geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey65") + geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), linetype = 2, color = "grey75") + geom_point(alpha = .75, size = 1.4) + ggrepel::geom_text_repel(data = labels, aes(label = ifelse(is.na(name), feature_id, name)), size = 2.7, max.overlaps = 20, show.legend = FALSE) + scale_color_manual(values = direction_palette[c(direction_control, direction_case, "not_sig")], breaks = c(direction_control, direction_case, "not_sig"), drop = FALSE) + scale_shape_manual(values = c(robust = 16, exploratory = 17, not_sig = 1), labels = c(robust = "BH + FC supported", exploratory = "Nominal-P candidate", not_sig = "Not selected"), drop = FALSE) + theme_classic(base_size = 10) + labs(x = paste0("log2 fold change (", case_group, " - ", control_group, ")"), y = "-log10 limma P", color = "Direction", shape = "Evidence")
    save_plot(volcano, file.path(dirs$plot, paste0(comparison_prefix, "_volcano")), 7, 5.5)
    top <- head(res |> filter(classification != "not_sig") |> arrange(desc(candidate_score)) |> pull(feature_id), 20); if (length(top)) { long <- as.data.frame(p$log2[top, , drop = FALSE], check.names = FALSE) |> tibble::rownames_to_column("feature_id") |> pivot_longer(-feature_id, names_to = "sample", values_to = "log2_intensity") |> left_join(si, by = "sample") |> filter(group != qc_group); long$group <- factor(long$group, levels = c(control_group, case_group)); bx <- ggplot(long, aes(group, log2_intensity, fill = group)) + geom_boxplot(outlier.shape = NA, linewidth = .2) + geom_point(position = position_jitter(width = .1, seed = 42), size = 1, alpha = .85) + facet_wrap(~feature_id, scales = "free_y", ncol = 4) + theme_classic(base_size = 8) + scale_fill_manual(values = palette) + labs(x = NULL, y = "log2 normalized intensity", fill = NULL); save_plot(bx, file.path(dirs$plot, paste0(comparison_prefix, "_top_feature_boxplots")), 10, 8); hm <- p$scaled[top, si_bio$sample, drop = FALSE]; ann_col <- data.frame(group = si_bio$group); rownames(ann_col) <- si_bio$sample; save_heatmap(hm, file.path(dirs$plot, paste0(comparison_prefix, "_top_feature_heatmap")), 9,
      annotation_col = ann_col, annotation_colors = list(group = palette), border_color = NA, show_colnames = TRUE) }
    invisible(TRUE)
  }
  draw_pathway <- function(vals) {
    dirs$plot <- file.path(result_dir, "06_pathway", "plots")
    top <- head(arrange(vals, p), 20)
    if (nrow(top)) {
      pp <- ggplot(top, aes(-log10(pmax(p, .Machine$double.xmin)), reorder(pathway, -log10(pmax(p, .Machine$double.xmin))), fill = log2FC)) +
        geom_col() + scale_fill_gradient2(low = pfos_palette("diverging")[["low"]], mid = pfos_palette("diverging")[["mid"]], high = pfos_palette("diverging")[["high"]], midpoint = 0) +
        labs(x = "-log10 P", y = NULL, fill = paste0("log2FC\n", case_group, "-", control_group))
      save_plot(pp, file.path(dirs$plot, "pathway_quantitative_statistics"), 8, 6)
    }
  }
  differential <- function() {
    p <- readRDS(file.path(dirs$qc_objects, "processed_data.rds")); ann <- readr::read_csv(file.path(dirs$annotation, "annotation_full.csv"), show_col_types = FALSE); si <- p$sample_info; bio <- si$group %in% c(case_group, control_group); si_bio <- si[bio, , drop = FALSE]; mat <- p$log2[, si_bio$sample, drop = FALSE]; grp <- factor(si_bio$group, levels = c(control_group, case_group)); design <- model.matrix(~0 + grp); colnames(design) <- make.names(levels(grp)); if (run_order_adjustment == "on") { design <- cbind(design, run_order_scaled = as.numeric(scale(si_bio$run_order))) }
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
      cpath <- resolve_path(cfg_get("fallback_statistics_table", "")); if (file.exists(cpath)) { cdat <- readr::read_csv(cpath, show_col_types = FALSE); cid <- intersect(c("name", "metabolite", "peak_name"), names(cdat))[1]; if (!is.na(cid) && all(c("p_values", "fold_changes") %in% names(cdat))) { cs <- cdat |> transmute(feature_id = as.character(.data[[cid]]), company_p_value = num(p_values, NA_real_), company_fold_change_raw = num(fold_changes, NA_real_)); res <- left_join(res, cs, by = "feature_id"); res$company_direction_convention <- as.character(cfg_get("company_comparison_direction", "unverified")) } }
    }
    res <- arrange(res, limma_p); robust_out <- res |> filter(classification == "robust"); exploratory_out <- res |> filter(classification == "exploratory"); write_csv2(res, file.path(dirs$stat, paste0(comparison_prefix, "_differential_features.csv"))); write_csv2(robust_out, file.path(dirs$stat, paste0(comparison_prefix, "_robust_differential_features.csv"))); write_csv2(exploratory_out, file.path(dirs$stat, paste0(comparison_prefix, "_exploratory_candidates.csv"))); write_csv2(res |> count(classification, direction, name = "n_features"), file.path(dirs$stat, paste0(comparison_prefix, "_summary.csv")))
    company_cols <- intersect(c("feature_id", "name", "polarity", "log2FC", "fold_change", "limma_p", "fdr", "classification", "direction", "company_p_value", "company_fold_change_raw", "company_direction_convention"), names(res))
    company_match <- if (company_comparison == "on" && "company_p_value" %in% names(res)) res |> filter(!is.na(company_p_value)) |> select(all_of(company_cols)) else res[0, company_cols, drop = FALSE]
    write_csv2(company_match, file.path(dirs$stat, paste0(comparison_prefix, "_company_comparison.csv")))
    writeLines(if (company_comparison == "on") paste0("Company statistics were joined by the supplied identifier without direction reversal. Matched rows: ", nrow(company_match), ".") else "Company comparison disabled by configuration.", file.path(dirs$stat_logs, "company_comparison_status.txt"))
    draw_differential(res, p)
    invisible(TRUE)
  }

  pathway <- function() {
    qpath <- resolve_path(cfg_get("fallback_pathway_table", "")); epath <- resolve_path(cfg_get("pathway_enrichment_file", "")); if (nzchar(epath) && file.exists(epath)) write_csv2(readr::read_csv(epath, show_col_types = FALSE), file.path(dirs$pathway, "external_metdna_pathway_enrichment.csv"))
    if (!nzchar(qpath) || !file.exists(qpath)) { writeLines("No external pathway quantitative table was available.", file.path(dirs$pathway_logs, "pathway_status.txt")); return(invisible(TRUE)) }
    dat <- readr::read_csv(qpath, show_col_types = FALSE); pc <- intersect(c("pathway_name", "pathway"), names(dat))[1]; samples <- intersect(names(dat), c(paste0(case_group, 1:5), paste0(control_group, 1:5))); if (is.na(pc) || length(samples) < 4) return(invisible(TRUE)); vals <- lapply(seq_len(nrow(dat)), function(i) { x <- num(dat[i, paste0(case_group, 1:5)], NA_real_); y <- num(dat[i, paste0(control_group, 1:5)], NA_real_); z <- c(x, y); shift <- min(z[z > 0], na.rm = TRUE) / 2; data.frame(pathway = dat[[pc]][i], log2FC = mean(log2(x + shift), na.rm = TRUE) - mean(log2(y + shift), na.rm = TRUE), p = tryCatch(stats::t.test(log2(x + shift), log2(y + shift))$p.value, error = function(e) NA_real_), stringsAsFactors = FALSE) }) |> bind_rows(); vals$fdr <- p.adjust(vals$p, "BH"); write_csv2(arrange(vals, p), file.path(dirs$pathway, "pathway_quantitative_statistics.csv")); draw_pathway(vals)
    membership <- resolve_path(cfg_get("pathway_membership_file", "")); if (!nzchar(membership) || !file.exists(membership)) writeLines("KEGG pathway membership file not supplied; external MetDNA enrichment is reported separately and no new ORA is claimed.", file.path(dirs$pathway_logs, "pathway_ora_status.txt"))
    invisible(TRUE)
  }

  report <- function() {
    proc <- readRDS(file.path(dirs$qc_objects, "processed_data.rds"))
    res <- readr::read_csv(file.path(dirs$stat, paste0(comparison_prefix, "_differential_features.csv")), show_col_types = FALSE)
    inv <- if (file.exists(file.path(dirs$audit, "inventory_summary.csv"))) readr::read_csv(file.path(dirs$audit, "inventory_summary.csv"), show_col_types = FALSE) else data.frame()
    spectra <- if (file.exists(file.path(dirs$audit, "mzml_spectra_summary.csv"))) readr::read_csv(file.path(dirs$audit, "mzml_spectra_summary.csv"), show_col_types = FALSE) else data.frame()
    ms1_total <- if (nrow(spectra)) sum(spectra$ms1_spectra, na.rm = TRUE) else NA_real_
    ms2_total <- if (nrow(spectra)) sum(spectra$ms2_spectra, na.rm = TRUE) else NA_real_
    robust <- sum(res$classification == "robust", na.rm = TRUE)
    exploratory <- sum(res$classification == "exploratory", na.rm = TRUE)
    annotated <- sum(!is.na(res$name) & nzchar(as.character(res$name)), na.rm = TRUE)
    qc_bad <- sum(proc$feature_qc$qc_rsd > qc_rsd_cutoff, na.rm = TRUE)
    order_text <- paste(proc$sample_info |> arrange(run_order) |> pull(sample), collapse = " → ")
    peak_source_text <- if (nrow(proc$source)) as.character(proc$source$source[[1]]) else "modern_xcms"
    peak_detection_note <- if (grepl("legacy", peak_source_text, ignore.case = TRUE)) "- 注意：当前缓存峰表的填峰前/后矩阵相同，不能把填峰值检出率解释为原始检测率；请用 `--force` 重建严格 XCMS 矩阵。" else "- 检出率基于填峰前矩阵计算，随后才进行填峰值缺失处理。"
    sample_tab <- paste(capture.output(knitr::kable(proc$sample_info |> select(any_of(c("sample", "group", "run_order", "acquired_at"))))), collapse = "\n")
    top <- res |> filter(classification != "not_sig") |> arrange(desc(candidate_score)) |> head(15)
    top_tab <- if (nrow(top)) paste(capture.output(knitr::kable(top |> select(any_of(c("feature_id", "name", "polarity", "mz", "rt", "log2FC", "fold_change", "limma_p", "fdr", "classification", "direction"))))), collapse = "\n") else "无符合当前策略的候选特征。"
    report_path <- file.path(dirs$audit_logs, "PFOS_metabolomics_full_report.md")
    lines <- c("# PFOS vs WT 小鼠肝脏代谢组分析报告", "", paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "", "## 1. 研究设计与流程", "", "- 研究对象：F1 小鼠肝脏；比较方向固定为 `PFOS - WT`。", paste0("- QC 表示质量控制样本；本批次 QC 来源尚未确认，因此当前默认按 `qc_role=", qc_role, "` 处理。"), paste0("- 有效配置：`", file.path(dirs$audit_logs, "effective_config.yml"), "`。"), paste0("- 原始进样顺序：", order_text), "", sample_tab, "", "## 2. 原始数据与 XCMS", "", paste0("- WIFF 文件：", if (nrow(inv)) inv$value[inv$metric == "wiff_files"] else "NA", "；mzML 文件：", if (nrow(inv)) inv$value[inv$metric == "mzml_files"] else "NA", "。"), paste0("- MS1 谱图数：", ifelse(is.finite(ms1_total), format(ms1_total, scientific = FALSE), "NA"), "；MS2 谱图数：", ifelse(is.finite(ms2_total), format(ms2_total, scientific = FALSE), "NA"), "。"), paste0("- 峰表特征数：", nrow(proc$log2), "；填峰前检出矩阵和填峰后矩阵分别保存。"), "- POS/NEG 分开处理，保留 MS2 注释证据。", paste0("- 峰表来源：", peak_source_text, "。"), peak_detection_note, "", "## 3. QC、过滤与归一化", "", paste0("- 生物组检出率阈值：", detection_cutoff, "；QC RSD 阈值：", qc_rsd_cutoff, "。"), paste0("- QC 硬过滤：", qc_filter, "；漂移校正：", drift_correction, "；归一化：", normalization, "。"), paste0("- 过滤后特征数：", nrow(proc$log2), "；QC RSD 超阈值特征数：", qc_bad, "。"), "- 输出 TIC/BPC 峰面积代理、峰数、缺失率、RLA、D-ratio、真实进样顺序、分极性 PCA 和相关性诊断。", "", "## 4. PFOS vs WT 差异结果", "", paste0("- 稳健差异特征：", robust, "（FDR < ", fdr_cutoff, " 且 |FC| ≥ ", fc_cutoff, "）。"), paste0("- 探索候选特征：", exploratory, "（名义 P < ", p_cutoff, " 且 |FC| ≥ ", fc_cutoff, "，不表述为 FDR 显著）。"), paste0("- 有注释名称的特征：", annotated, "；未注释条目仅称代谢特征。"), "", top_tab, "", "## 5. 注释与通路", "", "- MetDNA 先按特征 ID 匹配，再按极性、m/z 和 RT 容差匹配；报告匹配误差和歧义。", "- 不使用相关系数自动反转公司结果方向；公司统计仅作独立对照。", "- 通路定量和外部 MetDNA 富集结果分别保存于 `result/06_pathway/`；未提供本地 KEGG 成员表时，不声称重新完成 ORA。", "", "## 6. 结果解释边界", "", "- 小样本条件下，名义候选用于后续验证，不等同于多重校正显著代谢物。", "- 注释 level 2/3.x 代表相应证据等级，不等同于 level 1 标准品确认。", "- 异常样本只做标记，不由 PCA 自动删除。", "", "## 7. 输出与复现", "", "- QC：`result/03_qc/`；统计：`result/04_statistics/`；注释：`result/05_annotation/`；通路：`result/06_pathway/`。", "- 审计：`result/01_audit/`，包括有效配置、软件版本、R sessionInfo、输入文件哈希和真实进样顺序。")
    writeLines(lines, report_path)
    figures <- list.files(result_dir, "[.]png$", recursive = TRUE)
    page <- c('<!doctype html><html lang="zh"><meta charset="utf-8"><style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:30px}img{max-width:100%}</style>',
      paste0("<p>", lines, "</p>"), paste0('<figure><h3>', basename(figures), '</h3><img src="', figures, '"><figcaption>原分析完整图件；解释边界与统计见上文。</figcaption></figure>'), "</html>")
    writeLines(page, file.path(result_dir, "index.html"))
    invisible(report_path)
  }

  if (plots_only) {
    needed <- c(file.path(dirs$qc_objects, "processed_data.rds"),
      file.path(dirs$stat, paste0(comparison_prefix, "_differential_features.csv")))
    pathway_input <- file.path(dirs$pathway, "pathway_quantitative_statistics.csv")
    if (file.exists(pathway_input)) needed <- c(needed, pathway_input)
    missing <- needed[!file.exists(needed)]
    if (length(missing)) stop("缺少冻结的绘图来源：", paste(missing, collapse = ", "), call. = FALSE)
    p <- readRDS(needed[1])
    if (anyDuplicated(p$sample_info$sample) || !identical(colnames(p$log2), as.character(p$sample_info$sample))) stop("绘图矩阵与样本身份不一致")
    set.seed(20261012)
    plot_qc()
    draw_differential(readr::read_csv(needed[2], show_col_types = FALSE), p)
    if (length(needed) > 2) draw_pathway(readr::read_csv(needed[3], show_col_types = FALSE))
    message("Completed plots-only; 未调用统计检验或上游流程")
    return(invisible(TRUE))
  }
  module <- opts$module %||% "all"
  if (isTRUE(opts$dry_run)) { message("dry-run module=", module, " result_dir=", result_dir); return(invisible(TRUE)) }
  if (module %in% c("audit", "inventory", "all")) { inv <- inventory(); tool_audit(); effective_config(); manifest(c(inv$raw$raw_file, inv$mzml$file, sample_info_file, config_file)) }
  if (module %in% c("convert", "upstream", "all")) { inv <- if (exists("inv")) inv else inventory(); convert_wiff(inv) }
  if (module %in% c("xcms", "upstream", "all")) {
    if (!file.exists(file.path(dirs$peak, "abundance_matrix_detected.csv")) || !file.exists(file.path(dirs$peak, "abundance_matrix_filled.csv"))) {
      if (!isTRUE(opts$force) && bool(cfg_get("xcms.skip_existing", TRUE), TRUE) && file.exists(file.path(peak_dir, "xcms_peak_table_combined.csv"))) hydrate_cached_peak_table() else run_xcms_modern()
    } else if (isTRUE(opts$force)) run_xcms_modern()
  }
  if (module %in% c("preprocess", "qc", "downstream", "all")) preprocess()
  if (module %in% c("plots", "qc", "downstream", "all")) plot_qc()
  if (module %in% c("annotate", "annotation", "downstream", "all")) annotate()
  if (module %in% c("stats", "differential", "downstream", "all")) differential()
  if (module %in% c("pathway", "downstream", "all")) pathway()
  if (module %in% c("report", "downstream", "all")) report()
  message("Completed module: ", module)
  invisible(TRUE)
}
