# Shared publication figure contract for mouse WGBS projects.
# Plot layers and statistical results remain owned by each analysis module;
# this file only standardizes typography, physical size, and paired export.

wgbs_figure_width_mm <- 183
wgbs_figure_min_height_mm <- 70
wgbs_figure_max_height_mm <- 240
wgbs_figure_dpi <- 300
wgbs_figure_font_family <- "Liberation Sans"
wgbs_figure_base_size <- 8
wgbs_figure_tick_size <- 7

wgbs_figure_dimensions <- function(width = 8, height = 6) {
  width <- suppressWarnings(as.numeric(width)[1])
  height <- suppressWarnings(as.numeric(height)[1])
  ratio <- height / width
  if(!is.finite(ratio) || ratio <= 0) ratio <- 0.75
  height_mm <- min(
    wgbs_figure_max_height_mm,
    max(wgbs_figure_min_height_mm, wgbs_figure_width_mm * ratio)
  )
  # Cairo PDF records integer PostScript points. Round explicitly so every
  # file has the same physical page box.
  pdf_width_mm <- round(wgbs_figure_width_mm / 25.4 * 72) * 25.4 / 72
  pdf_height_mm <- round(height_mm / 25.4 * 72) * 25.4 / 72
  list(
    width_mm = wgbs_figure_width_mm,
    height_mm = height_mm,
    pdf_width_mm = pdf_width_mm,
    pdf_height_mm = pdf_height_mm,
    width_in = pdf_width_mm / 25.4,
    height_in = pdf_height_mm / 25.4
  )
}

wgbs_assert_figure_runtime <- function() {
  if(!isTRUE(capabilities("cairo"))) {
    stop("WGBS publication export requires Cairo support in R.", call. = FALSE)
  }
  invisible(TRUE)
}

wgbs_png_path <- function(pdf_file) {
  if(grepl("\\.pdf$", pdf_file, ignore.case = TRUE)) {
    sub("\\.pdf$", ".png", pdf_file, ignore.case = TRUE)
  } else {
    paste0(pdf_file, ".png")
  }
}

wgbs_finalize_rgb_png <- function(path, dpi = wgbs_figure_dpi) {
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

wgbs_normalize_png_widths <- function(result_dir) {
  if(!dir.exists(result_dir)) return(invisible(character()))
  if(!requireNamespace("magick", quietly = TRUE)) {
    stop("The magick package is required to normalize legacy WGBS PNG previews.", call. = FALSE)
  }

  pdf_files <- list.files(result_dir, pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE)
  for(pdf_file in pdf_files) {
    pdf_stem <- sub("[.]pdf$", "", pdf_file, ignore.case = TRUE)
    has_preview <- file.exists(paste0(pdf_stem, ".png")) ||
      length(Sys.glob(paste0(pdf_stem, "_page_*.png"))) > 0
    if(!has_preview) wgbs_normalize_existing_pdf(pdf_file)
  }

  target_width <- as.integer(round(wgbs_figure_width_mm / 25.4 * wgbs_figure_dpi))
  file_bin <- Sys.which("file")
  if(!nzchar(file_bin)) stop("The file utility is required to inspect WGBS PNG dimensions.", call. = FALSE)
  png_files <- list.files(result_dir, pattern = "[.]png$", recursive = TRUE, full.names = TRUE)
  for(path in png_files) {
    file_info <- suppressWarnings(system2(file_bin, c("-b", shQuote(path)), stdout = TRUE))
    width_match <- regexec("PNG image data, ([0-9]+) x", file_info)
    width_parts <- regmatches(file_info, width_match)[[1]]
    if(length(width_parts) < 2L) stop("Unable to read PNG dimensions: ", path, call. = FALSE)
    if(identical(as.integer(width_parts[[2]]), target_width)) next
    image <- magick::image_read(path)
    image <- magick::image_scale(image, paste0(target_width, "x"))
    image <- magick::image_background(image, "white", flatten = TRUE)
    magick::image_write(image, path = path, format = "png", depth = 8)
    magick::image_destroy(image)
    wgbs_finalize_rgb_png(path)
  }
  invisible(png_files)
}

wgbs_publication_typography <- function(coordinate_system = c("axes", "none")) {
  coordinate_system <- match.arg(coordinate_system)
  out <- ggplot2::theme(
    text = ggplot2::element_text(
      family = wgbs_figure_font_family,
      size = wgbs_figure_base_size,
      color = "#333333"
    ),
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = wgbs_figure_base_size, color = "#333333"),
    axis.text = ggplot2::element_text(size = wgbs_figure_tick_size, color = "#333333"),
    strip.text = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
    legend.title = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
    legend.text = ggplot2::element_text(size = wgbs_figure_tick_size, color = "#333333"),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
  )
  if(identical(coordinate_system, "none")) {
    out <- out + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  }
  out
}

wgbs_text_layers <- function(plot) {
  if(inherits(plot,'patchwork'))plot$patches$plots<-lapply(plot$patches$plots,wgbs_text_layers)
  if(inherits(plot,'ggplot'))for(i in seq_along(plot$layers)) {
    layer<-plot$layers[[i]]
    if(any(class(layer$geom)%in%c('GeomText','GeomTextRepel','GeomLabel','GeomLabelRepel')) && is.null(layer$mapping$size) && is.null(plot$mapping$size)) {
      layer$aes_params$size<-7/ggplot2::.pt;layer$aes_params$family<-'Liberation Sans';plot$layers[[i]]<-layer
    }
  }
  plot
}

wgbs_standardize_plot <- function(plot, coordinate_system = c("axes", "none")) {
  plot <- wgbs_text_layers(plot)
  coordinate_system <- match.arg(coordinate_system)
  typography <- wgbs_publication_typography(coordinate_system)
  if(inherits(plot, "patchwork")) return(plot & typography)
  if(inherits(plot, "gglist")) {
    plot[] <- lapply(plot, function(panel) {
      if(inherits(panel, "ggplot")) panel + typography else panel
    })
    return(plot)
  }
  plot + typography
}

wgbs_publication_theme <- function(legend = "right") {
  ggplot2::theme_classic(
    base_size = wgbs_figure_base_size,
    base_family = wgbs_figure_font_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = wgbs_figure_font_family, size = wgbs_figure_base_size, color = "#333333"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = wgbs_figure_base_size, color = "#333333"),
      axis.text = ggplot2::element_text(size = wgbs_figure_tick_size, color = "#333333"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#333333", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "#333333", linewidth = 0.35),
      strip.background = ggplot2::element_rect(fill = "#F3F4F5", color = NA),
      strip.text = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
      legend.position = legend,
      legend.title = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
      legend.text = ggplot2::element_text(size = wgbs_figure_tick_size, color = "#333333"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
    )
}

wgbs_publication_void_theme <- function(legend = "right") {
  ggplot2::theme_void(
    base_size = wgbs_figure_base_size,
    base_family = wgbs_figure_font_family
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = wgbs_figure_font_family, size = wgbs_figure_base_size, color = "#333333"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F3F4F5", color = NA),
      strip.text = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
      legend.position = legend,
      legend.title = ggplot2::element_text(size = wgbs_figure_base_size, face = "bold", color = "#333333"),
      legend.text = ggplot2::element_text(size = wgbs_figure_tick_size, color = "#333333"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(4, 6, 4, 4, unit = "mm")
    )
}

wgbs_ggsave <- function(filename, plot = ggplot2::last_plot(), device = grDevices::cairo_pdf,
                         path = NULL, scale = 1, width = 8, height = 6,
                         units = "in", dpi = wgbs_figure_dpi, limitsize = TRUE,
                         bg = "white", create.dir = FALSE, ..., coordinate_system = "axes") {
  if(!is.null(path)) filename <- file.path(path, filename)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  dims <- wgbs_figure_dimensions(width, height)
  plot <- wgbs_standardize_plot(plot, coordinate_system = coordinate_system)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = dims$pdf_width_mm,
    height = dims$pdf_height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    family = wgbs_figure_font_family,
    bg = "white",
    limitsize = FALSE,
    ...
  )
  if(grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    png_file <- wgbs_png_path(filename)
    grDevices::png(
      filename = png_file,
      width = dims$width_mm,
      height = dims$height_mm,
      units = "mm",
      res = wgbs_figure_dpi,
      pointsize = wgbs_figure_base_size,
      type = "cairo-png",
      bg = "white"
    )
    tryCatch(print(plot), finally = grDevices::dev.off())
    wgbs_finalize_rgb_png(png_file)
  }
  invisible(filename)
}

.wgbs_pdf_devices <- new.env(parent = emptyenv())
.wgbs_pdf_devices$files <- list()

wgbs_cairo_pdf <- function(filename, width = 8, height = 6, ...) {
  wgbs_assert_figure_runtime()
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  dims <- wgbs_figure_dimensions(width, height)
  dots <- list(...)
  args <- list(
    filename = filename,
    width = dims$width_in,
    height = dims$height_in,
    pointsize = wgbs_figure_base_size,
    family = wgbs_figure_font_family,
    bg = "white",
    onefile = TRUE
  )
  if(length(dots) > 0) args[names(dots)] <- dots
  do.call(grDevices::cairo_pdf, args)
  graphics::par(cex.axis=7/8,cex.main=1,cex.lab=1,family=wgbs_figure_font_family)
  .wgbs_pdf_devices$files[[as.character(grDevices::dev.cur())]] <- list(
    file = filename,
    width_px = as.integer(round(dims$width_mm / 25.4 * wgbs_figure_dpi)),
    height_px = as.integer(round(dims$height_mm / 25.4 * wgbs_figure_dpi))
  )
  invisible(grDevices::dev.cur())
}

wgbs_render_pdf_png <- function(pdf_file, width_px, height_px, dpi = wgbs_figure_dpi) {
  if(!file.exists(pdf_file)) return(invisible(""))
  gs_bin <- Sys.which("gs")
  if(!nzchar(gs_bin)) stop("Ghostscript is required for WGBS base-graphics previews.", call. = FALSE)
  png_file <- wgbs_png_path(pdf_file)
  status <- suppressWarnings(system2(
    gs_bin,
    args = c(
      "-dSAFER", "-dBATCH", "-dNOPAUSE",
      "-dPDFFitPage",
      "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
      "-dFirstPage=1", "-dLastPage=1",
      "-sDEVICE=png16m",
      paste0("-r", dpi),
      paste0("-g", width_px, "x", height_px),
      shQuote(paste0("-sOutputFile=", png_file)),
      shQuote(pdf_file)
    ),
    stdout = FALSE,
    stderr = FALSE
  ))
  if(!identical(status, 0L) || !file.exists(png_file)) {
    stop("Failed to render PNG preview: ", pdf_file, call. = FALSE)
  }
  wgbs_finalize_rgb_png(png_file, dpi = dpi)
  invisible(png_file)
}

wgbs_normalize_existing_pdf <- function(pdf_file, width = 12, height = 8) {
  if(!file.exists(pdf_file)) return(invisible(""))
  gs_bin <- Sys.which("gs")
  if(!nzchar(gs_bin)) stop("Ghostscript is required to normalize package-generated PDFs.", call. = FALSE)
  dims <- wgbs_figure_dimensions(width, height)
  width_pt <- as.integer(round(dims$pdf_width_mm / 25.4 * 72))
  height_pt <- as.integer(round(dims$pdf_height_mm / 25.4 * 72))
  normalized_pdf <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(pdf_file)), "_normalized_"),
    tmpdir = dirname(pdf_file),
    fileext = ".pdf"
  )
  on.exit(unlink(normalized_pdf, force = TRUE), add = TRUE)
  status <- suppressWarnings(system2(
    gs_bin,
    args = c(
      "-dSAFER", "-dBATCH", "-dNOPAUSE",
      "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.4",
      "-dEmbedAllFonts=true", "-dSubsetFonts=true",
      "-dFIXEDMEDIA", "-dPDFFitPage",
      paste0("-dDEVICEWIDTHPOINTS=", width_pt),
      paste0("-dDEVICEHEIGHTPOINTS=", height_pt),
      shQuote(paste0("-sOutputFile=", normalized_pdf)),
      shQuote(pdf_file)
    ),
    stdout = FALSE,
    stderr = FALSE
  ))
  if(!identical(status, 0L) || !file.exists(normalized_pdf)) {
    stop("Failed to normalize package-generated PDF: ", pdf_file, call. = FALSE)
  }
  if(!file.copy(normalized_pdf, pdf_file, overwrite = TRUE)) {
    stop("Failed to replace package-generated PDF: ", pdf_file, call. = FALSE)
  }
  wgbs_render_pdf_png(
    pdf_file,
    width_px = as.integer(round(dims$width_mm / 25.4 * wgbs_figure_dpi)),
    height_px = as.integer(round(dims$height_mm / 25.4 * wgbs_figure_dpi))
  )
  invisible(pdf_file)
}

wgbs_dev_off <- function(which = grDevices::dev.cur()) {
  dev_key <- as.character(as.integer(which))
  export_spec <- .wgbs_pdf_devices$files[[dev_key]]
  result <- grDevices::dev.off(which = which)
  if(!is.null(export_spec)) {
    .wgbs_pdf_devices$files[[dev_key]] <- NULL
    wgbs_render_pdf_png(export_spec$file, export_spec$width_px, export_spec$height_px)
  }
  invisible(result)
}

wgbs_capture_base_pdf <- function(file, width = 8, height = 6, draw) {
  if(!is.function(draw)) stop("draw must be a function", call. = FALSE)
  device_id <- as.integer(wgbs_cairo_pdf(file, width = width, height = height))
  closed <- FALSE
  on.exit({
    active_devices <- grDevices::dev.list()
    if(!closed && !is.null(active_devices) && device_id %in% active_devices) {
      wgbs_dev_off(device_id)
    }
  }, add = TRUE)
  value <- draw()
  wgbs_dev_off(device_id)
  closed <- TRUE
  value
}

wgbs_draw_grid_dual <- function(file, width = 8, height = 6, draw) {
  draw_sub <- substitute(draw)
  draw_env <- parent.frame()
  dims <- wgbs_figure_dimensions(width, height)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    file,
    width = dims$width_in,
    height = dims$height_in,
    pointsize = wgbs_figure_base_size,
    family = wgbs_figure_font_family,
    bg = "white"
  )
  tryCatch(eval(draw_sub, envir = draw_env), finally = grDevices::dev.off())
  png_file <- wgbs_png_path(file)
  grDevices::png(
    filename = png_file,
    width = dims$width_mm,
    height = dims$height_mm,
    units = "mm",
    res = wgbs_figure_dpi,
    pointsize = wgbs_figure_base_size,
    type = "cairo-png",
    bg = "white"
  )
  tryCatch(eval(draw_sub, envir = draw_env), finally = grDevices::dev.off())
  wgbs_finalize_rgb_png(png_file)
  invisible(file)
}

# 中心规范的本地冻结副本；只统一出图，不改变统计变量。
.wgbs_source_files <- vapply(sys.frames(), function(f) if(is.null(f$ofile)) "" else as.character(f$ofile), character(1))
.wgbs_source_files <- .wgbs_source_files[nzchar(.wgbs_source_files)]
.wgbs_style_file <- file.path(dirname(tail(.wgbs_source_files, 1)), "figure_style.R")
if(!file.exists(.wgbs_style_file)) stop("缺少本地冻结绘图规范: ", .wgbs_style_file)
source(.wgbs_style_file, local = globalenv())
wgbs_ggsave <- function(filename, plot = ggplot2::last_plot(), device = grDevices::cairo_pdf,
                        path = NULL, scale = 1, width = 8, height = 6,
                        units = "in", dpi = 300, limitsize = TRUE, bg = "white",
                        create.dir = FALSE, ..., coordinate_system = "axes") {
  if(!is.null(path))filename<-file.path(path,filename)
  pfos_save(wgbs_standardize_plot(plot,coordinate_system),filename,
            height_mm=max(80,183*height/width))
  invisible(filename)
}

# 原 methylKit 相关图只修改固定标签字号；相关系数字号仍按 r 映射。
wgbs_plot_correlation <- function(meth) {
  mat<-as.data.frame(methylKit::percMethylation(meth)/100)
  corr<-methylKit:::.plotCorrelation
  shrink<-if(ncol(mat)>2).66 else .83
  styled_pairs<-function(...)graphics::pairs(...,cex.labels=1/shrink)
  env<-new.env(parent=environment(corr));env$pairs<-styled_pairs;environment(corr)<-env
  old<-par(cex.axis=7/(8*shrink));on.exit(par(old))
  corr(mat,title='CpG base pearson cor.',method='pearson')
}

# 原生 methylKit 分布图仅替换直方图填色，保留分箱、频数与百分比标签。
wgbs_plot_sample_stats <- function(object, kind = c("methylation", "coverage")) {
  kind <- match.arg(kind)
  group <- sample_table$group[match(object@sample.id, sample_table$sample)]
  if(length(group) != 1L || is.na(group) || !group %in% names(viz_pal$group)) stop("直方图文库身份未确认")
  color <- unname(viz_pal$group[group])
  method <- if(kind == "methylation") "getMethylationStats" else "getCoverageStats"
  draw <- methods::selectMethod(method, class(object)[1])
  replacements <- 0L
  replace_hist_color <- function(e) {
    if(!is.call(e)) return(e)
    if(identical(e[[1]], as.name("hist")) && !is.null(e$col)) {
      e$col <- color
      replacements <<- replacements + 1L
    }
    for(j in seq_along(e)[-1]) if(!is.null(e[[j]]) && !identical(e[[j]], quote(expr=))) e[[j]] <- replace_hist_color(e[[j]])
    e
  }
  body(draw) <- replace_hist_color(body(draw))
  if(replacements == 0L) stop("原生直方图颜色入口不存在")
  draw(object, plot = TRUE, both.strands = FALSE)
}

# 与原 methylKit 默认方法相同：完整 strict-unite CpG、比例矩阵、SD > 50%分位。
wgbs_methylkit_display_matrix <- function(meth) {
  source <- methylKit::getData(meth)
  source <- source[rowSums(is.na(source)) == 0, , drop = FALSE]
  matrix <- source[, meth@numCs.index, drop = FALSE] /
    (source[, meth@numCs.index, drop = FALSE] + source[, meth@numTs.index, drop = FALSE])
  names(matrix) <- meth@sample.ids
  # 中位数并列边界须与原包逐位一致，使用 methylKit 自带 rowSds。
  deviations <- methylKit:::rowSds(as.matrix(matrix))
  cutoff <- stats::quantile(deviations, .5)
  keep <- deviations > cutoff
  list(matrix = matrix[keep, , drop = FALSE], cutoff = unname(cutoff),
       common_cpg = nrow(source), selected_cpg = sum(keep),
       selected_identity = digest::digest(source[keep, 1:4, drop = FALSE], algo = "sha256"))
}

wgbs_plot_cluster <- function(meth) {
  source <- wgbs_methylkit_display_matrix(meth)
  draw <- methylKit:::.cluster
  # 替换原显示函数的rainbow表达式；距离和hclust代码原样保留。
  replace_palette <- function(expr) {
    if(!is.call(expr)) return(expr)
    if(identical(expr[[1]], as.name("rainbow"))) {
      return(as.call(list(as.name("rep"), unname(viz_pal$group[c(config$group1_name, config$group2_name)]), length.out = expr[[2]])))
    }
    for(j in seq_along(expr)[-1]) if(!identical(expr[[j]], quote(expr=))) expr[[j]] <- replace_palette(expr[[j]])
    expr
  }
  body(draw) <- replace_palette(body(draw))
  actual <- draw(source$matrix, dist.method = "correlation", hclust.method = "ward", plot = TRUE,
    treatment = meth@treatment, sample.ids = meth@sample.ids, context = meth@context)
  reference <- methylKit::clusterSamples(meth, dist = "correlation", method = "ward", plot = FALSE)
  stopifnot(identical(actual$merge, reference$merge), isTRUE(all.equal(actual$height, reference$height, tolerance = 1e-10)))
  audit <- file.path("result", "00_audit", "tables", "cluster_display_state.json")
  dir.create(dirname(audit), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(source = "methylBase strict-unite numerical counts", distance = "correlation", linkage = "ward (original methylKit)",
    sd_function = "methylKit:::rowSds",
    common_cpg = source$common_cpg, selected_cpg = source$selected_cpg, selected_identity_sha256 = source$selected_identity,
    sample_order = meth@sample.ids, merge = actual$merge, height = actual$height, original_merge_identical = identical(actual$merge, reference$merge),
    original_height_max_abs_error = max(abs(actual$height - reference$height)), group_colors = as.list(viz_pal$group)),
    audit, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  invisible(actual)
}

wgbs_plot_pca <- function(meth, out_pdf, screeplot = FALSE) {
  source <- wgbs_methylkit_display_matrix(meth)
  actual <- stats::prcomp(t(source$matrix), scale. = TRUE, center = TRUE)
  grDevices::pdf(NULL)
  reference <- tryCatch(methylKit::PCASamples(meth, screeplot = FALSE, obj.return = TRUE), finally = grDevices::dev.off())
  score_error <- max(abs(actual$x - reference$x))
  if(!is.finite(score_error) || score_error > 1e-9) stop("PCA显示坐标与原methylKit不同")
  proportions <- actual$sdev^2 / sum(actual$sdev^2) * 100
  coordinates <- data.frame(Sample = meth@sample.ids, Group = sample_group_factor(meth@sample.ids, config),
    PC1 = actual$x[, 1], PC2 = actual$x[, 2], row.names = NULL)
  if(screeplot) {
    data <- data.frame(Component = factor(paste0("PC", seq_along(proportions)), levels = paste0("PC", seq_along(proportions))), Explained = proportions)
    p <- ggplot2::ggplot(data, ggplot2::aes(Component, Explained)) +
      ggplot2::geom_col(fill = unname(viz_pal$group[config$group1_name]), width = .7) +
      ggplot2::labs(x = "Principal component", y = "Variance explained (%)") + theme_wgbs_classic()
  } else {
    p <- ggplot2::ggplot(coordinates, ggplot2::aes(PC1, PC2, color = Group, shape = Group)) +
      ggplot2::geom_point(size = 2.4) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = 7 / ggplot2::.pt, family = "Liberation Sans",
        seed = 1, show.legend = FALSE, max.overlaps = Inf, box.padding = .4, point.padding = .35, min.segment.length = 0) +
      ggplot2::scale_color_manual(values = viz_pal$group, breaks = c(config$group1_name, config$group2_name)) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17), c(config$group1_name, config$group2_name))) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .18)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = .18)) +
      ggplot2::labs(x = sprintf("PC1 (%.1f%%)", proportions[1]), y = sprintf("PC2 (%.1f%%)", proportions[2]),
        caption = "Library-level display; biological independence unknown") + theme_wgbs_classic() +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, hjust = 0))
  }
  ggsave(out_pdf, p, width = 8, height = 6)
  audit <- file.path(sub("/02_qc/.*$", "", out_pdf), "00_audit", "tables", "pca_display_state.json")
  dir.create(dirname(audit), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(source = "methylBase strict-unite numerical counts", sd_function = "methylKit:::rowSds", common_cpg = source$common_cpg,
    filter = "row SD > median SD; identical to methylKit default", cutoff = source$cutoff, selected_cpg = source$selected_cpg,
    selected_identity_sha256 = source$selected_identity, sample_order = meth@sample.ids,
    biological_unit = "library; independent animal/sire/litter relationship unknown", coordinates = coordinates,
    variance_percent = proportions, original_methylkit_score_max_abs_error = score_error), audit, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  invisible(actual)
}

# 复用原 MethylSeekR 画图函数；只固定文字，星号 CpG 点的 cex 保持原值。
wgbs_plot_final_segmentation <- function(...) {
  colors <- stats::setNames(unname(pfos_palette('accent')[c('teal','warm','purple')]),c('UMR','LMR','PMD'))
  fixed_labels <- function(e) {
    if(!is.call(e)) return(e)
    head <- paste(deparse(e[[1]]),collapse='')
    if(head %in% c('<-', '=') && identical(e[[2]],as.name('cols.seg'))) e[[3]] <- unname(colors[c('UMR','LMR')])
    if(head %in% c('<-', '=') && identical(e[[2]],as.name('col.PMD'))) e[[3]] <- unname(colors['PMD'])
    if(head=='par' && !is.null(e$mar) && is.call(e$mar) && length(e$mar)==5L) {e$mar[[3]]<-6;e$mar[[5]]<-5;e$mgp<-quote(c(2.5,1,0))}
    if(head=='axis') e$cex.axis <- quote(7/(8*graphics::par('cex')))
    if(head=='plot') {e$cex.lab <- quote(8/(8*graphics::par('cex')));e$cex.main <- quote(8/(8*graphics::par('cex')))}
    if(head=='text') e$cex <- quote(7/(8*graphics::par('cex')))
    for(j in seq_along(e)[-1]) if(!is.null(e[[j]]) && !identical(e[[j]],quote(expr=))) e[[j]] <- fixed_labels(e[[j]])
    if(head=='plot' && identical(e$xlab,'')) {
      key <- substitute(if(gg==1L) graphics::legend(
        x=graphics::par('usr')[2],y=graphics::par('usr')[4]+.14*diff(graphics::par('usr')[3:4]),
        xjust=1,yjust=0,legend=c('UMR','LMR','PMD'),pch=c(22,24,15),
        pt.bg=COLORS,col=c('black','black',PMD_COLOR),horiz=TRUE,bty='n',xpd=NA,
        cex=7/(8*graphics::par('cex'))),list(COLORS=unname(colors),PMD_COLOR=unname(colors['PMD'])))
      return(as.call(list(as.name('{'),e,key)))
    }
    e
  }
  draw <- MethylSeekR::plotFinalSegmentation
  body(draw) <- fixed_labels(body(draw))
  draw(...)
}
