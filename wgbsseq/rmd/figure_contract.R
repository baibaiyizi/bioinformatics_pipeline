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

wgbs_standardize_plot <- function(plot, coordinate_system = c("axes", "none")) {
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
      paste0("-sOutputFile=", png_file),
      pdf_file
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
      paste0("-sOutputFile=", normalized_pdf),
      pdf_file
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
