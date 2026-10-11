# 图形公共契约：保留现有图型、布局和 pastel 配色，统一最终印刷字号及导出。
# Academic Figure Skill / Bizard：继承物理尺寸、字体嵌入及视觉验收原则。
# 用户指定 Liberation Sans、183 mm、PDF 单次绘制后生成 PNG，优先于技能默认模板。

FIGURE_WIDTH_MM <- 183
FIGURE_FONT <- "Liberation Sans"
FIGURE_DPI <- 300L

figure_measure <- function(fn, width = 183 / 25.4, height = 7) {
  # grob construction must use Cairo font metrics too, before the final page opens.
  # Otherwise R silently creates Rplots.pdf and measures Liberation Sans as Helvetica.
  scratch <- tempfile("figure_metrics_", fileext = ".pdf")
  grDevices::cairo_pdf(scratch, width = width, height = height,
    family = FIGURE_FONT, pointsize = 8)
  on.exit({ grDevices::dev.off(); unlink(scratch) }, add = TRUE)
  fn()
}

figure_theme <- function() {
  ggplot2::theme(
    text = ggplot2::element_text(family = FIGURE_FONT, size = 8),
    axis.title = ggplot2::element_text(family = FIGURE_FONT, size = 8),
    axis.text = ggplot2::element_text(family = FIGURE_FONT, size = 7),
    legend.title = ggplot2::element_text(family = FIGURE_FONT, size = 8),
    legend.text = ggplot2::element_text(family = FIGURE_FONT, size = 7),
    # Horizontal legends use the entire figure width, including long y labels.
    legend.location = "plot",
    strip.text = ggplot2::element_text(family = FIGURE_FONT, size = 8),
    plot.title = ggplot2::element_text(family = FIGURE_FONT, size = 9),
    plot.subtitle = ggplot2::element_text(family = FIGURE_FONT, size = 8),
    plot.caption = ggplot2::element_text(family = FIGURE_FONT, size = 7),
    plot.tag = ggplot2::element_text(family = FIGURE_FONT, size = 10, face = "bold"),
    panel.grid = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )
}

figure_wrap_label <- function(label, width_mm) {
  if (!is.character(label) || length(label) != 1L) return(label)
  lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
  # A panel title starts after its y-axis allocation. Reserve that space before
  # wrapping, and measure real Cairo glyph widths instead of average characters.
  target_mm <- max(14, width_mm - 18)
  figure_measure(function() paste(vapply(lines, function(line) {
    chars <- max(8L, floor(target_mm / 1.6))
    wrap_line <- function(n) {
      tokens <- strsplit(line, "[[:space:]]+")[[1L]]
      tokens <- vapply(tokens, function(token) {
        starts <- seq.int(1L, max(1L, nchar(token)), by = n)
        paste(substring(token, starts, starts + n - 1L), collapse = " ")
      }, character(1))
      strwrap(paste(tokens, collapse = " "), width = n)
    }
    wrapped <- wrap_line(chars)
    while (length(wrapped) && chars > 1L &&
           max(graphics::strwidth(wrapped, units = "inches", cex = 9 / 8, font = 2)) * 25.4 > target_mm) {
      chars <- chars - 1L
      wrapped <- wrap_line(chars)
    }
    paste(wrapped, collapse = "\n")
  }, character(1)), collapse = "\n"))
}

figure_normalize_plot <- function(p, available_width_mm = FIGURE_WIDTH_MM) {
  if (!inherits(p, "ggplot")) return(p)
  # A theme_void axis must stay absent; changing only its font must not revive it.
  original_theme <- ggplot2::theme_get() + p$theme
  blank_names <- c("axis.title", "axis.title.x", "axis.title.y", "axis.text", "axis.text.x", "axis.text.y", "legend.title")
  blanks <- stats::setNames(lapply(blank_names, function(nm) ggplot2::calc_element(nm, original_theme)), blank_names)
  blanks <- blanks[vapply(blanks, inherits, logical(1), "element_blank")]
  # Polar reference rings carry quantitative information and remain in radar plots.
  polar_grid <- if (inherits(p$coordinates, c("CoordPolar", "CoordRadar"))) {
    elements <- c("panel.grid", "panel.grid.major", "panel.grid.minor",
      "panel.grid.major.x", "panel.grid.major.y", "panel.grid.minor.x", "panel.grid.minor.y")
    stats::setNames(lapply(elements, function(nm) ggplot2::calc_element(nm, ggplot2::theme_get() + p$theme)), elements)
  } else NULL
  # Clone layers before changing their defaults: a layer can be shared by plots.
  p$layers <- lapply(p$layers, function(layer) {
    if (inherits(layer$geom, c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel"))) {
      # ggproto keeps its parent expression lazy. Reassigning `layer` here would
      # make the clone inherit from itself and overflow on its first text layer.
      copy <- ggplot2::ggproto(NULL, layer)
      copy$aes_params <- as.list(layer$aes_params)
      copy$aes_params$family <- FIGURE_FONT
      copy$aes_params$size <- 8 / ggplot2::.pt
      return(copy)
    }
    layer
  })
  if (inherits(p, "patchwork")) {
    ncol <- p$patches$layout$ncol
    if (is.null(ncol) || !length(ncol) || !is.finite(ncol)) ncol <- ceiling(sqrt(length(p$patches$plots) + 1L))
    child_width <- available_width_mm / ncol
    p$patches$plots <- lapply(p$patches$plots, figure_normalize_plot, available_width_mm = child_width)
    # Children were normalized independently, including their polar reference
    # rings. Apply this theme only to the current plot instead of overwriting them.
    p <- p + figure_theme()
    p$patches$annotation$theme <- p$patches$annotation$theme + figure_theme()
    for (nm in c("title", "subtitle", "caption")) {
      p$patches$annotation[[nm]] <- figure_wrap_label(p$patches$annotation[[nm]], available_width_mm)
    }
    available_width_mm <- child_width
  } else {
    p <- p + figure_theme()
  }
  if (length(blanks)) p <- p + do.call(ggplot2::theme, blanks)
  if (length(polar_grid)) p <- p + do.call(ggplot2::theme, polar_grid)
  for (nm in c("title", "subtitle", "caption")) p$labels[[nm]] <- figure_wrap_label(p$labels[[nm]], available_width_mm)
  p
}

figure_plot_metric <- function(p) {
  if (!inherits(p, "ggplot")) return("see_source_data")
  labels <- p$labels[intersect(c("y", "fill", "x"), names(p$labels))]
  values <- vapply(labels, function(x) {
    if (is.null(x)) return("")
    paste(if (is.character(x)) x else deparse(x), collapse = " ")
  }, character(1))
  values <- paste0(names(values)[nzchar(values)], ": ", values[nzchar(values)])
  if (inherits(p, "patchwork")) {
    values <- c(values, vapply(p$patches$plots, figure_plot_metric, character(1)))
  }
  values <- unique(values[nzchar(values) & values != "see_source_data"])
  if (length(values)) paste(values, collapse = " | ") else "see_source_data"
}

figure_context <- function(metric, denominator, applicability) {
  lookup <- function(name, default = "unspecified") {
    value <- get0(name, envir = parent.frame(2), inherits = TRUE, ifnotfound = NULL)
    if (is.null(value) || !length(value) || !nzchar(as.character(value[[1]]))) value <- Sys.getenv(name, unset = default)
    if (!length(value) || !nzchar(as.character(value[[1]]))) default else as.character(value[[1]])
  }
  project <- lookup("PANDORA_PROJECT_NAME", "")
  if (!nzchar(project)) project <- basename(lookup("ROOT_DIR"))
  # Only an explicit axis/legend denominator is inferred; within-class and
  # transformed measures must not be silently labelled clean-read abundance.
  if (is.null(denominator)) {
    metric_labels <- strsplit(metric, " | ", fixed = TRUE)[[1]]
    y_labels <- metric_labels[startsWith(metric_labels, "y:")]
    if (!length(y_labels)) y_labels <- metric_labels
    clean <- grepl("clean[ _-]?reads", y_labels, ignore.case = TRUE)
    denominator <- if (all(clean)) "clean_reads" else if (any(clean)) "mixed; see_source_data" else "unspecified"
  }
  data.frame(Project = project, Species = lookup("PANDORA_SPECIES"),
    Mode = lookup("PANDORA_MODE"), AnalysisStage = lookup("PANDORA_ANALYSIS_STAGE"),
    AnalysisScope = lookup("PANDORA_ANALYSIS_SCOPE"), Module = lookup("PANDORA_CURRENT_MODULE"),
    Metric = metric, Denominator = denominator,
    Applicability = if (is.null(applicability)) lookup("PANDORA_FIGURE_APPLICABILITY") else applicability,
    stringsAsFactors = FALSE)
}

figure_normalize_grob <- function(g, path = "") {
  if (!grid::is.grob(g)) return(g)
  here <- paste(path, g$name, sep = "/")
  if (inherits(g, "text")) {
    size <- if (grepl("tag", here)) 10 else
      if (grepl("guide|legend", here) && grepl("/title", here)) 8 else
      if (grepl("(plot.title|/title|/main)", here)) 9 else
      if (grepl("axis|guide|legend|row_names|col_names", here)) 7 else 8
    if (length(g$label) > 30L && !grepl("axis|guide|legend", here)) size <- 6
    gp <- if (is.null(g$gp)) list() else as.list(g$gp)
    gp$fontfamily <- FIGURE_FONT
    gp$fontsize <- size
    # cex multiplies fontsize in grid; avoid legacy size multiplying the contract.
    gp$cex <- 1
    g$gp <- do.call(grid::gpar, gp)
  }
  if (!is.null(g$grobs)) {
    g$grobs <- lapply(seq_along(g$grobs), function(i) {
      label <- if (!is.null(g$layout$name)) g$layout$name[[i]] else as.character(i)
      figure_normalize_grob(g$grobs[[i]], paste(here, label, sep = "/"))
    })
  }
  if (!is.null(g$children) && length(g$children)) {
    g$children <- do.call(grid::gList, lapply(g$children, figure_normalize_grob, path = here))
  }
  g
}

figure_content_height <- function(g, width, height) {
  if (length(width) != 1L || length(height) != 1L || !is.finite(width) ||
      !is.finite(height) || width <= 0 || height <= 0) stop("Figure dimensions must be positive.")
  needed <- max(80, FIGURE_WIDTH_MM * height / width)
  inspect <- function(x) {
    if (!grid::is.grob(x)) return(invisible(NULL))
    if (inherits(x, "gtable")) {
      rows <- x$layout$t[grepl("^panel", x$layout$name)]
      if (length(rows)) needed <<- max(needed, 30 + 42 * length(unique(rows)))
      ri <- which(x$layout$name == "row_names")
      for (i in ri) {
        n <- length(x$grobs[[i]]$label)
        needed <<- max(needed, 40 + 2.8 * n)
      }
      lapply(x$grobs, inspect)
    }
    if (!is.null(x$children)) lapply(x$children, inspect)
    invisible(NULL)
  }
  inspect(g)
  needed
}

figure_plot_data <- function(p, prefix = "plot") {
  out <- list()
  if (!inherits(p, "ggplot")) return(out)
  if (is.data.frame(p$data)) out[[paste0(prefix, "_data")]] <- p$data
  for (i in seq_along(p$layers)) {
    if (is.data.frame(p$layers[[i]]$data)) out[[paste0(prefix, "_layer", i)]] <- p$layers[[i]]$data
  }
  if (inherits(p, "patchwork")) {
    for (i in seq_along(p$patches$plots)) {
      out <- c(out, figure_plot_data(p$patches$plots[[i]], paste0(prefix, "_panel", i)))
    }
  }
  out
}

figure_write_data <- function(data, stem, status) {
  if (is.null(data) || (is.list(data) && !is.data.frame(data) && !length(data))) {
    return(list(status = "not_supplied", files = "", rows = 0L))
  }
  data_dir <- paste0(stem, ".source_data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  # RDS preserves list columns, matrix names and exact values; CSV is a readable companion.
  rds <- file.path(data_dir, "source_data.rds")
  saveRDS(data, rds)
  tables <- if (is.data.frame(data) || is.matrix(data)) list(data = data) else data
  index <- list()
  for (i in seq_along(tables)) {
    x <- tables[[i]]
    if (!is.data.frame(x) && !is.matrix(x)) next
    if (is.matrix(x)) x <- data.frame(Feature = rownames(x), x, check.names = FALSE)
    if (!ncol(x)) next
    nm <- if (!is.null(names(tables)) && nzchar(names(tables)[[i]])) names(tables)[[i]] else paste0("data", i)
    name <- paste0(sprintf("%02d_", i), gsub("[^A-Za-z0-9_.-]", "_", nm), ".csv")
    list_cols <- vapply(x, is.list, logical(1))
    for (j in which(list_cols)) x[[j]] <- vapply(x[[j]], toString, character(1))
    file <- file.path(data_dir, name)
    utils::write.csv(x, file, row.names = FALSE, na = "")
    index[[length(index) + 1L]] <- data.frame(Source = nm, File = file, Rows = nrow(x),
      Columns = ncol(x), ListColumnsSerialized = any(list_cols), stringsAsFactors = FALSE)
  }
  if (length(index)) {
    index <- do.call(rbind, index)
    utils::write.csv(index, file.path(data_dir, "source_data_manifest.csv"), row.names = FALSE)
    return(list(status = status, files = paste(c(rds, index$File), collapse = ";"), rows = sum(index$Rows)))
  }
  list(status = paste0(status, "_rds_only"), files = rds, rows = NA_integer_)
}

figure_export <- function(file_path, draw_fn, width, height, data = NULL,
                          data_status = "explicit", figure_id = NULL, kind = "base",
                          metric = "see_source_data", denominator = NULL, applicability = NULL) {
  if (!grepl("[.]pdf$", file_path, ignore.case = TRUE)) stop("Figure master must have a .pdf suffix: ", file_path)
  if (!isTRUE(capabilities("cairo"))) stop("Cairo support is required for embedded-font PDF output.")
  gs <- Sys.which("gs")
  if (!nzchar(gs)) stop("Ghostscript (gs) is required to render the PDF preview.")
  if (!requireNamespace("systemfonts", quietly = TRUE) ||
      !FIGURE_FONT %in% systemfonts::system_fonts()$family) stop("Install Liberation Sans before rendering figures.")
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  file_path <- normalizePath(file_path, winslash = "/", mustWork = FALSE)
  stem <- sub("[.]pdf$", "", file_path, ignore.case = TRUE)
  png_path <- paste0(stem, ".png")
  grDevices::cairo_pdf(file_path, width = FIGURE_WIDTH_MM / 25.4,
    height = height / 25.4, family = FIGURE_FONT, pointsize = 8, bg = "white", onefile = TRUE)
  tryCatch({
    graphics::par(family = FIGURE_FONT, ps = 8, cex = 1, cex.axis = 7 / 8, bty = "l",
      cex.lab = 1, cex.main = 9 / 8, cex.sub = 7 / 8)
    draw_fn()
  }, finally = grDevices::dev.off())
  gs_output <- system2(gs, c("-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
    paste0("-r", FIGURE_DPI), "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
    shQuote(paste0("-sOutputFile=", png_path)), shQuote(file_path)), stdout = TRUE, stderr = TRUE)
  if ((!is.null(attr(gs_output, "status")) && attr(gs_output, "status") != 0L) || !file.exists(png_path)) {
    stop("Ghostscript preview failed: ", paste(gs_output, collapse = "\n"))
  }
  source <- figure_write_data(data, stem, data_status)
  manifest <- data.frame(FigureID = if (is.null(figure_id)) basename(stem) else figure_id,
    PDF = file_path, PNG = png_path, PlotKind = kind,
    WidthMM = FIGURE_WIDTH_MM, HeightMM = height, DPI = FIGURE_DPI, Font = FIGURE_FONT,
    SourceDataStatus = source$status, SourceDataFiles = source$files, SourceDataRows = source$rows,
    PDFMD5 = unname(tools::md5sum(file_path)), PNGMD5 = unname(tools::md5sum(png_path)),
    RenderedUTC = format(Sys.time(), tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE)
  manifest <- cbind(manifest, figure_context(metric, denominator, applicability))
  utils::write.csv(manifest, paste0(stem, ".figure_manifest.csv"), row.names = FALSE)
  invisible(file_path)
}

save_ggplot_file <- function(plot_obj, file_path, width = 11, height = 7,
                             dpi = 300, data = NULL, figure_id = NULL,
                             metric = NULL, denominator = NULL, applicability = NULL) {
  if (is.null(metric)) metric <- figure_plot_metric(plot_obj)
  data_status <- if (is.null(data)) "extracted_plot_and_layer_data" else "explicit"
  if (is.null(data)) data <- figure_plot_data(plot_obj)
  p <- figure_normalize_plot(plot_obj)
  kind <- if (inherits(p, "patchwork")) "patchwork" else if (inherits(p, "ggplot")) "ggplot" else "grid"
  g <- figure_measure(function() {
    if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else
      if (inherits(p, "ggplot")) ggplot2::ggplotGrob(p) else p
  })
  if (!grid::is.grob(g)) stop("Expected a ggplot, patchwork or grid grob.")
  g <- figure_normalize_grob(g)
  height_mm <- figure_content_height(g, width, height)
  figure_export(file_path, function() { grid::grid.newpage(); grid::grid.draw(g) },
    width, height_mm, data, data_status, figure_id, kind, metric, denominator, applicability)
}

save_plot_device <- function(file_path, draw_fn, width = 11, height = 7,
                              dpi = 300, data = NULL, figure_id = NULL,
                              metric = "see_source_data", denominator = NULL, applicability = NULL) {
  is_grid <- grepl("grid::|grid[.]draw|grid[.]newpage", paste(deparse(body(draw_fn)), collapse = " "))
  if (is_grid) {
    # Capture the existing grid drawing, including manually added pheatmap blocks.
    scratch <- tempfile("figure_capture_", fileext = ".pdf")
    on.exit(unlink(scratch), add = TRUE)
    # grid.grabExpr restores its previous device with dev.set(); a null previous
    # device would create an unwanted Rplots.pdf. Supply a real measuring device.
    g <- figure_measure(function() grid::grid.grabExpr(draw_fn(),
      width = FIGURE_WIDTH_MM / 25.4,
      height = max(80, FIGURE_WIDTH_MM * height / width) / 25.4, wrap = FALSE,
      device = function(width, height) grDevices::cairo_pdf(scratch,
        width = width, height = height, family = FIGURE_FONT, pointsize = 8)))
    g <- figure_normalize_grob(g)
    height_mm <- figure_content_height(g, width, height)
    return(figure_export(file_path, function() { grid::grid.newpage(); grid::grid.draw(g) },
      width, height_mm, data, "explicit", figure_id, "grid_callback", metric, denominator, applicability))
  }
  figure_export(file_path, draw_fn, width, figure_content_height(NULL, width, height),
    data, "explicit", figure_id, "base_callback", metric, denominator, applicability)
}

figure_pheatmap <- function(..., width = 11, height = 7) {
  args <- list(...)
  mat <- if (!is.null(args$mat)) args$mat else args[[1L]]
  if (nrow(mat)<2L) args$cluster_rows <- FALSE
  if (ncol(mat)<2L) args$cluster_cols <- FALSE
  finite_values <- as.numeric(mat[is.finite(mat)])
  if (is.null(args$breaks) && length(finite_values) && diff(range(finite_values))==0) {
    center <- finite_values[[1L]]
    delta <- max(abs(center)*0.01,0.5)
    color_n <- if (is.null(args$color)) 100L else length(args$color)
    args$breaks <- seq(center-delta,center+delta,length.out=color_n+1L)
  }
  args$filename <- NA_character_
  args$silent <- TRUE
  args$fontsize <- 8
  args$fontsize_row <- if (nrow(mat) > 30) 6 else 7
  args$fontsize_col <- 7
  args$fontfamily <- FIGURE_FONT
  # pheatmap opens pdf(NULL) internally even with silent=TRUE. Override that
  # device in a local function environment, never in the package namespace,
  # so its layout measures the same embedded font as the final Cairo page.
  metric_file <- tempfile("pheatmap_metrics_", fileext = ".pdf")
  on.exit(unlink(metric_file), add = TRUE)
  pheatmap_fn <- pheatmap::pheatmap
  environment(pheatmap_fn) <- list2env(list(pdf = function(file = NULL, ...) {
    grDevices::cairo_pdf(metric_file, width = FIGURE_WIDTH_MM / 25.4,
      height = max(80, FIGURE_WIDTH_MM * height / width) / 25.4,
      family = FIGURE_FONT, pointsize = 8)
  }), parent = environment(pheatmap_fn))
  heatmap <- figure_measure(function() do.call(pheatmap_fn, args))
  heatmap$gtable <- figure_normalize_grob(heatmap$gtable)
  heatmap
}

save_pheatmap_file <- function(..., filename, width, height, data = NULL, figure_id = NULL,
                              metric = "see_source_data", denominator = NULL, applicability = NULL) {
  args <- list(...)
  mat <- if (!is.null(args$mat)) args$mat else args[[1L]]
  if (is.null(data)) data <- list(matrix = mat, annotation_row = args$annotation_row, annotation_col = args$annotation_col)
  heatmap <- do.call(figure_pheatmap, c(args, list(width = width, height = height)))
  g <- heatmap$gtable
  figure_export(filename, function() { grid::grid.newpage(); grid::grid.draw(g) }, width,
    figure_content_height(g, width, height), data, "explicit_matrix_and_annotations", figure_id, "pheatmap",
    metric, denominator, applicability)
}
