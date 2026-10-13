# PFOS / 中心模板共同成图规范。项目运行其本地冻结副本。
PFOS_PROFILE_RELATIVE <- "../config/figure_style.json"
pfos_profile <- NULL
pfos_figure_style <- NULL

pfos_set_profile <- function(profile_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("统一绘图需要已有依赖 jsonlite")
  path <- normalizePath(profile_path, mustWork = TRUE)
  profile <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  required <- c("group", "de", "neutral", "diverging", "sequential", "binary", "quadrant")
  if (!all(required %in% names(profile$palettes))) stop("绘图profile缺少必需语义调色板")
  style <- profile$typography
  style$colors <- c(unlist(profile$palettes$group), NS = profile$palettes$de$ns)
  style$diverging <- unname(unlist(profile$palettes$diverging))
  style$profile_id <- profile$profile_id
  style$profile_path <- path
  pfos_profile <<- profile
  pfos_figure_style <<- style
  invisible(profile)
}

pfos_palette <- function(name) {
  value <- pfos_profile$palettes[[name]]
  if (is.null(value)) stop("未注册的语义调色板: ", name)
  unlist(value, use.names = TRUE)
}

pfos_identity_colors <- function(keys) {
  # 对无方向含义的 term/cluster 身份编码；同一名称不受子集或出现顺序影响。
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  rule <- pfos_profile$category_identity
  if (!identical(rule$method, "murmur3_hcl")) stop("缺少固定分类身份配色规则")
  hashes <- as.double(digest::digest2int(keys, seed = as.integer(rule$seed))) + 2147483648
  stats::setNames(grDevices::hcl(h = hashes / 4294967296 * 360,
    c = rule$chroma, l = rule$luminance, fixup = TRUE), keys)
}

pfos_profile_sources <- vapply(sys.frames(), function(frame) {
  value <- frame$ofile
  if (is.character(value) && length(value) == 1L) value else ""
}, character(1))
pfos_profile_source <- tail(pfos_profile_sources[basename(pfos_profile_sources) == "figure_style.R"], 1)
if (length(pfos_profile_source) != 1L) stop("请使用source加载项目本地figure_style.R")
pfos_set_profile(file.path(dirname(pfos_profile_source), PFOS_PROFILE_RELATIVE))

pfos_theme <- function() {
  s <- pfos_figure_style
  ggplot2::theme_classic(base_size = s$text_pt, base_family = s$family) +
    ggplot2::theme(text = ggplot2::element_text(family = s$family, size = s$text_pt),
      axis.text = ggplot2::element_text(size = s$tick_pt, color = "#333333"),
      axis.title = ggplot2::element_text(size = s$text_pt),
      legend.text = ggplot2::element_text(size = s$tick_pt),
      legend.title = ggplot2::element_text(size = s$text_pt),
      strip.text = ggplot2::element_text(size = s$text_pt),
      plot.title = ggplot2::element_text(size = s$text_pt),
      axis.line = ggplot2::element_line(linewidth = s$line_mm),
      axis.ticks = ggplot2::element_line(linewidth = s$line_mm),
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(8, 12, 8, 8, unit = "pt"))
}

pfos_style_plot <- function(plot) {
  # 分别处理子图，保留各自的空轴、隐藏标题、旋转和布局。
  s <- pfos_figure_style
  if (inherits(plot,"patchwork")) plot$patches$plots <- lapply(plot$patches$plots,pfos_style_plot)
  elements <- list(text=ggplot2::element_text(family=s$family,size=s$text_pt),
    panel.grid=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),
    panel.background=ggplot2::element_rect(fill="white",color=NA),
    legend.background=ggplot2::element_rect(fill="white",color=NA),
    legend.key=ggplot2::element_rect(fill="white",color=NA),
    strip.background=ggplot2::element_rect(fill="white",color=NA),
    plot.background=ggplot2::element_rect(fill="white",color=NA))
  # 普通笛卡尔图补回缺失轴线；保留显式无轴布局、既有自定义轴及矩阵/流图。
  ordinary <- inherits(plot$coordinates, "CoordCartesian") &&
    !(length(plot$layers) && inherits(plot$layers[[1]]$geom, c("GeomTile", "GeomRaster"))) &&
    !any(vapply(plot$layers, function(layer)
      inherits(layer$geom, c("GeomAlluvium", "GeomStratum", "GeomFlow")), logical(1)))
  if (ordinary) {
    full <- ggplot2::theme_get() + plot$theme
    line <- ggplot2::element_line(colour="#333333", linewidth=s$line_mm, inherit.blank=TRUE)
    repaired <- character()
    for (axis in c("x", "y")) {
      side <- if (axis == "x") "bottom" else "left"
      suffix <- paste(axis, side, sep=".")
      visible <- vapply(paste0(c("axis.text.", "axis.ticks."), suffix), function(n) {
        value <- ggplot2::calc_element(n, full)
        !is.null(value) && !inherits(value, "element_blank")
      }, logical(1))
      title <- plot$labels[[axis]]
      title_visible <- !is.null(title) && length(title)>0 && !identical(title,"") &&
        !inherits(ggplot2::calc_element(paste0("axis.title.",axis),full),"element_blank")
      if ((!any(visible) && !title_visible) || !inherits(ggplot2::calc_element(paste0("axis.line.",suffix),full),"element_blank")) next
      for (prefix in c("axis.line", "axis.ticks")) {
        for (n in c(prefix, paste(prefix,axis,sep="."), paste(prefix,suffix,sep=".")))
          if (inherits(ggplot2::calc_element(n,full),"element_blank")) elements[[n]] <- line
      }
      elements[[paste0("axis.ticks.",suffix)]] <- line
      repaired <- c(repaired,suffix)
    }
    if (length(repaired)) for (suffix in setdiff(c("x.bottom","x.top","y.left","y.right"),repaired)) {
      for (prefix in c("axis.line.","axis.ticks.")) {
        n <- paste0(prefix,suffix)
        if (inherits(ggplot2::calc_element(n,full),"element_blank")) elements[[n]] <- ggplot2::element_blank()
      }
    }
  }
  sizes <- c(axis.text=s$tick_pt,axis.text.x=s$tick_pt,axis.text.y=s$tick_pt,
    axis.title=s$text_pt,axis.title.x=s$text_pt,axis.title.y=s$text_pt,
    legend.text=s$tick_pt,legend.title=s$text_pt,strip.text=s$text_pt,plot.title=s$text_pt)
  for(n in names(sizes)) {
    element<-ggplot2::calc_element(n,plot$theme)
    if(!inherits(element,"element_blank"))elements[[n]]<-ggplot2::element_text(family=s$family,size=unname(sizes[n]))
  }
  plot+do.call(ggplot2::theme,elements)
}

pfos_save <- function(plot, pdf, height_mm = 140) {
  stopifnot(is.finite(height_mm), height_mm > 0, grepl("\\.pdf$", pdf))
  Sys.setenv(MAGICK_THREAD_LIMIT="1")
  if (!requireNamespace("magick", quietly = TRUE)) stop("成图需要已安装的 magick，以核验并写出 RGB PNG")
  s <- pfos_figure_style
  dir.create(dirname(pdf), recursive = TRUE, showWarnings = FALSE)
  png <- sub("\\.pdf$", ".png", pdf)
  if (inherits(plot, "ggplot")) plot <- pfos_style_plot(plot)
  draw <- function() {
    if (is.function(plot)) plot() else if (inherits(plot, c("gtable", "grob"))) {
      grid::grid.newpage(); grid::grid.draw(plot)
    } else print(plot)
  }
  grDevices::cairo_pdf(pdf, width = s$width_mm / 25.4, height = height_mm / 25.4,
                      family = s$family, pointsize = s$text_pt, bg = "white")
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(png, width = s$width_mm, height = height_mm, units = "mm", res = s$dpi,
                type = "cairo-png", bg = "white", pointsize = s$text_pt)
  tryCatch(draw(), finally = grDevices::dev.off())
  im <- magick::image_read(png)
  magick::image_write(magick::image_background(im, "white", flatten = TRUE), png,
                      format = "png", density = paste0(s$dpi,"x",s$dpi),
                      defines = c("png:color-type" = "2"))
  invisible(c(pdf = pdf, png = png))
}
