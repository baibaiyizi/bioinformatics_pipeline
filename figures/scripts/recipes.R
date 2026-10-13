# 源表已完成模型拟合；本文件只实现原生R图形。
recipe_file <- tail(unlist(lapply(sys.frames(), function(frame) frame$ofile)), 1)
source(file.path(dirname(recipe_file), "figure_style.R"))

plot_violin <- function(data, output_stem, group_colors, value_label) {
  stopifnot(all(c("category", "group", "value") %in% names(data)),
            is.numeric(data$value), !is.null(names(group_colors)),
            !anyNA(data$category), !anyNA(data$group),
            all(as.character(data$group) %in% names(group_colors)))
  data$category <- factor(data$category, levels = sort(unique(as.character(data$category))))
  data$group <- factor(data$group, levels = names(group_colors))
  summary <- function(z) {
    z <- z[is.finite(z)]
    data.frame(y = median(z), ymin = unname(quantile(z, .25, type = 7)),
               ymax = unname(quantile(z, .75, type = 7)))
  }
  p <- ggplot2::ggplot(data, ggplot2::aes(category, value, fill = group)) +
    ggplot2::geom_violin(colour = NA, alpha = .30, trim = FALSE, position = ggplot2::position_dodge(.8)) +
    ggplot2::stat_summary(fun.data = summary, geom = "linerange", linewidth = .25,
      colour = pfos_palette("neutral")[["dark"]], position = ggplot2::position_dodge(.8), show.legend = FALSE) +
    ggplot2::stat_summary(fun = median, geom = "point", size = 1.1,
      colour = pfos_palette("neutral")[["dark"]], position = ggplot2::position_dodge(.8), show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(x = NULL, y = value_label, fill = "Group")
  pfos_save(p, paste0(output_stem, ".pdf"), height_mm = 137.25)
  invisible(p)
}

plot_ora <- function(data, output_stem) {
  stopifnot(all(c("term", "count", "query_size", "pvalue") %in% names(data)),
            !anyDuplicated(data$term), all(data$query_size > 0),
            all(data$count >= 0 & data$count <= data$query_size),
            all(data$pvalue >= 0 & data$pvalue <= 1))
  data$gene_ratio <- data$count / data$query_size
  data$minus_log10_p <- -log10(pmax(data$pvalue, .Machine$double.xmin))
  data$term <- factor(data$term, levels = rev(data$term))
  p <- ggplot2::ggplot(data, ggplot2::aes(gene_ratio, term, size = count, color = minus_log10_p)) +
    ggplot2::geom_point() + ggplot2::scale_size_area(max_size = 7, name = "Overlap genes") +
    ggplot2::scale_color_viridis_c(limits = unlist(pfos_profile$scales$ora_minus_log10_p),
      oob = scales::squish, breaks = c(0,5,10), labels = c("0", "5", "≥10"), name = "−log10(nominal P)") +
    ggplot2::labs(x = "Gene ratio (Count / tested query)", y = NULL)
  pfos_save(p, paste0(output_stem, ".pdf"), height_mm = 180)
}

plot_effect <- function(data, output_stem, value_label) {
  stopifnot(all(c("feature", "effect") %in% names(data)), !anyDuplicated(data$feature))
  data$direction <- ifelse(data$effect >= 0, "up", "down")
  data$feature <- factor(data$feature, levels = rev(data$feature))
  p <- ggplot2::ggplot(data, ggplot2::aes(effect, feature, color = direction)) +
    ggplot2::geom_vline(xintercept = 0, color = pfos_palette("neutral")[["dark"]], linewidth = .35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect, yend = feature), linewidth = .4) +
    ggplot2::geom_point(size = 2) + ggplot2::scale_color_manual(values = pfos_palette("de")[c("up", "down")]) +
    ggplot2::labs(x = value_label, y = NULL, color = "Signed direction")
  if (all(c("lower", "upper") %in% names(data))) p <- p + ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = feature), color = "#374151")
  pfos_save(p, paste0(output_stem, ".pdf"), height_mm = 150)
}

plot_marker <- function(data, output_stem) {
  stopifnot(all(c("gene", "group", "fraction", "value") %in% names(data)),
            !anyDuplicated(paste(data$gene, data$group)), all(data$fraction >= -1e-12 & data$fraction <= 1 + 1e-12))
  data$fraction <- pmin(1, pmax(0, data$fraction))
  groups <- names(pfos_palette("clusters"))
  stopifnot(all(as.character(data$group) %in% groups))
  data$group <- factor(data$group, levels = groups[groups %in% as.character(data$group)])
  p <- ggplot2::ggplot(data, ggplot2::aes(gene, group, size = fraction, color = value)) +
    ggplot2::geom_point() + ggplot2::scale_size_area(max_size = 6, breaks = c(.25,.5,1), labels = c("25%","50%","100%"), name = "Detected cells") +
    ggplot2::scale_color_gradient2(low = pfos_palette("diverging")[["low"]], mid = pfos_palette("diverging")[["mid"]], high = pfos_palette("diverging")[["high"]], midpoint = 0, limits = unlist(pfos_profile$scales$marker_zscore), oob = scales::squish, breaks = c(-2.5,0,2.5), labels = c("≤−2.5","0","≥2.5"), name = "Within-gene\ncluster z-score") +
    ggplot2::labs(x = "Marker gene", y = "Cluster") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  pfos_save(p, paste0(output_stem, ".pdf"), height_mm = 150)
}

plot_program_matrix <- function(data, captured, output_stem) {
  stopifnot(all(c("feature","library_id","condition","cell_type","value") %in% names(data)),
            !anyDuplicated(paste(data$feature, data$library_id, data$cell_type)))
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("需要已安装的ComplexHeatmap")
  data$sample <- paste(data$library_id, data$cell_type, sep = " / ")
  annotations <- unique(data[c("sample", "library_id", "condition", "cell_type")])
  annotations <- annotations[order(match(annotations$condition, names(pfos_palette("group"))), match(annotations$cell_type, names(pfos_palette("cell_type")))), ]
  matrix <- tidyr::pivot_wider(data[c("feature", "sample", "value")], names_from = "sample", values_from = "value")
  rn <- matrix$feature;matrix <- as.matrix(matrix[annotations$sample]);rownames(matrix) <- rn
  ncol_name <- intersect(c("cells", "n_cells", "count", "captured_cells"), names(captured))[1]
  stopifnot(!is.na(ncol_name))
  ix <- match(paste(annotations$library_id, annotations$cell_type), paste(captured$library_id, captured$cell_type))
  n <- captured[[ncol_name]][ix]
  stopifnot(!anyNA(n))
  draw <- function() {
    annotation <- ComplexHeatmap::HeatmapAnnotation(
      condition = annotations$condition, cell_type = annotations$cell_type,
      captured_cells = ComplexHeatmap::anno_barplot(n, border = FALSE, gp = grid::gpar(fill = pfos_palette("group")[annotations$condition]), axis_param = list(gp = grid::gpar(fontsize = 7))),
      col = list(condition = pfos_palette("group"), cell_type = pfos_palette("cell_type")),
      annotation_name_gp = grid::gpar(fontsize = 7), annotation_legend_param = list(condition = list(title = "Library condition", at = c("WT", "PFOS"))))
    heatmap <- ComplexHeatmap::Heatmap(matrix,
      name = "Module-wise z", col = circlize::colorRamp2(c(-2.5,0,2.5), unname(pfos_palette("diverging"))),
      na_col = pfos_profile$scales$missing_color, top_annotation = annotation,
      column_split = factor(annotations$condition, levels = names(pfos_palette("group"))),
      cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 7),
      column_names_rot = 65, column_title_gp = grid::gpar(fontsize = 8),
      column_title = "Descriptive programs: one library per condition; no group inference",
      heatmap_legend_param = list(at = c(-2.5,0,2.5), labels = c("≤−2.5","0","≥2.5")),
      rect_gp = grid::gpar(col = "white", lwd = .5))
    ComplexHeatmap::draw(heatmap, merge_legends = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  }
  pfos_save(draw, paste0(output_stem, ".pdf"), height_mm = 190)
}
