# 保留原热图矩阵、注释和聚类，仅统一文字与成图参数。
pfos_heatmap_labels <- function(labels) {
  vapply(as.character(labels),function(label) {
    if(nchar(label)<=35L)return(label)
    stringr::str_wrap(gsub("_"," ",label,fixed=TRUE),width=35)
  },character(1),USE.NAMES=FALSE)
}

pfos_heatmap <- function(...) {
  args<-list(...)
  args$row_names_gp<-grid::gpar(fontsize=7,fontfamily="Liberation Sans")
  args$column_names_gp<-grid::gpar(fontsize=7,fontfamily="Liberation Sans")
  args$row_title_gp<-grid::gpar(fontsize=8,fontfamily="Liberation Sans")
  args$column_title_gp<-grid::gpar(fontsize=8,fontfamily="Liberation Sans")
  args$heatmap_legend_param<-modifyList(args$heatmap_legend_param %||% list(),list(
    labels_gp=grid::gpar(fontsize=7,fontfamily="Liberation Sans"),
    title_gp=grid::gpar(fontsize=8,fontfamily="Liberation Sans")))
  matrix<-if(!is.null(args$matrix))args$matrix else args[[1]]
  args$na_col <- pfos_palette("neutral")[["missing"]]
  if (tolower(args$name %||% "") %in% c("z-score", "row z-score", "column z-score")) {
    bounds <- unlist(pfos_profile$scales$row_zscore)
    args$col <- heatmap_diverging_col_fun(breaks=c(bounds[1],0,bounds[2]))
    args$heatmap_legend_param <- modifyList(args$heatmap_legend_param, list(
      at=c(bounds[1],0,bounds[2]), labels=c(paste0("<=",bounds[1]),"0",paste0(">=",bounds[2]))))
  }
  row_labels<-args$row_labels %||% rownames(matrix)
  if(length(row_labels)) {
    args$row_labels<-pfos_heatmap_labels(row_labels)
    args$row_names_max_width<-ComplexHeatmap::max_text_width(args$row_labels,gp=args$row_names_gp)+grid::unit(5,"mm")
  }
  if(is.null(args$width))args$width<-grid::unit(if(any(nchar(row_labels)>35))75 else 100,"mm")
  column_labels<-args$column_labels %||% colnames(matrix)
  if(length(column_labels))args$column_labels<-pfos_heatmap_labels(column_labels)
  do.call(ComplexHeatmap::Heatmap,args)
}

pfos_annotation <- function(...,which="column") {
  args<-list(...)
  args$annotation_name_gp<-grid::gpar(fontsize=8,fontfamily="Liberation Sans")
  args$annotation_legend_param<-modifyList(args$annotation_legend_param %||% list(),list(
    labels_gp=grid::gpar(fontsize=7,fontfamily="Liberation Sans"),
    title_gp=grid::gpar(fontsize=8,fontfamily="Liberation Sans")))
  args$which<-which
  do.call(ComplexHeatmap::HeatmapAnnotation,args)
}
pfos_row_annotation <- function(...)pfos_annotation(...,which="row")

pfos_draw_heatmap <- function(object,...) {
  args<-list(...)
  if(is.null(args$newpage))args$newpage<-FALSE
  if(is.null(args$heatmap_legend_side))args$heatmap_legend_side<-"bottom"
  if(is.null(args$annotation_legend_side))args$annotation_legend_side<-"bottom"
  if(is.null(args$merge_legend))args$merge_legend<-TRUE
  if(is.null(args$padding))args$padding<-grid::unit(c(5,8,5,5),"mm")
  do.call(ComplexHeatmap::draw,c(list(object),args))
}

pfos_pheatmap <- function(...) {
  args<-list(...)
  args$fontsize<-8;args$fontsize_row<-7;args$fontsize_col<-7
  args$fontfamily<-"Liberation Sans"
  # 该入口用于 Top-DE 行标准化表达，沿用原矩阵与聚类。
  bounds <- unlist(pfos_profile$scales$row_zscore)
  args$color <- grDevices::colorRampPalette(unname(pfos_palette("diverging")))(100)
  args$breaks <- seq(bounds[1],bounds[2],length.out=101)
  args$legend_breaks <- c(bounds[1],0,bounds[2])
  args$legend_labels <- c(paste0("<=",bounds[1]),"0",paste0(">=",bounds[2]))
  args$na_col <- pfos_palette("neutral")[["missing"]]
  do.call(pheatmap::pheatmap,args)
}

heatmap_group_annotation <- function(groups,palette=group_colors,annotation_name="Group",which=c("column","row"),...) {
  which<-match.arg(which)
  args<-c(setNames(list(as.character(groups)),annotation_name),list(
    col=setNames(list(heatmap_group_palette(groups,palette)),annotation_name),
    annotation_name_gp=grid::gpar(fontsize=8,fontfamily="Liberation Sans"),
    annotation_legend_param=list(labels_gp=grid::gpar(fontsize=7,fontfamily="Liberation Sans"),
                                 title_gp=grid::gpar(fontsize=8,fontfamily="Liberation Sans"))),list(...))
  do.call(if(which=="row")ComplexHeatmap::rowAnnotation else ComplexHeatmap::HeatmapAnnotation,args)
}
