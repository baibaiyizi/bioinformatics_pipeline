#!/usr/bin/env Rscript
# 复用成熟跨情境绘图块；冻结 overlap / ORA 数值，不运行其检验与覆盖构建块。
a<-commandArgs(TRUE);stopifnot(length(a)%in%c(2L,3L));scope<-if(length(a)==3L)a[3] else "all";stopifnot(scope%in%c("all","loci","heatmaps","correlation"))
root<-normalizePath(a[1],mustWork=TRUE);pair<-a[2]
Sys.setenv(CROSS_WGBS_ROOT_DIR=file.path(root,"scripts/cross_wgbs",pair))
source(file.path(root,"scripts/joint_plot_helpers.R"))
source(file.path(root,"scripts/figure_style.R"))
j<-joint_init(root,file.path("result/cross_context/wgbs",pair),paste0("wgbs_",pair))
lines<-readLines(file.path(root,"scripts/cross_wgbs",pair,"com.Rmd"),warn=FALSE)
chunks<-list();inside<-FALSE;buf<-character()
for(line in lines) {
 if(grepl('^```\\{r',line)){inside<-TRUE;buf<-character();next}
 if(inside&&grepl('^```',line)){chunks[[length(chunks)+1L]]<-paste(buf,collapse="\n");inside<-FALSE;next}
 if(inside)buf<-c(buf,line)
}
e<-new.env(parent=globalenv())
# 配置、纯函数和成熟依赖定义。临时 out_dir 不含旧图片。
eval(parse(text=chunks[[1L]]),e)
route<-function(file)joint_figure_path(root,paste0("wgbs_",pair),tools::file_path_sans_ext(basename(file)),"pdf")
e$ggsave<-function(filename,plot=ggplot2::last_plot(),width=8,height=6,...) {
  if(scope=="correlation"&&basename(filename)!="03_Cross_Context_Correlation.pdf")return(invisible(NULL))
  for(i in seq_along(plot$layers))if(inherits(plot$layers[[i]]$geom,c("GeomText","GeomLabel")))plot$layers[[i]]$aes_params$size<-pfos_figure_style$text_pt/ggplot2::.pt
  pfos_save(plot+ggplot2::theme(plot.title=ggplot2::element_blank(),plot.subtitle=ggplot2::element_blank()),route(filename),height_mm=min(240,max(110,183*height/width)))
}
e$cairo_pdf<-function(filename,width=8,height=6,...)wgbs_cairo_pdf(route(filename),width=width,height=height,...)
e$write_csv<-function(x,file,...)invisible(joint_cached(j,tools::file_path_sans_ext(basename(file)),"csv"))
e$merged_dmr<-readr::read_csv(joint_cached(j,"01_overlap_details","csv"),show_col_types=FALSE)
e$draw_pheatmap_pdf<-function(mat,out_pdf,annotation_row=NULL,annotation_colors=NULL,width=8,height=10,main="",scale="row",show_rownames=TRUE,cluster_cols=FALSE,cluster_rows=nrow(mat)>1,fontsize_row=NULL){
 display_label<-function(x) {labels<-c(Larger_in_source2="Larger in #2",Smaller_in_source2="Smaller in #2",Opposite_direction="Opposite",Common_Hyper="Shared hyper",Common_Hypo="Shared hypo",EffectMagnitude="Effect size",EffectState="Effect size",Upstream.Promoter="Upstream");y<-unname(labels[x]);gsub("_"," ",ifelse(is.na(y),x,y),fixed=TRUE)}
 if(!is.null(annotation_row)) {for(n in names(annotation_row))if(is.character(annotation_row[[n]])||is.factor(annotation_row[[n]])) {old<-as.character(annotation_row[[n]]);if(!is.null(annotation_colors[[n]]))annotation_colors[[n]]<-annotation_colors[[n]][intersect(names(annotation_colors[[n]]),unique(old))];annotation_row[[n]]<-display_label(old);if(!is.null(annotation_colors[[n]]))names(annotation_colors[[n]])<-display_label(names(annotation_colors[[n]]))}}
 if(!is.null(annotation_colors))names(annotation_colors)<-display_label(names(annotation_colors))
 h<-pheatmap::pheatmap(mat,scale=scale,show_rownames=show_rownames,cluster_cols=cluster_cols,cluster_rows=cluster_rows,
   annotation_row=if(is.null(annotation_row))NULL else {a<-annotation_row; names(a)<-display_label(names(a));a},annotation_colors=annotation_colors,color=pfos_figure_style$diverging,
   fontsize=8,fontsize_row=7,fontsize_col=7,fontfamily=pfos_figure_style$family,silent=TRUE,main="",
   labels_row=stringr::str_wrap(rownames(mat),24))
 pfos_save(h$gtable,route(out_pdf),height_mm=min(240,max(150,nrow(mat)*3+55)))
}
e$run_enrichment<-function(entrez_ids,title_prefix,out_prefix){
 for(collection in c("GO_BP","KEGG")) {
   z<-readr::read_csv(joint_cached(j,paste0("06_Table_",collection,"_",out_prefix),"csv"),show_col_types=FALSE)
   top<-z %>% dplyr::filter(Count>0,pvalue<e$integration_config$p_cutoff) %>% dplyr::slice_head(n=15L)
   if(nrow(top)) {
     p<-ggplot2::ggplot(top,ggplot2::aes(FoldEnrichment,reorder(stringr::str_wrap(Description,45),FoldEnrichment),size=Count,color=-log10(pvalue)))+
       ggplot2::geom_point()+pfos_theme()+ggplot2::labs(x="Fold enrichment",y=NULL,color="-log10(P)")
     pfos_save(p,joint_figure_path(root,paste0("wgbs_",pair),paste0("06_Enrich_",collection,"_",out_prefix),"pdf"),height_mm=180)
   }
 }
}
for(chunk in chunks[-1L]) {
 if(scope=="loci"&&!grepl("Common-direction 位点浏览图",chunk,fixed=TRUE))next
 if(scope=="heatmaps"&&!grepl("draw_pheatmap_pdf",chunk,fixed=TRUE))next
 if(scope=="correlation"&&!grepl("03_Cross_Context_Correlation.pdf",chunk,fixed=TRUE))next
 if(grepl('source1_df <- read_dmr|result_report|write_csv\\(tibble\\(analysis="cross_context_motif"|j\\$analysis<-paste0',chunk))next
 # Gviz 固定 7/8 pt，不把字号随旧画布缩小。
 chunk<-gsub('cex = 0.7','cex = 1',chunk,fixed=TRUE)
 chunk<-gsub('cex.axis = 0.7','cex.axis = 1',chunk,fixed=TRUE)
 chunk<-gsub('cex.title = 0.55','cex.title = 1',chunk,fixed=TRUE)
 chunk<-gsub('cex.main = 0.8','cex.main = 1',chunk,fixed=TRUE)
 chunk<-gsub('name = locus_row$locus_label','name = "Genes"',chunk,fixed=TRUE)
 eval(parse(text=chunk),e)
}
cat("Frozen cross-context redraw complete:",pair,"\n")
