#!/usr/bin/env Rscript
# 只读已完成数值；使用现有绘图函数和原绘图片段，不 knit 分析 Rmd。
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1]], mustWork = TRUE)
output <- normalizePath(args[[2]], mustWork = TRUE)
rmd <- file.path(root, 'rmd')
Sys.setenv(WGBS_ROOT_DIR = root, WGBS_RMD_DIR = rmd)
setup <- tempfile(fileext = '.R')
knitr::purl(file.path(rmd, '01_setup.Rmd'), output = setup, quiet = TRUE)
source(setup, local = globalenv()); unlink(setup)
setwd(dirname(output))
# Center template has its helper inside scripts/.
canonical <- function(path) file.path(output,sub('^result/','',path))
# All analysis write calls in selected display fragments are suppressed. Existing
# tables are signed and validated by the caller; no tests or tables are rerun.
write_csv <- function(x, file, ...) invisible(file)
write_tsv <- function(x, file, ...) invisible(file)
write.csv <- function(x, file, ...) invisible(file)
read_stage <- function(stage, name) {
  hits <- list.files(file.path(output, stage), recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == name]
  if(length(hits) != 1L) stop('MISSING/AMBIGUOUS frozen source: ', stage, '/', name)
  readr::read_csv(hits, show_col_types = FALSE, name_repair = 'minimal')
}
source_file <- function(stage, name) {
  hits <- list.files(file.path(output, stage), recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == name]
  if(length(hits) != 1L) stop('MISSING/AMBIGUOUS frozen object: ', name)
  hits[[1]]
}
ggsave <- function(filename, plot = ggplot2::last_plot(), width = 8, height = 6, ..., coordinate_system = 'axes') {
  target <- if(startsWith(filename, output)) filename else canonical(filename)
  pfos_save(wgbs_standardize_plot(plot, coordinate_system), target, height_mm = max(80, 183 * height / width))
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
old_standardize<-wgbs_standardize_plot
wgbs_standardize_plot<-function(plot,coordinate_system='axes')old_standardize(wgbs_text_layers(plot),coordinate_system)
# 本地复制原函数；仅固定 pheatmap 的字号，聚类和矩阵保持原算法。
style_heatmap_call <- function(e) {
 if(is.call(e)) {
  if(grepl('pheatmap$',paste(deparse(e[[1]]),collapse='')))e$fontsize<-7
  for(j in seq_along(e)[-1])if(!identical(e[[j]],quote(expr=)))e[[j]]<-style_heatmap_call(e[[j]])
 }
 e
}
body(draw_pheatmap_pdf)<-style_heatmap_call(body(draw_pheatmap_pdf))
pheatmap <- function(...) {
  dots<-list(...);dots$fontsize<-7;dots$fontsize_row<-7;dots$fontsize_col<-7;dots$fontsize_number<-7
  do.call(pheatmap::pheatmap,dots)
}
plotTracks <- function(...) {
  dots<-list(...)
  dots$fontsize<-8;dots$fontsize.group<-7;dots$cex.group<-1;dots$fontfamily.group<-'Liberation Sans';dots$cex<-1;dots$cex.axis<-7/8;dots$cex.title<-1;dots$cex.main<-1;dots$fontfamily<-'Liberation Sans'
  dots$title.width<-2.2
  dots$margin<-6/25.4 * (grDevices::dev.size('px')[1]/grDevices::dev.size('in')[1])
  if(length(dots) && is.list(dots[[1]]))for(j in seq_along(dots[[1]]))Gviz::displayPars(dots[[1]][[j]])<-list(rotation.title=if(methods::is(dots[[1]][[j]],'DataTrack'))90 else 0)
  do.call(Gviz::plotTracks,dots)
}
old_cairo_pdf <- cairo_pdf
cairo_pdf <- function(filename, width = 8, height = 6, ...) {
  old_cairo_pdf(if(startsWith(filename, output)) filename else canonical(filename), width = width, height = height, ...)
  par(cex.axis=7/8,cex.main=1,cex.lab=1,family='Liberation Sans')
}
# Do not let helpers resolve inputs against old result addresses.
dmr_anno_df <- standardize_dmr_result_columns(read_stage('03_dmr','dmr_annotated.csv'))
dmrs <- standardize_dmr_result_columns(read_stage('03_dmr','dmr_raw.csv'))
if(!identical(unique(dmr_anno_df$group1_label), config$group1_name) ||
   !identical(unique(dmr_anno_df$group2_label), config$group2_name)) stop('保存 DMR 对比方向与 reads_info 不一致')
dmr_gr <- GRanges(dmr_anno_df$seqnames, IRanges(dmr_anno_df$start, dmr_anno_df$end))
mcols(dmr_gr) <- DataFrame(dmr_anno_df[,setdiff(names(dmr_anno_df),c('seqnames','start','end','width','strand'))])
gene_methylation_long <- read_stage('08_gene_methylation','gene_methylation_long.csv')
ensure_dmr_outputs_loaded <- function() invisible(TRUE)
ensure_gene_methylation_loaded <- function() invisible(TRUE)
ensure_dmr_cgi_context <- function(required = FALSE, persist = FALSE) {
  resource <- get_cgi_context_resources(config, target_gr = dmr_gr, required = required)
  if(!is.null(resource)) assign('gr_all_cgi_context',resource$gr_all_cgi_context,envir=.GlobalEnv)
  'CGI_context' %in% names(dmr_anno_df)
}
# Existing data sample identities must be checked before the permitted display
# preparation. assign_sample_ids_to_obj itself validates counts.
# 本样式补绘只读取已冻结展示矩阵/表和轨道，不恢复大型原始对象。
wgbs_original_ensure_obj_cpg_loaded <- ensure_obj_cpg_loaded
ensure_obj_cpg_loaded <- function() {
wgbs_original_ensure_obj_cpg_loaded()
expected <- read_stage('02_qc', 'sample_cpg_methylation_support.csv')
actual_records <- vapply(obj.cpg, function(x) nrow(methylKit::getData(x)), integer(1))
expected<-expected[match(sample_ids,expected$Sample),,drop=FALSE]
if(!identical(as.character(expected$Sample),as.character(sample_ids)) ||
   !identical(as.integer(expected$n_cpg_records),unname(actual_records))) stop('CpG 样本或记录数与原图来源不一致')

  invisible(TRUE)
}

# Existing source snippets contain functions with inference in other sections.
# Read only explicit display ranges and reject statistical/upstream calls.
code <- function(name) {
  x <- readLines(file.path(rmd,name),warn=FALSE)
  x <- x[!grepl('^```',x)]
  if(name=='04_dmroverview.Rmd')x<-gsub('cowplot::plot_grid(p_chr_count, p_chr_density, p_chr_ratio', 'cowplot::plot_grid(pfos_style_plot(p_chr_count), pfos_style_plot(p_chr_density), pfos_style_plot(p_chr_ratio)',x,fixed=TRUE)
  if(name=='09_browser.Rmd') {
    if(!any(grepl('rotation.title = 0',x,fixed=TRUE)))x<-gsub('collapseTranscripts = "meta"','collapseTranscripts = "meta", rotation.title = 0',x,fixed=TRUE)
    x<-gsub('name = paste0(gene_row$SYMBOL, " | abs(delta)="','name = paste0(gene_row$SYMBOL, "\\nabs(delta)="',x,fixed=TRUE)
    x<-gsub('cairo_pdf(combined_pdf, width = 12, height = 8)','cairo_pdf(combined_pdf, width = 12, height = 12)',x,fixed=TRUE)
  }
  if(name=='04_dmroverview.Rmd') {
    i<-grep('  ensure_obj_cpg_loaded()',x,fixed=TRUE)[1]
    j<-grep('  group1_col <-',x,fixed=TRUE)[1]
    stopifnot(!is.na(i),!is.na(j),j>i)
    x<-c(x[seq_len(i-1)],'  chr_meth_diff <- read_stage("04_dmroverview","chromosome_methylation_difference.csv")',x[j:length(x)])
  }
  x <- vapply(x, function(line) {
    hit <- gregexpr('["\']result/[^"\']+["\']',line,perl=TRUE)
    z <- regmatches(line,hit)[[1]]
    if(length(z))regmatches(line,hit)<-list(vapply(z,function(v)paste0(substr(v,1,1),canonical(substr(v,2,nchar(v)-1)),substr(v,1,1)),character(1)))
    line
  },character(1))
  x
}
fragment <- function(name, start, end = NULL) {
  x <- code(name); i <- grep(start,x,fixed=TRUE)[1]
  if(is.na(i)) stop('绘图源码锚点缺失: ',name,' ',start)
  j <- if(is.null(end)) length(x)+1L else { z <- which(seq_along(x)>i & grepl(end,x,fixed=TRUE));if(!length(z))stop('绘图源码结束锚点缺失');z[1] }
  expr <- parse(text=x[i:(j-1L)])
  forbidden <- c('DMLtest','callDML','callDMR','enrichGO','enrichKEGG','gseGO','gseKEGG','fgsea','wilcox.test','fisher.test','t.test','calculateDiffMeth','segmentPMDs','segmentUMRsLMRs','BSmooth','export.bw')
  inspect <- function(e) {
    if(is.call(e)) {
      head <- as.character(e[[1]])
      if(any(head %in% forbidden))stop('禁止在纯绘图片段中调用统计: ',head)
      lapply(as.list(e)[-1],inspect)
    }
  }
  lapply(expr,inspect);eval(expr,envir=.GlobalEnv)
}
selected <- Sys.getenv('PFOS_PLOT_STEPS')
steps <- if(!nzchar(selected)) c(2:10,13L,14L) else unique(unlist(lapply(strsplit(gsub(',',' ',selected),' +')[[1]],function(x) {
  z<-as.integer(strsplit(x,'-',fixed=TRUE)[[1]]);if(length(z)==2)seq(z[1],z[2])else z
})))
for(step in steps) {
  cat('PLOTS_ONLY stage',step,'\n');flush.console()
  if(step == 2L) {
    # Only the old QC inference block is excluded. No whole Rmd is rendered.
    x <- code('02_qc.Rmd')
    start <- grep('  if(length(unique(global_meth$Group)) == 2)',x,fixed=TRUE)[1]
    end <- grep('  # 7 染色体级别平均甲基化热图',x,fixed=TRUE)[1]
    if(is.na(start)||is.na(end))stop('QC inference boundary not found')
    x <- x[-seq(start,end-1L)]
    if(any(grepl('wilcox.test',x,fixed=TRUE)))stop('QC inference remains')
    # Result reads must consume frozen copied tables, never a changed upstream QC.
    x <- gsub('file.path(config$work_dir, "01rawdata/qc_results_GRCm39/qc_read_alignment_summary.csv")','source_file("02_qc", "qc_read_alignment_summary.csv")',x,fixed=TRUE)
    x <- gsub('file.path(config$work_dir, "04methykit/cpg_coverage_threshold_summary.csv")','file.path(root, "__not_used_frozen_qc__")',x,fixed=TRUE)
    # Existing relative output paths remain read-only; resolve saved numeric sources.
    x <- gsub('qc_read_alignment_file <- file.path(qc_table_dir, "qc_read_alignment_summary.csv")','qc_read_alignment_file <- source_file("02_qc", "qc_read_alignment_summary.csv")',x,fixed=TRUE)
    x <- gsub('qc_coverage_file <- file.path(qc_table_dir, "cpg_coverage_threshold_summary.csv")','qc_coverage_file <- tryCatch(source_file("02_qc", "cpg_coverage_threshold_summary.csv"), error=function(e) file.path(root,"__absent__"))',x,fixed=TRUE)
    eval(parse(text=x),globalenv())
    # Persist only new display coordinates/matrices, never replace scientific tables.
    if(exists('meth')) {
      matrix <- methylKit::percMethylation(meth)
      out <- canonical('result/02_qc/01_group_comparison/tables/display_strict_common_cpg_matrix.rds')
      dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE);saveRDS(matrix,out)
    }
  }
  if(step==14L) {
    ora<-read_stage('05_enrichment','enrichment_ora_long.csv')
    tab<-ora[which(ora$is_selected_significant & ora$gene_set_source=='GO' & ora$ontology=='BP'),,drop=FALSE]
    tab<-tab[order(tab$p.adjust),,drop=FALSE];tab<-head(tab,15)
    if(nrow(tab)) {
      sets<-strsplit(tab$geneID,'/');names(sets)<-ifelse(nchar(tab$Description)>30,paste0(substr(tab$Description,1,30),'...'),tab$Description)
      u<-UpSetR::upset(UpSetR::fromList(sets),nsets=length(sets),nintersects=20,order.by='freq',sets.bar.color=viz_pal$upset_sets,main.bar.color=viz_pal$upset_main,number.colors=viz_pal$neutral_dark,text.scale=c(8/8.3,1,8/8.3,1,1,7/(2.2*ggplot2::.pt)))
      u$Sizes$widths[1]<-u$Sizes$widths[1]+grid::unit(12,'pt')
      pfos_save(function(){print(u);invisible(NULL)},canonical('result/05_enrichment/01_group_comparison/plots/go_bp_upset.pdf'),150)
    }
  }
  if(step==13L) {
    qc_plot_dir<-'result/02_qc/01_group_comparison/plots'
    qc_read_alignment<-read_stage('02_qc','qc_read_alignment_summary.csv')
    fragment('02_qc.Rmd','  if(exists("qc_read_alignment") && "bismark_duplicate_rate_pct"')
    mat<-as.data.frame(readRDS(source_file('02_qc','display_strict_common_cpg_matrix.rds'))/100)
    corr<-methylKit:::.plotCorrelation
    styled_pairs<-function(...)graphics::pairs(...,cex.labels=1/.66)
    env<-new.env(parent=environment(corr));env$pairs<-styled_pairs;environment(corr)<-env
    cairo_pdf(file.path(qc_plot_dir,'sample_correlation.pdf'),width=8,height=6)
    par(cex.axis=7/(8*.66));corr(mat,title='CpG base pearson cor.',method='pearson');dev.off()
  }
  if(step==12L && file.exists(file.path(rmd,'10_sperm.Rmd'))) {
    sperm_plot_dir<-'result/10_sperm/01_group_comparison/plots'
    te_sample_mean<-read_stage('10_sperm','sperm_te_sample_methylation.csv')
    fragment('10_sperm.Rmd','save_sperm_tile_heatmap <-','ensure_dmr_outputs_loaded()')
    fragment('10_sperm.Rmd','      if(nrow(te_sample_mean) > 0) {','      # 主要转座子类别甲基化比较')
  }
  if(step==11L) {
    # 使用原 methylKit 展示函数，原 strict-unite 矩阵的样本顺序及 50% SD 过滤不变。
    mat<-as.data.frame(readRDS(source_file('02_qc','display_strict_common_cpg_matrix.rds'))/100)
    stopifnot(identical(names(mat),as.character(sample_ids)))
    treat<-ifelse(sample_ids %in% sample_ids[config$group2_idx],1,0)
    sd<-methylKit:::rowSds(as.matrix(mat));filtered<-mat[sd>quantile(sd,.5),,drop=FALSE]
    qp<-'result/02_qc/01_group_comparison/plots';qc_plot_dir<-qp;qc_table_dir<-'result/02_qc/01_group_comparison/tables'
    cairo_pdf(file.path(qp,'sample_correlation.pdf'),width=8,height=6)
    methylKit:::.plotCorrelation(meth.mat=mat,title='CpG base pearson cor.',method='pearson');dev.off()
    ensure_meth_loaded()
    cairo_pdf(file.path(qp,'sample_clustering.pdf'),width=8,height=6)
    wgbs_plot_cluster(meth);dev.off()
    for(scree in c(TRUE,FALSE)) wgbs_plot_pca(meth,canonical(file.path(qp,if(scree)'pca_screeplot.pdf'else'pca_sample_layout.pdf')),screeplot=scree)
    final_df_meta<-read_stage('06_genomecontext','metagene_profile_bin_support.csv');feature_order_meta<-profile_segment_levels
    fragment('06_genomecontext.Rmd','feature_bg_meta <-','# ==============================================================================')
    cgi_final_df<-read_stage('06_genomecontext','cgi_profile_bin_support.csv');region_order_cgi<-c('N_Shelf','N_Shore','Island','S_Shore','S_Shelf')
    fragment('06_genomecontext.Rmd','  feature_bg_cgi <-','  # 4. 包含 Others')
    stat_df<-read_stage('06_genomecontext','global_cgi_stats.csv');methyl_group_cols<-c(config$group1_name,config$group2_name)
    fragment('06_genomecontext.Rmd','  # 绘制柱状图','  # 5. 为后续分析')
    pairwise_cor_df<-read_stage('02_qc','sample_pairwise_methylation_correlation.csv')
    fragment('02_qc.Rmd','    pairwise_plot_df <-','  }')
    # 原文库对象只用于原分布图及逐样本总量，绝不 unite 或建模。
    load(file.path(root,config$input_rdata))
    obj.cpg<-assign_sample_ids_to_obj(obj.cpg,sample_ids)
    expected<-read_stage('02_qc','sample_cpg_methylation_support.csv');expected<-expected[match(sample_ids,expected$Sample),]
    stopifnot(identical(as.integer(expected$n_cpg_records),unname(vapply(obj.cpg,function(x)nrow(methylKit::getData(x)),integer(1)))))
    fragment('02_qc.Rmd','  for(i in seq_along(obj.cpg)) {','  # 5 相关性与 PCA')
    fragment('02_qc.Rmd','  global_meth <- data.frame(','  if(length(unique(global_meth$Group)) == 2)')
    chr_meth<-dplyr::bind_rows(lapply(seq_along(obj.cpg),function(i) {
      d<-methylKit::getData(obj.cpg[[i]])
      d<-data.frame(chr_clean=strip_chr_prefix(d$chr),Methylation=100*d$numCs/(d$numCs+d$numTs))
      d%>%filter(chr_clean %in% c(as.character(1:19),'X','Y'))%>%group_by(chr_clean)%>%summarise(mean_meth=mean(Methylation,na.rm=TRUE),.groups='drop')%>%mutate(Sample=sample_ids[i])
    }))
    fragment('02_qc.Rmd','  if(nrow(chr_meth) > 0) {','  # 8 CpG 密度分层甲基化比较')
    fragment('02_qc.Rmd','  # 9 BS 方向性检查可视化')
  }
  if(step == 3L) {
    dml_df <- readRDS(source_file('03_dmr','dml_results.rds'))
    dml_df <- standardize_dmr_result_columns(dml_df)
    dml_pval_col <- first_existing_col(dml_df,c('pval','pvalue'))
    dmr_plot_dir <- 'result/03_dmr/01_calling/plots'
    fragment('03_dmr.Rmd','# DML 火山图')
  }
  if(step == 4L) {
    dmr_overview_plot_dir <- 'result/04_dmroverview/01_group_comparison/plots'
    dmr_overview_table_dir <- 'result/04_dmroverview/01_group_comparison/tables'
    fragment('04_dmroverview.Rmd','# 1 甲基化水平散点图','promoter_dmrs <-')
  }
  if(step == 5L) {
    ora <- read_stage('05_enrichment','enrichment_ora_long.csv')
    ora <- ora[which(ora$is_selected_significant),,drop=FALSE]
    gene_df_final <- dmr_anno_df %>% group_by(gene_id_clean) %>% slice_max(abs(diffMethy),n=1,with_ties=FALSE) %>% ungroup()
    geneList <- gene_df_final$diffMethy;names(geneList)<-gene_df_final$SYMBOL;geneList<-sort(geneList,decreasing=TRUE)
    make_ora <- function(tab,ontology) {
      if(!nrow(tab))return(NULL)
      tab<-as.data.frame(tab);rownames(tab)<-tab$ID
      obj<-methods::new('enrichResult',result=as.data.frame(tab),pvalueCutoff=1,pAdjustMethod='BH',qvalueCutoff=1,organism='Mus musculus',ontology=ontology,keytype='SYMBOL',gene=unique(unlist(strsplit(tab$geneID,'/'))),universe=character(),geneSets=setNames(strsplit(tab$geneID,'/'),tab$ID),readable=TRUE)
      if(ontology!='KEGG')obj<-enrichplot::pairwise_termsim(obj)
      obj
    }
    for(ont in c('KEGG','BP','MF','CC')) {
      tab<-if(ont=='KEGG')ora[ora$gene_set_source=='KEGG',]else ora[ora$gene_set_source=='GO'&ora$ontology==ont,]
      obj<-make_ora(tab,ont);if(is.null(obj))next
      prefix<-if(ont=='KEGG')'kegg'else paste0('go_',tolower(ont))
      pp<-paste0('result/05_enrichment/01_group_comparison/plots/',prefix)
      save_ora_barplot(as.data.frame(obj),paste0(pp,'_barplot.pdf'),ont,top_n=15)
      save_ora_dotplot(as.data.frame(obj),paste0(pp,'_dotplot.pdf'),ont,top_n=15)
      save_enrichment_cnetplot(obj,geneList,paste0(pp,'_cnetplot.pdf'),ont,showCategory=5)
      if(ont!='KEGG')save_enrichment_treeplot(obj,paste0(pp,'_treeplot.pdf'),ont,showCategory=15)
      if(ont=='BP') {
        gp<-build_goplot_inputs(as.data.frame(obj),geneList,category_label='BP',top_n=10)
        save_go_goplot_bundle(gp,paste0(pp,'_chord.pdf'),paste0(pp,'_cluster.pdf'),title_suffix='BP')
        tab<-as.data.frame(obj)%>%arrange(p.adjust)%>%head(15)
        sets<-setNames(strsplit(tab$geneID,'/'),ifelse(nchar(tab$Description)>30,paste0(substr(tab$Description,1,30),'...'),tab$Description))
        u<-UpSetR::upset(UpSetR::fromList(sets),nsets=length(sets),nintersects=20,order.by='freq',sets.bar.color=viz_pal$upset_sets,main.bar.color=viz_pal$upset_main,number.colors=viz_pal$neutral_dark,text.scale=c(8/8.3,1,8/8.3,1,1,7/(2.2*ggplot2::.pt)))
        u$Sizes$widths[1]<-u$Sizes$widths[1]+grid::unit(12,'pt')
        pfos_save(function(){print(u);invisible(NULL)},canonical(paste0(pp,'_upset.pdf')),height_mm=150)
      }
    }
  }
  if(step == 6L) {
    final_df_meta <- read_stage('06_genomecontext','metagene_profile_bin_support.csv');feature_order_meta<-profile_segment_levels
    fragment('06_genomecontext.Rmd','feature_bg_meta <-','# ==============================================================================')
    cgi_final_df<-read_stage('06_genomecontext','cgi_profile_bin_support.csv');region_order_cgi<-c('N_Shelf','N_Shore','Island','S_Shore','S_Shelf')
    fragment('06_genomecontext.Rmd','  feature_bg_cgi <-','  # 4. 包含 Others')
    stat_df<-read_stage('06_genomecontext','global_cgi_stats.csv');methyl_group_cols<-c(config$group1_name,config$group2_name)
    fragment('06_genomecontext.Rmd','  # 绘制柱状图','  # 5. 为后续分析')
    # Violin data are deterministic annotation of the same frozen CpGs.
    ensure_meth_loaded()
    if(!exists('meth_raw_dt'))meth_raw_dt<-build_sample_methylation_dt(obj.cpg,sample_ids,config)
    meth_gr_sites<-GRanges(meth_raw_dt$chr,IRanges(meth_raw_dt$start,meth_raw_dt$end));mcols(meth_gr_sites)<-DataFrame(meth_raw_dt[,c('Sample','Group','Methylation')])
    gr_all_cgi_context<-get_cgi_context_resources(config,target_gr=meth_gr_sites,required=TRUE)$gr_all_cgi_context
    fragment('06_genomecontext.Rmd','  # 6. Violin + jitter','  message("✅ CGI 分析完成。")')
    if(any(basename(list.files(output,recursive=TRUE))=='segmentation_summary.csv')) {
      seg_summary_df<-read_stage('06_genomecontext','segmentation_summary.csv')
      fragment('06_genomecontext.Rmd','      p_seg_summary <-','    }')
      # Model diagnostic PDFs require unfrozen fitted models and are preserved.
      # Final segmentation is redrawn from the exact saved segment coordinates.
      meth_raw_df<-methylKit::getData(meth)
      for(grp in c(config$group1_name,config$group2_name)) {
        seg<-read_stage('06_genomecontext',paste0(tolower(grp),'_umr_lmr_segments.csv'))
        pmd<-read_stage('06_genomecontext',paste0(tolower(grp),'_pmd_segments.csv'))
        seggr<-makeGRangesFromDataFrame(as.data.frame(seg),keep.extra.columns=TRUE)
        pmdgr<-makeGRangesFromDataFrame(as.data.frame(pmd),keep.extra.columns=TRUE)
        ix<-if(grp==config$group1_name)config$group1_idx else config$group2_idx
        ccols<-grep('^coverage',names(meth_raw_df),value=TRUE);mcols0<-grep('^numCs',names(meth_raw_df),value=TRUE)
        seek<-GRanges(meth_raw_df$chr,IRanges(meth_raw_df$start,meth_raw_df$end),T=rowSums(meth_raw_df[,ccols[ix],drop=FALSE]),M=rowSums(meth_raw_df[,mcols0[ix],drop=FALSE]))
        seek<-filter_gr_to_main_chr(seek);seqlevels(seek)<-standardize_chr_style(seqlevels(seek),TRUE)
        pfos_save(function(){set.seed(1);n<-if(identical(names(grDevices::dev.cur()),'png'))1L else 3L;for(i in seq_len(n))wgbs_plot_final_segmentation(m=seek,segs=seggr,PMDs=pmdgr,meth.cutoff=.5,numRegions=1,pdfFilename=NULL)},canonical(paste0('result/06_genomecontext/02_segmentation/plots/',tolower(grp),'_final_segmentation.pdf')),height_mm=130)
      }
    }
  }
  if(step == 7L) {
    has_dmr_cgi_context<-TRUE
    fragment('07_regulatory.Rmd','# 1 DMR 的 CpG density context 分布与效应','# 2 显著差异')
    motif_enrich_df<-read_stage('07_regulatory','motif_enrichment_top500_matched_background.csv')
    fragment('07_regulatory.Rmd','    p_motif_enrich <-','  }')
    motif_df<-read_stage('07_regulatory','motif_top20.csv')
    fragment('07_regulatory.Rmd','  # 绘图','}')
    fragment('07_regulatory.Rmd','# 4 印记基因','# 5 RepeatMasker')
    repeat_class_summary<-read_stage('07_regulatory','repeatmasker_class_summary.csv');repeat_family_summary<-read_stage('07_regulatory','repeatmasker_family_top20.csv')
    fragment('07_regulatory.Rmd','      p_repeat_class <-','    }')
  }
  if(step == 8L) {
    p<-ggplot(gene_methylation_long,aes(Layer,diffMethy))+geom_boxplot(outlier.shape=NA)+geom_hline(yintercept=0,linetype='dashed')+theme_wgbs_classic()+labs(x='Genomic layer',y=paste('Methylation proportion difference:',group_diff_label))
    ggsave('result/08_gene_methylation/01_group_comparison/plots/gene_methylation_layers.pdf',p,width=8,height=6)
  }
  if(step == 9L) {
    ensure_dmr_cgi_context(required=FALSE,persist=FALSE)
    # Drawing-only browser chunk; track export and BSmooth chunks are excluded.
    fragment('09_browser.Rmd','quiet_gviz_deprecations <-','# 2 BSseq 平滑')
  }
  if(step == 10L && file.exists(file.path(rmd,'10_sperm.Rmd'))) {
    # Sperm topic contains only deterministic region annotation/aggregation and
    # drawing. It regenerates display summaries from the verified original CpGs;
    # original scientific table writes remain suppressed.
    x<-code('10_sperm.Rmd')
    if(any(grepl('test\\(|segmentPMDs|segmentUMRsLMRs|BSmooth',x)))stop('精子显示模块出现统计调用')
    eval(parse(text=x),globalenv())
  }
  # GSEA dotplots consume only the saved full statistics, never fgsea.
  if(step %in% c(4L,5L)) {
    stage<-if(step==4L)'04_dmroverview'else'05_enrichment';prefix<-if(step==4L)'promoter'else'methylation'
    d<-read_stage(stage,paste0(prefix,'_gsea_long.csv'))
    for(src in c('GO','KEGG')) {
      tab<-d[d$gene_set_source==src & d$is_selected_significant,,drop=FALSE]
      if(!nrow(tab))next
      if(src=='GO')tab<-tab[tab$ontology=='BP',,drop=FALSE]
      if(!nrow(tab))next
      save_gsea_dotplot(tab,paste0('result/',stage,'/01_group_comparison/plots/',prefix,'_',if(src=='GO')'GO_BP'else'KEGG','_GSEA_Dotplot.pdf'),prefix)
    }
    if(step==4L) {
      # Reconstruct only the deterministic running-score display. No permutations
      # or enrichment P values are calculated. Its ES and leading-edge members
      # must match the frozen result before the original curve can be replaced.
      top <- d[d$gene_set_source=='GO' & d$ontology=='BP' & d$is_selected_significant,,drop=FALSE]
      if(nrow(top)) {
        top<-top[1,,drop=FALSE]
        promoter_dmrs <- dmr_anno_df %>% filter(!is.na(gene_id_clean),grepl('Promoter',annotation),!is.na(diffMethy)) %>% group_by(gene_id_clean) %>% summarise(diffMethy=mean(diffMethy,na.rm=TRUE),.groups='drop')
        map <- clusterProfiler::bitr(promoter_dmrs$gene_id_clean,fromType='ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb=org.Mm.eg.db)
        mapped <- inner_join(promoter_dmrs,map,by=c('gene_id_clean'='ENSEMBL')) %>% filter(!is.na(ENTREZID)) %>% group_by(ENTREZID) %>% summarise(diffMethy=mean(diffMethy,na.rm=TRUE),.groups='drop')
        ranking <- mapped$diffMethy;names(ranking)<-mapped$ENTREZID;ranking<-sort(ranking,decreasing=TRUE)
        membership <- AnnotationDbi::select(org.Mm.eg.db,keys=names(ranking),keytype='ENTREZID',columns='GOALL')
        members <- unique(membership$ENTREZID[membership$GOALL==top$ID & !is.na(membership$GOALL)])
        score <- DOSE:::gseaScores(geneList=ranking,geneSet=members,exponent=1,fortify=FALSE)
        if(is.finite(score$ES) && abs(score$ES-top$enrichmentScore)<1e-7) {
          obj <- methods::new('gseaResult',result=structure(as.data.frame(top),row.names=top$ID),geneSets=setNames(list(members),top$ID),geneList=ranking,params=list(exponent=1),organism='Mus musculus',setType='GO',keytype='ENTREZID')
          p <- enrichplot::gseaplot2(obj,geneSetID=1,title=top$Description)
          ggsave('result/04_dmroverview/01_group_comparison/plots/promoter_gsea_go_top1.pdf',p,width=8,height=6)
        } else message('UNVERIFIED_DISPLAY_SOURCE: original GO running score does not match frozen ES; preserving top1 curve')
      }
    }
  }
}
cat('纯绘图完成；缺模型/完整成员来源的历史图件保留并由 SHA 清单明确标识。\n')
