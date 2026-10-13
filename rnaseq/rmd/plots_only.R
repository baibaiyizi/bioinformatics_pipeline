#!/usr/bin/env Rscript
# 从冻结的表、原表达矩阵和对象重绘；不运行完整 Rmd 或统计模块。
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("usage: plots_only.R MODULE_ROOT STAGING_RESULT [MODULES]")
working_dir <- normalizePath(args[1], mustWork = TRUE)
result_root_dir <- normalizePath(args[2], mustWork = TRUE)
rmd_dir <- file.path(working_dir, "rmd")
source(file.path(rmd_dir, "figure_style.R"))
setwd(working_dir)
Sys.setenv(RNASEQ_ROOT_DIR = working_dir, RNASEQ_RMD_DIR = rmd_dir,
           RNASEQ_CLEAN_MODULE_RESULT = "false", RNASEQ_CLEAN_MODULE_CACHE = "false")
packages <- c("ggplot2", "ggrepel", "dplyr", "tibble", "tidyr", "readr", "stringr",
              "purrr", "reshape2", "patchwork", "pheatmap", "scales", "circlize",
              "ComplexHeatmap", "clusterProfiler", "enrichplot", "org.Mm.eg.db",
              "AnnotationDbi", "GOplot", "UpSetR", "igraph", "edgeR")
for (pkg in packages) {
  message("[plots-only] dependency ",pkg)
  if (!requireNamespace(pkg, quietly = TRUE)) stop("缺少已声明绘图依赖: ", pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (n in getNamespaceExports("dplyr")) if (n %in% c("filter", "select", "rename", "count", "mutate", "arrange", "summarise", "summarize", "slice")) assign(n, getExportedValue("dplyr", n))

# 只载入定义，绝不执行 setup 中的输入生产、安装或清理。
setup <- parse(file.path(rmd_dir, "01_main_setup.R"))
assigned_name <- function(x) if (is.call(x) && is.symbol(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=") && is.symbol(x[[2]])) as.character(x[[2]]) else ""
is_function_definition <- function(x) nzchar(assigned_name(x)) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))
for (x in setup) if (is_function_definition(x)) eval(x)
source(file.path(rmd_dir,"plot_heatmaps.R"))
# 原setup默认将result/绑定现役目录；纯绘图强制绑定staging。
normalize_project_output_path <- function(path) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)) return(path)
  current <- file.path(working_dir,"result")
  if(grepl("^result(/|$)",path)) path<-file.path(result_root_dir,sub("^result/?","",path))
  else if(startsWith(path,paste0(current,"/"))) path<-file.path(result_root_dir,substring(path,nchar(current)+2L))
  else if(grepl("^\\.cache(/|$)",path)) path<-file.path(working_dir,path)
  normalizePath(path,winslash="/",mustWork=FALSE)
}
contrast_source <- utils::read.csv(file.path(result_root_dir,"05_de/tables/5_0_contrast_manifest.csv"))
contrast_groups <- strsplit(contrast_source$analysis_contrast[1]," - ",fixed=TRUE)[[1]]
if(length(contrast_groups)!=2L)stop("冻结比较方向无效")
group2_label <- contrast_groups[1];group1_label <- contrast_groups[2]
safe_globals <- c("figure_significance_metric", "rnaseq_random_seed", "rnaseq_figure_width_mm",
  "rnaseq_figure_min_height_mm", "rnaseq_figure_max_height_mm", "rnaseq_figure_dpi",
  "rnaseq_figure_font_family", "rnaseq_figure_base_size", "rnaseq_figure_tick_size",
  "rnaseq_figure_text_size", "rnaseq_figure_tick_text_size", "DIRECTION_LEVELS", "DIRECTION_LABELS",
  "PLOT_FILE_EXTENSIONS", "RESULT_LAYOUT_REGISTRY", "PRIMARY_CHROM_REGEX", "FIG_PNG_DPI")
for (x in setup) if (assigned_name(x) %in% safe_globals) eval(x)
rnaseq_species <- "mouse"; annotation_orgdb <- org.Mm.eg.db
orgdb_package <- "org.Mm.eg.db";kegg_organism_code<-"mmu";reactome_organism<-"mouse"
frozen_context<-utils::read.csv(file.path(result_root_dir,"01_setup/tables/01_analysis_context.csv"),check.names=FALSE)
context_values<-setNames(frozen_context$value,frozen_context$metric)
de_pvalue_cutoff <- as.numeric(context_values[["de_pvalue_cutoff"]])
de_log2fc_cutoff <- as.numeric(context_values[["de_log2fc_cutoff"]])
de_significance_column <- context_values[["de_significance_column"]]
de_use_padj <- identical(de_significance_column,"padj")
de_significance_label <- if(de_use_padj)"adjusted p-value" else "nominal p-value"
if(any(!is.finite(c(de_pvalue_cutoff,de_log2fc_cutoff))))stop("冻结DE阈值缺失")
for (x in setup) {
  n <- assigned_name(x)
  if (n %in% c("rna_group_palette", "group1_color", "group2_color", "rna_de_palette", "rna_sig_palette",
      "rna_neutral_palette", "rna_accent_palette", "rna_red_yellow_blue_palette", "rna_red_yellow_blue_binary_palette",
      "rna_blue_red_mid", "rna_sashimi_color_string", "viz_pal", "group_colors", "viz_style")) eval(x)
}
for (n in names(RESULT_LAYOUT_REGISTRY)) assign(paste0("result_", substr(n,1,2), "_dir"), file.path(result_root_dir,n))
cache_root_dir <- file.path(working_dir, ".cache/rmd")
cache_names <- c(cache_06_tables_path="06_main_enrichment_tables.rds",cache_06_objects_path="06_main_enrichment_objects.rds",
  cache_06_human_tables_path="06_human_ortholog_enrichment_tables.rds",cache_09_gsva_state_path="09_gsva_state.rds",
  cache_09_human_gsva_state_path="09_human_ortholog_gsva_state.rds",cache_10_tf_state_path="10_tf_regulator_state.rds",
  cache_11_ppi_state_path="11_ppi_state.rds",cache_12_wgcna_path="12_wgcna_state.rds",
  cache_15_splicing_state_path="15_splicing_state.rds",cache_16_immune_state_path="16_immune_state.rds",
  cache_mouse_native_gene_sets_path="mouse_native_gene_sets.rds",cache_human_ortholog_gene_sets_path="human_ortholog_gene_sets.rds")
for (n in names(cache_names)) assign(n, file.path(cache_root_dir, cache_names[[n]]))
aliases <- c(step1_input_dir="02",step1_overview_dir="02",step2_qc_dir="03",step3_pca_dir="04",step4_de_dir="05",
  step5_function_dir="06",step6_go_dir="07",step7_gsea_dir="08",step8_wgcna_dir="12",step9_circle_dir="18",
  step10_root_dir="15",step11_dir="13",step12_dir="16",step13_dir="14",step14_dir="17",step15_dir="10",step16_dir="09",step17_dir="11")
for (n in names(aliases)) assign(n,get(paste0("result_",aliases[[n]],"_dir")))
cleaned_module_result_steps <- new.env(parent = emptyenv()); MODULE_CACHE_REGISTRY <- list()
plot_audit_records <- list(); rendered <- list(); drawing_sources <- character()

# 所有统计写入和缓存写入都禁用。原数值结果完整保留于 staging。
write.csv <- write.table <- saveRDS <- write_csv_safe <- write_csv_compat <- function(...) invisible(NULL)
record_plot_status <- flush_plot_audit <- function(...) invisible(NULL)
install_and_load <- function(pkgs, ...) {
  for (p in pkgs) {if (!requireNamespace(p,quietly=TRUE)) stop("缺少依赖: ",p); suppressPackageStartupMessages(library(p,character.only=TRUE))}
}
try_install_and_load <- function(pkgs, ...) {install_and_load(pkgs);list(ok=TRUE,loaded=pkgs,missing=character())}
for (n in c("run_limma_contrast_table", "run_ssgsea_matrix", "run_fgsea_bundle", "run_gsva_score_matrix")) assign(n,function(...) stop("纯绘图禁止统计重算"))
for (n in c("fisher.test","enrichGO","enrichKEGG","gseGO","gseKEGG","gsva","DESeq")) assign(n,function(...) stop("纯绘图禁止统计重算"))

read_source_csv <- function(relative) {
  p <- file.path(result_root_dir, relative)
  if (!base::file.exists(p)) stop("缺少冻结数值源: ", relative)
  utils::read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
}
read_source_matrix <- function(relative) {
  p <- file.path(result_root_dir,relative)
  if (!base::file.exists(p)) stop("缺少冻结矩阵: ",relative)
  as.matrix(utils::read.delim(p,row.names=1,check.names=FALSE))
}
de_result <- read_source_csv("05_de/tables/5_1_annotated_DE_results.csv")
restore_main_de_result <- function(...) de_result
sig_de <- de_result %>% filter(direction != "ns")
source_manifest <- utils::read.delim(file.path(rmd_dir,"plot_inputs.tsv"),stringsAsFactors=FALSE)
inputs <- setNames(source_manifest$path,source_manifest$key)
read_input_matrix <- function(key) read_matrix(file.path(working_dir,inputs[[key]]))
selected_samples_tbl <- read_source_csv("02_selection/tables/sample_qc_overview.csv")
reads_info <- read_reads_info(file.path(working_dir,"rawdata/reads_info.txt"))
sample_info <- reads_info[selected_samples_tbl$sample,c("group","stage"),drop=FALSE]
gene_exp <- read_input_matrix("selected_expression")
if (!identical(colnames(gene_exp),selected_samples_tbl$sample)) stop("原表达矩阵与原图样本身份不一致")
for(spec in list(c("sampleA","baseMeanA"),c("sampleB","baseMeanB"))) {
  group<-unique(de_result[[spec[1]]]);if(length(group)!=1L)stop("DE分组来源不唯一")
  samples<-rownames(sample_info)[sample_info$group==group]
  observed<-rowMeans(gene_exp[de_result$gene_id,samples,drop=FALSE])
  if(any(abs(observed-de_result[[spec[2]]])>1e-6,na.rm=TRUE))stop("原表达矩阵与DE组均值不一致")
}
gene_count_full <- read_input_matrix("full_counts"); gene_exp_full <- read_input_matrix("full_expression")
gene_count_mat <- gene_count_full[,colnames(gene_exp),drop=FALSE]
gene_logcpm_full <- read_input_matrix("full_logcpm")
gene_logcpm <- gene_logcpm_full[,colnames(gene_exp),drop=FALSE]
full_sample_info <- reads_info[colnames(gene_exp_full),c("group","stage"),drop=FALSE]
sample_name_info <- list(group1=rownames(sample_info)[sample_info$group==group1_label],group2=rownames(sample_info)[sample_info$group==group2_label])
rawdata_dir <- file.path(working_dir,"rawdata")
reference_dir <- normalizePath(file.path(working_dir,"../../../../../reference/GRCm39"))
post_tx_root_dir <- as_root_dir <- file.path(working_dir,"5.post_transcriptional_regulation")
splicing_dir <- file.path(as_root_dir,"splicing"); rmats_dir <- file.path(splicing_dir,"rmats")
exon_usage_dir <- file.path(splicing_dir,"exon_usage"); isoform_results_dir <- file.path(as_root_dir,"isoform_switch/results")
apa_dir <- file.path(as_root_dir,"apa"); der_dir <- file.path(as_root_dir,"der"); integrated_as_dir <- file.path(as_root_dir,"integrated")
selected_samples_path <- file.path(rmd_dir,"plot_samples.tsv")
auto_selection_summary <- utils::read.delim(file.path(working_dir,"result_pca/auto_selection_summary.tsv"),check.names=FALSE)

# 每张图记录它实际使用的表/对象；Rmd 中画图代码仍是唯一实现。
register_draw <- function(primary, data=NULL) {
  if(!startsWith(primary,paste0(result_root_dir,"/")))stop("拒绝在staging之外成图: ",primary)
  rel <- substring(primary,nchar(result_root_dir)+2L)
  src <- ""
  if (is.data.frame(data) && nrow(data)) {
    src <- sub("/plots/","/tables/",sub("\\.pdf$","_plot_source.tsv",rel))
    dir.create(dirname(file.path(result_root_dir,src)),recursive=TRUE,showWarnings=FALSE)
    for(n in names(data)) if(is.list(data[[n]]))data[[n]]<-vapply(data[[n]],function(x)paste(as.character(x),collapse=";"),character(1))
    utils::write.table(data,file.path(result_root_dir,src),sep="\t",row.names=FALSE,quote=TRUE)
  }
  rendered[[length(rendered)+1L]] <<- data.frame(figure_id=sub("\\.pdf$","",rel),pdf=rel,png=sub("\\.pdf$",".png",rel),source_table=src,
    source_files=paste(drawing_sources,collapse=";"),status="COMPLETED",scientific_status=figure_scientific_status(primary),stringsAsFactors=FALSE)
}
save_plot_versions <- function(plot,primary,width=8,height=6,coordinate_system=c("axes","none"),...) {
  coordinate_system<-match.arg(coordinate_system)
  prepared<-prepare_labelled_barplot(plot,primary,width,height);plot<-prepared$plot;height<-prepared$height
  primary <- route_output_path(primary,preferred="plot")
  if(!startsWith(primary,paste0(result_root_dir,"/")))stop("拒绝在staging之外成图")
  height_mm <- max(75,183*height/width)
  if (inherits(plot,"ggplot") && !inherits(plot,"patchwork")) {
    if(coordinate_system=="none")plot<-plot+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
    plot <- plot + theme(plot.margin=margin(8,12,8,8,unit="pt"))
  }
  pfos_save(frozen_adapter_display(plot,primary,height_mm),primary,height_mm)
  register_draw(primary,if(inherits(plot,"ggplot")) plot$data else NULL)
}
save_base_pdf_versions <- function(primary,width=8,height=6,expr) {
  expr_sub <- substitute(expr); expr_env <- parent.frame()
  primary <- route_output_path(primary,preferred="plot")
  if(!startsWith(primary,paste0(result_root_dir,"/")))stop("拒绝在staging之外成图")
  draw <- function() {
    if(grepl("upset",basename(primary),ignore.case=TRUE))grid::pushViewport(grid::viewport(width=.9,height=.94))
    graphics::par(family="Liberation Sans",cex.axis=7/8,cex.lab=1,bty="l")
    eval(expr_sub,envir=expr_env)
  }
  height_mm<-max(75,183*height/width)
  if(grepl("09_10_|10_8_TF_pathway|17_8_",basename(primary)))height_mm<-max(height_mm,260)
  pfos_save(frozen_adapter_display(draw,primary,height_mm),primary,height_mm);register_draw(primary)
}
save_pheatmap_versions <- function(primary,width=8,height=6,expr) {
  p <- eval(substitute(expr),envir=parent.frame()); g<-if(is.list(p)&&!is.null(p$gtable))p$gtable else p
  primary<-route_output_path(primary,preferred="plot")
  if(!startsWith(primary,paste0(result_root_dir,"/")))stop("拒绝在staging之外成图")
  pfos_save(g,primary,max(100,183*height/width));register_draw(primary)
}
chunk_code <- function(module,number) {
  p <- list.files(rmd_dir,pattern=paste0("^",sprintf("%02d",module),"_.*\\.Rmd$"),full.names=TRUE)
  if(length(p)!=1L)stop("模块文件不唯一: ",module)
  lines<-readLines(p,warn=FALSE); starts<-grep("^```\\{r",lines);ends<-grep("^```$",lines)
  if(number>length(starts))stop("绘图chunk不存在")
  begin<-starts[number];end<-ends[ends>begin][1]
  paste(lines[(begin+1L):(end-1L)],collapse="\n")
}
draw_chunk <- function(module,number) {
  message("[plots-only] module=",module," chunk=",number)
  code_text <- chunk_code(module,number)
  if(module==2L) {
    code_text <- gsub("Full cohort samples by group","Original plot cohort samples by group",code_text,fixed=TRUE)
    code_text <- gsub("Full-cohort PCA; no PCA-based sample selection","Original full-pool PCA; frozen selected-sample identity",code_text,fixed=TRUE)
  }
  code <- parse(text=code_text)
  # 依旧使用原绘图片段，但用冻结表替换片段中的统计赋值。
  replace_values <- list()
  if(module==11L && number==2L) replace_values <- list(
    community=quote(NULL), ppi_topology=quote(frozen_ppi$ppi_topology),
    ppi_centrality=quote(frozen_ppi$ppi_centrality))
  if(module==10L && number==2L) replace_values <- list(tf_target_edges=quote(frozen_tf$tf_target_edges))
  if(module==15L && number==4L) replace_values <- list(
    dsfa_event_table=quote(read_source_csv("15_splicing/01_event_level/tables/15_8b_DSFA_event_table.csv")),
    dsfa_gene_burden=quote(read_source_csv("15_splicing/01_event_level/tables/15_8d_DSFA_gene_burden.csv")))
  if(module==15L && number==7L) replace_values <- list(
    das_go=quote(read_source_csv("15_splicing/02_gene_level/tables/15_12_DAS_GO_enrichment.csv")),
    das_kegg=quote(read_source_csv("15_splicing/02_gene_level/tables/15_14_DAS_KEGG_enrichment.csv")))
  replace_ast <- function(x) {
    if(is.call(x) && identical(x[[1]],quote(ComplexHeatmap::Heatmap)))x[[1]]<-as.name("pfos_heatmap")
    if(is.call(x) && identical(x[[1]],quote(pheatmap::pheatmap)))x[[1]]<-as.name("pfos_pheatmap")
    n<-assigned_name(x)
    if(n %in% names(replace_values)) {x[3]<-list(replace_values[[n]]); return(x)}
    if(is.call(x) || is.expression(x)) for(i in seq_along(x)) if(!is.symbol(x[[i]])) x[i]<-list(replace_ast(x[[i]]))
    x
  }
  for(x in code) eval(replace_ast(x),envir=.GlobalEnv)
}
hydrate_list <- function(path) {x<-base::readRDS(path);for(n in names(x))assign(n,x[[n]],envir=.GlobalEnv);invisible(x)}
enrichment_tables<-hydrate_list(cache_06_tables_path);enrichment_objects<-hydrate_list(cache_06_objects_path)
human_enrichment_tables<-base::readRDS(cache_06_human_tables_path)
for(n in c("gsea_human_hallmark_df","gsea_human_c2_df","gsea_human_c5_df"))assign(n,human_enrichment_tables[[n]])
geneList_symbol_fc <- build_gene_foldchange_symbol(de_result)
human_plot_dir <- file.path(result_06_dir,"human_ortholog/plots")
human_table_dir <- file.path(result_06_dir,"human_ortholog/tables")

modules <- if(length(args)>=3L)as.integer(strsplit(args[3],"[, ]+")[[1]])else c(2:11,13:18)
for(module in modules) {
  drawing_sources <- c("rmd/plot_inputs.tsv","result/05_de/tables/5_1_annotated_DE_results.csv",
                       paste0(".cache/rmd/",unname(cache_names)))
  if(module %in% 2:4) for(k in seq_along(grep("^```\\{r",readLines(list.files(rmd_dir,pattern=paste0("^",sprintf("%02d",module),"_.*\\.Rmd$"),full.names=TRUE),warn=FALSE))))draw_chunk(module,k)
  if(module==2L) {
    top_selection_tbl <- auto_selection_summary %>% arrange(desc(score),max_outlier_distance,desc(between_group_distance),sample_key) %>% slice_head(n=15)
    p <- top_selection_tbl %>% mutate(combination_label=paste0(combination_type,": ",group1_samples," vs ",group2_samples),combination_label=factor(combination_label,levels=rev(combination_label))) %>%
      ggplot(aes(combination_label,score,fill=selected))+geom_col(width=.75)+coord_flip()+
      scale_fill_manual(values=viz_rekey_colors(viz_pal$sig_gradient,c(yes="high",no="low")))+
      labs(title="Historical sample-combination diagnostic",subtitle="Frozen original scores; current analysis uses the full cohort",x=NULL,y="Original selection score")+theme_rnaseq_publication()
    save_plot_versions(p,"result/02_selection/2_5_auto_selection_scores.pdf",width=11,height=max(5,.45*nrow(top_selection_tbl)))
  }
  if(module==5L) {
    draw_chunk(5,2);draw_chunk(5,3)
    source(file.path(rmd_dir,"plot_sensitivity.R"),local=.GlobalEnv)
  }
  if(module==6L)draw_chunk(6,2)
  if(module %in% c(7L,8L,18L))draw_chunk(module,1)
  if(module==9L) {
    mouse_gsva_table_dir<-file.path(result_09_dir,"tables");mouse_gsva_plot_dir<-file.path(result_09_dir,"plots")
    human_gsva_table_dir<-file.path(result_09_dir,"human_ortholog/tables");human_gsva_plot_dir<-file.path(result_09_dir,"human_ortholog/plots")
    mouse_gsva<-setNames(lapply(c("09_1_Mouse_MH","09_2_Mouse_KEGG","09_3_Mouse_Reactome","09_4_Mouse_M5"),function(n)read_source_matrix(paste0("09_gsva/tables/",n,"_GSVA_matrix.tsv"))),c("Mouse_MH","Mouse_KEGG","Mouse_Reactome","Mouse_M5"))
    human_gsva<-setNames(lapply(c("09_1_Human_Hallmark","09_2_Human_KEGG_Legacy","09_3_Human_Reactome","09_4_Human_C5"),function(n)read_source_matrix(paste0("09_gsva/human_ortholog/tables/",n,"_GSVA_matrix.tsv"))),c("Human_Hallmark","Human_KEGG_Legacy","Human_Reactome","Human_C5"))
    if(any(!vapply(c(mouse_gsva,human_gsva),function(x)identical(colnames(x),rownames(sample_info)),logical(1))))stop("GSVA缓存样本身份不符")
    mouse_gsva_leading_summary<-read_source_csv("09_gsva/tables/09_9_leading_GSVA_pathway_summary.csv")
    human_gsva_leading_summary<-read_source_csv("09_gsva/human_ortholog/tables/09_9_leading_GSVA_pathway_summary.csv")
    draw_chunk(9,2)
  }
  if(module==10L){frozen_tf<-hydrate_list(cache_10_tf_state_path);tf_regulon<-data.frame();draw_chunk(10,2);draw_chunk(10,3)}
  if(module==11L){
    frozen_ppi<-hydrate_list(cache_11_ppi_state_path)
    sig_de_symbols<-sig_de %>% filter(!is.na(gene_label),nzchar(gene_label)) %>% distinct(gene_label,.keep_all=TRUE)
    string_db<-NULL;string_mapped<-data.frame();string_species_id<-10090L
    for(k in 2:10)draw_chunk(11,k)
  }
  if(module %in% c(13L,14L)) {
    for(n in c("apa_all_events","apa_sig_events","apa_gene_summary_tbl","der_all_regions","der_sig_regions","der_gene_summary_tbl"))if(!exists(n))assign(n,data.frame())
    # 输入标准化只读既有APA/DER结果，不拟合事件检验。
    helper_code<-parse(text=chunk_code(15,1));for(x in helper_code)if(is_function_definition(x))eval(x)
    if(module==13L) {
      apa_all_events<-read_source_csv("13_apa/tables/13_1_APA_all_events.csv")
      apa_sig_events<-apa_all_events%>%filter(significant_apa)
      apa_gene_summary_tbl<-utils::read.csv(file.path(apa_dir,"apa_gene_summary.csv"),check.names=FALSE)
      apa_de_scatter_tbl<-read_source_csv("13_apa/tables/13_APA_gene_evidence.csv")
      apa_summary<-read_source_csv("13_apa/tables/13_APA_summary.csv")
      apa_sig_ids<-unique(apa_sig_events$gene_id_clean)
      apa_overlap_summary<-make_overlap_summary(apa_sig_ids,unique(sig_de$ensembl_id),"APA","DE")
    } else {
      der_all_regions<-read_source_csv("14_der/tables/14_1_all_DERs.csv")
      der_sig_regions<-der_all_regions%>%filter(significant_der)
      der_gene_summary_tbl<-utils::read.csv(file.path(der_dir,"der_gene_summary.csv"),check.names=FALSE)
      der_de_scatter_tbl<-read_source_csv("14_der/tables/14_DER_gene_evidence.csv")
      der_summary<-read_source_csv("14_der/tables/14_DER_summary.csv")
      der_sig_ids<-unique(der_sig_regions$gene_id_clean)
      der_overlap_summary<-make_overlap_summary(der_sig_ids,unique(sig_de$ensembl_id),"DER","DE")
    }
    draw_chunk(module,2)
  }
  if(module==15L) {
    hydrate_list(cache_12_wgcna_path)
    for(k in 1:10) {
      draw_chunk(15,k)
      if(k==1L) {
        frozen_rmats<-read_source_csv("15_splicing/01_event_level/tables/15_1_rMATS_all_events.csv")
        key<-paste(rmats_combined$event_key,rmats_combined$count_mode)
        frozen_key<-paste(frozen_rmats$event_key,frozen_rmats$count_mode)
        if(!identical(key,frozen_key))stop("rMATS事件或顺序与冻结结果不一致")
        for(n in c("delta_psi","event_pvalue","event_fdr"))if(any(abs(rmats_combined[[n]]-frozen_rmats[[n]])>1e-10,na.rm=TRUE))stop("rMATS冻结数值不一致: ",n)
      }
    }
  }
  if(module==16L) {
    immune_ssgsea_mat<-read_source_matrix("16_immune/tables/16_1_ssgsea_score_matrix.tsv")
    mmcp_mat<-read_source_matrix("16_immune/tables/16_2_mmcpcounter_score_matrix.tsv")
    immune_diff<-read_source_csv("16_immune/tables/16_immune_differential_scores.csv")
    immune_ssgsea_diff<-immune_diff%>%filter(method=="ssGSEA");immune_mmcp_diff<-immune_diff%>%filter(method=="mMCPcounter")
    immune_method_concordance<-read_source_csv("16_immune/tables/16_5_method_concordance.csv")
    draw_chunk(16,2)
  }
  if(module==17L) {
    hydrate_list(cache_12_wgcna_path); hydrate_list(cache_15_splicing_state_path);hydrate_list(cache_16_immune_state_path)
    mouse_activity_table_dir<-file.path(result_17_dir,"tables");mouse_activity_plot_dir<-file.path(result_17_dir,"plots")
    human_activity_table_dir<-file.path(result_17_dir,"human_ortholog/tables");human_activity_plot_dir<-file.path(result_17_dir,"human_ortholog/plots")
    pathway_activity_differential<-read_source_csv("17_pathway_activity/tables/17_pathway_activity_differential.csv")
    human_pathway_activity_differential<-read_source_csv("17_pathway_activity/human_ortholog/tables/17_pathway_activity_differential.csv")
    for(spec in list(c("hallmark","Mouse_MH", "1"),c("kegg","Mouse_KEGG","2"),c("reactome","Mouse_Reactome","3"))) {
      assign(paste0(spec[1],"_activity"),read_source_matrix(paste0("17_pathway_activity/tables/17_",spec[3],"_",spec[2],"_activity_matrix.tsv")))
      assign(paste0(spec[1],"_activity_diff"),pathway_activity_differential %>% filter(collection==spec[2]))
    }
    for(spec in list(c("human_hallmark","Human_Hallmark", "1"),c("human_kegg","Human_KEGG_Legacy","2"),c("human_reactome","Human_Reactome","3"))) {
      assign(paste0(spec[1],"_activity"),read_source_matrix(paste0("17_pathway_activity/human_ortholog/tables/17_",spec[3],"_",spec[2],"_activity_matrix.tsv")))
      assign(paste0(spec[1],"_activity_diff"),human_pathway_activity_differential %>% filter(collection==spec[2]))
    }
    draw_chunk(17,2);draw_chunk(17,3)
  }
}
manifest<-bind_rows(rendered)
dir.create(file.path(result_root_dir,"00_audit/tables"),recursive=TRUE,showWarnings=FALSE)
prior_path<-file.path(result_root_dir,"00_audit/tables/figure_manifest.tsv")
if(base::file.exists(prior_path)) {
  prior<-utils::read.delim(prior_path,check.names=FALSE)
  manifest<-bind_rows(prior[!prior$figure_id %in% manifest$figure_id,,drop=FALSE],manifest)
}
manifest$scientific_status <- vapply(manifest$pdf, figure_scientific_status, character(1))
utils::write.table(manifest,file.path(result_root_dir,"00_audit/tables/figure_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
message("COMPLETED_PLOTS ",nrow(manifest))
