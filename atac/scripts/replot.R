#!/usr/bin/env Rscript
# 现役冻结源表和绘图对象重放；不执行统计 Rmd、DiffBind、chromVAR 或富集检验。
options(device=function(...)grDevices::cairo_pdf(tempfile(fileext='.pdf'),family='Liberation Sans'))
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2L)
root<-normalizePath(a[1],mustWork=TRUE);output<-normalizePath(a[2],mustWork=TRUE);rmd<-file.path(root,'rmd')
Sys.setenv(ATAC_ROOT_DIR=root,ATAC_RMD_DIR=rmd,ATAC_RUN_PROFILE='report',ATAC_REUSE_STEP_TABLES='FALSE',ATAC_SAVE_PNG_PREVIEWS='TRUE')
source(file.path(rmd,'01_main_setup.R'),local=globalenv())
result_dir<-output
save_step_table<-function(x,step,file_name,...)invisible(file_name)
write.csv<-function(x,file,...)invisible(file)
save_csv_df<-function(x,file)invisible(file)
write_skipped_state<-function(...)invisible(NULL)
remove_output_with_preview<-function(...)invisible(NULL)
# Only temporary staged artifacts may be touched. All old numerical tables are
# checked byte-for-byte after drawing by replot.py.
unlink<-function(...)invisible(0L)
save_state<-function(...)invisible(NULL)
readtab<-function(stage,name) {
  p<-file.path(output,stage,'tables',name);if(!file.exists(p))stop('MISSING: ',p)
  read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
}
# Existing plot functions are shared with ordinary analysis. Use exactly their
# plot expression, extracted from the parsed R syntax tree without evaluating
# any containing statistical branch.
parsed<-new.env(parent=emptyenv())
expressions<-function(file) {
  if(!exists(file,parsed,inherits=FALSE)) {
    f<-tempfile(fileext='.R');knitr::purl(file.path(rmd,file),output=f,quiet=TRUE)
    e<-parse(f);base::unlink(f);assign(file,e,parsed)
  }
  get(file,parsed,inherits=FALSE)
}
lookup<-function(file,name) {
  hit<-list()
  scan<-function(e) {
    if(!is.call(e)&&!is.expression(e)&&!is.pairlist(e))return()
    if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])==name)hit[[length(hit)+1L]]<<-e
    for(j in seq_along(e)) {
      if(identical(e[[j]],quote(expr=)))next
      x<-e[[j]];if(is.call(x)||is.expression(x)||is.pairlist(x))scan(x)
    }
  }
  scan(expressions(file));if(!length(hit))stop('绘图定义缺失: ',file,':',name)
  hit[[1]]
}
expr<-function(file,name)eval(lookup(file,name),globalenv())
save<-function(p,file,height=140)pfos_save(atac_standardize_plot(p),file,height)
atac_save_object_dual<-function(plot_obj,file,width=10,height=8,standardize=FALSE,bg='white',coordinate_system=c('axes','none')) {
  plot_obj<-atac_text_layers(plot_obj)
  if(standardize)plot_obj<-atac_standardize_plot(plot_obj,match.arg(coordinate_system))
  pfos_save(plot_obj,file,max(80,183*height/width));invisible(file)
}
# pheatmap original defaults are normalized once; data and clustering are retained.
old_pheatmap<-pheatmap::pheatmap
pheatmap<-function(...,fontsize=8,fontsize_row=7,fontsize_col=7,fontfamily='Liberation Sans')old_pheatmap(...,fontsize=fontsize,fontsize_row=fontsize_row,fontsize_col=fontsize_col,fontfamily=fontfamily)
# Frozen differential tables are authoritative, independent of cache age.
all_res_df<-standardize_diff_peak_df(readtab('05_differential_accessibility','Differential_Peaks_AllSites.csv'))
res_df<-standardize_diff_peak_df(readtab('05_differential_accessibility','Differential_Peaks.csv'))
audit<-readtab('05_differential_accessibility','Differential_Peak_Source_Audit.csv')
if(!all(audit$AnalysisSource=='native_diffbind'))stop('原差异来源不是 native_diffbind')
if(any(audit$ReferenceGroup!=baseline_group)||any(audit$ComparisonGroup!=treatment_group))stop('对比方向不同')
if(any(audit$SampleIDs!=paste(samples,collapse=';')))stop('原图样本与当前 reads_info 不同')
load_frozen_cache<-function(name) {
  f<-file.path(root,'.cache/rmd',name);if(!file.exists(f))stop('MISSING frozen drawing object: ',f)
  x<-readRDS(f)
  # Source configuration must agree; path-sensitive historical signatures remain
  # historical and are not rewritten. Numerical tables remain authoritative.
  cfg<-x$.__analysis_signature_config__
  if(!is.null(cfg))for(k in intersect(names(cfg),names(analysis_signature_config))) {
    if(k=='diffbind_version')next
    if(!identical(as.character(cfg[[k]]),as.character(analysis_signature_config[[k]])))stop('缓存配置不一致: ',k)
  }
  x$.__analysis_signature__<-NULL;x$.__analysis_signature_config__<-NULL
  list2env(x,globalenv());invisible(x)
}
load_frozen_cache('04_peak_context.rds')
sel<-Sys.getenv('PFOS_PLOT_STEPS');steps<-if(!nzchar(sel))3:13 else unique(unlist(lapply(strsplit(gsub(',',' ',sel),' +')[[1]],function(x){z<-as.integer(strsplit(x,'-',fixed=TRUE)[[1]]);if(length(z)==2)seq(z[1],z[2])else z})))
for(step in steps) {
  cat('PLOTS_ONLY stage',step,'\n');flush.console()
  if(step==3) {
    e<-expressions('03_atac_qc.Rmd')
    bad<-c('t.test','wilcox.test','fisher.test','dba.analyze','enrichGO','enrichKEGG','system2')
    check<-function(e){if(is.call(e)){if(any(as.character(e[[1]])%in%bad))stop('QC 禁止统计调用');lapply(as.list(e)[-1],check)}}
    lapply(e,check);eval(e,globalenv())
  }
  if(step==4) {
    file<-'04_peak_catalog_annotation.Rmd';stage<-'04_peak_catalog_annotation'
    annotation_summary<-readtab(stage,'Peak_Annotation_Summary.csv')
    feature_cols<-resolve_named_colors(feature_levels_main,viz_pal$feature,fallback_colors=viz_pal$feature_accent,default_color=viz_pal$neutral_mid)
    expr(file,'p_anno');save(p_anno,step_plot_file('Step4','Peak_Annotation_Barplot.pdf'),110)
    peak_widths_by_group<-readtab(stage,'Peak_Width_By_Group.csv');expr(file,'p_width');save(p_width,step_plot_file('Step4','Peak_Width_Distribution.pdf'),110)
    peak_widths_df<-readtab(stage,'Peak_Length_By_Sample.csv');sample_peak_feature_summary<-readtab(stage,'Peak_Feature_Composition_By_Sample.csv');macs2_peak_signal_df<-readtab(stage,'MACS2_Peak_Signal_By_Sample.csv')
    for(item in list(c('p_sample_peak_length','Peak_Length_Distribution_By_Sample',130),c('p_sample_peak_feature','Peak_Feature_Composition_By_Sample',120),c('p_macs2_peak_signal','MACS2_Peak_Signal_Distribution_By_Sample',140))) {
      expr(file,item[1]);save(get(item[1]),step_plot_file('Step4',paste0(item[2],'.pdf')),as.numeric(item[3]))
    }
    overlap_by_group_df<-readtab(stage,'Peak_Overlap_Membership.csv')
    u<-atac_plot_peak_upset(overlap_by_group_df)
    pfos_save(function(){print(u);invisible(NULL)},step_plot_file('Step4','Peak_UpSet.pdf'),140)
    dist_to_tss_df<-readtab(stage,'Peak_DistToTSS_Data.csv');expr(file,'p_dist_tss');save(p_dist_tss,step_plot_file('Step4','Peak_DistToTSS.pdf'),110)
    group_colors<-viz_pal$group;expr(file,'peaks_by_group');expr(file,'p_cov');save(p_cov,step_plot_file('Step4','Peak_Chrom_Distribution.pdf'),200)
    tss_profile_df<-readtab(stage,'Peak_TSS_Profile_Data.csv');expr(file,'p_tss_profile');save(p_tss_profile,step_plot_file('Step4','Peak_TSS_Profile.pdf'),110)
  }
  if(step==5) {
    file<-'05_differential_accessibility.Rmd';stage<-'05_differential_accessibility'
    expr(file,'diff_sample_count_columns');counts_norm<-readtab(stage,'DiffBind_Normalized_Counts.csv');sample_cols<-diff_sample_count_columns(counts_norm)
    if(!setequal(sample_cols,samples))stop('冻结计数矩阵样本不一致')
    log_counts<-log2(as.matrix(counts_norm[,sample_cols])+1)
    pca_df<-readtab(stage,'Sample_PCA_Scores.csv');var_explained<-round(prcomp(t(log_counts),scale.=TRUE)$sdev^2/sum(prcomp(t(log_counts),scale.=TRUE)$sdev^2)*100,1)
    pca_colors<-c(viz_pal$group,Unknown=viz_pal$neutral_mid);expr(file,'p_pca');save(p_pca,step_plot_file('Step5','Sample_PCA.pdf'),160)
    cor0<-readtab(stage,'Sample_Correlation_Matrix.csv');sample_cor<-as.matrix(cor0[,sample_cols]);rownames(sample_cor)<-cor0$SampleID
    ann<-data.frame(Group=factor(sample_group(sample_cols),levels=group_levels),row.names=sample_cols)
    ht<-pheatmap(sample_cor,cluster_rows=FALSE,cluster_cols=FALSE,annotation_col=ann,annotation_row=ann,annotation_colors=list(Group=viz_pal$group),color=viz_pal$heat_meth,breaks=seq(-1,1,length.out=length(viz_pal$heat_meth)+1),legend_breaks=c(-1,0,1),legend_labels=c("-1","0","1"),na_col=unname(pfos_palette("neutral")["missing"]),silent=TRUE)
    ht<-atac_label_heatmap_legend(ht,"Pearson r")
    pfos_save(ht$gtable,step_plot_file('Step5','Sample_Correlation_Heatmap.pdf'),150)
    volcano_df<-all_res_df%>%mutate(PlotP=SelectedSignificanceValue,logP=-log10(pmax(PlotP,1e-300)))
    expr(file,'p_volcano');save(p_volcano,step_plot_file('Step5','Volcano_Plot.pdf'),140)
    baseline_cols<-sample_cols[sample_group(sample_cols)==baseline_group];treatment_cols<-sample_cols[sample_group(sample_cols)==treatment_group]
    counts_norm$Mean_Baseline<-rowMeans(counts_norm[,baseline_cols]);counts_norm$Mean_Treatment<-rowMeans(counts_norm[,treatment_cols])
    scatter_data<-counts_norm%>%select(peak_id_join,Mean_Baseline,Mean_Treatment)%>%left_join(all_res_df%>%select(peak_id_join,pvalue,FDR,Fold,PrimaryDiffPeak),by='peak_id_join')%>%mutate(Significance=factor(case_when(PrimaryDiffPeak&Fold>0~'primaryUp',PrimaryDiffPeak&Fold<0~'primaryDown',TRUE~'nonSig'),levels=c('nonSig','primaryDown','primaryUp')))
    expr(file,'p_scatter');save(p_scatter,step_plot_file('Step5','Scatter_Plot_DiffPeaks.pdf'),150)
    hm<-readtab(stage,'Top_DiffPeaks_Sample_Matrix.csv');heatmap_mat<-t(scale(t(log2(as.matrix(hm[,sample_cols])+1))));rownames(heatmap_mat)<-hm$peak_id;cluster_display_mat<-heatmap_mat;cluster_display_mat[!is.finite(cluster_display_mat)]<-0;heatmap_row_cluster<-stats::hclust(stats::dist(cluster_display_mat))
    ht<-pheatmap(heatmap_mat,cluster_rows=heatmap_row_cluster,cluster_cols=FALSE,annotation_col=ann,annotation_colors=list(Group=viz_pal$group),color=viz_pal$heat_meth,breaks=seq(-2.5,2.5,length.out=length(viz_pal$heat_meth)+1),legend_breaks=c(-2.5,0,2.5),legend_labels=c("≤-2.5","0","≥2.5"),na_col=unname(pfos_palette("neutral")["missing"]),show_rownames=FALSE,silent=TRUE)
    ht<-atac_label_heatmap_legend(ht,"Row Z-score")
    pfos_save(ht$gtable,step_plot_file('Step5','Top_DiffPeaks_Heatmap.pdf'),210)
    # dba.plotVenn is a display method on the frozen analyzed DBA object, not dba.analyze.
    dcache<-load_frozen_cache('05_differential_accessibility_state.rds')
    if(exists('atac_dba')&&!is.null(atac_dba))pfos_save(function(){par(cex.main=1,cex.axis=7/8,cex.lab=1);atac_plot_direction_venn(atac_dba,contrast=1,bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE,th=if(atac_significance_mode=='pvalue')atac_primary_p_cutoff else atac_diff_fdr_cutoff,bUsePval=atac_significance_mode=='pvalue',method=DiffBind::DBA_DESEQ2)},step_plot_file('Step5','DiffPeaks_Venn.pdf'),140)
    direction_context_df<-readtab(stage,'Directional_Chromatin_Context.csv');direction_feature_summary<-readtab(stage,'Directional_Feature_Summary.csv')
    direction_peak_colors<-setNames(c(viz_pal$de[reference_open_label],viz_pal$de[comparison_open_label],viz_pal$de['Neutral']),c(reference_open_label,comparison_open_label,'Neutral'))
    for(item in list(c('p_direction_context','Directional_Chromatin_Context'),c('p_direction_tss','Directional_DistToTSS'),c('p_direction_width','Directional_Peak_Width'))) {
      expr(file,item[1]);save(get(item[1]),step_plot_file('Step5',paste0(item[2],'.pdf')),135)
    }
    save(p_direction_context,step_plot_file('Step5','Directional_Annotation_Composition.pdf'),135)
    save(p_direction_tss,step_plot_file('Step5','Directional_TSS_Density.pdf'),135)
    save(p_direction_width,step_plot_file('Step5','Directional_Peak_Width_Ridge.pdf'),155)
    # Cache display objects cannot override the authoritative current differential table.
    all_res_df<-standardize_diff_peak_df(readtab(stage,'Differential_Peaks_AllSites.csv'));res_df<-standardize_diff_peak_df(readtab(stage,'Differential_Peaks.csv'))
  }
  if(step==6) {
    # Original Gviz display code reads the frozen RPGC tracks. No alignment/counting.
    eval(expressions('06_genome_tracks.Rmd'),globalenv())
  }
  if(step==7) {
    load_frozen_cache('07_peak_gene_enrichment_state.rds')
    # Confirm cached result tables match the frozen on-disk numerical source.
    for(key in names(go_results)) {
      tab<-readtab('07_peak_gene_enrichment',paste0('GO_',key,'.csv'));obj<-go_results[[key]]
      if(!setequal(tab$ID,as.data.frame(obj)$ID))stop('GO cache/source ID mismatch: ',key)
      plot_enrichment_suite(obj,paste0('GO_',key),key,step_plot_dir('Step7'),show_n=15)
    }
    for(key in names(kegg_results)) {
      tab<-readtab('07_peak_gene_enrichment',paste0('KEGG_',key,'.csv'));obj<-kegg_results[[key]]
      if(!setequal(tab$ID,as.data.frame(obj)$ID))stop('KEGG cache/source ID mismatch: ',key)
      plot_enrichment_suite(obj,paste0('KEGG_',key),key,step_plot_dir('Step7'),show_n=15)
    }
    plot_compare_cluster_suite(ck_go_raw,'GO_Compare','GO comparison',step_plot_dir('Step7'),show_n=15)
    plot_compare_cluster_suite(ck_kegg_raw,'KEGG_Compare','KEGG comparison',step_plot_dir('Step7'),show_n=15)
  }
  if(step==8) {
    file<-'08_motif_activity.Rmd';stage<-'08_motif_activity'
    load_frozen_cache('08_motif_state.rds')
    motif_res<-readtab(stage,'Motif_Enrichment.csv');significant_motifs<-motif_res%>%filter(is.finite(pvalue),pvalue>=0,pvalue<.05)%>%arrange(pvalue,desc(count))%>%head(20)
    # Original bar plot variable name differs from its full motif table.
    top_motifs_plot<-significant_motifs%>%mutate(PlotScore=-log10(pmax(pvalue,1e-300)))
    expr(file,'p_motif');save(p_motif,step_plot_file('Step8','Motif_Enrichment_Barplot.pdf'),180)
    directional_motif_df<-readtab(stage,'Directional_Motif_Enrichment.csv');if(!'PlotScore'%in%names(directional_motif_df))directional_motif_df$PlotScore<--log10(pmax(directional_motif_df$pvalue,1e-300))
    expr(file,'p_directional_motif');save(p_directional_motif,step_plot_file('Step8','Directional_Motif_Enrichment.pdf'),160)
    z<-readtab(stage,'ChromVAR_Deviation_ZScores.csv');id<-names(z)[1];m<-as.matrix(z[,samples,drop=FALSE]);rownames(m)<-z[[id]]
    chromvar_z_extent <- max(abs(m),na.rm=TRUE);if(!is.finite(chromvar_z_extent)||chromvar_z_extent==0)chromvar_z_extent<-1
    summary<-readtab(stage,'Motif_Activity_Summary.csv');top<-summary%>%filter(is.finite(pvalue),pvalue>=0,pvalue<.05)%>%head(30)
    m<-m[intersect(top$motif_name,rownames(m)),,drop=FALSE]
    ht<-pheatmap(m,cluster_rows=nrow(m)>1,cluster_cols=TRUE,annotation_col=data.frame(Group=factor(conditions,levels=group_levels),row.names=samples),annotation_colors=list(Group=viz_pal$group),color=viz_pal$heat_meth,breaks=seq(-chromvar_z_extent,chromvar_z_extent,length.out=length(viz_pal$heat_meth)+1),na_col=unname(pfos_palette("neutral")["missing"]),silent=TRUE)
    ht<-atac_label_heatmap_legend(ht,"chromVAR Z")
    pfos_save(ht$gtable,step_plot_file('Step8','ChromVAR_Motif_Activity_Heatmap.pdf'),220)
    motif_activity_enrichment_tbl<-readtab(stage,'Motif_Activity_Enrichment_Combined.csv');top_motif_evidence<-motif_activity_enrichment_tbl%>%filter(!is.na(motif_name),motif_name!='')%>%head(25)
    expr(file,'motif_heat_df');expr(file,'p_motif_dir_heat');save(p_motif_dir_heat,step_plot_file('Step8','Motif_Direction_Heatmap.pdf'),200)
    expr(file,'p_motif_activity_vs_enrichment');save(p_motif_activity_vs_enrichment,step_plot_file('Step8','Motif_Activity_vs_Enrichment.pdf'),160)
  }
  if(step==9) {
    file<-'09_footprinting_tf_evidence.Rmd';stage<-'09_footprinting_tf_evidence'
    profiles<-readtab(stage,'Footprint_Profile_Group_Values.csv')
    # The sequence logo is retrieved by motif identity from the same JASPAR release.
    opts<-list(species=10090,all_versions=FALSE);pfm<-TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,opts)
    for(motif_name in unique(profiles$motif_name)) {
      profile_data<-profiles[profiles$motif_name==motif_name,,drop=FALSE]
      ix<-which(vapply(pfm,function(x)TFBSTools::name(x)==motif_name,logical(1)))[1]
      if(is.na(ix))stop('MISSING motif PWM: ',motif_name)
      pwm_mat<-as.matrix(TFBSTools::toPWM(pfm[[ix]],type='prob'));motif_len<-ncol(pwm_mat);window_size<-100
      expr(file,'p_logo');profile_data$x_idx<-profile_data$distance+window_size+1;center_idx<-window_size+1;y_max<-max(profile_data$insertions,na.rm=TRUE)*1.05
      expr(file,'poly_df');expr(file,'p_fp');expr(file,'p_combined')
      save(p_combined,step_plot_file('Step9',paste0('Footprint_',gsub('::','_',motif_name,fixed=TRUE),'.pdf')),150)
    }
  }
  if(step%in%c(10,12)) {
    stage<-if(step==10)'10_tf_gene_network'else'12_machine_learning'
    name<-if(step==10)'TF_Target_Network_Skipped_State.csv'else'Machine_Learning_Classification_Skipped_State.csv'
    write_module_status_plot(paste0('Step',step),if(step==10)'Module 10 TF-gene network status'else'Module 12: supervised prediction not performed')
  }
  if(step==11) {
    file<-'11_sequence_conservation_gc.Rmd';conservation_df<-readtab('11_sequence_conservation_gc','Peak_Conservation_Analysis.csv')
    expr(file,'gc_hist');save(gc_hist,step_plot_file('Step11','GC_Content_Distribution.pdf'),130)
    expr(file,'gc_vs_fc');save(gc_vs_fc,step_plot_file('Step11','GC_vs_FoldChange.pdf'),140)
  }
  if(step==13) {
    file<-'13_keygenes_locus_reports.Rmd';stage<-'13_keygenes_locus_reports'
    step15_plot_dir<-step_plot_dir('Step13');gene_candidates<-readtab(stage,'ATAC_Gene_Candidates.csv');expr(file,'p');save(p,step_plot_file('Step13','key_genes_accessibility.pdf'),150)
    expr(file,'plot_locus_report')
    loci<-readtab(stage,'Representative_Loci.csv');for(g in head(loci$geneSymbol,8))plot_locus_report(g)
  }
}
cat('纯绘图完成；统计方法与科学数值保持原样。\n')
