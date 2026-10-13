#!/usr/bin/env Rscript
# 从当前完整RNA cohort及直接CpG记录恢复WGBS-RNA图型；不重跑DSS。
args <- commandArgs(TRUE)
if(length(args)!=2L) stop("usage: rebuild_wgbs_rna.R INTEGRATION_ROOT CONTEXT")
root <- normalizePath(args[[1]],mustWork=TRUE)
context <- args[[2L]]
cfg <- yaml::read_yaml(file.path(root,"config/project.yml"))
if(!context %in% names(cfg$joint_rna_context)) stop("context必须在joint_rna_context显式声明: ",context)
source(file.path(root,"scripts/joint_plot_helpers.R"))
message("WGBS-RNA当前cohort: ",context)
j <- joint_init(root,file.path("result/multiomics/wgbs_rna",context),paste0("wgbs_rna_",context))
sets <- joint_gene_sets(j)
w <- joint_input(j,paste0("wgbs_",context)); r <- joint_input(j,paste0("rna_",j$cfg$joint_rna_context[[context]]))
dmr <- joint_dmr(j,w);rna <- joint_rna(j,r)
risk <- if(context=="2cell") "2cell sensitivity only: low direct coverage" else if(context=="sperm") "Sperm WGBS / testis RNA: cross-cell-type; testis F0" else "Library-level exploratory; biological pairing unknown"
# 区域表完整保留所有映射；gene-centric展示的极值汇总另外列出冲突与区域数。
all_links <- left_join(dmr %>% dplyr::select(-any_of(c("SYMBOL","ENTREZID","GENENAME"))),rna,by="gene_id_clean") %>%
  mutate(MethylationDirection=ifelse(diffMethy>0,"Hyper","Hypo"))
joint_table(j,"all_region_rna_links",all_links)
region <- all_links %>% filter(is.finite(log2FoldChange),is.finite(pvalue)) %>%
  mutate(Quadrant=joint_quadrant(diffMethy,log2FoldChange))
all_region_table <- joint_table(j,"all_region_rna",region)
p <- ggplot(region,aes(diffMethy,log2FoldChange,color=Quadrant)) + geom_point(alpha=.6,size=1) +
  geom_hline(yintercept=c(-j$cfg$lfc_cutoff,j$cfg$lfc_cutoff),linetype=2) + geom_vline(xintercept=c(-.1,.1),linetype=2) +
  scale_color_manual(values=joint_colors) + labs(title=paste(context,"DMR-RNA quadrants"),subtitle=risk,x=paste0("Methylation fraction delta (",j$contrast,")"),y=paste0("RNA log2FC (",j$contrast,")"))
joint_plot(j,"region_quadrant_scatter","区域甲基化与RNA效应四象限",p,all_region_table)
layers <- bind_rows(region,region %>% filter(!is.na(CGI_context),CGI_context!="Open_Sea") %>% mutate(Layer="CGI"))
genes <- joint_gene_summary(layers) %>% mutate(Quadrant=joint_quadrant(diffMethy,log2FoldChange))
gene_table <- joint_table(j,"gene_centric_long",genes)
p <- ggplot(genes,aes(diffMethy,log2FoldChange,color=Quadrant)) + geom_point(alpha=.6,size=1) +
  facet_wrap(~Layer,scales="free") + scale_color_manual(values=joint_colors) +
  labs(title=paste(context,"gene-centric effects"),subtitle=risk,x=paste0("Methylation fraction delta (",j$contrast,")"),y=paste0("RNA log2FC (",j$contrast,")"))
joint_plot(j,"gene_centric_scatter","Promoter/GeneBody/CGI相关区域/Distal分层效应散点",p,gene_table,width=9,height=6.5)
p <- ggplot(genes,aes(Layer,log2FoldChange,fill=MethylationDirection)) +
  geom_violin(alpha=.3,trim=FALSE,position=position_dodge(.8),colour=NA) +
  stat_summary(fun.data=function(z) { z <- z[is.finite(z)]; data.frame(y=median(z),ymin=unname(quantile(z,.25,type=7)),ymax=unname(quantile(z,.75,type=7))) },geom="linerange",position=position_dodge(.8),linewidth=.25,colour=pfos_palette("neutral")[["dark"]],show.legend=FALSE) +
  stat_summary(fun=median,geom="point",position=position_dodge(.8),size=1.1,colour=pfos_palette("neutral")[["dark"]],show.legend=FALSE) +
  scale_fill_manual(values=c(Hyper=unname(pfos_palette("de")[["up"]]),Hypo=unname(pfos_palette("de")[["down"]]))) + labs(title=paste(context,"RNA distributions across methylation layers"),subtitle=risk,x=NULL,y="RNA log2FC")
joint_plot(j,"gene_centric_violin","甲基化分层中的RNA效应分布",p,gene_table)
composition <- genes %>% count(Layer,Quadrant,name="GeneCount")
p <- ggplot(composition,aes(Layer,GeneCount,fill=Quadrant)) + geom_col(position="fill") + scale_fill_manual(values=joint_colors) +
  labs(title=paste(context,"quadrant composition"),x=NULL,y="Fraction of annotated genes")
joint_plot(j,"gene_centric_composition","分层四象限组成",p,composition)
stats <- genes %>% filter(is.finite(diffMethy),is.finite(log2FoldChange)) %>% group_by(Layer) %>%
  summarise(n_genes=n(),spearman_rho=if(n()>2)cor(diffMethy,log2FoldChange,method="spearman") else NA_real_,
    conflicting_genes=sum(DirectionConflict),.groups="drop")
stats$interpretation <- "Across-feature descriptive rho; no biological treatment correlation P"
joint_table(j,"feature_effect_rho",stats)
flow <- region %>% filter(Quadrant!="NS",rna_selected) %>% count(Layer,MethylationDirection,Quadrant,name="n")
if(nrow(flow)) {
  p <- ggplot(flow,aes(axis1=Layer,axis2=MethylationDirection,axis3=Quadrant,y=n)) +
    ggalluvial::geom_alluvium(aes(fill=Quadrant),alpha=.7,width=.2) + ggalluvial::geom_stratum(width=.2,fill=unname(pfos_palette("neutral")[["missing"]])) +
    geom_text(stat="stratum",aes(label=after_stat(stratum)),size=8/ggplot2::.pt) + scale_fill_manual(values=joint_colors) +
    scale_x_discrete(limits=c("Feature","Methylation","RNA quadrant"),expand=c(.05,.05)) +
    labs(title=paste(context,"region feature / methylation / RNA"),subtitle="RNA raw P < 0.05; all regional conflicts retained",x=NULL,y="DMR count")
  joint_plot(j,"three_layer_alluvial","区域类别→甲基化方向→RNA四象限",p,flow,width=9,height=6)
} else joint_skip(j,"three_layer_alluvial","区域类别→甲基化方向→RNA四象限","无满足已冻结RNA阈值的四象限区域","EMPTY",all_region_table)
support <- read_csv(joint_path(j,w$sample_methylation),show_col_types=FALSE)
top <- all_links %>% mutate(Quadrant=ifelse(is.finite(log2FoldChange)&is.finite(pvalue),joint_quadrant(diffMethy,log2FoldChange),"RNA_unmeasured")) %>% distinct(DMR_id,.keep_all=TRUE) %>% arrange(desc(abs(diffMethy))) %>% slice_head(n=80L)
heat <- support %>% filter(RegionID %in% top$DMR_id) %>%
  left_join(top %>% dplyr::select(RegionID=DMR_id,SYMBOL,log2FoldChange,Quadrant,Layer,CGI_context),by="RegionID")
joint_heat(j,"top_dmr_sample_heatmap","Top DMR直接覆盖甲基化样本热图及RNA注释",heat,
  conclusion=paste(risk,"；NA表示该样本无直接覆盖，不以DSS平滑值补齐。"))
annotation <- top %>% transmute(RegionID=DMR_id,Sample="RNA log2FC",log2FoldChange)
joint_heat(j,"top_dmr_rna_annotation","Top DMR对应RNA效应注释",annotation,value="log2FoldChange")
# 完整gene/promoter可测量区域构成背景；不是已选DMR基因并集。
coverage <- joint_coverage(j,w,"measured_regions")
eligible_rna <- rna %>% filter(is.finite(pvalue),is.finite(log2FoldChange))
coverage_genes <- coverage$means %>% filter(eligible) %>% inner_join(eligible_rna,by="gene_id_clean")
eligible_table <- joint_table(j,"joint_eligible_gene_universe",coverage_genes)
deg <- coverage_genes %>% filter(rna_selected)
deg_table <- joint_table(j,"deg_focused_region_effects",deg)
p <- ggplot(deg,aes(log2FoldChange,diffMethy,color=Layer)) + geom_point(alpha=.6,size=1.2) + facet_wrap(~Layer) +
  geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2) +
  labs(title=paste(context,"DEG-focused measured-region effects"),subtitle=risk,x=paste0("RNA log2FC (",j$contrast,")"),y=paste0("Direct methylation delta (",j$contrast,")"))
joint_plot(j,"deg_focused_methylation_expression_scatter","DEG启动子及gene-body直接覆盖甲基化效应",p,deg_table,width=9)
for(layer in c("Promoter","Gene Body")) {
  top_degs <- deg %>% filter(Layer==layer) %>% arrange(desc(abs(log2FoldChange))) %>% slice_head(n=60L)
  long <- coverage$support %>% filter(Layer==layer,RegionID %in% top_degs$gene_id_clean) %>%
    left_join(top_degs %>% dplyr::select(RegionID=gene_id_clean,SYMBOL,log2FoldChange),by="RegionID")
  joint_heat(j,paste0("deg_",gsub(" ","_",tolower(layer)),"_sample_heatmap"),paste(context,"DEG",layer,"direct CpG methylation"),long)
}
cgi <- coverage$cgi_means %>% filter(eligible) %>% inner_join(eligible_rna,by="gene_id_clean")
cgi_deg <- cgi %>% filter(rna_selected)
p <- ggplot(cgi_deg,aes(CGI_context,log2FoldChange,fill=CGI_context)) + geom_boxplot(outlier.alpha=.2) +
  labs(title=paste(context,"DEGs at covered CGI contexts"),subtitle=risk,x=NULL,y="RNA log2FC")
joint_plot(j,"deg_cgi_boxplot","有直接覆盖的CGI背景中DEG效应分布",p,cgi_deg)
p <- ggplot(cgi_deg,aes(diffMethy,log2FoldChange,color=CGI_context)) + geom_point(alpha=.6,size=1) + facet_wrap(~CGI_context) +
  labs(title=paste(context,"CGI-related-region / DEG effects"),subtitle=risk,x=paste0("Direct methylation delta (",j$contrast,")"),y="RNA log2FC")
joint_plot(j,"deg_cgi_scatter","CGI区域甲基化与DEG效应",p,cgi_deg,width=9,height=6)
enrichment <- list()
for(layer in c("Promoter","Gene Body","CGI","Distal")) {
  universe <- if(layer=="CGI") unique(cgi$gene_id_clean) else unique(coverage_genes$gene_id_clean[coverage_genes$Layer==layer])
  for(quad in c("Hyper-Up","Hyper-Down","Hypo-Up","Hypo-Down")) {
    selected <- genes %>% filter(Layer==layer,Quadrant==quad,rna_selected)
    for(collection in names(sets)) {
      prefix <- paste(gsub(" ","_",tolower(layer)),gsub("-","_",tolower(quad)),tolower(collection),sep="_")
      title <- paste(context,layer,quad,collection)
      if(layer=="Distal") {
        for(family in if(collection %in% c("GO_MF","GO_CC")) c("dotplot","barplot","tree") else c("dotplot","barplot","tree","cnet","chord","cluster")) joint_skip(j,paste0(prefix,"_",family),title,"未构建与所有远端最近基因映射规则一致的直接覆盖背景；不使用DMR并集")
        next
      }
      z <- joint_ora_figures(j,prefix,selected$gene_id_clean,universe,sets[[collection]],title)
      if(nrow(z)) {z$Layer<-layer;z$Quadrant<-quad;z$Collection<-collection;enrichment[[prefix]]<-z}
    }
  }
}
enrichment <- bind_rows(enrichment)
if(nrow(enrichment)) {
  joint_table(j,"quadrant_enrichment_full",enrichment)
  for(collection in c("GO_BP","KEGG")) {
    z <- enrichment %>% filter(Collection==collection,pvalue<j$cfg$p_cutoff,Count>0) %>% arrange(pvalue) %>%
      group_by(term) %>% mutate(best=min(pvalue)) %>% ungroup() %>% arrange(best) %>%
      filter(term %in% head(unique(term),25L)) %>% mutate(Context=paste(Layer,Quadrant),Score=-log10(pmax(pvalue,1e-300)))
    if(nrow(z)) {
      p <- ggplot(z,aes(Context,Description,fill=Score)) + geom_tile() + scale_fill_gradientn(colors=unname(pfos_palette("sequential"))) +
        labs(title=paste(context,collection,"quadrant enrichment"),subtitle="Exploratory raw P; BH retained in source",x=NULL,y=NULL) + theme(axis.text.x=element_text(angle=60,hjust=1),axis.text.y=element_text(size=7))
      joint_plot(j,paste0(tolower(collection),"_quadrant_heatmap"),paste(collection,"四象限富集热图"),p,z,width=9,height=8)
    } else joint_skip(j,paste0(tolower(collection),"_quadrant_heatmap"),paste(collection,"四象限富集热图"),"没有raw-P合格条目","EMPTY")
  }
}
cgi_ora <- list()
for(cg in unique(cgi$CGI_context)) for(collection in names(sets)) {
  z <- joint_ora_figures(j,paste0("deg_",tolower(cg),"_",tolower(collection)),
    cgi_deg$gene_id_clean[cgi_deg$CGI_context==cg],cgi$gene_id_clean[cgi$CGI_context==cg],sets[[collection]],paste(context,"DEG",cg,collection))
  if(nrow(z)){z$CGI_context<-cg;z$Collection<-collection;cgi_ora[[paste(cg,collection)]]<-z}
}
cgi_ora <- bind_rows(cgi_ora)
joint_table(j,"deg_cgi_enrichment_full",cgi_ora)
for(collection in c("GO_BP","KEGG")) {
  z <- if(nrow(cgi_ora)) cgi_ora %>% filter(Collection==collection,pvalue<j$cfg$p_cutoff,Count>0) %>% arrange(pvalue) %>% filter(term %in% head(unique(term),25L)) else tibble()
  if(nrow(z)) {
    p <- ggplot(z,aes(CGI_context,stringr::str_wrap(Description,45),fill=-log10(pmax(pvalue,1e-300)))) + geom_tile() +
      scale_fill_gradientn(colors=unname(pfos_palette("sequential"))) + labs(title=paste(context,collection,"DEG CGI-related regions"),subtitle="Exploratory raw P; BH retained; context-specific directly covered background",x=NULL,y=NULL)
    joint_plot(j,paste0("deg_cgi_",tolower(collection),"_context_heatmap"),paste(collection,"DEG CGI相关区域富集热图"),p,z,height=8)
  } else joint_skip(j,paste0("deg_cgi_",tolower(collection),"_context_heatmap"),paste(collection,"DEG CGI相关区域富集热图"),"没有按预设raw-P阈值合格的条目","EMPTY")
}
# 相同coverage可评估背景下独立两层的富集比较，避免旧RNA enrichment缓存。
joint_universe <- unique(coverage_genes$gene_id_clean)
pathway_compare <- read_tsv(joint_cached(j,"wgbs_rna_kegg_scatter"),show_col_types=FALSE)
if(nrow(pathway_compare)) {
  p <- ggplot(pathway_compare,aes(WGBS_FoldEnrichment,RNA_FoldEnrichment,size=-log10(pmax(WGBS_p,1e-300)),color=-log10(pmax(RNA_p,1e-300)))) +
    geom_point(alpha=.65) + scale_color_gradientn(colors=unname(pfos_palette("sequential"))) + labs(title=paste(context,"WGBS / current RNA KEGG enrichment"),subtitle="Same current directly covered / RNA-testable universe; no causal claim",x="WGBS fold enrichment",y="RNA fold enrichment")
  joint_plot(j,"wgbs_rna_kegg_scatter","当前覆盖背景WGBS与RNA KEGG富集比较",p,pathway_compare)
  top_pathway <- pathway_compare %>% mutate(best=pmin(WGBS_p,RNA_p)) %>% arrange(best) %>% slice_head(n=20L) %>%
    pivot_longer(c(WGBS_FoldEnrichment,RNA_FoldEnrichment),names_to="Omics",values_to="FoldEnrichment")
  p <- ggplot(top_pathway,aes(Omics,stringr::str_wrap(Description,45),size=FoldEnrichment,color=FoldEnrichment)) + geom_point() +
    scale_color_gradientn(colors=unname(pfos_palette("sequential"))) + labs(title=paste(context,"KEGG enrichment comparison"),x=NULL,y=NULL)
  joint_plot(j,"wgbs_rna_kegg_dotplot_compare","当前WGBS-RNA KEGG富集对照点图",p,top_pathway,height=8)
} else for(id in c("wgbs_rna_kegg_scatter","wgbs_rna_kegg_dotplot_compare")) joint_skip(j,id,"当前WGBS与RNA KEGG富集比较","共同可评估背景下没有可展示条目","EMPTY")
# CGI相关区域全集同样可提供完整排序；仅极值汇总展示，保留所有区域及方向冲突表。
cgi_rank <- cgi %>% group_by(gene_id_clean) %>% arrange(desc(abs(diffMethy)),RegionID,.by_group=TRUE) %>%
  mutate(RegionCount=n(),DirectionConflict=any(diffMethy>0)&any(diffMethy<0)) %>% slice_head(n=1L) %>% ungroup()
joint_table(j,"complete_cgi_related_gene_ranking",cgi_rank)
for(collection in c("GO_BP","KEGG")) {
  joint_gsea(j,paste0("cgi_related_",tolower(collection),"_methylation"),setNames(cgi_rank$diffMethy,cgi_rank$gene_id_clean),sets[[collection]],paste(context,"CGI-related regions",collection,"direct methylation rank"))
  joint_gsea(j,paste0("cgi_related_",tolower(collection),"_effect_product"),setNames(-cgi_rank$diffMethy*cgi_rank$log2FoldChange,cgi_rank$gene_id_clean),sets[[collection]],paste(context,"CGI-related regions",collection,"predefined exploratory effect-product rank"))
  joint_skip(j,paste0("distal_",tolower(collection),"_gsea"),paste(collection,"远端甲基化完整排序"),"没有与所有远端最近基因映射规则一致的完整可评估coverage排序，不能用已选DMR代替")
}
# 仅完整覆盖支持的gene/promoter排序可执行GSEA；RNA统计与甲基化值保留来源。
for(layer in c("Promoter","Gene Body")) {
  z <- coverage_genes %>% filter(Layer==layer)
  methyl_rank <- setNames(z$diffMethy,z$gene_id_clean)
  joint_rank <- setNames(if(layer=="Promoter") -z$diffMethy*z$log2FoldChange else z$diffMethy*z$log2FoldChange,z$gene_id_clean)
  for(collection in c("GO_BP","KEGG")) {
    prefix <- paste0(gsub(" ","_",tolower(layer)),"_",tolower(collection))
    joint_gsea(j,paste0(prefix,"_methylation"),methyl_rank,sets[[collection]],paste(context,layer,collection,"direct methylation rank"))
    joint_gsea(j,paste0(prefix,"_coupling"),joint_rank,sets[[collection]],paste(context,layer,collection,"predefined exploratory effect-product rank"))
  }
}
joint_skip(j,"sample_level_crossomics_correlation","跨组学样本相关","组学间动物/父本/窝及样本配对关系unknown，不能按库名假设匹配")
joint_finish(j)
cat("WGBS-RNA completed context: ",context,"\n")
