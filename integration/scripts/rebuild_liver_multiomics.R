#!/usr/bin/env Rscript
# 复用原 ATAC 09/11 图型，当前RNA完整cohort；所有区域链接是注释/坐标描述。
args <- commandArgs(TRUE)
if(length(args)!=1L) stop("usage: rebuild_liver_multiomics.R INTEGRATION_ROOT")
root <- normalizePath(args[[1]],mustWork=TRUE)
source(file.path(root,"scripts/joint_plot_helpers.R"))
j <- joint_init(root,"result/multiomics/liver/figures","liver_atac_rna_wgbs")
r <- joint_input(j,"rna_liver");a <- joint_input(j,"atac_liver");w <- joint_input(j,"wgbs_liver")
rna <- joint_rna(j,r);dmr <- joint_dmr(j,w)
atac <- read_csv(joint_path(j,a$path),show_col_types=FALSE) %>% dplyr::rename(atac_pvalue=pvalue) %>% mutate(gene_id_clean=sub("\\..*$","",ENSEMBL),Layer=joint_layer(annotation))
atac$ATACDirection <- ifelse(atac$Fold>0,"Open","Closed")
all_peak_links <- atac %>% dplyr::select(-any_of(c("SYMBOL","ENTREZID","GENENAME"))) %>% left_join(rna,by="gene_id_clean")
joint_table(j,"all_atac_region_rna_links",all_peak_links)
peak <- all_peak_links %>% filter(is.finite(log2FoldChange),is.finite(pvalue)) %>%
  mutate(Quadrant=case_when(Fold>0 & log2FoldChange>0 ~ "Open-Up",Fold>0 & log2FoldChange<0 ~ "Open-Down",Fold<0 & log2FoldChange>0 ~ "Closed-Up",TRUE ~ "Closed-Down"))
cols <- c("Open-Up"=unname(pfos_palette("de")[["up"]]),"Open-Down"=unname(pfos_palette("accent")[["purple"]]),"Closed-Up"=unname(pfos_palette("accent")[["teal"]]),"Closed-Down"=unname(pfos_palette("de")[["down"]]))
peak_table <- joint_table(j,"rna_testable_atac_region_links",peak)
p <- ggplot(peak,aes(Fold,log2FoldChange,color=Quadrant)) + geom_point(alpha=.55,size=1) +
  geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2) + scale_color_manual(values=cols) +
  labs(title="Liver ATAC-RNA region effects",subtitle=paste0("",j$contrast,"; region-to-nearest-gene annotation; not sample correlations"),x="ATAC log2FC",y="RNA log2FC")
joint_plot(j,"atac_rna_peak_quadrants","肝脏开放区域与RNA效应四象限",p,peak_table)
gene <- joint_gene_summary(peak,effect="Fold",region="peak_id")
gene_table <- joint_table(j,"atac_rna_gene_centric",gene)
p <- ggplot(gene,aes(Fold,log2FoldChange,color=Quadrant)) + geom_point(alpha=.6,size=1.2) + facet_wrap(~Layer) +
  scale_color_manual(values=cols) + labs(title="Liver ATAC-RNA gene-centric effects",subtitle="Maximum-absolute region shown; all mapped regions and conflicts retained",x="ATAC log2FC",y="RNA log2FC")
joint_plot(j,"atac_rna_gene_scatter","肝脏ATAC-RNA基因层分区域散点",p,gene_table,width=9)
comp <- gene %>% count(Layer,Quadrant,name="GeneCount")
p <- ggplot(comp,aes(Layer,GeneCount,fill=Quadrant)) + geom_col(position="fill") + scale_fill_manual(values=cols) +
  labs(title="Liver ATAC-RNA quadrant composition",x=NULL,y="Fraction")
joint_plot(j,"atac_rna_quadrant_composition","肝脏ATAC-RNA四象限组成",p,comp)
p <- ggplot(gene,aes(Quadrant,log2FoldChange,fill=Quadrant)) + geom_violin(alpha=.3,trim=FALSE,colour=NA) +
  stat_summary(fun.data=function(z) { z <- z[is.finite(z)]; data.frame(y=median(z),ymin=unname(quantile(z,.25,type=7)),ymax=unname(quantile(z,.75,type=7))) },geom="linerange",linewidth=.25,colour=pfos_palette("neutral")[["dark"]],show.legend=FALSE) +
  stat_summary(fun=median,geom="point",size=1.1,colour=pfos_palette("neutral")[["dark"]],show.legend=FALSE) +
  facet_wrap(~Layer) + scale_fill_manual(values=cols) + labs(title="RNA effect distribution by accessibility quadrant",x=NULL,y="RNA log2FC") + theme(axis.text.x=element_text(angle=45,hjust=1))
joint_plot(j,"atac_rna_gene_violin","各开放-表达象限RNA效应分布",p,gene_table,width=9,height=6)
distal <- peak %>% filter(Layer=="Distal")
p <- ggplot(distal,aes(ATACDirection,log2FoldChange,fill=ATACDirection)) + geom_boxplot(outlier.alpha=.2) +
  scale_fill_manual(values=c(Open=unname(pfos_palette("de")[["up"]]),Closed=unname(pfos_palette("de")[["down"]]))) + labs(title="RNA effects at annotated distal ATAC regions",x=NULL,y="RNA log2FC")
joint_plot(j,"distal_expression_boxplot","远端开放区域最近基因的RNA效应分布",p,distal)
distal_top <- distal %>% arrange(desc(abs(Fold)*abs(log2FoldChange))) %>% slice_head(n=40L) %>%
  mutate(LinkLabel=paste0(coalesce(SYMBOL,gene_id_clean)," | ",peak_id),DistanceKb=distanceToTSS/1000)
p <- ggplot(distal_top,aes(DistanceKb,reorder(LinkLabel,DistanceKb),color=ATACDirection,size=abs(Fold))) +
  geom_segment(aes(x=0,xend=DistanceKb,yend=LinkLabel),linewidth=.4,color=unname(pfos_palette("neutral")[["missing"]])) + geom_point() +
  scale_color_manual(values=c(Open=unname(pfos_palette("de")[["up"]]),Closed=unname(pfos_palette("de")[["down"]]))) + labs(title="Distal peak / nearest-gene annotation distances",subtitle="Distance annotation is not a validated enhancer-promoter link",x="Signed distance to annotated TSS (kb)",y=NULL) + theme(axis.text.y=element_text(size=7))
joint_plot(j,"distal_peak_gene_links","远端peak到最近基因的注释距离",p,distal_top,height=9)
primary <- gene %>% group_by(gene_id_clean) %>% arrange(desc(abs(Fold)),peak_id,.by_group=TRUE) %>% slice_head(n=1L) %>% ungroup()
wm <- joint_gene_summary(bind_rows(dmr,dmr %>% filter(!is.na(CGI_context),CGI_context!="Open_Sea") %>% mutate(Layer="CGI")))
wm <- wm %>% dplyr::select(gene_id_clean,WGBSLayer=Layer,WGBSRegion=DMR_id,WGBSDelta=diffMethy,WGBSRegionCount=RegionCount,WGBSDirectionConflict=DirectionConflict)
joint <- primary %>% left_join(wm,by="gene_id_clean")
joint$MeasuredLayer <- ifelse(is.finite(joint$WGBSDelta),"ATAC + RNA + annotated DMR","ATAC + RNA")
joint$DescriptiveEffectMagnitude <- abs(joint$Fold)+abs(joint$log2FoldChange)
joint_table(j,"multiomics_gene_evidence_full",joint)
top <- joint %>% arrange(desc(DescriptiveEffectMagnitude)) %>% distinct(gene_id_clean,.keep_all=TRUE) %>% slice_head(n=40L)
top$Label <- make.unique(ifelse(is.na(top$SYMBOL),top$gene_id_clean,top$SYMBOL))
p <- ggplot(top,aes(log2FoldChange,reorder(Label,DescriptiveEffectMagnitude),size=abs(Fold),color=WGBSDelta,shape=MeasuredLayer)) +
  geom_point(alpha=.85) + scale_color_gradientn(colors=unname(pfos_palette("diverging")),na.value=unname(pfos_palette("neutral")[["missing"]])) +
  labs(title="Liver multiomics candidate effect overview",subtitle="No missing WGCNA/splicing/TF-binding score is imputed",x="RNA log2FC",y=NULL,color="Methylation fraction delta",size="|ATAC log2FC|")
joint_plot(j,"multiomics_candidate_bubble","肝脏候选基因实际多组学效应气泡图",p,top,height=9)
expr <- read_tsv(joint_path(j,r$expression),show_col_types=FALSE) %>% mutate(gene_id_clean=sub("\\..*$","",gene_id))
samples <- read_tsv(joint_path(j,r$samples),show_col_types=FALSE)
if(!all(c("sample","group") %in% names(samples)) || anyDuplicated(samples$sample) ||
   !setequal(setdiff(names(expr),c("gene_id","gene_id_clean")),samples$sample) ||
   !setequal(unique(samples$group),c(j$cfg$reference,j$cfg$comparison))) stop("RNA显式样本表与表达样本/配置组别不一致")
expr_long <- expr %>% filter(gene_id_clean %in% head(top$gene_id_clean,35L)) %>% dplyr::select(-gene_id) %>%
  pivot_longer(-gene_id_clean,names_to="Sample",values_to="NormalizedExpression") %>%
  left_join(samples %>% dplyr::select(Sample=sample,Group=group),by="Sample") %>%
  mutate(LogExpression=log2(NormalizedExpression+1)) %>% group_by(gene_id_clean) %>%
  mutate(RowZ=if(sd(LogExpression)>0)(LogExpression-mean(LogExpression))/sd(LogExpression) else 0) %>% ungroup()
joint_heat(j,"top_candidates_rna_sample_heatmap","肝脏候选完整cohort归一化RNA样本热图",expr_long,row="gene_id_clean",value="RowZ")
evidence <- top %>% dplyr::select(gene_id_clean,Label,Fold,log2FoldChange,WGBSDelta) %>%
  pivot_longer(c(Fold,log2FoldChange,WGBSDelta),names_to="Omics",values_to="Effect") %>%
  mutate(Label=factor(Label,levels=rev(top$Label)),Direction=case_when(
    !is.finite(Effect) ~ "Missing", Effect==0 ~ "Zero",
    Omics=="Fold" & Effect>0 ~ "Open",Omics=="Fold" & Effect<0 ~ "Closed",
    Omics=="log2FoldChange" & Effect>0 ~ "Up",Omics=="log2FoldChange" & Effect<0 ~ "Down",
    Omics=="WGBSDelta" & Effect>0 ~ "Hyper",TRUE ~ "Hypo"))
p <- ggplot(evidence,aes(Omics,Label,fill=Direction)) + geom_tile(color="white") + geom_text(aes(label=Direction),size=8/ggplot2::.pt) +
  scale_fill_manual(values=c(Open=unname(pfos_palette("de")[["up"]]),Closed=unname(pfos_palette("de")[["down"]]),Up=unname(pfos_palette("de")[["up"]]),Down=unname(pfos_palette("de")[["down"]]),Hyper=unname(pfos_palette("de")[["up"]]),Hypo=unname(pfos_palette("de")[["down"]]),Missing=unname(pfos_palette("neutral")[["missing"]]),Zero=unname(pfos_palette("neutral")[["light"]]))) +
  labs(title="Measured evidence directions",subtitle="Direction only; continuous effects and distinct units are retained in the source table",x=NULL,y=NULL)
joint_plot(j,"multiomics_evidence_matrix","实际ATAC/RNA/WGBS证据矩阵",p,evidence,height=9)
links <- read_tsv(joint_cached(j,"atac_wgbs_coordinate_overlap"),show_col_types=FALSE)
link_table <- joint_table(j,"atac_wgbs_coordinate_overlap",links)
if(nrow(links)) {
  p <- ggplot(links,aes(WGBSDelta,ATACFold,color=SameAnnotatedGene,size=OverlapBp)) + geom_point(alpha=.7) +
    geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2) +
    labs(title="Coordinate-overlapping ATAC regions and DMRs",subtitle="All overlap pairs retained; annotation mismatches shown",x=paste0("Methylation fraction delta (",j$contrast,")"),y="ATAC log2FC")
  joint_plot(j,"atac_wgbs_coordinate_scatter","坐标重叠ATAC-DMR效应散点",p,link_table)
} else joint_skip(j,"atac_wgbs_coordinate_scatter","坐标重叠ATAC-DMR效应散点","所选区域无坐标重叠","EMPTY",link_table)
p <- ggplot(joint %>% filter(is.finite(WGBSDelta)),aes(WGBSDelta,Fold,color=log2FoldChange)) + geom_point(alpha=.65,size=1.3) +
  facet_wrap(~WGBSLayer) + scale_color_gradientn(colors=unname(pfos_palette("diverging"))) +
  labs(title="Annotated gene-level ATAC / methylation effects",subtitle="Nearest-gene correspondence does not require coordinate overlap",x="Methylation fraction delta",y="ATAC log2FC",color="RNA log2FC")
joint_plot(j,"atac_wgbs_gene_centric_scatter","肝脏ATAC-WGBS同注释基因效应",p,joint %>% filter(is.finite(WGBSDelta)),width=9)
# 富集背景是共识peak可关联基因与当前RNA可评估基因的交集。
sets <- joint_gene_sets(j)
bg <- read_csv(joint_path(j,a$universe),show_col_types=FALSE)
db <- get(j$cfg$orgdb,envir=asNamespace(j$cfg$orgdb))
mapping <- AnnotationDbi::select(db,keys=intersect(as.character(bg$ENTREZID),AnnotationDbi::keys(db,keytype="ENTREZID")),keytype="ENTREZID",columns="ENSEMBL")
universe <- intersect(mapping$ENSEMBL,rna$gene_id_clean[is.finite(rna$pvalue)&is.finite(rna$log2FoldChange)])
joint_table(j,"atac_rna_enrichment_universe",data.frame(gene_id_clean=universe))
enrich <- list()
for(quad in names(cols)) for(collection in c("GO_BP","KEGG")) {
  selected <- gene %>% filter(Quadrant==quad,rna_selected)
  z <- joint_ora_figures(j,paste0(tolower(gsub("-","_",quad)),"_",tolower(collection)),selected$gene_id_clean,universe,sets[[collection]],paste("Liver ATAC-RNA",quad,collection))
  if(nrow(z)){z$Quadrant<-quad;z$Collection<-collection;enrich[[paste(quad,collection)]]<-z}
}
enrich <- bind_rows(enrich);joint_table(j,"atac_rna_quadrant_enrichment_full",enrich)
if(nrow(enrich)) {
  z <- enrich %>% filter(pvalue<j$cfg$p_cutoff,Count>0) %>% arrange(pvalue) %>% filter(term %in% head(unique(term),25L))
  if(nrow(z)) {
    p <- ggplot(z,aes(Quadrant,Description,fill=-log10(pmax(pvalue,1e-300)))) + geom_tile() +
      scale_fill_gradientn(colors=unname(pfos_palette("sequential"))) + labs(title="Liver ATAC-RNA quadrant pathways",subtitle="Exploratory raw-P; BH retained",x=NULL,y=NULL)
    joint_plot(j,"atac_rna_quadrant_pathway_heatmap","肝脏ATAC-RNA四象限通路热图",p,z,height=8)
  }
}
# 原生共识interval全集提供完整ATAC排序；重新注释，不把所选差异区域作为GSEA全集。
all_sites <- read_csv(joint_path(j,a$all_sites),show_col_types=FALSE)
if(!all(c("seqnames","start","end","Fold","peak_id") %in% names(all_sites))) stop("完整ATAC源表缺少坐标或效应字段")
full_key <- paste(sub("^chr","",all_sites$seqnames),all_sites$start,all_sites$end,sep=":")
selected_key <- paste(sub("^chr","",atac$seqnames),atac$start,atac$end,sep=":")
selected_match <- match(selected_key,full_key)
if(anyDuplicated(full_key) || anyNA(selected_match) || any(!is.finite(all_sites$Fold[selected_match])) ||
   any(abs(atac$Fold-all_sites$Fold[selected_match])>1e-8)) stop("AllSites必须覆盖全部现役显著ATAC坐标，且Fold与已核定方向的现役表一致")
joint_table(j,"complete_atac_direction_source_check",data.frame(peak_id=atac$peak_id,seqnames=atac$seqnames,start=atac$start,end=atac$end,selected_fold=atac$Fold,all_sites_fold=all_sites$Fold[selected_match],comparison=j$cfg$comparison,reference=j$cfg$reference))
full_rank <- read_tsv(joint_cached(j,"complete_atac_gene_ranking_source"),show_col_types=FALSE)
joint_table(j,"complete_atac_gene_ranking_source",full_rank)
atac_rank <- setNames(full_rank$Fold,full_rank$gene_id_clean)
rna_rank <- setNames(rna$stat,rna$gene_id_clean)
for(collection in c("GO_BP","KEGG")) {
  za <- joint_gsea(j,paste0("complete_atac_",tolower(collection)),atac_rank,sets[[collection]],paste("Liver complete ATAC",collection,"effect rank"))
  zr <- joint_gsea(j,paste0("complete_rna_",tolower(collection)),rna_rank,sets[[collection]],paste("Liver current full-cohort RNA",collection,"statistic rank"))
  if(is.data.frame(za) && is.data.frame(zr)) {
    comp <- inner_join(za %>% dplyr::select(pathway,ATAC_NES=NES,ATAC_p=pval,ATAC_fdr=padj),
                       zr %>% dplyr::select(pathway,RNA_NES=NES,RNA_p=pval,RNA_fdr=padj),by="pathway")
    comp <- comp %>% mutate(Concordance=ifelse(sign(ATAC_NES)==sign(RNA_NES),"Same NES direction","Different NES direction"))
    p <- ggplot(comp,aes(ATAC_NES,RNA_NES,color=Concordance)) + geom_point(alpha=.6,size=1.2) +
      geom_hline(yintercept=0,linetype=2) + geom_vline(xintercept=0,linetype=2) +
      labs(title=paste("Liver",collection,"complete-ranking NES comparison"),subtitle="Across-pathway descriptive comparison; no matched-sample or activity claim",x="ATAC normalized enrichment score",y="RNA normalized enrichment score")
    joint_plot(j,paste0(tolower(collection),"_pathway_nes_scatter"),paste(collection,"完整排序ATAC与RNA通路NES比较"),p,comp)
    top_comp <- comp %>% mutate(best=pmin(ATAC_fdr,RNA_fdr)) %>% arrange(best) %>% slice_head(n=20L) %>%
      pivot_longer(c(ATAC_NES,RNA_NES),names_to="Omics",values_to="NES")
    p <- ggplot(top_comp,aes(Omics,stringr::str_wrap(gsub("_"," ",pathway),45),color=NES,size=abs(NES))) + geom_point() +
      scale_color_gradientn(colors=unname(pfos_palette("diverging"))) +
      labs(title=paste("Liver",collection,"ATAC / RNA pathway ranks"),x=NULL,y=NULL)
    joint_plot(j,paste0(tolower(collection),"_pathway_nes_dotplot"),paste(collection,"ATAC与RNA完整排序通路点图"),p,top_comp,height=8)
  } else joint_skip(j,paste0(tolower(collection),"_pathway_nes_scatter"),paste(collection,"完整排序ATAC-RNA NES比较"),"完整排序或符合大小阈值的基因集不足","NOT_RUN")
}
joint_skip(j,"binding_tf_network","结合TF→gene网络","BINDetect结果UNVERIFIED，motif命中及原始Tn5描述不能升级为TF结合证据")
joint_skip(j,"wgcna_splicing_priority","WGCNA/剪接候选评分","当前无对应有效产物；缺失保持NOT_RUN，不填0评分")
joint_skip(j,"sample_level_crossomics_correlation","跨组学样本相关","动物/父本/窝及组学间配对unknown，禁止按样本名拼接相关性")
# 真实坐标区域+RNA样本+直接CpG样本组成位点图，不伪造连续轨道。
sup <- read_csv(joint_path(j,w$sample_methylation),show_col_types=FALSE)
for(i in seq_len(min(12L,nrow(top)))) {
  target <- top$gene_id_clean[i]; label <- top$Label[i]
  ap <- atac %>% filter(gene_id_clean==target) %>% transmute(Region=peak_id,Start=start,End=end,Omics="ATAC",Effect=Fold)
  wp <- dmr %>% filter(gene_id_clean==target) %>% transmute(Region=DMR_id,Start=start,End=end,Omics="WGBS DMR",Effect=diffMethy)
  pos <- bind_rows(ap,wp)
  p1 <- ggplot(pos,aes(x=Start,xend=End,y=Region,yend=Region,color=Effect)) + geom_segment(linewidth=2.5) + facet_grid(Omics~.,scales="free_y",space="free_y") +
    scale_color_gradientn(colors=unname(pfos_palette("diverging"))) + labs(title=paste(label,"annotated genomic regions"),x="GRCm39 coordinate (bp)",y=NULL) + pfos_theme()
  ex <- expr_long %>% filter(gene_id_clean==target)
  p2 <- ggplot(ex,aes(Sample,NormalizedExpression,fill=Group)) + geom_col() + scale_fill_manual(values=pfos_figure_style$colors[c("WT","PFOS")]) + labs(title="Current full-cohort normalized RNA",x=NULL,y="Normalized expression") +
    pfos_theme() + theme(axis.text.x=element_text(angle=60,hjust=1),legend.position="none")
  ds <- sup %>% filter(RegionID %in% wp$Region)
  p3 <- ggplot(ds,aes(Sample,RegionID,fill=Methylation)) + geom_tile() + scale_fill_gradientn(colors=unname(pfos_palette("sequential")),na.value=unname(pfos_palette("neutral")[["missing"]]),limits=unlist(pfos_profile$scales$methylation_percent)) +
    labs(title="Direct-read CpG methylation; grey = no coverage",x=NULL,y=NULL) + pfos_theme() + theme(axis.text.x=element_text(angle=60,hjust=1))
  composite <- cowplot::plot_grid(p1,p2,p3,ncol=1,rel_heights=c(1.5,1,1),labels=c("A","B","C"),label_size=10,label_fontfamily=pfos_figure_style$family)
  source <- bind_rows(pos %>% mutate(panel="region"),ex %>% mutate(panel="RNA"),ds %>% mutate(panel="CpG"))
  joint_plot(j,paste0("locus_",target),paste(label,"注释区域-RNA-直接CpG复合位点图"),composite,source,width=9,height=9)
}
joint_finish(j)
cat("Liver ATAC-RNA-WGBS figure restoration completed\n")
