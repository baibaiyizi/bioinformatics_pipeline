#!/usr/bin/env Rscript
# 当前完整 cohort 的跨情境图；基因集显式冻结，背景限于实际有效检验基因。
suppressPackageStartupMessages({library(ggplot2); library(ggrastr); library(yaml); library(UpSetR); library(patchwork)})
root <- normalizePath(commandArgs(TRUE)[1], mustWork=TRUE)
cfg <- read_yaml(file.path(root,"config/project.yml"))
source(file.path(root,"scripts/joint_plot_helpers.R"))
source(file.path(root,"scripts/figure_style.R"))
j <- joint_init(root,"result/cross_context/rnaseq","cross_context_rna")
out <- joint_result_root(root)
plotfile <- function(id,extension)joint_figure_path(root,"cross_context_rna",id,extension)
read_input <- function(p) read.csv(file.path(root,p),check.names=FALSE,stringsAsFactors=FALSE)
specs <- Filter(function(x) x$kind=="rna",cfg$inputs)
contexts <- vapply(specs,`[[`,"", "context")
contrast <- paste(cfg$comparison,cfg$reference,sep=" - ")
order <- cfg$rna_context_order
if (is.null(order)) order <- contexts
order <- order[order %in% contexts]
dat <- do.call(rbind,lapply(specs,function(s){d<-read_input(s$path);d$context<-s$context;d}))
dat$gene_id <- sub("\\.[0-9]+$","",dat$ensembl_id)
dat$nominal <- is.finite(dat$pvalue)&is.finite(dat$log2FoldChange)&dat$pvalue<cfg$p_cutoff&abs(dat$log2FoldChange)>=cfg$lfc_cutoff
dat$bh <- is.finite(dat$padj)&is.finite(dat$log2FoldChange)&dat$padj<cfg$fdr_cutoff&abs(dat$log2FoldChange)>=cfg$lfc_cutoff
manifest <- list()
rel <- function(p) joint_rel(j,p)
table_out <- function(d,name) rel(joint_cached(j,name,"csv"))
theme_set(pfos_theme())
record <- function(id,source,title,meaning,distribution,conclusion,status="COMPLETED",reason="") {
 manifest[[length(manifest)+1]] <<- data.frame(analysis="cross_context_rna",figure_id=id,png=if(status=="COMPLETED")rel(plotfile(id,"png"))else"",pdf=if(status=="COMPLETED")rel(plotfile(id,"pdf"))else"",source_table=source,status=status,reason=reason,title=title,biological_meaning=meaning,distribution=distribution,conclusion=conclusion)
}
save_plot <- function(p,id,source,title,meaning,distribution,conclusion,w=7.2,h=4.6){
 height<-max(110,183*h/w)
 if(grepl("_all$",id))height<-240
 if(grepl("heatmap|dotplot",id))height<-max(height,180)
 height<-min(height,240)
 p<-p+theme(plot.title=element_blank(),plot.subtitle=element_blank())
 pfos_save(p,plotfile(id,"pdf"),height_mm=height)
 record(id,source,title,meaning,distribution,conclusion)
}
upset <- function(sets,id,source,title,meaning){
 if(length(sets)>1 && length(unique(unlist(sets)))>0){
  draw <- function(){print(UpSetR::upset(UpSetR::fromList(sets),sets=names(sets),keep.order=TRUE,nsets=length(sets),nintersects=20,order.by="freq",text.scale=c(8/8.3,1,8/8.3,1,1,7/(2.2*ggplot2::.pt)),mainbar.y.label="Intersection size",sets.x.label="Set size"),newpage=FALSE)}
  pfos_save(draw,plotfile(id,"pdf"),height_mm=160)
  record(id,source,title,meaning,"柱高为实际集合交集成员数；只展示最多20个非空交集。","集合交汇描述共性，不能证明组织间传递或配对效应。")
 } else if(length(sets)==1 && length(sets[[1]])>0) {
  d<-data.frame(context=names(sets),members=length(sets[[1]]))
  p<-ggplot(d,aes(context,members))+geom_col(fill=unname(pfos_palette("assay")[["rna"]]))+labs(x=NULL,y="Set size",title="Single configured context")
  save_plot(p,id,source,title,meaning,"配置只有一个情境，柱高为实际集合成员数。","单一情境不能计算跨情境交集。")
 } else {
  p<-ggplot()+annotate("text",x=0,y=0,label="No non-empty eligible intersections")+theme_void()
  save_plot(p,id,source,title,meaning,"当前阈值没有非空交集。","如实展示空结果，不降低阈值。")
 }
}
sets <- setNames(lapply(order,function(k)unique(dat$gene_id[dat$context==k & dat$nominal])),order)
members <- do.call(rbind,lapply(names(sets),function(k)data.frame(context=rep(k,length(sets[[k]])),gene_id=sets[[k]])))
src <- table_out(members,"deg_set_members")
upset(sets,"deg_upset",src,"RNA探索候选集合交汇","比较相同阈值下各情境候选的共享与特异。")
wide <- reshape(dat[,c("gene_id","context","log2FoldChange")],idvar="gene_id",timevar="context",direction="wide")
pairs <- if(length(order)>1) combn(order,2,simplify=FALSE) else list()
allpoints <- list()
for(pair in pairs){
 d<-data.frame(gene_id=wide$gene_id,x=wide[[paste0("log2FoldChange.",pair[1])]],y=wide[[paste0("log2FoldChange.",pair[2])]])
 d<-d[is.finite(d$x)&is.finite(d$y),];d$comparison<-paste(pair,collapse=" / ")
 rho<-if(nrow(d)>2)cor(d$x,d$y,method="spearman") else NA_real_
 d$spearman_rho<-rho;allpoints[[length(allpoints)+1]]<-d
 src<-table_out(d,paste0("effect_scatter_",paste(pair,collapse="__")))
 p<-ggplot(d,aes(x,y))+ggrastr::geom_point_rast(size=.45,alpha=.25,color=unname(pfos_palette("assay")[["rna"]]),raster.dpi=300)+geom_hline(yintercept=0,linewidth=.2,color=unname(pfos_palette("neutral")[["mid"]]))+geom_vline(xintercept=0,linewidth=.2,color=unname(pfos_palette("neutral")[["mid"]]))+labs(x=paste(pair[1],"RNA log2FC"),y=paste(pair[2],"RNA log2FC"),title=sprintf("%s; genes = %s; descriptive rho = %.2f",contrast,nrow(d),rho))
 save_plot(p,paste0("effect_scatter_",paste(pair,collapse="__")),src,paste(pair,collapse="／"),"对照两个独立情境的基因效应方向。",paste0("仅保留两侧效应均有限的",nrow(d),"个基因；rho是基因层描述量。"),"点不代表动物；不报告以基因为独立实验重复的显著性。")
}
if(length(allpoints)){
 d<-do.call(rbind,allpoints);src<-table_out(d,"effect_scatter_all")
 p<-ggplot(d,aes(x,y))+ggrastr::geom_point_rast(size=.12,alpha=.15,color=unname(pfos_palette("assay")[["rna"]]),raster.dpi=300)+geom_hline(yintercept=0,linewidth=.15,color=unname(pfos_palette("neutral")[["mid"]]))+geom_vline(xintercept=0,linewidth=.15,color=unname(pfos_palette("neutral")[["mid"]]))+facet_wrap(~comparison,ncol=4,scales="free")+labs(x="Left context RNA log2FC",y="Right context RNA log2FC",title="Current full-cohort RNA effect comparisons")
 save_plot(p,"effect_scatter_all",src,"全部情境两两效应总览","比较共享效应及情境差异。","每个小面板对应一对情境；各面板使用自身有效基因。","跨情境一致性属于描述，不是同一个体随访或因果传递。",w=11,h=12)
}
bundle <- readRDS(file.path(root,cfg$gene_sets))
pathways <- bundle$libraries$Mouse_KEGG
stopifnot(length(pathways)>0)
ora<-read.csv(joint_cached(j,"kegg_ora","csv"),check.names=FALSE)
eligible<-list(read.csv(joint_cached(j,"kegg_eligible_gene_universe","csv"),check.names=FALSE))
src<-table_out(ora,"kegg_ora")
table_out(do.call(rbind,eligible),"kegg_eligible_gene_universe")
# 原跨情境散点比较的是KEGG富集因子，不能用基因RNA效应散点替代。
kegg_points <- list()
for(pair in pairs){
 a<-ora[ora$context==pair[1],];b<-ora[ora$context==pair[2],]
 d<-merge(a,b,by="pathway",suffixes=c("_left","_right"))
 d<-d[is.finite(d$fold_enrichment_left)&is.finite(d$fold_enrichment_right)&is.finite(d$p_left)&is.finite(d$p_right),]
 d$comparison<-paste(pair,collapse=" / ")
 d$support<-ifelse(is.finite(d$fdr_left)&is.finite(d$fdr_right)&d$fdr_left<cfg$fdr_cutoff&d$fdr_right<cfg$fdr_cutoff,"Both BH",
                   ifelse(d$p_left<cfg$p_cutoff&d$p_right<cfg$p_cutoff,"Both raw P","Other tested"))
 d$support<-factor(d$support,levels=c("Other tested","Both raw P","Both BH"))
 kegg_points[[length(kegg_points)+1]]<-d
 pair_src<-table_out(d,paste0("kegg_enrichment_scatter_",paste(pair,collapse="__")))
 p<-ggplot(d,aes(fold_enrichment_left,fold_enrichment_right,color=support))+geom_point(alpha=.65,size=1.3)+
    scale_color_manual(values=c("Other tested"=unname(pfos_palette("neutral")[["mid"]]),"Both raw P"=unname(pfos_palette("accent")[["warm"]]),"Both BH"=unname(pfos_palette("assay")[["rna"]])),drop=FALSE)+
    labs(x=paste(pair[1],"ORA fold enrichment"),y=paste(pair[2],"ORA fold enrichment"),color="Support",title=paste("Current eligible KEGG ORA:",paste(pair,collapse=" / ")))
 save_plot(p,paste0("kegg_enrichment_scatter_",paste(pair,collapse="__")),pair_src,paste(pair,collapse="／"),
           "比较两个情境的KEGG过度代表幅度。","每点为两侧均可检验的同一KEGG条目；颜色区分双侧原始P、双侧BH与其他条目。",
           "富集因子分别基于各情境实际检验背景；原始P交汇属于探索，通路不是独立生物学重复。")
}
if(length(kegg_points)){
 d<-do.call(rbind,kegg_points);pair_src<-table_out(d,"kegg_enrichment_scatter_all")
 p<-ggplot(d,aes(fold_enrichment_left,fold_enrichment_right,color=support))+geom_point(alpha=.5,size=.55)+
    scale_color_manual(values=c("Other tested"=unname(pfos_palette("neutral")[["mid"]]),"Both raw P"=unname(pfos_palette("accent")[["warm"]]),"Both BH"=unname(pfos_palette("assay")[["rna"]])),drop=FALSE)+
    facet_wrap(~comparison,ncol=4,scales="free")+labs(x="Left context ORA fold enrichment",y="Right context ORA fold enrichment",color="Support",title="KEGG ORA comparisons across current cohorts")
 save_plot(p,"kegg_enrichment_scatter_all",pair_src,"全部KEGG富集因子两两比较","浏览当前各情境的功能条目共现。",
           "每个小面板为一对情境；背景、命中数、原始P和BH均保留源表。","共享条目不等于通路活性、跨组织传递或样本层相关。",w=11,h=12)
}
common_counts<-table(ora$pathway[is.finite(ora$p)&ora$p<cfg$p_cutoff])
common<-names(common_counts)[common_counts==length(order)]
common_table<-ora[ora$pathway %in% common,];common_src<-table_out(common_table,"kegg_common_all_contexts")
if(nrow(common_table)){
 score_common<-aggregate(p~pathway,common_table,min)
 keep<-head(score_common$pathway[order(score_common$p)],15)
 d<-common_table[common_table$pathway %in% keep,];d$context<-factor(d$context,levels=order)
 p<-ggplot(d,aes(context,sub("^mmu[0-9]+_","",pathway),size=fold_enrichment,color=-log10(pmax(p,1e-300))))+
    geom_point()+scale_color_gradientn(colors=unname(pfos_palette("sequential")),limits=unlist(pfos_profile$scales$ora_minus_log10_p),oob=scales::squish,breaks=c(0,5,10),labels=c("0","5","≥10"))+labs(x="Independent context cohorts",y=NULL,size="Fold enrichment",color="-log10 raw P",title="KEGG terms with nominal ORA support in every configured context")+
    theme(axis.text.x=element_text(angle=35,hjust=1))
 save_plot(p,"kegg_common_context_dotplot",common_src,"全部情境共同KEGG条目比较","比较所有配置情境均命中的探索性功能条目。",
           "面积为富集因子，颜色为原始P；源表同时保留BH。","所有情境共同原始P支持仍是探索，不等于所有情境BH显著或通路激活。",h=6)
} else {
 record("kegg_common_context_dotplot",common_src,"全部情境共同KEGG条目比较","比较跨情境探索性条目。",
        "没有条目在全部配置情境均达到原始P阈值。","保留空结果，不删减情境或放宽阈值。",status="EMPTY",reason="所有配置情境的nominal KEGG交集为空")
}
ksets<-setNames(lapply(order,function(k)ora$pathway[ora$context==k & is.finite(ora$fdr)&ora$fdr<cfg$fdr_cutoff]),order)
upset(ksets,"kegg_upset",src,"KEGG富集条目交汇","比较当前探索候选在可检验背景中的KEGG富集。")
score<-aggregate(fdr~pathway,ora,min,na.rm=TRUE);top<-head(score$pathway[order(score$fdr)],20)
d<-ora[ora$pathway %in% top,];d$context<-factor(d$context,levels=order);d$pathway_label<-sub("^mmu[0-9]+_","",d$pathway);d$minus_log10_fdr<- -log10(pmax(d$fdr,1e-300));d$ratio<-d$overlap/pmax(d$selected,1)
p<-ggplot(d,aes(context,pathway_label))+geom_point(aes(size=overlap,color=minus_log10_fdr))+scale_color_gradientn(colors=unname(pfos_palette("sequential")))+scale_size_area(max_size=7)+labs(x=NULL,y=NULL,color="-log10 BH",size="Overlap",title="KEGG ORA: top 20 by minimum BH across contexts")+theme(axis.text.x=element_text(angle=35,hjust=1))
save_plot(p,"kegg_dotplot",src,"KEGG富集点图","描述各情境候选涉及的功能条目。","面积为命中基因数，颜色为BH；显示按最小BH排序的20条。","背景为该情境实际有效检验且有SYMBOL的基因；无显著条目仍保留原BH。",h=8)
p<-ggplot(d,aes(context,pathway_label,fill=minus_log10_fdr))+geom_tile()+scale_fill_gradientn(colors=unname(pfos_palette("sequential")),na.value=unname(pfos_palette("neutral")[["missing"]]))+labs(x=NULL,y=NULL,fill="-log10 BH",title="KEGG ORA, eligible tested-gene universe")+theme(axis.text.x=element_text(angle=35,hjust=1))
save_plot(p,"kegg_heatmap",src,"KEGG富集热图","比较功能条目的情境依赖。","色深为同一ORA表的负对数BH；非显著值也显示。","通路成员重叠且不是独立发现；不能将ORA解释成通路活性。",h=8)
panel<-cfg$rna_gene_panel
if(length(panel)){
 d<-dat[dat$SYMBOL %in% panel,c("context","gene_id","SYMBOL","log2FoldChange","pvalue","padj")]
 d$context<-factor(d$context,levels=order);d$SYMBOL<-factor(d$SYMBOL,levels=rev(panel));src<-table_out(d,"methylation_regulator_rna_panel")
 lim<-max(abs(d$log2FoldChange[is.finite(d$log2FoldChange)]),1)
 p<-ggplot(d,aes(context,SYMBOL,fill=log2FoldChange))+geom_tile(color="white")+scale_fill_gradientn(colors=unname(pfos_palette("diverging")),limits=c(-lim,lim),na.value=unname(pfos_palette("neutral")[["missing"]]))+labs(x=NULL,y=NULL,fill="RNA log2FC",title="Methylation regulators: RNA expression effects")
 save_plot(p,"methylation_regulator_rna_heatmap",src,"甲基化调节基因的RNA效应热图","观察预先指定调节基因的转录响应。","颜色来自RNA log2FC；缺失灰色保留，不填0。","这里测量的是表达，不能称为基因甲基化或酶活性。")
 embryo<-cfg$embryo_contexts
 e<-d[as.character(d$context)%in%embryo,];e$context<-factor(e$context,levels=embryo)
 if(nrow(e)){
 p<-ggplot(e,aes(context,log2FoldChange,group=SYMBOL,color=SYMBOL))+geom_hline(yintercept=0,color=unname(pfos_palette("neutral")[["mid"]]))+geom_line(na.rm=TRUE)+geom_point(na.rm=TRUE)+scale_color_manual(values=pfos_identity_colors(unique(e$SYMBOL)))+facet_wrap(~SYMBOL,ncol=3)+guides(color="none")+labs(x="Independent embryo stage cohorts",y="RNA log2FC",title="Stage-ordered RNA effects; not longitudinal animals")
 save_plot(p,"methylation_regulator_embryo_lines",src,"早胚阶段调节基因RNA效应","按真实早胚阶段顺序比较表达响应。","只连接配置中的早胚阶段，各阶段为独立样本。","连线辅助阅读，不能作为同一个体的纵向变化。",h=6)
 }
 a<-d[!as.character(d$context)%in%embryo,]
 if(nrow(a)){
 p<-ggplot(a,aes(context,log2FoldChange,color=context))+geom_hline(yintercept=0,color=unname(pfos_palette("neutral")[["mid"]]))+geom_point(size=2,na.rm=TRUE)+scale_color_manual(values=pfos_palette("context"))+facet_wrap(~SYMBOL,ncol=3)+guides(color="none")+labs(x="Separate tissue cohorts",y="RNA log2FC",title="Tissue-specific RNA effects")+theme(axis.text.x=element_text(angle=35,hjust=1))
 save_plot(p,"methylation_regulator_tissue_points",src,"成体组织调节基因RNA效应","并列比较各组织转录响应。","组织仅用散点并列，无发育阶段连线；睾丸F0。","组织差异不构成发育轨迹或甲基化遗传证据。",h=6)
 }
}
write.table(do.call(rbind,manifest),file.path(out,"00_audit/render_cross_context_rna.tsv"),sep="\t",quote=TRUE,row.names=FALSE,na="")
writeLines(capture.output(sessionInfo()),file.path(out,"logs/rna_session.txt"))
cat("RNA figure pairs:",length(manifest),"\n")
