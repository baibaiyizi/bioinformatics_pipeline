# 联合分析图件共用函数；单组学成熟区域汇总逻辑来自 WGBS build_region_sample_methylation。
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(ggalluvial)
  library(GenomicRanges); library(GenomeInfoDb)
})
joint_helper_source <- sys.frame(1)$ofile
if(is.null(joint_helper_source)) stop("请使用source加载joint_plot_helpers.R")
source(file.path(dirname(joint_helper_source),"figure_style.R"))
joint_plot_mode <- function() identical(Sys.getenv("PFOS_PLOTS_ONLY"), "1")
joint_result_root <- function(root) {
  value <- Sys.getenv("FROZEN_RESULT_DIR")
  if(nzchar(value)) normalizePath(value,mustWork=TRUE) else file.path(root,"result")
}
.joint_frozen_cache <- new.env(parent=emptyenv())
joint_catalog <- function(root) {
 if(!exists("catalog",.joint_frozen_cache,inherits=FALSE)) .joint_frozen_cache$catalog<-readr::read_tsv(file.path(joint_result_root(root),"00_audit/tables/source_catalog.tsv"),show_col_types=FALSE)
 .joint_frozen_cache$catalog
}
joint_cached <- function(j,name,extension="tsv") {
  key <- paste0(j$key,"/",name,".",extension)
  cat <- joint_catalog(j$root);hit <- cat$path[cat$source_key==key]
  if(length(hit)!=1L) stop("缺少唯一冻结数值源: ",key)
  p<-file.path(joint_result_root(j$root),hit)
  if(!file.exists(p)) stop("MISSING冻结源: ",key)
  p
}
joint_figure_path <- function(root,analysis,id,extension) {
  if(!exists("manifest",.joint_frozen_cache,inherits=FALSE)) .joint_frozen_cache$manifest<-readr::read_tsv(file.path(joint_result_root(root),"00_audit/tables/figure_manifest.tsv"),show_col_types=FALSE)
  m<-.joint_frozen_cache$manifest
  hit<-m[m$analysis==analysis & m$figure_id==id & m$status=="COMPLETED",,drop=FALSE]
  if(nrow(hit)!=1L) stop("现役图身份非唯一: ",analysis,"/",id)
  p<-file.path(joint_result_root(root),sub("^result/","",hit[[extension]]))
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);p
}
joint_init <- function(root, relative, analysis) {
  root <- normalizePath(root, mustWork = TRUE)
  if(!joint_plot_mode()) stop("联合图件仅通过 --plots-only 或 --integrate 从冻结结果绘制")
  e <- new.env(parent=globalenv());e$root<-root;e$cfg<-yaml::read_yaml(file.path(root,"config/project.yml"))
  e$key<-relative;e$out<-joint_result_root(root);e$analysis<-analysis;e$rows<-list()
  e$contrast<-paste(e$cfg$comparison,"-",e$cfg$reference);e
}
joint_path <- function(j, value, required = TRUE) {
  if(is.null(value) || !length(value) || !nzchar(value)) {
    if(required) stop("缺少显式输入配置")
    return(NA_character_)
  }
  p <- if(grepl("^/", value)) value else file.path(j$root, value)
  normalizePath(p, mustWork = required, winslash = "/")
}
joint_input <- function(j, id) {
  hits <- Filter(function(x) identical(x$id, id), j$cfg$inputs)
  if(length(hits) != 1L) stop("输入ID必须唯一: ", id)
  hits[[1L]]
}
joint_rel <- function(j,p) {
  prefix<-paste0(joint_result_root(j$root),"/")
  if(startsWith(p,prefix)) paste0("result/",substring(p,nchar(prefix)+1L)) else p
}
joint_table <- function(j,name,x) {
  # x 仅来自现有绘图函数的显示整理；既有科学表格不重写。
  joint_cached(j,name)
}
joint_record <- function(j, id, title, table = "", status = "COMPLETED", reason = "", meaning = "", distribution = "", conclusion = "", png = "", pdf = "") {
  if(!nzchar(meaning)) meaning <- paste0(title,"；描述实际可测量的区域/基因效应或基因集成员关系。")
  if(!nzchar(distribution)) {
    z <- if(nzchar(table) && file.exists(table)) if(grepl("\\.csv$",table)) readr::read_csv(table,col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE) else readr::read_tsv(table,col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE) else data.frame()
    distribution <- paste0("源表保留",nrow(z),"条记录；缺失覆盖保留NA。")
    for(n in intersect(c("diffMethy","WGBSDelta","Fold","log2FoldChange","pvalue","padj","Methylation"),names(z))) {
      v<-as.numeric(z[[n]]);v<-v[is.finite(v)]
      if(length(v)) distribution<-paste0(distribution,n,":",format(min(v),digits=3),"至",format(max(v),digits=3),"；")
    }
  }
  if(!nzchar(conclusion)) conclusion <- "探索性描述；未配对组学不支持样本相关、调控因果或甲基化传递。"
  risk <- if(grepl("2cell",j$analysis)) "高风险低覆盖2cell：仅敏感性展示；" else if(grepl("sperm",j$analysis) && grepl("rna",j$analysis)) "sperm WGBS–testis RNA跨细胞类型；睾丸F0；" else ""
  conclusion <- paste0(risk,conclusion)
  j$rows[[length(j$rows) + 1L]] <- data.frame(analysis = j$analysis, figure_id = id,
    png = joint_rel(j, png), pdf = joint_rel(j, pdf), source_table = joint_rel(j, table),
    status = status, reason = reason, title = title, biological_meaning = meaning,
    distribution = distribution, conclusion = conclusion, stringsAsFactors = FALSE)
}
joint_plot <- function(j,id,title,plot,data,width=7.2,height=5.4,...) {
  table<-if(is.character(data)&&length(data)==1L)data else joint_table(j,id,data)
  pdf<-joint_figure_path(j$root,j$analysis,id,"pdf");png<-sub("\\.pdf$",".png",pdf)
  if(!exists("pfos_save",mode="function"))source(file.path(j$root,"scripts/figure_style.R"))
  # 原图型、分面、组合及展示成员保留；更长标签通过高度及留白容纳。
  height_mm<-max(100,183*height/width)
  if(inherits(plot,"ggplot"))plot<-plot+theme(plot.title=element_blank(),plot.subtitle=element_blank())
  if(grepl("heatmap|annotation|cluster|cnet|tree|chord|locus|bubble|links",id)) height_mm<-max(height_mm,180)
  if(grepl("_tree$",id))height_mm<-max(height_mm,220)
  if(grepl("_(cluster|cnet)$",id))height_mm<-max(height_mm,200)
  if(grepl("_chord$",id))height_mm<-max(height_mm,183)
  height_mm<-min(height_mm,240)
  pfos_save(plot,pdf,height_mm=height_mm)
  joint_record(j,id,title,table,png=png,pdf=pdf,...);invisible(table)
}
joint_skip <- function(j, id, title, reason, status = "NOT_RUN", table = "") {
  joint_record(j, id, title, table, status, reason)
}
joint_finish <- function(j) {
  z<-bind_rows(j$rows)
  dir<-file.path(joint_result_root(j$root),"00_audit")
  write_tsv(z,file.path(dir,paste0("render_",j$analysis,".tsv")),na="")
  writeLines(capture.output(sessionInfo()),file.path(joint_result_root(j$root),"logs",paste0(j$analysis,"_session.txt")))
  invisible(z)
}
joint_rna <- function(j, item) {
  x <- read_csv(joint_path(j, item$path), show_col_types = FALSE)
  if(any(!is.na(x$sampleA) & x$sampleA != j$cfg$comparison) || any(!is.na(x$sampleB) & x$sampleB != j$cfg$reference)) stop("RNA方向不匹配")
  x$gene_id_clean <- sub("\\..*$", "", x$ensembl_id)
  x$rna_selected <- is.finite(x$pvalue) & x$pvalue < j$cfg$p_cutoff & is.finite(x$log2FoldChange) & abs(x$log2FoldChange) >= j$cfg$lfc_cutoff
  x$rna_fdr <- is.finite(x$padj) & x$padj < j$cfg$fdr_cutoff & abs(x$log2FoldChange) >= j$cfg$lfc_cutoff
  x %>% distinct(gene_id_clean, .keep_all = TRUE)
}
joint_layer <- function(annotation) {
  case_when(grepl("^Promoter", annotation) ~ "Promoter", grepl("Intron|Exon|UTR", annotation) ~ "Gene Body", TRUE ~ "Distal")
}
joint_dmr <- function(j, item) {
  x <- read_csv(joint_path(j, item$path), show_col_types = FALSE)
  if(any(x$group1_label != j$cfg$reference) || any(x$group2_label != j$cfg$comparison)) stop("WGBS方向不匹配")
  if(any(abs((x$meanMethy_group2 - x$meanMethy_group1) - x$diffMethy) > 1e-7, na.rm = TRUE)) stop("DMR效应与组均值不一致")
  x$gene_id_clean <- sub("\\..*$", "", x$gene_id_clean)
  x$Layer <- joint_layer(x$annotation)
  x
}
joint_quadrant <- function(m, r, delta = .1, lfc = .58) {
  case_when(m > delta & r > lfc ~ "Hyper-Up", m > delta & r < -lfc ~ "Hyper-Down",
    m < -delta & r > lfc ~ "Hypo-Up", m < -delta & r < -lfc ~ "Hypo-Down", TRUE ~ "NS")
}
joint_colors <- pfos_palette("quadrant")
joint_gene_summary <- function(x, effect = "diffMethy", region = "DMR_id") {
  x %>% group_by(gene_id_clean, Layer) %>%
    mutate(RegionCount = n_distinct(.data[[region]]), DirectionConflict = any(.data[[effect]] > 0) & any(.data[[effect]] < 0)) %>%
    arrange(desc(abs(.data[[effect]])), .data[[region]], .by_group = TRUE) %>% slice_head(n = 1L) %>% ungroup()
}
joint_heat <- function(j, id, title, long, row = "RegionID", sample = "Sample", value = "Methylation", ...) {
  if(n_distinct(long[[row]]) < 2L) return(joint_skip(j, id, title, "少于2个有支持的区域/基因", "EMPTY", joint_table(j,id,long)))
  long[[row]] <- factor(long[[row]], levels = rev(unique(long[[row]])))
  p <- ggplot(long, aes(x = .data[[sample]], y = .data[[row]], fill = .data[[value]])) + geom_tile()
  if (value == "RowZ") {
    p <- p + scale_fill_gradientn(colours=unname(pfos_palette("diverging")),limits=unlist(pfos_profile$scales$row_zscore),oob=scales::squish,breaks=c(-2.5,0,2.5),labels=c("≤-2.5","0","≥2.5"),na.value=pfos_profile$scales$missing_color,name="row z-score")
  } else if (value == "Methylation") {
    p <- p + scale_fill_gradientn(colours=unname(pfos_palette("sequential")),limits=unlist(pfos_profile$scales$methylation_percent),oob=scales::squish,na.value=pfos_profile$scales$missing_color,name="Methylation (%)")
  } else {
    limit <- max(abs(long[[value]]),na.rm=TRUE); if(!is.finite(limit)||limit==0) limit<-1
    p <- p + scale_fill_gradientn(colours=unname(pfos_palette("diverging")),limits=c(-limit,limit),oob=scales::squish,na.value=pfos_profile$scales$missing_color,name=value)
  }
  p <- p + labs(title = title, x = NULL, y = NULL) + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 50, hjust = 1))
  joint_plot(j,id,title,p,long,height=max(4.5,(45+3*n_distinct(long[[row]]))/25.4), ...)
}
# 直接由覆盖记录构建区域样本汇总，不从DSS平滑值或已选DMR推断覆盖。
joint_region_support <- function(raw_list, samples, regions, ids) {
  if(length(regions) != length(ids)) stop("区域ID长度不符")
  lapply(seq_along(raw_list), function(i) {
    raw <- methylKit::getData(raw_list[[i]])
    raw <- raw[is.finite(raw$coverage) & raw$coverage > 0, c("chr","start","end","coverage","numCs")]
    gr <- GRanges(sub("^chr", "", as.character(raw$chr)), IRanges(raw$start, raw$end))
    ov <- findOverlaps(regions, gr, ignore.strand = TRUE)
    z <- data.table::data.table(RegionID = ids[queryHits(ov)], Coverage = raw$coverage[subjectHits(ov)], NumCs = raw$numCs[subjectHits(ov)])
    z <- z[, .(Coverage = sum(Coverage), NumCs = sum(NumCs), n_cpg_records = .N), by = RegionID]
    z <- as.data.frame(z); z$Sample <- rep(samples$sample[i],nrow(z)); z$Group <- rep(samples$group[i],nrow(z))
    z$Methylation <- 100 * z$NumCs/z$Coverage
    z
  }) %>% bind_rows()
}
joint_coverage <- function(j,item,out_prefix,include_cgi=TRUE) {
  support<-read_tsv(joint_cached(j,paste0(out_prefix,"_coverage_support")),show_col_types=FALSE)
  means<-read_tsv(joint_cached(j,paste0(out_prefix,"_coverage_eligibility")),show_col_types=FALSE)
  result<-list(support=support,means=means,table=joint_cached(j,paste0(out_prefix,"_coverage_support")))
  if(include_cgi) {
    result$cgi_support<-read_tsv(joint_cached(j,paste0(out_prefix,"_cgi_coverage_support")),show_col_types=FALSE)
    result$cgi_means<-read_tsv(joint_cached(j,paste0(out_prefix,"_cgi_coverage_eligibility")),show_col_types=FALSE)
  }
  result
}
joint_gene_sets <- function(j) {
  x <- readRDS(joint_path(j,j$cfg$gene_sets))
  if(!is.list(x$libraries) || !all(c("Mouse_KEGG","Mouse_M5") %in% names(x$libraries))) stop("静态基因集结构不符合已核验协议")
  libs <- list(KEGG=x$libraries$Mouse_KEGG,GO_BP=x$libraries$Mouse_M5[grepl("^GOBP_",names(x$libraries$Mouse_M5))],
    GO_MF=x$libraries$Mouse_M5[grepl("^GOMF_",names(x$libraries$Mouse_M5))],GO_CC=x$libraries$Mouse_M5[grepl("^GOCC_",names(x$libraries$Mouse_M5))])
  if(!length(libs$GO_BP)) stop("静态Mouse_M5未提供GO BP集合")
  values <- unique(unlist(libs,use.names=FALSE))
  if(all(grepl("^ENSMUSG",values))) {
    convert <- function(v) sub("\\..*$","",as.character(v))
  } else {
    keytype <- if(all(grepl("^[0-9]+$",values))) "ENTREZID" else "SYMBOL"
    db <- get(j$cfg$orgdb,envir=asNamespace(j$cfg$orgdb))
    valid <- intersect(values,AnnotationDbi::keys(db,keytype=keytype))
    map <- AnnotationDbi::select(db,keys=valid,keytype=keytype,columns="ENSEMBL")
    map <- map[!is.na(map$ENSEMBL),,drop=FALSE]
    lookup <- split(map$ENSEMBL,map[[keytype]])
    convert <- function(v) unique(unlist(lookup[as.character(v)],use.names=FALSE))
  }
  lapply(libs,function(lib) lapply(lib,convert))
}
joint_ora <- function(...) stop("纯绘图禁止运行富集检验；必须读取已冻结 ORA 表")
joint_ora_figures <- function(j,prefix,genes,universe,sets,title) {
  z <- read_tsv(joint_cached(j,paste0(prefix,"_ora")),show_col_types=FALSE)
  table <- joint_table(j,paste0(prefix,"_ora"),z)
  top <- z %>% filter(Count>0,pvalue<j$cfg$p_cutoff) %>% slice_head(n=15L)
  families <- if(grepl("go_(mf|cc)(_|$)",prefix)) c("dotplot","barplot","tree") else c("dotplot","barplot","tree","cnet","chord","cluster")
  if(!nrow(top)) {
    for(f in families) joint_skip(j,paste0(prefix,"_",f),paste(title,f),"按预设raw-P探索阈值没有可展示条目；完整检验与BH仍保留","EMPTY",table)
    return(z)
  }
  top$Label <- factor(stringr::str_wrap(top$Description,45),levels=rev(stringr::str_wrap(top$Description,45)))
  top$GeneRatio <- top$Count / pmax(top$gene_n, 1)
  p <- ggplot(top,aes(GeneRatio,Label,size=Count,color=-log10(pmax(pvalue,1e-300)))) + geom_point() +
    scale_color_gradientn(colours=unname(pfos_palette("sequential")),limits=unlist(pfos_profile$scales$ora_minus_log10_p),oob=scales::squish,breaks=c(0,5,10),labels=c("0","5","≥10")) + scale_size_area(max_size=7) +
    labs(title=paste(title,"exploratory raw P < 0.05"),x="Gene ratio",y=NULL,color="-log10(raw P)",size="Count")
  joint_plot(j,paste0(prefix,"_dotplot"),paste(title,"dotplot"),p,table,height=6.5)
  p <- ggplot(top,aes(Count,Label,fill=-log10(pmax(pvalue,1e-300)))) + geom_col() +
    scale_fill_gradientn(colours=unname(pfos_palette("sequential")),limits=unlist(pfos_profile$scales$ora_minus_log10_p),oob=scales::squish,breaks=c(0,5,10),labels=c("0","5","≥10")) + labs(title=title,x="Observed genes",y=NULL,fill="-log10(raw P)")
  joint_plot(j,paste0(prefix,"_barplot"),paste(title,"barplot"),p,table,height=6.5)
  member <- lapply(top$geneID,function(x) strsplit(x,"/",fixed=TRUE)[[1]])
  names(member) <- as.character(top$term)
  hits <- unique(unlist(member)); mat <- matrix(vapply(member,function(x) as.numeric(hits %in% x),numeric(length(hits))),nrow=length(hits),ncol=length(member))
  dimnames(mat) <- list(hits,top$Description)
  if(ncol(mat)>1L) {
    sim <- outer(seq_along(member),seq_along(member),Vectorize(function(a,b)length(intersect(member[[a]],member[[b]]))/max(1,length(union(member[[a]],member[[b]])))))
    hc <- hclust(as.dist(1-sim),method="average")
    hc$labels <- stringr::str_wrap(top$Description,30)
    if(any(sim[lower.tri(sim)]<1)) {
      leaf_lines<-max(lengths(strsplit(hc$labels,"\n",fixed=TRUE)))
      tree_height_mm<-max(220,50+(length(hc$labels)-1)*(leaf_lines*3.4+1))
      joint_plot(j,paste0(prefix,"_tree"),paste(title,"gene-overlap tree"),function(){par(mar=c(4,2,2,23),family=pfos_figure_style$family);plot(as.dendrogram(hc),horiz=TRUE,main="",xlab="1 - gene membership Jaccard",cex=7/8)},table,height=tree_height_mm*7.2/183)
    } else joint_skip(j,paste0(prefix,"_tree"),paste(title,"gene-overlap tree"),"所有展示术语的当前命中基因完全相同，Jaccard距离全为0，不能形成有意义的层次","EMPTY",table)
    if("cluster" %in% families) {
    long <- expand_grid(Term=colnames(mat),Gene=rownames(mat)); long$Membership <- factor(as.vector(mat),levels=c(0,1))
    long <- long %>% mutate(Term=factor(stringr::str_wrap(as.character(Term),30),levels=stringr::str_wrap(colnames(mat)[hc$order],30)))
    p <- ggplot(long,aes(Term,Gene,fill=Membership)) + geom_tile() +
      scale_fill_manual(values=c("0"=pfos_palette("binary")[["zero"]],"1"=pfos_palette("binary")[["one"]])) +
      labs(title=paste(title,"membership clusters"),x=NULL,y=NULL) + theme(axis.text.x=element_text(angle=70,hjust=1),axis.text.y=element_blank(),plot.margin=margin(12,40,12,40,unit="pt"))
    joint_plot(j,paste0(prefix,"_cluster"),paste(title,"membership clusters"),p,joint_table(j,paste0(prefix,"_membership"),long),height=6.5)
    }
  } else {
    for(f in intersect(c("tree","cluster"),families)) joint_skip(j,paste0(prefix,"_",f),paste(title,f),"只有一个合格条目，不能聚类","EMPTY",table)
  }
  if(!any(c("cnet","chord") %in% families)) return(z)
  edges <- bind_rows(lapply(seq_along(member),function(i)data.frame(from=top$Description[i],to=member[[i]])))
  net <- igraph::graph_from_data_frame(edges,directed=FALSE)
  # 确定性布局只展示实际术语-基因成员关系；不把成员图解释为调控网络。
  set.seed(20261012); xy <- igraph::layout_with_fr(net)
  vertex_colors <- pfos_identity_colors(igraph::V(net)$name)
  joint_plot(j,paste0(prefix,"_cnet"),paste(title,"term-gene membership network"),function(){oldpar<-par(cex.main=1,font.main=1);on.exit(par(oldpar));plot(net,layout=xy,xlim=c(-1.6,1.6),ylim=c(-1.6,1.6),vertex.size=ifelse(igraph::V(net)$name %in% top$Description,5,2),vertex.label=ifelse(igraph::V(net)$name %in% top$Description,stringr::str_wrap(igraph::V(net)$name,25),NA),vertex.label.cex=7/8,vertex.label.family=pfos_figure_style$family,vertex.label.color=pfos_palette("neutral")[["dark"]],vertex.color=unname(vertex_colors[igraph::V(net)$name]),edge.color=pfos_palette("neutral")[["mid"]],main="Term-gene membership",cex.main=1,font.main=1)},joint_table(j,paste0(prefix,"_network_edges"),edges),height=7)
  chord_terms <- head(member,6L); chord_genes <- head(sort(table(unlist(chord_terms)),decreasing=TRUE),15L)
  chord <- bind_rows(lapply(seq_along(chord_terms),function(i){g<-intersect(chord_terms[[i]],names(chord_genes));data.frame(Term=rep(top$Description[i],length(g)),Gene=g,value=rep(1,length(g)))}))
  if(nrow(chord)>0) joint_plot(j,paste0(prefix,"_chord"),paste(title,"membership chord"),function(){
    oldpar<-par(mar=c(1,1,3,1),family=pfos_figure_style$family,cex=1)
    circlize::circos.clear();on.exit({circlize::circos.clear();par(oldpar)})
    sectors<-unique(c(chord$Term,chord$Gene));labels<-stringr::str_wrap(sectors,26)
    sector_colors<-pfos_identity_colors(sectors)
    label_height_mm<-max(strwidth(labels,units="inches",cex=7/8))*25.4+5
    label_lines<-lengths(strsplit(labels,"\n",fixed=TRUE))
    circlize::circos.par(start.degree=90,gap.after=setNames(2+2*pmax(label_lines-1,0),sectors),track.margin=c(.001,.001),canvas.xlim=c(-1.08,1.08),canvas.ylim=c(-1.08,1.08))
    drawn<-circlize::chordDiagram(chord,annotationTrack="grid",transparency=.4,grid.col=sector_colors,
      preAllocateTracks=list(track.height=circlize::mm_h(label_height_mm)))
    stopifnot(nrow(drawn)==nrow(chord),identical(sort(paste(drawn$rn,drawn$cn)),sort(paste(chord$Term,chord$Gene))),all(drawn$value1==chord$value),all(drawn$value2==chord$value))
    circlize::circos.trackPlotRegion(track.index=1,panel.fun=function(x,y){
      label<-circlize::get.cell.meta.data("sector.index")
      circlize::circos.text(circlize::get.cell.meta.data("xcenter"),circlize::get.cell.meta.data("ylim")[1],stringr::str_wrap(label,26),facing="clockwise",niceFacing=TRUE,adj=c(0,0.5),cex=7/8)
    },bg.border=NA)
    title("Term-gene membership; top six terms / fifteen genes",line=0,cex.main=1,font.main=1)
  },joint_table(j,paste0(prefix,"_chord_edges"),chord),height=7)
  z
}
joint_gsea <- function(j,prefix,rank,sets,title) {
  rank <- rank[is.finite(rank) & !is.na(names(rank)) & nzchar(names(rank))]
  rank <- sort(rank[!duplicated(names(rank))],decreasing=TRUE)
  table <- joint_table(j,paste0(prefix,"_ranking"),data.frame(gene_id_clean=names(rank),rank=as.numeric(rank)))
  if(length(rank)<100L || length(unique(rank))<10L) return(joint_skip(j,paste0(prefix,"_dotplot"),title,"完整可评估排序不足100基因或缺少变化","NOT_RUN",table))
  z<-read_tsv(joint_cached(j,paste0(prefix,"_gsea")),show_col_types=FALSE)
  table<-joint_cached(j,paste0(prefix,"_gsea"))
  top <- head(z,20L)
  if(!nrow(top)) return(joint_skip(j,paste0(prefix,"_dotplot"),title,"无符合集合大小的条目","EMPTY",table))
  p <- ggplot(top,aes(NES,reorder(stringr::str_wrap(gsub("_"," ",pathway),45),NES),size=size,color=-log10(pmax(padj,1e-300)))) + geom_point() +
    scale_color_gradientn(colours=unname(pfos_palette("sequential"))) + scale_x_continuous(limits=unlist(pfos_profile$scales$nes),oob=scales::squish,breaks=c(-3,0,3),labels=c("≤-3","0","≥3")) + labs(title=title,x="Normalized enrichment score",y=NULL,color="-log10(BH)",size="Set size")
  joint_plot(j,paste0(prefix,"_dotplot"),title,p,table,height=7)
  invisible(z)
}
