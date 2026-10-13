#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(GSVA); library(fgsea); library(BiocParallel)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3L) stop("用法: deepening_analysis.R <module_root> <stage_result> <analyse|plot|plot-library>")
root <- normalizePath(args[[1]], mustWork=TRUE)
stage <- normalizePath(args[[2]], mustWork=TRUE)
mode <- match.arg(args[[3]], c("analyse","plot","plot-library"))
cfg <- yaml::read_yaml(file.path(root,"config/project.yml"))
source(file.path(root,"scripts/figure_style.R"))
set.seed(cfg$deepening$seed)

read_any <- function(path) data.table::fread(path, data.table=FALSE, showProgress=FALSE)
write_tsv <- function(x, rel) {
  path <- file.path(stage,rel); dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(as.data.frame(x),path,sep="\t",quote=FALSE,na="")
  paste0("result/",rel)
}
clean_gene <- function(x) sub("\\.[0-9]+$","",as.character(x))
truth <- function(x) as.character(x) %in% c("TRUE","True","true","1")
resolve <- function(rel) normalizePath(file.path(root,rel),mustWork=TRUE)
template_path <- function(key,context) resolve(gsub("\\{context\\}",context,cfg$deepening$posttranscriptional[[key]]))
rna_spec <- function(context) cfg$inputs[[which(vapply(cfg$inputs,function(x) identical(x$kind,"rna")&&identical(x$context,context),logical(1)))[1]]]
spec <- function(kind,context) cfg$inputs[[which(vapply(cfg$inputs,function(x) identical(x$kind,kind)&&identical(x$context,context),logical(1)))[1]]]
sha256 <- function(path) as.character(digest::digest(file=path,algo="sha256",serialize=FALSE))
scales_fixed <- pfos_profile$scales
div_cols <- unname(pfos_palette("diverging")); seq_cols <- unname(pfos_palette("sequential"))

manifest_row <- function(analysis,id,source,title,meaning,distribution,conclusion,kind) data.frame(
  analysis=analysis,figure_id=id,
  png=paste0("result/",analysis,"/plots/",id,".png"),
  pdf=paste0("result/",analysis,"/plots/",id,".pdf"),source_table=source,status="COMPLETED",reason="",
  title=title,biological_meaning=meaning,distribution=distribution,conclusion=conclusion,render_kind=kind,
  asset_class="standalone",figure_group="PFOS_integration_deepening",panel_id="",stringsAsFactors=FALSE)

analyse_programs <- function() {
  gsdb <- readRDS(resolve(cfg$gene_sets)); collections <- cfg$deepening$program_collections
  score_rows <- list(); summary_rows <- list(); gsea_rows <- list(); cohort_rows <- list()
  for (context in cfg$rna_context_order) {
    s <- rna_spec(context); de_path <- resolve(s$path); de <- read_any(de_path)
    expression_path <- file.path(dirname(de_path),"normalized_expression.tsv")
    samples_path <- file.path(dirname(de_path),"samples.tsv")
    expr <- read_any(expression_path); samples <- read_any(samples_path)
    expr$gene_id_clean <- clean_gene(expr[[1]]); de$gene_id_clean <- clean_gene(de$ensembl_id)
    map <- de %>% filter(!is.na(SYMBOL),SYMBOL!="") %>% arrange(desc(baseMean),gene_id_clean) %>% distinct(gene_id_clean,.keep_all=TRUE) %>% select(gene_id_clean,SYMBOL)
    x <- inner_join(expr,map,by="gene_id_clean") %>% arrange(SYMBOL,desc(rowMeans(across(all_of(samples$sample)),na.rm=TRUE))) %>% distinct(SYMBOL,.keep_all=TRUE)
    mat <- as.matrix(x[,samples$sample,drop=FALSE]); storage.mode(mat)<-"double"; rownames(mat)<-x$SYMBOL; mat<-log2(mat+1)
    hallmark <- gsdb$libraries$Mouse_MH
    param <- GSVA::ssgseaParam(mat,hallmark,minSize=10,maxSize=500,normalize=TRUE,checkNA="auto",verbose=FALSE)
    scores <- GSVA::gsva(param,verbose=FALSE,BPPARAM=BiocParallel::SerialParam())
    long <- as.data.frame(scores) %>% tibble::rownames_to_column("term") %>% pivot_longer(-term,names_to="sample",values_to="score") %>% left_join(samples[,c("sample","group")],by="sample") %>% mutate(context=context)
    score_rows[[context]] <- long
    summary_rows[[context]] <- long %>% group_by(context,term,group) %>% summarise(median_score=median(score),.groups="drop") %>% pivot_wider(names_from=group,values_from=median_score) %>% mutate(pfos_minus_wt=PFOS-WT)
    cohort_rows[[context]] <- data.frame(context=context,analysis="ssGSEA_full_cohort",samples=paste(sort(samples$sample),collapse=";"),n=nrow(samples),sample_sha256=sha256(samples_path),expression_sha256=sha256(expression_path),de_sha256=sha256(de_path))
    rank <- de %>% filter(is.finite(stat),!is.na(SYMBOL),SYMBOL!="") %>% arrange(SYMBOL,desc(abs(stat)),gene_id_clean) %>% distinct(SYMBOL,.keep_all=TRUE) %>% select(SYMBOL,stat)
    rv <- rank$stat; names(rv)<-rank$SYMBOL; rv<-sort(rv,decreasing=TRUE)
    for (collection in collections) {
      set.seed(cfg$deepening$seed)
      fg <- fgsea::fgseaMultilevel(gsdb$libraries[[collection]],rv,minSize=10,maxSize=500,eps=0,nproc=1,BPPARAM=BiocParallel::SerialParam())
      fg <- as.data.frame(fg); fg$leadingEdge<-vapply(fg$leadingEdge,paste,collapse=";",character(1)); fg$context<-context; fg$collection<-collection
      fg$padj_context_collection <- p.adjust(fg$pval,"BH")
      gsea_rows[[paste(context,collection)]] <- fg
    }
  }
  scores <- bind_rows(score_rows); summaries<-bind_rows(summary_rows); gsea<-bind_rows(gsea_rows); cohorts<-bind_rows(cohort_rows)
  write_tsv(scores,"02_program_integration/03_context_programs/tables/hallmark_ssgsea_library_scores.tsv")
  write_tsv(summaries,"02_program_integration/03_context_programs/tables/hallmark_ssgsea_context_summary.tsv")
  write_tsv(gsea,"02_program_integration/03_context_programs/tables/full_statistic_fgsea.tsv")
  write_tsv(cohorts,"02_program_integration/03_context_programs/tables/program_input_cohorts.tsv")
  le <- gsea %>% filter(collection=="Mouse_MH",padj_context_collection<.05,leadingEdge!="") %>% select(context,collection,pathway,NES,padj_context_collection,leadingEdge) %>% separate_rows(leadingEdge,sep=";") %>% rename(gene=leadingEdge) %>% mutate(member=1L)
  write_tsv(le,"02_program_integration/03_context_programs/tables/hallmark_leading_edge_membership.tsv")
}

gene_annotation <- function() {
  command <- sprintf("awk -F '\\t' '$3==\"gene\"' %s",shQuote(resolve(cfg$joint_inputs$gtf)))
  gtf <- data.table::fread(cmd=command,sep="\t",header=FALSE,quote="",data.table=FALSE,showProgress=FALSE)
  names(gtf)<-c("chr","source","feature","start","end","score","strand","frame","attribute")
  gtf <- gtf[gtf$feature=="gene",]; gtf$gene_id<-sub('.*gene_id "([^"]+)".*','\\1',gtf$attribute); gtf$gene_id<-clean_gene(gtf$gene_id)
  gtf$chr<-paste0("chr",sub("^chr","",gtf$chr))
  gtf$tss<-ifelse(gtf$strand=="-",gtf$end,gtf$start)
  gtf <- gtf[,c("gene_id","chr","start","end","strand","tss")]
  names(gtf)[names(gtf)=="start"] <- "gene_start"; names(gtf)[names(gtf)=="end"] <- "gene_end"
  gtf
}

analyse_calibration_loci <- function() {
  de_path <- resolve(rna_spec("liver")$path); rna <- read_any(de_path); rna$gene_id<-clean_gene(rna$ensembl_id)
  atac <- read_any(file.path(stage,"03_liver_multiomics/02_pairwise_evidence/05_complete_rank/tables/complete_atac_gene_ranking_source.tsv")); atac$gene_id<-clean_gene(atac$gene_id_clean)
  genes <- gene_annotation(); universe <- inner_join(rna,atac,by="gene_id",suffix=c("_rna","_atac")) %>% left_join(genes,by="gene_id") %>% filter(is.finite(log2FoldChange),is.finite(pvalue_rna),is.finite(Fold),is.finite(pvalue_atac))
  midpoint <- (universe$start+universe$end)/2
  universe$region_category <- ifelse(is.finite(universe$tss)&abs(midpoint-universe$tss)<=3000,"promoter",ifelse(is.finite(universe$gene_start)&midpoint>=universe$gene_start&midpoint<=universe$gene_end,"gene_body","distal"))
  universe$expression_bin <- paste0("Q",dplyr::ntile(log10(pmax(universe$baseMean,0)+1),5))
  universe$peak_count_bin <- ifelse(universe$RegionCount>=3,"3+",as.character(universe$RegionCount))
  universe$stratum <- paste(universe$expression_bin,universe$peak_count_bin,universe$region_category,sep="|")
  universe$rna_state <- ifelse(is.finite(universe$padj)&universe$padj<.05&universe$log2FoldChange>=.58,"Up",ifelse(is.finite(universe$padj)&universe$padj<.05&universe$log2FoldChange<=-.58,"Down","Other"))
  universe$atac_state <- ifelse(is.finite(universe$pvalue_atac)&universe$pvalue_atac<.05&universe$Fold>0,"Open",ifelse(is.finite(universe$pvalue_atac)&universe$pvalue_atac<.05&universe$Fold<0,"Closed","Other"))
  universe$joint_state <- paste(universe$atac_state,universe$rna_state,sep="-")
  keep <- c("gene_id","SYMBOL","baseMean","log2FoldChange","pvalue_rna","padj","Fold","pvalue_atac","FDR","RegionCount","region_category","expression_bin","peak_count_bin","stratum","rna_state","atac_state","joint_state")
  write_tsv(universe[,keep],"03_liver_multiomics/03_background_calibration/tables/rna_atac_complete_tested_universe.tsv")
  strata <- universe %>% count(stratum,expression_bin,peak_count_bin,region_category,name="genes")
  write_tsv(strata,"03_liver_multiomics/03_background_calibration/tables/permutation_strata.tsv")
  observed <- c(`Open-Up`=sum(universe$joint_state=="Open-Up"),`Closed-Down`=sum(universe$joint_state=="Closed-Down")); observed<-c(observed,Concordant=sum(observed))
  B<-cfg$deepening$permutations
  null <- matrix(0L,nrow=B,ncol=3,dimnames=list(NULL,names(observed)))
  set.seed(cfg$deepening$seed)
  # 分层内随机置换等价于固定RNA/ATAC边际的多元超几何抽样；向量化实现保留10,000次真实置换分布。
  for (idx in split(seq_len(nrow(universe)),universe$stratum)) {
    n<-length(idx); u<-sum(universe$rna_state[idx]=="Up"); down<-sum(universe$rna_state[idx]=="Down")
    open<-sum(universe$atac_state[idx]=="Open"); closed<-sum(universe$atac_state[idx]=="Closed"); other<-n-open-closed
    u_open<-rhyper(B,open,n-open,u)
    u_closed<-rhyper(B,closed,other,u-u_open)
    u_other<-u-u_open-u_closed
    rem_open<-open-u_open; rem_closed<-closed-u_closed; rem_other<-other-u_other
    d_open<-rhyper(B,rem_open,rem_closed+rem_other,down)
    d_closed<-rhyper(B,rem_closed,rem_other,down-d_open)
    null[,1]<-null[,1]+u_open; null[,2]<-null[,2]+d_closed
  }
  null[,3]<-null[,1]+null[,2]
  null_long <- as.data.frame(null,check.names=FALSE) %>% mutate(permutation=seq_len(B)) %>% pivot_longer(-permutation,names_to="test",values_to="count")
  tests <- data.frame(test=names(observed),observed=as.integer(observed),empirical_p=vapply(seq_along(observed),function(i)(1+sum(null[,i]>=observed[i]))/(B+1),numeric(1)),B=B,seed=cfg$deepening$seed)
  tests$bh_three_tests<-p.adjust(tests$empirical_p,"BH")
  write_tsv(null_long,"03_liver_multiomics/03_background_calibration/tables/stratified_permutation_null.tsv")
  write_tsv(tests,"03_liver_multiomics/03_background_calibration/tables/stratified_permutation_tests.tsv")
  elig <- data.frame(layer=c("RNA","ATAC","WGBS"),eligible=c(TRUE,TRUE,FALSE),basis=c("full DE tested universe","complete annotated ATAC gene ranking","selected DMR/coverage eligibility only"),inference=c("stratified permutation","stratified permutation","no three-layer P"))
  write_tsv(elig,"03_liver_multiomics/03_background_calibration/tables/layer_inference_eligibility.tsv")

  candidates <- read_any(file.path(stage,"03_liver_multiomics/01_joint_discovery/tables/liver_gene_layer_support.tsv")) %>% filter(truth(three_layer_measurable_support),truth(rna_fdr_supported)) %>% pull(gene_id)
  atac_links <- read_any(file.path(stage,"03_liver_multiomics/02_pairwise_evidence/01_atac_rna/tables/all_atac_region_rna_links.tsv"))
  wgbs <- read_any(resolve(spec("wgbs","liver")$path)); locus_rows<-list()
  for (g in candidates) {
    rr<-rna[!is.na(rna$gene_id)&rna$gene_id==g,]; atac_gene<-clean_gene(atac_links$gene_id_clean); aa<-atac_links[!is.na(atac_gene)&atac_gene==g,]; wgbs_gene<-clean_gene(wgbs$gene_id_clean); ww<-wgbs[!is.na(wgbs_gene)&wgbs_gene==g,]; gg<-genes[!is.na(genes$gene_id)&genes$gene_id==g,]
    symbol<-if(nrow(rr)) as.character(rr$SYMBOL[1]) else g
    z<-bind_rows(
      data.frame(gene_id=g,symbol=symbol,assay="RNA",chr=gg$chr[1],start=gg$tss[1],end=gg$tss[1],effect=rr$log2FoldChange[1],unit="RNA log2FC",p=rr$pvalue[1],fdr=rr$padj[1],evidence_scope="gene-level full-cohort DE"),
      aa %>% transmute(gene_id=g,symbol=symbol,assay="ATAC",chr=seqnames,start=start,end=end,effect=Fold,unit="ATAC log2FC",p=atac_pvalue,fdr=FDR,evidence_scope="selected differential peak; nearest-gene annotation"),
      ww %>% transmute(gene_id=g,symbol=symbol,assay="WGBS",chr=paste0("chr",sub("^chr","",seqnames)),start=start,end=end,effect=diffMethy,unit="methylation fraction delta",p=NA_real_,fdr=NA_real_,evidence_scope="selected DMR; no region P exported here"))
    z$direction<-ifelse(z$effect>0,"positive",ifelse(z$effect<0,"negative","zero")); locus_rows[[g]]<-z
  }
  write_tsv(bind_rows(locus_rows),"03_liver_multiomics/04_shared_loci/tables/five_candidate_locus_evidence.tsv")
}

method_frame <- function(context,method,df,gene_col,support_col,bh_support_col=NA_character_,p_col=NA_character_,fdr_col=NA_character_,effect_col=NA_character_,scope,eligible=TRUE) {
  out<-data.frame(context=context,method=method,gene_id=clean_gene(df[[gene_col]]),tested_universe=eligible,supported=if(eligible) truth(df[[support_col]]) else TRUE,
                  evidence_state=if(eligible) "tested" else "reported_significant_only",scope=scope,stringsAsFactors=FALSE)
  out$candidate_bh<-if(eligible&&!is.na(bh_support_col)&&bh_support_col%in%names(df)) truth(df[[bh_support_col]]) else NA
  out$p<-if(!is.na(p_col)&&p_col%in%names(df)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_
  out$fdr<-if(!is.na(fdr_col)&&fdr_col%in%names(df)) suppressWarnings(as.numeric(df[[fdr_col]])) else NA_real_
  out$effect<-if(!is.na(effect_col)&&effect_col%in%names(df)) suppressWarnings(as.numeric(df[[effect_col]])) else NA_real_
  out[!is.na(out$gene_id)&out$gene_id!="",]
}

analyse_posttx <- function() {
  rows<-list(); cohorts<-list(); method_cohorts<-list(); method_summary<-list(); ora<-list(); gsdb<-readRDS(resolve(cfg$gene_sets)); hallmark<-gsdb$libraries$Mouse_MH
  for(context in cfg$rna_context_order) {
    selected_samples<-read_any(template_path("selected_samples",context)); full_path<-file.path(dirname(resolve(rna_spec(context)$path)),"samples.tsv"); full_samples<-read_any(full_path)
    same<-identical(sort(as.character(selected_samples$sample)),sort(as.character(full_samples$sample)))
    cohorts[[context]]<-data.frame(context=context,selected_n=nrow(selected_samples),full_n=nrow(full_samples),selected_samples=paste(sort(selected_samples$sample),collapse=";"),full_samples=paste(sort(full_samples$sample),collapse=";"),same_sample_signature=same,overlap_allowed=same,selected_sha256=sha256(template_path("selected_samples",context)),full_sha256=sha256(full_path))
    selected_pairs<-sort(paste(selected_samples$sample,selected_samples$group,sep="="))
    for(method in c("DAS","DEU","APA","DER","IS")) {
      key<-c(DAS="das_metadata",DEU="deu_metadata",APA="apa_metadata",DER="der_metadata",IS="is_metadata")[[method]]
      metadata_path<-template_path(key,context)
      if(method=="DAS") {
        lines<-readLines(metadata_path,warn=FALSE); kv<-strsplit(lines,"=",fixed=TRUE); names(kv)<-vapply(kv,`[`,character(1),1)
        wt<-strsplit(kv$group1_samples[2],",",fixed=TRUE)[[1]]; pfos<-strsplit(kv$group2_samples[2],",",fixed=TRUE)[[1]]
        metadata_pairs<-sort(c(paste(wt,kv$group1[2],sep="="),paste(pfos,kv$group2[2],sep="=")))
      } else {
        metadata<-read_any(metadata_path); id_col<-intersect(c("sample","sampleID"),names(metadata))[1]; group_col<-intersect(c("group","condition"),names(metadata))[1]
        if(length(id_col)!=1L||length(group_col)!=1L)stop("方法队列metadata缺少sample/group列: ",metadata_path)
        metadata_pairs<-sort(paste(metadata[[id_col]],metadata[[group_col]],sep="="))
      }
      same_method<-identical(selected_pairs,metadata_pairs)
      method_cohorts[[paste(context,method)]]<-data.frame(context=context,method=method,metadata_path=gsub("\\{context\\}",context,cfg$deepening$posttranscriptional[[key]]),metadata_sha256=sha256(metadata_path),metadata_n=length(metadata_pairs),selected_n=nrow(selected_samples),metadata_sample_group=paste(metadata_pairs,collapse=";"),selected_sample_group=paste(selected_pairs,collapse=";"),same_selected_signature=same_method,status=ifelse(same_method,"PASS","FAIL"))
    }
    de<-read_any(template_path("selected_de",context)); de$selected<-as.character(de$direction)!="ns"; de$bh_selected<-de$selected&is.finite(de$padj)&de$padj<.05
    das<-read_any(template_path("das_events",context)) %>% group_by(gene_id_clean) %>% summarise(supported=any(truth(significant_as)),candidate_bh=any(truth(significant_as)&is.finite(event_fdr)&event_fdr<.05),p=min(event_pvalue,na.rm=TRUE),fdr=min(event_fdr,na.rm=TRUE),effect=max(abs(delta_psi),na.rm=TRUE),.groups="drop")
    deu<-read_any(template_path("deu_gene_tests",context)); deu$bh_selected<-truth(deu$significant_deu)&is.finite(deu$FDR)&deu$FDR<.05
    apa<-read_any(template_path("apa_events",context)) %>% group_by(gene_id_clean) %>% summarise(supported=any(truth(significant_apa)),candidate_bh=any(truth(significant_apa)&is.finite(apa_padj)&apa_padj<.05),p=min(apa_pvalue,na.rm=TRUE),fdr=min(apa_padj,na.rm=TRUE),effect=max(abs(delta_PAU),na.rm=TRUE),.groups="drop")
    iso<-read_any(template_path("isoform_selected",context)); der<-read_any(template_path("der_regions",context)) %>% group_by(gene_id_clean) %>% summarise(supported=any(truth(significant_der)),candidate_bh=any(truth(significant_der)&is.finite(adj.P.Val)&adj.P.Val<.05),p=min(P.Value,na.rm=TRUE),fdr=min(adj.P.Val,na.rm=TRUE),effect=max(abs(der_effect),na.rm=TRUE),.groups="drop")
    frames<-list(
      method_frame(context,"DE_selected6",de,"ensembl_id","selected","bh_selected","pvalue","padj","log2FoldChange","upstream selected-cohort DE candidate rule; gene P/BH retained",TRUE),
      method_frame(context,"DAS",das,"gene_id_clean","supported","candidate_bh","p","fdr","effect","exploratory nominal event P < 0.05 + |delta PSI| >= 0.1; gene row stores minimum event P/FDR as a summary, not gene P/FDR",TRUE),
      method_frame(context,"DEU",deu,"gene_id_clean","significant_deu","bh_selected","P.Value","FDR",NA,"exploratory nominal gene-level Simes P < 0.05; gene FDR retained",TRUE),
      method_frame(context,"APA",apa,"gene_id_clean","supported","candidate_bh","p","fdr","effect","exploratory nominal event P < 0.05 + |delta PAU| >= 0.1; gene row stores minimum event P/FDR as a summary, not gene P/FDR",TRUE),
      method_frame(context,"DER",der,"gene_id_clean","supported","candidate_bh","p","fdr","effect","exploratory nominal region P < 0.05; gene row stores minimum region P/FDR as a summary, not gene P/FDR",TRUE),
      method_frame(context,"IS",iso,"gene_id_clean","significant_switch_gene",NA,"gene_switch_q_value","gene_switch_q_value","max_abs_dIF","significant-only isoform-switch export; absent genes unreported and tested universe unknown",FALSE))
    z<-bind_rows(frames); rows[[context]]<-z
    method_summary[[context]]<-z %>% group_by(context,method,evidence_state,scope) %>% summarise(
      tested_genes=ifelse(first(tested_universe),n_distinct(gene_id),NA_integer_),candidate_genes=n_distinct(gene_id[supported]),
      candidate_bh_genes=ifelse(first(tested_universe),n_distinct(gene_id[candidate_bh%in%TRUE]),NA_integer_),
      candidate_min_reported_fdr=if(any(supported&is.finite(fdr)))min(fdr[supported&is.finite(fdr)])else NA_real_,
      candidate_max_reported_fdr=if(any(supported&is.finite(fdr)))max(fdr[supported&is.finite(fdr)])else NA_real_,
      reported_fdr_scope=ifelse(first(method)%in%c("DAS","APA","DER"),"event/region FDR summary; not gene FDR",ifelse(first(method)=="IS","unknown tested universe","gene-level FDR")),.groups="drop")
    symbol_map<-de %>% filter(!is.na(SYMBOL),SYMBOL!="") %>% distinct(ensembl_id,.keep_all=TRUE) %>% transmute(gene_id=clean_gene(ensembl_id),SYMBOL)
    for(m in c("DE_selected6","DAS","DEU","APA","DER")) {
      zz<-z[z$method==m,]; universe<-unique(inner_join(zz,symbol_map,by="gene_id")$SYMBOL); selected<-unique(inner_join(zz[zz$supported,],symbol_map,by="gene_id")$SYMBOL)
      for(term in names(hallmark)) {
        set<-intersect(hallmark[[term]],universe); a<-length(intersect(selected,set)); b<-length(setdiff(selected,set)); c0<-length(setdiff(set,selected)); d<-length(setdiff(universe,union(selected,set)))
        pv<-if(length(universe)>0&&length(selected)>0&&length(set)>=10) fisher.test(matrix(c(a,b,c0,d),2),alternative="greater")$p.value else NA_real_
        ora[[length(ora)+1L]]<-data.frame(context=context,method=m,term=term,overlap=a,selected_n=length(selected),set_n=length(set),universe_n=length(universe),gene_ratio=if(length(selected))a/length(selected) else NA_real_,p=pv,scope=paste0("exploratory candidate-defined ORA; upstream rule: ",unique(zz$scope)[1],"; event/region methods have unequal detection opportunities"))
      }
    }
    ora[[length(ora)+1L]]<-data.frame(context=context,method="IS",term="NOT_RUN",overlap=NA,selected_n=nrow(iso),set_n=NA,universe_n=NA,gene_ratio=NA,p=NA,scope="NOT_RUN: significant-only export has no tested universe")
  }
  evidence<-bind_rows(rows); summaries<-bind_rows(method_summary); enrichment<-bind_rows(ora) %>% group_by(context,method) %>% mutate(bh=p.adjust(p,"BH")) %>% ungroup()
  write_tsv(evidence,"06_synthesis/02_posttranscriptional_patterns/tables/rebuilt_method_gene_evidence.tsv")
  write_tsv(summaries,"06_synthesis/02_posttranscriptional_patterns/tables/method_availability_summary.tsv")
  write_tsv(bind_rows(cohorts),"06_synthesis/02_posttranscriptional_patterns/tables/cohort_signature_audit.tsv")
  method_cohort_table<-bind_rows(method_cohorts)
  if(nrow(method_cohort_table)!=length(cfg$rna_context_order)*5L||any(method_cohort_table$status!="PASS"))stop("转录后方法队列与selected cohort签名不一致")
  write_tsv(method_cohort_table,"06_synthesis/02_posttranscriptional_patterns/tables/method_cohort_signature_audit.tsv")
  write_tsv(enrichment,"06_synthesis/02_posttranscriptional_patterns/tables/method_hallmark_ora.tsv")
  old <- bind_rows(lapply(cfg$rna_context_order,function(context){x<-read_any(template_path("integrated_evidence",context));data.frame(context=context,legacy_adapter_selected_genes=sum(truth(x$DER)),legacy_adapter_rows=nrow(x))}))
  new <- evidence %>% filter(method=="DER") %>% group_by(context) %>% summarise(reconstructed_nominal_candidate_genes=sum(supported),reconstructed_candidate_bh_genes=sum(candidate_bh%in%TRUE),tested_genes=n(),.groups="drop")
  write_tsv(left_join(old,new,by="context"),"06_synthesis/02_posttranscriptional_patterns/tables/der_adapter_reconstruction_audit.tsv")
}

build_manifest <- function() {
  rows<-list(
    manifest_row("02_program_integration/03_context_programs","hallmark_ssgsea_context_matrix","result/02_program_integration/03_context_programs/tables/hallmark_ssgsea_context_summary.tsv","七情境 Hallmark 程序差异","比较每个情境PFOS与WT文库的Hallmark程序得分中位数差。","色彩为ssGSEA分数组内中位数差，使用本数据对称范围；它不是NES。","这是逐文库描述，生物单位与配对关系unknown，不作跨情境重复测量推断。","deep_program_context"),
    manifest_row("02_program_integration/03_context_programs","hallmark_ssgsea_library_matrix","result/02_program_integration/03_context_programs/tables/hallmark_ssgsea_library_scores.tsv","Hallmark 逐文库程序谱","展示变异最大的Hallmark在全部真实文库中的相对活性。","每个Hallmark跨文库行标准化并固定在-2.5至2.5；原始ssGSEA分数保留在源表。","单文库分数只描述样本结构，不能代替差异检验。","deep_program_library"),
    manifest_row("02_program_integration/03_context_programs","full_statistic_fgsea_matrix","result/02_program_integration/03_context_programs/tables/full_statistic_fgsea.tsv","完整DE统计量的程序富集矩阵","以七情境完整DE统计量队列检验Hallmark、KEGG及Reactome。","色彩为NES并固定-3至3，圆点大小为同情境同集合BH支持，端点标示饱和。","Mouse_MH为主集合；KEGG与Reactome为补充。不同情境独立，不构成纵向轨迹。","deep_gsea"),
    manifest_row("02_program_integration/03_context_programs","hallmark_leading_edge_matrix","result/02_program_integration/03_context_programs/tables/hallmark_leading_edge_membership.tsv","Hallmark leading-edge成员矩阵","显示多情境富集由哪些leading-edge基因承载。","二元色仅表示成员关系；不编码表达幅度。","成员重叠提示可追踪程序支点，不证明共同调控因果。","deep_leading_edge"),
    manifest_row("03_liver_multiomics/03_background_calibration","rna_atac_full_universe","result/03_liver_multiomics/03_background_calibration/tables/rna_atac_complete_tested_universe.tsv","RNA–ATAC完整检验背景","在完整可匹配RNA和ATAC基因背景中显示效应与预设方向象限。","所有匹配基因均保留；颜色仅标识Open-Up、Closed-Down及其他状态。","象限是未配对证据共现，不是样本级相关或调控因果。","deep_full_universe"),
    manifest_row("03_liver_multiomics/03_background_calibration","stratified_permutation_null","result/03_liver_multiomics/03_background_calibration/tables/stratified_permutation_null.tsv","RNA–ATAC分层置换校准","检验Open-Up与Closed-Down共现是否超过表达量、链接数及区域类别匹配的随机背景。","密度为10000次固定分层置换，竖线为观察计数；BH只覆盖三个预设检验。","WGBS只有selected DMR/coverage资格，不进入三层P值。","deep_permutation"),
    manifest_row("03_liver_multiomics/04_shared_loci","five_shared_locus_tracks","result/03_liver_multiomics/04_shared_loci/tables/five_candidate_locus_evidence.tsv","五个三层RNA-BH候选共享位点","把Banp、Tob1、Pnrc1、Fbxo31和Dleu2的RNA、ATAC及WGBS证据放到各自共同基因组坐标。","线段位置是GRCm39坐标；色彩只编码PFOS-WT方向，效应单位保留在源表。","多层同位显示注释收敛，不能证明样本配对、共定位因果或调控方向。","deep_locus"),
    manifest_row("06_synthesis/02_posttranscriptional_patterns","posttx_method_availability","result/06_synthesis/02_posttranscriptional_patterns/tables/method_availability_summary.tsv","转录后方法可用性与候选量","审计DE、DAS、DEU、IS、APA及DER在七情境的测试域、上游候选量及候选中的BH支持量。","圆面积表示上游候选基因数；文字依次给出测试域和候选BH数。IS只有显著表，未报告基因保持unknown。","DAS、APA、DER候选遵循上游nominal P加效应或nominal P规则；事件/区域P与FDR摘要均不称为gene P/FDR。旧DER适配器全FALSE已失效，新入口从现役区域表重建。","deep_posttx_availability"),
    manifest_row("06_synthesis/02_posttranscriptional_patterns","posttx_evidence_matrix","result/06_synthesis/02_posttranscriptional_patterns/tables/rebuilt_method_gene_evidence.tsv","转录后多方法证据矩阵","展示跨情境多方法重复出现的上游候选基因。","点为方法内候选；实心表示候选中对应gene/event/region BH小于0.05，空心表示探索候选，叉号表示显著表没有测试域。完整证据表保留全部测试域。","selected6与full10按样本签名分队列；签名不同处不计算同队列交集。事件/区域min P与min FDR只是基因映射摘要，不是gene P/FDR。","deep_posttx_matrix")
  )
  for(method in c("DE_selected6","DAS","DEU","APA","DER")) {
    id<-paste0("posttx_hallmark_ora_",tolower(method))
    rows[[length(rows)+1L]]<-manifest_row("06_synthesis/02_posttranscriptional_patterns",id,"result/06_synthesis/02_posttranscriptional_patterns/tables/method_hallmark_ora.tsv",paste(method,"Hallmark程序卡"),paste0("在",method,"实际测试域内对上游候选基因作探索性Hallmark富集。"),"横轴GeneRatio、圆面积线性表示Count、颜色为-log10原始P并固定0至10；形状标识该情境×方法内ORA BH。","上游候选规则多为nominal P并保留效应门槛；事件和区域方法有不同检测机会。ORA BH不把上游候选改写为gene-level显著。IS因缺少测试域明确NOT_RUN。","deep_posttx_ora")
  }
  write_tsv(bind_rows(rows),"00_audit/tables/deepening_manifest.tsv")
}

save_plot <- function(p,analysis,id,height=145) {
  if(!is.finite(height)||height>240)stop("独立图高度必须不超过240 mm: ",id)
  pfos_save(p,file.path(stage,analysis,"plots",paste0(id,".pdf")),height_mm=height)
}
sat_labels <- function(lo,hi) c(paste0("≤",lo),"0",paste0("≥",hi))

plot_program_context <- function() {
  d<-read_any(file.path(stage,"02_program_integration/03_context_programs/tables/hallmark_ssgsea_context_summary.tsv")); lim<-max(abs(d$pfos_minus_wt),na.rm=TRUE)
  d$term<-factor(gsub("HALLMARK_|_"," ",d$term),levels=rev(gsub("HALLMARK_|_"," ",unique(d$term))))
  p<-ggplot(d,aes(context,term,fill=pfos_minus_wt))+geom_tile()+scale_fill_gradientn(colors=div_cols,limits=c(-lim,lim),oob=scales::squish,name="PFOS - WT\nmedian ssGSEA")+labs(x=NULL,y=NULL)+pfos_theme()+theme(axis.text.y=element_text(size=7))
  save_plot(p,"02_program_integration/03_context_programs","hallmark_ssgsea_context_matrix",240)
}
plot_program_library <- function() {
  d<-read_any(file.path(stage,"02_program_integration/03_context_programs/tables/hallmark_ssgsea_library_scores.tsv")); v<-d %>% group_by(term) %>% summarise(v=var(score),.groups="drop") %>% slice_max(v,n=20); d<-semi_join(d,v,by="term") %>% group_by(term) %>% mutate(row_z=as.numeric(scale(score))) %>% ungroup()
  d$term<-gsub("HALLMARK_|_"," ",d$term); p<-ggplot(d,aes(sample,term,fill=row_z))+geom_tile()+facet_grid(.~context,scales="free_x",space="free_x",labeller=labeller(context=c(blastocyst="blasto-\ncyst")))+scale_fill_gradientn(colors=div_cols,limits=unlist(scales_fixed$row_zscore),oob=scales::squish,breaks=c(-2.5,0,2.5),labels=sat_labels(-2.5,2.5),name="row z-score")+labs(x=NULL,y=NULL)+pfos_theme()+theme(axis.text.x=element_text(angle=70,hjust=1),axis.text.y=element_text(size=7),strip.text.x=element_text(lineheight=.9))
  save_plot(p,"02_program_integration/03_context_programs","hallmark_ssgsea_library_matrix",190)
}
plot_gsea <- function() {
  d<-read_any(file.path(stage,"02_program_integration/03_context_programs/tables/full_statistic_fgsea.tsv")); top<-d %>% group_by(collection,pathway) %>% summarise(q=min(padj_context_collection,na.rm=TRUE),.groups="drop") %>% group_by(collection) %>% slice_min(q,n=10,with_ties=FALSE) %>% ungroup(); d<-semi_join(d,top,by=c("collection","pathway")); d$support<-pmin(-log10(pmax(d$padj_context_collection,1e-10)),10)
  d$label<-gsub("^(HALLMARK_|KEGG_|REACTOME_)","",d$pathway); d$label<-stringr::str_wrap(gsub("_"," ",d$label),32)
  p<-ggplot(d,aes(context,label))+geom_tile(aes(fill=NES))+geom_point(data=d[is.finite(d$support),],aes(size=support),shape=21,fill=NA,color=pfos_palette("neutral")[["dark"]],stroke=.35)+facet_grid(collection~.,scales="free_y",space="free_y")+scale_fill_gradientn(colors=div_cols,limits=unlist(scales_fixed$nes),oob=scales::squish,breaks=c(-3,0,3),labels=sat_labels(-3,3),name="NES")+scale_size_continuous(range=c(0,4),name="-log10 BH\n(capped 10)")+labs(x=NULL,y=NULL)+pfos_theme()+theme(axis.text.y=element_text(size=7))
  save_plot(p,"02_program_integration/03_context_programs","full_statistic_fgsea_matrix",240)
}
plot_leading <- function() {
  d<-read_any(file.path(stage,"02_program_integration/03_context_programs/tables/hallmark_leading_edge_membership.tsv")); if(!nrow(d)) { d<-data.frame(context="none",pathway="no BH-supported Hallmark",gene="none",member=0) }
  terms<-d %>% count(pathway,sort=TRUE) %>% slice_head(n=10); genes<-semi_join(d,terms,by="pathway") %>% count(gene,sort=TRUE) %>% slice_head(n=35); d<-semi_join(semi_join(d,terms,by="pathway"),genes,by="gene")
  grid<-expand_grid(context=cfg$rna_context_order,pathway=unique(d$pathway),gene=unique(d$gene)) %>% left_join(d[,c("context","pathway","gene","member")],by=c("context","pathway","gene")) %>% mutate(member=factor(replace_na(member,0L)),pathway=stringr::str_wrap(gsub("HALLMARK_|_"," ",pathway),24))
  p<-ggplot(grid,aes(context,gene,fill=member))+geom_tile()+facet_wrap(~pathway,ncol=2)+scale_fill_manual(values=c("0"=pfos_palette("binary")[["zero"]],"1"=pfos_palette("binary")[["one"]]),name="Leading edge")+labs(x=NULL,y=NULL)+pfos_theme()+theme(axis.text.x=element_text(angle=45,hjust=1,size=7),axis.text.y=element_text(size=7))
  save_plot(p,"02_program_integration/03_context_programs","hallmark_leading_edge_matrix",240)
}
plot_calibration <- function() {
  d<-read_any(file.path(stage,"03_liver_multiomics/03_background_calibration/tables/rna_atac_complete_tested_universe.tsv")); cols<-c("Open-Up"=pfos_palette("de")[["up"]],"Closed-Down"=pfos_palette("de")[["down"]],"Open-Down"=pfos_palette("neutral")[["mid"]],"Closed-Up"=pfos_palette("neutral")[["mid"]],"Other-Other"=pfos_palette("neutral")[["light"]])
  d$display<-ifelse(d$joint_state%in%c("Open-Up","Closed-Down"),d$joint_state,"Other")
  p<-ggplot(d,aes(log2FoldChange,Fold,color=display))+geom_point(size=.7,alpha=.55)+geom_vline(xintercept=0,color=pfos_palette("neutral")[["mid"]],linewidth=.3)+geom_hline(yintercept=0,color=pfos_palette("neutral")[["mid"]],linewidth=.3)+scale_color_manual(values=c("Open-Up"=pfos_palette("de")[["up"]],"Closed-Down"=pfos_palette("de")[["down"]],"Other"=pfos_palette("neutral")[["mid"]]))+labs(x="RNA log2FC (PFOS - WT)",y="ATAC log2FC (PFOS - WT)",color=NULL)+pfos_theme()
  save_plot(p,"03_liver_multiomics/03_background_calibration","rna_atac_full_universe",145)
  n<-read_any(file.path(stage,"03_liver_multiomics/03_background_calibration/tables/stratified_permutation_null.tsv")); t<-read_any(file.path(stage,"03_liver_multiomics/03_background_calibration/tables/stratified_permutation_tests.tsv"))
  p2<-ggplot(n,aes(count,fill=test,color=test))+geom_density(alpha=.25,linewidth=.5)+geom_vline(data=t,aes(xintercept=observed,color=test),linewidth=.7)+facet_wrap(~test,scales="free",ncol=1)+scale_color_manual(values=c("Open-Up"=pfos_palette("de")[["up"]],"Closed-Down"=pfos_palette("de")[["down"]],"Concordant"=pfos_palette("accent")[["purple"]]),aesthetics=c("color","fill"))+labs(x="Stratified-permutation count",y="Density",color=NULL,fill=NULL)+pfos_theme()+theme(legend.position="none")
  save_plot(p2,"03_liver_multiomics/03_background_calibration","stratified_permutation_null",180)
}
plot_loci <- function() {
  d<-read_any(file.path(stage,"03_liver_multiomics/04_shared_loci/tables/five_candidate_locus_evidence.tsv")); d$mid<-(d$start+d$end)/2; d$assay<-factor(d$assay,levels=c("RNA","ATAC","WGBS")); d$label<-paste0(d$symbol," (",d$chr,")")
  p<-ggplot(d,aes(y=assay,color=direction))+geom_segment(aes(x=start,xend=end,yend=assay),linewidth=2.4,lineend="round")+geom_point(data=d[d$start==d$end,],aes(x=mid),size=2.6)+facet_wrap(~label,scales="free_x",ncol=1)+scale_color_manual(values=c("negative"=pfos_palette("de")[["down"]],"positive"=pfos_palette("de")[["up"]],"zero"=pfos_palette("de")[["ns"]]))+scale_x_continuous(labels=scales::label_number(scale_cut=scales::cut_short_scale()))+labs(x="GRCm39 coordinate",y=NULL,color="PFOS - WT direction")+pfos_theme()
  save_plot(p,"03_liver_multiomics/04_shared_loci","five_shared_locus_tracks",210)
}
plot_posttx <- function() {
  s<-read_any(file.path(stage,"06_synthesis/02_posttranscriptional_patterns/tables/method_availability_summary.tsv")); s$tested_label<-ifelse(is.na(s$tested_genes),"NA",paste0(s$tested_genes," / ",s$candidate_bh_genes)); s$method_color<-sub("_selected6$","",s$method)
  p<-ggplot(s,aes(context,method,size=candidate_genes,color=method_color,shape=evidence_state))+geom_point()+geom_text(aes(label=tested_label),nudge_y=.28,size=7/ggplot2::.pt,color=unname(pfos_palette("neutral")[["dark"]]))+scale_color_manual(values=pfos_palette("rna_method"),name="Method")+scale_shape_manual(values=c("tested"=16,"reported_significant_only"=17),labels=c("tested"="tested","reported_significant_only"="significant-only export"),name="Availability\n(label = tested / candidate BH;\nNA = test universe unavailable)")+scale_size_area(max_size=7,name="Candidate genes")+labs(x=NULL,y=NULL)+pfos_theme()
  save_plot(p,"06_synthesis/02_posttranscriptional_patterns","posttx_method_availability",145)
  d<-read_any(file.path(stage,"06_synthesis/02_posttranscriptional_patterns/tables/rebuilt_method_gene_evidence.tsv")); top<-d[d$supported,] %>% count(gene_id,sort=TRUE) %>% slice_head(n=50); z<-semi_join(d[d$supported,],top,by="gene_id") %>% group_by(context,method,gene_id) %>% summarise(candidate_bh=if(all(is.na(candidate_bh)))NA else any(candidate_bh%in%TRUE),.groups="drop")
  symbols<-bind_rows(lapply(cfg$rna_context_order,function(context){x<-read_any(template_path("selected_de",context));data.frame(gene_id=clean_gene(x$ensembl_id),symbol=x$SYMBOL)})) %>% filter(!is.na(symbol),symbol!="") %>% distinct(gene_id,.keep_all=TRUE)
  z<-left_join(z,symbols,by="gene_id") %>% mutate(gene_label=ifelse(is.na(symbol)|symbol=="",gene_id,symbol),method_color=sub("_selected6$","",method),bh_state=ifelse(is.na(candidate_bh),"tested universe unknown",ifelse(candidate_bh,"candidate + feature BH < 0.05","exploratory candidate only")))
  p2<-ggplot(z,aes(interaction(context,method,sep=" | "),gene_label,color=method_color,shape=bh_state))+geom_point(size=2,stroke=.45)+scale_color_manual(values=pfos_palette("rna_method"),name="Method")+scale_shape_manual(values=c("candidate + feature BH < 0.05"=16,"exploratory candidate only"=1,"tested universe unknown"=4),name="Candidate status")+labs(x=NULL,y=NULL)+pfos_theme()+theme(axis.text.x=element_text(angle=70,hjust=1,size=7),axis.text.y=element_text(size=7))
  save_plot(p2,"06_synthesis/02_posttranscriptional_patterns","posttx_evidence_matrix",200)
  o<-read_any(file.path(stage,"06_synthesis/02_posttranscriptional_patterns/tables/method_hallmark_ora.tsv")); o<-o[is.finite(o$p)&o$overlap>0,]; top<-o %>% group_by(context,method) %>% slice_min(p,n=3,with_ties=FALSE) %>% ungroup(); top$mlogp<-pmin(-log10(top$p),10); top$bh_support<-ifelse(is.finite(top$bh)&top$bh<.05,"ORA BH < 0.05","ORA BH ≥ 0.05"); top$label<-paste(top$context,stringr::str_wrap(gsub("HALLMARK_|_"," ",top$term),28),sep=" | ")
  for(method in c("DE_selected6","DAS","DEU","APA","DER")) {
    view<-top[top$method==method,]
    p3<-ggplot(view,aes(gene_ratio,reorder(label,gene_ratio),size=overlap,color=mlogp,shape=bh_support))+geom_point()+scale_color_gradientn(colors=seq_cols,limits=unlist(scales_fixed$ora_minus_log10_p),oob=scales::squish,breaks=c(0,5,10),labels=c("0","5","≥10"),name="-log10 raw P")+scale_size_area(max_size=6,name="Count")+labs(x="GeneRatio",y=NULL,shape=NULL)+pfos_theme()+theme(axis.text.y=element_text(size=7))
    save_plot(p3,"06_synthesis/02_posttranscriptional_patterns",paste0("posttx_hallmark_ora_",tolower(method)),220)
  }
}

if (mode=="analyse") {
  analyse_programs(); analyse_calibration_loci(); analyse_posttx(); build_manifest()
} else if (mode=="plot") {
  plot_program_context(); plot_program_library(); plot_gsea(); plot_leading(); plot_calibration(); plot_loci(); plot_posttx()
} else {
  plot_program_library()
}
