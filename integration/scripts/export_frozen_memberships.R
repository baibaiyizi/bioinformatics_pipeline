#!/usr/bin/env Rscript
# 只转换冻结成员与正式注释 ID；不进行富集检验。
a <- commandArgs(TRUE); stopifnot(length(a)==3)
bundle <- readRDS(a[1]); out <- a[2]; mapout <- a[3]
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Mm.eg.db);library(GO.db);library(data.table)})
libs <- list(KEGG=bundle$libraries$Mouse_KEGG, GO_BP=bundle$libraries$Mouse_M5[grepl('^GOBP_',names(bundle$libraries$Mouse_M5))])
values <- unique(unlist(libs,use.names=FALSE))
if(all(grepl('^ENSMUSG',values))) {
 lookup <- setNames(lapply(values,function(x) sub('\\..*$','',x)),values)
} else {
 kt <- if(all(grepl('^[0-9]+$',values))) 'ENTREZID' else 'SYMBOL'
 z <- AnnotationDbi::select(org.Mm.eg.db,keys=intersect(values,keys(org.Mm.eg.db,keytype=kt)),keytype=kt,columns='ENSEMBL')
 lookup <- split(z$ENSEMBL,z[[kt]])
}
go <- AnnotationDbi::select(GO.db,keys=keys(GO.db,keytype='GOID'),keytype='GOID',columns=c('TERM','ONTOLOGY'))
normalize <- function(x) gsub('[^A-Z0-9]+','_',toupper(x))
go$key <- normalize(go$TERM)
go <- go[go$ONTOLOGY=='BP',]
rows <- list(); maps <- list()
for(collection in names(libs)) for(term in names(libs[[collection]])) {
 if(collection=='KEGG') {
   id<-sub('_.*$','',term); description<-sub('^[^_]+_','',term); status<-'MATCHED_FROZEN_KEGG_ID'
 } else {
   hit <- unique(go$GOID[go$key==normalize(sub('^GOBP_','',term))])
   id<-if(length(hit)==1) hit else NA_character_
   description<-gsub('_',' ',sub('^GOBP_','',term))
   status<-if(length(hit)==1) 'MATCHED_GO_DB_EXACT_TERM' else if(length(hit)>1) 'AMBIGUOUS' else 'UNMAPPED'
 }
 members<-unique(na.omit(unlist(lookup[libs[[collection]][[term]]],use.names=FALSE)))
 rows[[length(rows)+1]]<-data.table(collection=collection,term_id=id,frozen_term=term,gene_id=members)
 maps[[length(maps)+1]]<-data.table(collection=collection,term_id=id,frozen_term=term,description=description,status=status)
}
fwrite(rbindlist(rows),out,sep='\t',na='NA');fwrite(rbindlist(maps),mapout,sep='\t',na='NA')
