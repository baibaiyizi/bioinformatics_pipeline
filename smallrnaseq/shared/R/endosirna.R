# endo 是共同计数构建的一部分。参考命中只作为证据；冲突类别保留 SPORTS 分配。
ENDO_ENABLED <- identical(Sys.getenv("RUN_ENDOSIRNA"), "1")
reverse_complement <- function(x) {
  vapply(strsplit(chartr("ACGTN", "TGCAN", toupper(x)), "", fixed = TRUE), function(z) paste(rev(z), collapse = ""), character(1), USE.NAMES = FALSE)
}
endo_reference_files <- function() {
  specs <- strsplit(Sys.getenv("ENDOSIRNA_REFS"), ";", fixed=TRUE)[[1]]
  if (!length(specs) || any(!nzchar(specs))) stop("Missing explicit ENDOSIRNA_REFS configuration.")
  tags <- vapply(specs,function(x)strsplit(x,"=",fixed=TRUE)[[1]][[1]],character(1))
  prefixes <- sub("^[^=]+=", "", specs)
  if (anyDuplicated(tags)) stop("Duplicate endo reference labels.")
  refs <- vapply(prefixes,function(prefix) {
    files <- paste0(prefix,c(".fasta",".fa"))
    found <- files[file.exists(files)]
    if (!length(found)) stop("Missing endo FASTA: ",prefix)
    found[[1]]
  },character(1))
  stats::setNames(refs,tags)
}
# 精确匹配不因同一序列的多条参考/两个方向重复计数。
load_endo_evidence <- function(sample_info) {
  refs <- endo_reference_files()
  if (any(!file.exists(refs))) stop("Missing endo reference: ", paste(refs[!file.exists(refs)],collapse=", "))
  ref_audit <- tibble(Reference=names(refs), Fasta=unname(refs), MD5=unname(tools::md5sum(refs)), Species=Sys.getenv("PANDORA_SPECIES"), MatchRule="Bowtie1 full-read ungapped v0; either orientation")
  evidence <- purrr::map_dfr(as.character(sample_info$Sample), function(sample_id) {
    purrr::map_dfr(names(refs), function(ref) {
      prefix <- file.path(SPORTS_RESULT_DIR,"05_endosirna","samples",sample_id,paste0(sample_id,".",ref))
      stats_path <- paste0(prefix,".stats.tsv")
      counts_path <- paste0(prefix,".sequence_counts.tsv")
      if (!file.exists(stats_path) || !file.exists(counts_path)) stop("Missing required endo evidence: ",prefix,"; run upstream for this project.")
      st <- read_tsv_flex(stats_path)
      if (nrow(st)!=1 || st$endosirna_mismatch[[1]]!=0 || st$sample[[1]]!=sample_id || st$ref_tag[[1]]!=ref) stop("Unverified/non-exact endo alignment: ",stats_path)
      ev <- read_tsv_flex(counts_path)
      if (nrow(ev)==0) return(tibble(Sample=character(),Sequence=character(),Reference=character()))
      if (!all(ev$sample==sample_id & ev$ref_tag==ref)) stop("Endo evidence sample/reference mismatch: ",counts_path)
      seqs <- unique(toupper(chartr("U","T",ev$sequence)))
      # Bowtie输出反向命中时会给出反向互补序列，因此两个方向均用于成员关系。
      tibble(Sample=sample_id, Sequence=unique(c(seqs,reverse_complement(seqs))),Reference=ref)
    })
  }) %>% distinct(Sample,Sequence,Reference) %>% group_by(Sample,Sequence) %>%
    summarise(EndoReferences=paste(sort(unique(Reference)),collapse=";"),ReferenceN=n_distinct(Reference),.groups="drop")
  # 用当前FASTA复核21-24nt全读长匹配，旧stats路径不作为参考身份的证明。
  if (!requireNamespace("Biostrings",quietly=TRUE)) stop("Biostrings is required to verify endo reference matches.")
  for (ref in names(refs)) {
    short <- evidence %>% filter(grepl(paste0("(^|;)",ref,"(;|$)"),EndoReferences),nchar(Sequence)%in%21:24) %>% pull(Sequence) %>% unique()
    seqs <- Biostrings::readDNAStringSet(refs[[ref]])
    subject <- Biostrings::DNAString(paste(as.character(seqs),collapse=strrep("N",25)))
    for (len in 21:24) {
      patterns <- short[nchar(short)==len]
      if (!length(patterns)) next
      hits <- Biostrings::countPDict(Biostrings::PDict(Biostrings::DNAStringSet(patterns)),subject) +
        Biostrings::countPDict(Biostrings::PDict(Biostrings::DNAStringSet(reverse_complement(patterns))),subject)
      if (any(hits==0)) stop("Endo evidence does not match current FASTA: ",ref," (",sum(hits==0)," sequences)")
    }
  }
  ref_audit$CurrentFastaMatchVerified21to24 <- TRUE
  list(evidence=evidence,references=ref_audit)
}
lncrna_transfer_allocations <- function(candidates, sample_info) {
  candidates <- candidates %>% filter(Annotation != "NO_Annotation")
  if (!nrow(candidates)) return(tibble(Sample=character(),Group=character(),Class=character(),Length=numeric(),Reads=numeric()))
  purrr::map_dfr(unique(candidates$Sample),function(sid) {
    ca <- candidates %>% filter(Sample==sid)
    si <- sample_info[as.character(sample_info$Sample)==sid,,drop=FALSE]
    processed <- file.path(dirname(si$OutputFile[[1]]),paste0(sid,"_processed"),paste0(sid,"_output_ensembl_match_genome"))
    if (!file.exists(processed)) stop("Cannot preserve lncRNA fractional allocation; missing: ",processed)
    al <- data.table::fread(processed,sep="\t",header=FALSE,fill=TRUE,quote="",select=c(1,4),showProgress=FALSE,data.table=FALSE)
    names(al)<-c("ID","ReferenceHeader")
    al <- as_tibble(al) %>% filter(ID %in% ca$ID) %>%
      mutate(Biotype=stringr::str_match(ReferenceHeader,"gene_biotype:([^ ]+)")[,2])
    if (any(!al$Biotype %in% c("lncRNA","lincRNA")) || !setequal(unique(al$ID),ca$ID)) stop("lncRNA allocation conflicts with SPORTS annotation: ",sid)
    al %>% add_count(ID,name="Alignments") %>% left_join(ca %>% select(ID,Reads,Length,Group),by="ID") %>%
      transmute(Sample=sid,Group,Class=paste0("ensembl-",Biotype,"_Match_Genome"),Length,Reads=Reads/Alignments) %>%
      group_by(Sample,Group,Class,Length) %>% summarise(Reads=sum(Reads),.groups="drop")
  })
}
integrate_endosirna <- function(summary, length, output, sample_info) {
  if (anyDuplicated(paste(output$Sample,output$Sequence))) stop("SPORTS output must contain one row per sample/sequence.")
  output$OriginalAnnotation <- output$Annotation
  output$EndoCandidate <- FALSE
  output$EndoReferenceOverlap <- FALSE
  output$EndoReferences <- ""
  output$ReferenceN <- 0L
  if (!ENDO_ENABLED) return(list(summary=summary,length=length,output=output,evidence=tibble(),audit=tibble(),references=tibble()))
  message("[INPUT] Reading audited endo evidence")
  reference_data <- load_endo_evidence(sample_info)
  output <- output %>% select(-EndoReferences,-ReferenceN) %>%
    left_join(reference_data$evidence,by=c("Sample","Sequence")) %>%
    mutate(EndoReferenceOverlap=!is.na(EndoReferences),EndoReferences=coalesce(EndoReferences,""),ReferenceN=coalesce(ReferenceN,0L))
  unique_annotations <- unique(output$Annotation)
  allowed <- unique_annotations[vapply(unique_annotations,function(a) {
    if (is.na(a)) return(FALSE)
    if (a=="NO_Annotation") return(TRUE)
    tokens <- strsplit(a,";",fixed=TRUE)[[1]]
    length(tokens)>0 && all(tokens %in% c("lncRNA","lincRNA"))
  },logical(1))]
  output$EndoCandidate <- output$EndoReferenceOverlap & output$Length>=21 & output$Length<=24 & output$MatchGenome=="Yes" & output$Annotation %in% allowed
  candidates <- output %>% filter(EndoCandidate)
  transfers <- bind_rows(
    candidates %>% filter(Annotation=="NO_Annotation") %>%
      transmute(Sample,Group,Class="Unannotated_Match_Genome",Length,Reads),
    lncrna_transfer_allocations(candidates,sample_info)
  ) %>% group_by(Sample,Group,Class,Length) %>% summarise(TransferredReads=sum(Reads),.groups="drop")
  root_transfers <- transfers %>% group_by(Sample,Class) %>% summarise(TransferredReads=sum(TransferredReads),.groups="drop")
  before_summary <- summary
  before_length <- length
  summary <- summary %>% left_join(root_transfers %>% mutate(SubClass="-"),by=c("Sample","Class","SubClass")) %>%
    mutate(Reads=Reads-coalesce(TransferredReads,0)) %>% select(-TransferredReads)
  length <- length %>% left_join(transfers %>% select(Sample,Class,Length,TransferredReads),by=c("Sample","Class","Length")) %>%
    mutate(Reads=Reads-coalesce(TransferredReads,0)) %>% select(-TransferredReads)
  if (any(summary$Reads < -0.51) || any(length$Reads < -0.011)) stop("Endo transfer exceeded original SPORTS allocation.")
  # 小于报告精度的四舍五入误差保持审计，不通过任意截断改写。
  endo_summary <- candidates %>% group_by(Sample,Group) %>% summarise(Reads=sum(Reads),.groups="drop") %>%
    left_join(sample_info %>% transmute(Sample=as.character(Sample),Replicate=as.character(Replicate)),by="Sample") %>%
    mutate(Class="endo-siRNA_Match_Genome",SubClass="-",Direction="sense",GenomeStatus="Match_Genome",SourceDB="endo reference",MainCategory="endo-siRNA")
  endo_length <- candidates %>% group_by(Sample,Group,Length) %>% summarise(Reads=sum(Reads),.groups="drop") %>%
    mutate(Class="endo-siRNA_Match_Genome",Direction="sense",GenomeStatus="Match_Genome",SourceDB="endo reference",MainCategory="endo-siRNA",LengthCategory="endo-siRNA")
  summary <- bind_rows(summary,endo_summary)
  length <- bind_rows(length,endo_length)
  audit <- candidates %>% group_by(Sample) %>% summarise(CandidateReads=sum(Reads),CandidateSequences=n(),.groups="drop") %>%
    right_join(sample_info %>% transmute(Sample=as.character(Sample)),by="Sample") %>%
    mutate(CandidateReads=coalesce(CandidateReads,0),CandidateSequences=coalesce(CandidateSequences,0L)) %>%
    left_join(transfers %>% group_by(Sample) %>% summarise(TransferredReads=sum(TransferredReads),.groups="drop"),by="Sample") %>%
    mutate(TransferredReads=coalesce(TransferredReads,0),Difference=CandidateReads-TransferredReads,Pass=abs(Difference)<0.00001)
  if (any(!audit$Pass)) stop("Endo transfer failed conservation.")
  evidence <- output %>% filter(EndoReferenceOverlap) %>%
    mutate(EvidenceStatus=case_when(EndoCandidate~"reference_supported_candidate",Length<21|Length>24~"outside_length_window",MatchGenome!="Yes"~"not_genome_matched",TRUE~"conflicting_or_unknown_annotation"))
  list(summary=summary,length=length,output=output,evidence=evidence,audit=audit,references=reference_data$references)
}
