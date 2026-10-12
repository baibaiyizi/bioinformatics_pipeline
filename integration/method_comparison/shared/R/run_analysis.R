#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name, default=NULL) {i<-match(name,args);if(is.na(i))default else args[[i+1L]]}
PANDORA_ANALYSIS_STAGE <- arg_value("--stage","describe")
stopifnot(PANDORA_ANALYSIS_STAGE %in% c("describe","differential","enrichment"))
root <- normalizePath(Sys.getenv("PANDORA_PROJECT_ROOT"),mustWork=TRUE)
shared <- normalizePath(Sys.getenv("PANDORA_SHARED_ROOT"),mustWork=TRUE)
species <- tolower(Sys.getenv("PANDORA_SPECIES"))
code_files <- c(file.path(dirname(shared),"run.sh"),file.path(dirname(shared),"scripts/check_inputs.py"),list.files(shared,pattern="[.](R|Rmd|sh|py)$",recursive=TRUE,full.names=TRUE))
code_signature <- tools::md5sum(code_files)
Sys.setenv(PANDORA_ANALYSIS_STAGE=PANDORA_ANALYSIS_STAGE)
source(file.path(shared,"R",paste0("setup_",species,".R")),local=globalenv())
source(file.path(shared,"R","endosirna.R"),local=globalenv())
source(file.path(shared,"R","runtime.R"),local=globalenv())
source(file.path(shared,"R","figure_contract.R"),local=globalenv())
dir.create(RESULT_DIR,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(RESULT_DIR,"reports"),recursive=TRUE,showWarnings=FALSE)
# 分析状态只在全部所选模块成功时完成。失败运行不冒充有效结果。
receipt_path <- file.path(RESULT_DIR,"run_status.tsv")
write_status <- function(state,detail="") write_tsv_safe(tibble(State=state,Project=Sys.getenv("PANDORA_PROJECT_NAME"),Species=species,Mode=if(ENDO_ENABLED)"with_endo"else"without_endo",Stage=PANDORA_ANALYSIS_STAGE,Time=as.character(Sys.time()),Detail=gsub("[\r\n\t]+"," | ",detail)),receipt_path)
write_status("running")
run_module <- function(file,analysis_scope=NULL,analysis_stage=PANDORA_ANALYSIS_STAGE) {
  module_env <- new.env(parent=globalenv())
  # 原模块helper依赖当前模块的数据与输出目录；为每次knit绑定独立作用域。
  for (nm in ls(globalenv(),all.names=TRUE)) {
    obj <- get(nm,envir=globalenv())
    if (is.function(obj) && identical(environment(obj),globalenv())) {
      environment(obj) <- module_env
      assign(nm,obj,envir=module_env)
    }
  }
  module_env$PANDORA_CURRENT_MODULE <- basename(file)
  module_env$PANDORA_ANALYSIS_STAGE <- analysis_stage
  module_env$PANDORA_ANALYSIS_SCOPE <- if(species=="mus") analysis_scope %||% "all" else analysis_scope
  module_env$PANDORA_DIFFERENTIAL_RESULT_DIR <- PANDORA_DIFFERENTIAL_RESULT_DIR
  dest <- file.path(RESULT_DIR,"reports",paste0(if(!is.null(analysis_scope))paste0(analysis_scope,"_"),tools::file_path_sans_ext(basename(file)),".md"))
  message("[MODULE] ",basename(file)," scope=",analysis_scope %||% "default", " stage=",analysis_stage)
  # 子队列重放保留独立chunk命名，避免与父文档的init等标签冲突。
  module_text <- readLines(file,warn=FALSE)
  label_prefix <- paste0(safe_dir_name(paste(basename(file),analysis_scope %||% "all",sep="_")),"__")
  named_chunk <- grepl("^```\\{r[[:space:]]+[^,[:space:]}]+",module_text)
  module_text[named_chunk] <- sub("^```\\{r[[:space:]]+",paste0("```{r ",label_prefix),module_text[named_chunk])
  knitr::knit(text=module_text,output=dest,envir=module_env,quiet=TRUE)
  invisible(dest)
}
tryCatch({
  ctx <- if(species=="mus")load_mii_context(REPORT_VIEW)else load_serum_context(REPORT_VIEW,"all")
  # 为全局 helper 保留当前作用域；模块局部仍用activate选择独立队列。
  list2env(ctx,envir=globalenv())
  # 联合算法需重建序列层特征；必须先与单组学标准类别计数逐值一致。
  standard <- read_tsv_flex(file.path(ROOT_DIR,"result/describe/1_setup/1_expression_matrix/class_reads.tsv"))
  standard_long <- standard %>% tidyr::pivot_longer(-MainCategory,names_to="Sample",values_to="StandardReads")
  regenerated <- ctx$sample_comp_view %>% transmute(Sample=as.character(Sample),MainCategory=as.character(MainCategory),Reads) %>%
    tidyr::pivot_wider(names_from=Sample,values_from=Reads,values_fill=0) %>%
    tidyr::pivot_longer(-MainCategory,names_to="Sample",values_to="Reads")
  boundary <- full_join(standard_long,regenerated,by=c("Sample","MainCategory"))
  if(anyNA(boundary$StandardReads) || anyNA(boundary$Reads) || any(abs(boundary$StandardReads-boundary$Reads)>1e-6)) stop("UNVERIFIED: 联合输入与单组学标准类别计数不一致")
  write_tsv_safe(boundary,file.path(RESULT_DIR,"single_assay_boundary_check.tsv"))
  if(species=="mus")write_mii_setup_expression_matrices()else write_serum_setup_expression_matrices()
  write_common_audit(ctx)
  files <- c(capture="09_differential_capture.Rmd",capture_bias="10_capture_bias_diagnostics.Rmd",similarity="11_method_similarity.Rmd",spikein="13_spikein_qc.Rmd")
  if(ENDO_ENABLED)files<-c(files,endo="09_endo_compare.Rmd")
  applicable <- if(PANDORA_ANALYSIS_STAGE=="describe")setdiff(names(files),c("capture","capture_bias")) else if(PANDORA_ANALYSIS_STAGE=="differential") {
    if(species=="mus")c("pca","sncrna","tsrna","rsrna","mirna","pirna",if(ENDO_ENABLED)"endo")else c("capture","capture_bias")
  }else if(species=="mus")"mirna"else character()
  if(nrow(ctx$sample_info)<3)applicable<-setdiff(applicable,"pca")
  if(PANDORA_DESIGN!="method_capture")applicable<-setdiff(applicable,c("capture","capture_bias","similarity","spikein"))
  requested <- arg_value("--modules")
  if(!is.null(requested)) {
    selected <- unlist(strsplit(requested,"[, ]+"))
    invalid <- setdiff(selected,names(files))
    if(length(invalid))stop("Unknown module: ",paste(invalid,collapse=", "))
    if(length(setdiff(selected,applicable)))stop("Module not applicable to this stage/design: ",paste(setdiff(selected,applicable),collapse=", "))
    if("capture_bias"%in%selected)selected<-unique(c("capture",selected))
    applicable<-intersect(names(files),selected)
  }
  statuses <- tibble(Module=names(files),File=unname(files),Status=if_else(names(files)%in%applicable,"pending","not_applicable_or_not_selected"))
  write_tsv_safe(statuses,file.path(RESULT_DIR,"module_status.tsv"))
  for(m in applicable) {
    module_file<-if(m=="endo")file.path(shared,"rmd",files[[m]])else file.path(SHARED_RMD_DIR,species,files[[m]])
    run_module(module_file)
    statuses$Status[statuses$Module==m]<-"completed"
    write_tsv_safe(statuses,file.path(RESULT_DIR,"module_status.tsv"))
  }
  inputs<-c(file.path(Sys.getenv("METHOD_COMPARISON_ROOT"),"config/project.sh"),METADATA_FILE,file.path(ROOT_DIR,"config/project.sh"),file.path(ROOT_DIR,"result/describe/input_signature.tsv"),file.path(ROOT_DIR,"result/describe/run_status.tsv"),ctx$sample_info$SummaryFile,ctx$sample_info$LengthFile,ctx$sample_info$OutputFile,code_files)
  if(ENDO_ENABLED)inputs<-c(inputs,list.files(file.path(SPORTS_RESULT_DIR,"05_endosirna"),pattern="(stats|sequence_counts)[.]tsv$",recursive=TRUE,full.names=TRUE),endo_reference_files())
  inputs<-unique(inputs[file.exists(inputs)])
  signature<-tibble(Path=inputs,Size=file.info(inputs)$size,MD5=unname(tools::md5sum(inputs)))
  write_tsv_safe(signature,file.path(RESULT_DIR,"input_signature.tsv"))
  if (!identical(unname(tools::md5sum(code_files)),unname(code_signature))) stop("STALE_OUTPUTS: analysis code changed during this run; rerun with frozen source files.")
  write_status("completed")
  source(file.path(shared,"R","report_index.R"),local=globalenv())
  write_report_index(ctx,statuses)
},interrupt=function(e){write_status("interrupted","Analysis interrupted before completion");stop("Analysis interrupted")},error=function(e){write_status("failed",conditionMessage(e));stop(e)})
