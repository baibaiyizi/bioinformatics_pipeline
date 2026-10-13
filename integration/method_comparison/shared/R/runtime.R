SHARED_RMD_DIR <- file.path(Sys.getenv("PANDORA_SHARED_ROOT"), "rmd")
PANDORA_DIFFERENTIAL_RESULT_DIR <- file.path(ROOT_DIR,"result","differential")
PANDORA_DESIGN <- Sys.getenv("PANDORA_DESIGN")
ENDOSIRNA_COMPARE_DIR <- file.path(RESULT_DIR,"9_endosirna")
ENDOSIRNA_COMPARE_FIG <- file.path(ENDOSIRNA_COMPARE_DIR,"plots")
ENDOSIRNA_COMPARE_TAB <- file.path(ENDOSIRNA_COMPARE_DIR,"tables")
PROJECT_STAGE <- Sys.getenv("PROJECT_STAGE")
sncrna_heatmap_gradient <- c("#A1C9F4","#FDFFB6","#FF9F9B")
sncrna_theme_minimal <- function(...) ggplot2::theme_minimal(...) + figure_theme() + theme(panel.grid=element_blank())
if (ENDO_ENABLED) {
  category_levels <- append(category_levels,"endo-siRNA",after=4)
  main_palette <- c(main_palette,`endo-siRNA`="#C7B9D9")
  sports_length_categories <- append(sports_length_categories,"endo-siRNA",after=2)
  radar_levels <- append(radar_levels,"endo-siRNA",after=match("piRNA",radar_levels))
  radar_labels <- c(radar_labels,`endo-siRNA`="endo-siRNA")
}
category_palette <- main_palette[category_levels]
sports_length_palette <- main_palette[sports_length_categories]
radar_sector_palette <- main_palette[radar_levels]
radar_label_palette <- make_sector_label_palette(radar_sector_palette,amount=.28)
# 含endo和不含endo都使用完整主类别；被归入other的组成不丢弃。
# 所有上下文复用同一次输入读取与endo归类；scope不同的serum队列保持隔离。
.context_cache <- new.env(parent=emptyenv())
cache_context <- function(loader) {
  force(loader)
  function(...) {
    values <- unlist(list(...))
    if (!length(values)) values <- REPORT_VIEW
    if (Sys.getenv("PANDORA_SPECIES")=="Homo" && length(values)==1) values <- c(values,"all")
    key <- paste(values,collapse="|")
    if (!exists(key,envir=.context_cache,inherits=FALSE)) {
      ctx <- loader(...)
      # output与独立主类别均必须恢复clean reads（SPORTS两位小数允许2 reads报告误差）。
      outcheck <- ctx$output_all %>% group_by(Sample) %>% summarise(OutputReads=sum(Reads),.groups="drop") %>%
        left_join(ctx$clean_reads_tbl,by="Sample") %>% mutate(Difference=OutputReads-CleanReads,Pass=abs(Difference)<1e-6)
      if (any(!outcheck$Pass)) stop("SPORTS sequence count does not recover clean reads.")
      ctx$output_integrity_check <- outcheck
      ld <- ctx$length_all %>% filter(!Class%in%c("Clean_Reads","Match_Genome"),!is_tsrna_nested_class(Class),!is_nested_length_class(Class))
      ctx$length_composition_integrity_check <- ld %>% group_by(Sample) %>% summarise(ExclusiveLengthReads=sum(Reads),.groups="drop") %>%
        left_join(ctx$clean_reads_tbl,by="Sample") %>% mutate(ReadDifference=ExclusiveLengthReads-CleanReads,Pass=abs(ReadDifference)<=2)
      if (any(!ctx$length_composition_integrity_check$Pass)) stop("Length classification failed read conservation.")
      assign(key,ctx,envir=.context_cache)
    }
    get(key,envir=.context_cache,inherits=FALSE)
  }
}
if (Sys.getenv("PANDORA_SPECIES")=="Mus") {
  load_mii_context <- cache_context(load_mii_context)
} else {
  load_serum_context <- cache_context(load_serum_context)
  initialize_profile_comparison_module <- function(env,default_fig_dir,default_tab_dir,replay_module,sections) {
    scope <- get0("PANDORA_ANALYSIS_SCOPE",envir=env,inherits=TRUE,ifnotfound=NULL)
    if (is.null(scope)) scope <- if (PANDORA_DESIGN=="method_capture") "method_comparison" else "all"
    assign("PANDORA_ANALYSIS_SCOPE",scope,envir=env)
    pair <- identical(scope,"spikein_pair")
    root <- file.path(SPIKEIN_DIR,"profiles",replay_module)
    fig <- if(pair) file.path(root,"plots") else default_fig_dir
    tab <- if(pair) file.path(root,"tables") else default_tab_dir
    assign("OUTPUT_FIG_DIR",fig,envir=env); assign("OUTPUT_TAB_DIR",tab,envir=env)
    prepare_module_outputs(fig,tab,sections)
    activate_serum_context(env,analysis_scope=scope)
    invisible(pair)
  }
}
write_common_audit <- function(ctx) {
  dir.create(SETUP_DIR,recursive=TRUE,showWarnings=FALSE)
  for (n in c("endo_evidence","endo_transfer_audit","endo_reference_audit","output_integrity_check","composition_integrity_check","length_composition_integrity_check")) {
    if (startsWith(n,"endo") && !ENDO_ENABLED) next
    write_tsv_safe(ctx[[n]],file.path(SETUP_DIR,paste0(n,".tsv")))
  }
  # 权威类别矩阵用于QC、PCA、DE和报告；不以首个注释替代SPORTS分配。
  cls <- ctx$sample_comp_view %>% transmute(Sample=as.character(Sample),MainCategory=as.character(MainCategory),Reads,RPM=1e6*Reads/CleanReads)
  for (metric in c("Reads","RPM")) {
    wide <- cls %>% select(MainCategory,Sample,all_of(metric)) %>% pivot_wider(names_from=Sample,values_from=all_of(metric),values_fill=0)
    write_tsv_safe(wide,file.path(EXPRESSION_MATRIX_DIR,paste0("class_",tolower(metric),".tsv")))
  }
  write_tsv_safe(ctx$sample_info,file.path(SETUP_DIR,"sample_metadata_used.tsv"))
}

# 下游只消费已完成且输入/源码签名仍有效的显式阶段产物。
validate_stage_artifact <- function(path,stage,module) {
  fail <- function(reason) list(valid=FALSE,reason=reason)
  stage_dir <- file.path(ROOT_DIR,"result",stage)
  required <- file.path(stage_dir,c("run_status.tsv","module_status.tsv","input_signature.tsv"))
  if (!file.exists(path) || any(!file.exists(required))) return(fail("not_run: missing explicit upstream artifact or receipt"))
  receipt <- read_tsv_flex(required[[1]])
  modules <- read_tsv_flex(required[[2]])
  if (nrow(receipt)!=1 || receipt$State[[1]]!="completed" ||
      !any(modules$Module==module & modules$Status=="completed")) return(fail("not_run: upstream stage/module is not completed"))
  signature <- read_tsv_flex(required[[3]])
  if (!nrow(signature) || any(!file.exists(signature$Path))) return(fail("STALE_OUTPUTS: upstream input is missing"))
  if (any(file.info(signature$Path)$size!=signature$Size) ||
      any(unname(tools::md5sum(signature$Path))!=signature$MD5)) return(fail("STALE_OUTPUTS: upstream inputs or analysis code changed"))
  list(valid=TRUE,reason="Completed upstream module with current input signature")
}
