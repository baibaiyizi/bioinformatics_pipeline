write_report_index <- function(ctx,statuses) {
  module_dirs <- list(qc=SINGLE_DIR,pca=PCA_DIR,sncrna=SNC_COMPARE_DIR,
    tsrna=TSRNA_COMPARE_DIR,rsrna=RSRNA_COMPARE_DIR,mirna=MIRNA_COMPARE_DIR,
    pirna=PIRNA_COMPARE_DIR,endo=ENDOSIRNA_COMPARE_DIR,
    capture=get0("CAPTURE_DIR"),capture_bias=get0("CAPTURE_BIAS_DIR"),
    similarity=file.path(RESULT_DIR,"11_method_similarity"),
    snrna=get0("SNRNA_COMPARE_DIR"),spikein=get0("SPIKEIN_DIR"))
  completed_dirs <- unlist(module_dirs[statuses$Module[statuses$Status=="completed"]],use.names=FALSE)
  pdfs <- unique(unlist(lapply(completed_dirs,function(d)
    list.files(d,pattern="[.]pdf$",recursive=TRUE,full.names=TRUE)),use.names=FALSE))
  relative<-substring(pdfs,nchar(RESULT_DIR)+2L)
  titles<-gsub("_"," ",tools::file_path_sans_ext(basename(pdfs)))
  figs<-tibble(Figure=relative,Title=titles,PNG=sub("[.]pdf$",".png",relative),SourceData=sub("[.]pdf$",".source_data",relative))
  write_tsv_safe(figs,file.path(RESULT_DIR,"figure_index.tsv"))
  cohort<-if(PANDORA_DESIGN=="method_capture")"serum 的方法比较来自供体2，过程配对来自供体3；分装库用于技术捕获评价，不能作为独立人群样本。"else if(nrow(ctx$sample_info)==1)"单库描述：可展示组成、长度和序列；不适用PCA、组间检验或差异推断。"else if(length(ctx$group_levels)==1)"单组描述：可展示库间差异和分布；不适用组间差异推断。"else {
    meta<-read.delim(METADATA_FILE,check.names=FALSE,stringsAsFactors=FALSE)
    names(meta)<-tolower(names(meta))
    biological_source<-if("biological_sample"%in%names(meta))trimws(tolower(as.character(meta$biological_sample)))else rep(NA_character_,nrow(meta))
    source_unknown<-any(is.na(biological_source)|biological_source%in%c("","unknown","unspecified","na","n/a"))
    source_note<-if(source_unknown) {
      if(grepl("mus|mouse",tolower(PROJECT_SPECIES)))"biological_sample未填写或unknown，来源母鼠/窝关系unknown，生物学独立性未核验。"else"biological_sample未填写或unknown，来源个体关系unknown，生物学独立性未核验。"
    }else"biological_sample已填写；来源标签本身不能证明生物学独立性，仍需核验来源个体及采样关系。"
    paste0("比较方向：",Sys.getenv("CONTRAST_CASE")," - ",Sys.getenv("CONTRAST_BASELINE"),"。按pool/library标签进行探索性比较（pool/library-level exploratory）；",source_note)
  }
  text<-c(paste0("# ",Sys.getenv("PANDORA_PROJECT_NAME")," sncRNA 分析"),"",paste0("阶段：",PANDORA_ANALYSIS_STAGE,"；流程：",if(ENDO_ENABLED)"含 endo-siRNA"else"不含 endo-siRNA","；运行状态：completed。"),"",cohort,"","所有原始 reads 均保留。全库主类别组成和全库 RPM 使用同一 clean reads 分母；类别内部比例的分母见轴标和源表；组级描述以样本等权汇总。SPORTS 母类/子类的分数分配保留，子集不重复加到总量。", "",if(ENDO_ENABLED)"endo-siRNA 为参考支持的候选类别：21–24 nt、基因组匹配、全读长零错配参考命中，原注释为 NO_Annotation 或全部为 lncRNA/lincRNA。其他 sncRNA 冲突只保留证据；不据此证明 Dicer/AGO2 依赖或功能。"else"本流程不加载或展示 endo-siRNA 证据。", "","图片统一183 mm宽，Liberation Sans，PDF与300dpi PNG配对导出；完整图目录及源数据状态见 figure_index.tsv 和每图 .figure_manifest.csv。", "","| 模块 | 状态 |","|---|---|",paste0("|",statuses$Module,"|",statuses$Status,"|"),"","| 图片 | PDF | PNG |","|---|---|---|",paste0("|",titles,"|[PDF](",relative,")|[PNG](",figs$PNG,")|"))
  writeLines(text,file.path(RESULT_DIR,"README.md"))
  esc<-function(x)gsub("&","&amp;",gsub("<","&lt;",x,fixed=TRUE),fixed=TRUE)
  cards<-paste0('<article data-name="',esc(titles),'"><h3>',esc(titles),'</h3><a href="',relative,'"><img loading="lazy" src="',figs$PNG,'" alt="',esc(titles),'"></a><p><a href="',relative,'">PDF</a> · <a href="',figs$PNG,'">PNG</a></p></article>')
  html<-c('<!doctype html><html lang="zh"><meta charset="utf-8"><title>sncRNA 结果</title><style>body{font:16px sans-serif;max-width:1300px;margin:auto;padding:25px;color:#263445;background:#fff}input{padding:10px;width:60%}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:22px}article{border:1px solid #ddd;padding:12px}img{width:100%;height:auto}h3{font-size:14px}</style>',paste0('<h1>',Sys.getenv("PANDORA_PROJECT_NAME"),' sncRNA ',PANDORA_ANALYSIS_STAGE,'</h1><p>',esc(cohort),'</p><p>原始 reads 全部保留；点击图片打开 PDF。<a href="README.md">分析说明</a></p><input id="q" placeholder="搜索图名或类型"><main>'),cards,'</main><script>document.querySelector("#q").oninput=e=>document.querySelectorAll("article").forEach(a=>a.hidden=!a.dataset.name.toLowerCase().includes(e.target.value.toLowerCase()))</script></html>')
  writeLines(html,file.path(RESULT_DIR,"index.html"))
}
