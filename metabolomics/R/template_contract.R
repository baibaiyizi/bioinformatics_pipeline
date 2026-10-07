# 输入和重用检查在创建结果目录之前执行。
metabo_validate <- function(config_file, root, result_override=NULL) {
  if (!requireNamespace("yaml",quietly=TRUE)) stop("需要 yaml 包")
  cfg <- yaml::read_yaml(config_file)
  resolve <- function(x) if (startsWith(x,"/")) x else file.path(root,x)
  required <- c("project_name","sample_info_file","case_group","control_group","input","result_dir")
  if (!all(required %in% names(cfg))) stop("配置缺少: ",paste(setdiff(required,names(cfg)),collapse=","))
  if (identical(cfg$case_group,cfg$control_group) || any(!nzchar(c(cfg$case_group,cfg$control_group)))) stop("case/control 必须明确且不同")
  if (length(unique(make.names(c(cfg$case_group,cfg$control_group)))) != 2) stop("组名在 R 设计矩阵中冲突")
  if (!cfg$input$mode %in% c("mzml","peak_table")) stop("input.mode 只能是 mzml 或 peak_table")
  sample_file <- resolve(cfg$sample_info_file)
  si <- read.delim(sample_file,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("","NA"))
  cols <- c("sample","group","biological_sample_id")
  if (!all(cols %in% names(si)) || !nrow(si) || anyNA(si[,cols]) || anyDuplicated(si$sample)) stop("样本表需要唯一 sample、group 和 biological_sample_id")
  if ("include" %in% names(si)) si <- si[tolower(as.character(si$include)) %in% c("true","1","yes"),,drop=FALSE]
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",si$sample))) stop("sample 仅允许字母数字下划线点横线")
  qc <- cfg$qc_group %||% "QC"
  if (!all(si$group %in% c(cfg$case_group,cfg$control_group,qc))) stop("发现未配置的组别")
  bio <- si[si$group != qc,,drop=FALSE]
  if (any(tolower(bio$biological_sample_id) %in% c("unknown","na")) || anyDuplicated(bio$biological_sample_id)) stop("本模板要求一行一个独立生物学样本；重复进样请预先审计汇总")
  if (any(table(factor(bio$group,levels=c(cfg$control_group,cfg$case_group))) < 3)) stop("两组均至少需要 3 个独立样本")
  if (any(c("subject_id","pair_id") %in% names(si))) stop("首版不支持配对或重复测量设计")
  files <- c(config_file,sample_file)
  if (cfg$input$mode == "peak_table") {
    for (name in c("filled_table","feature_table")) if (is.null(cfg$input[[name]]) || !nzchar(cfg$input[[name]])) stop("input 缺少 ",name)
    files <- c(files,resolve(cfg$input$filled_table),resolve(cfg$input$feature_table))
    if (!is.null(cfg$input$detected_table) && nzchar(cfg$input$detected_table)) files <- c(files,resolve(cfg$input$detected_table))
    peak <- read.csv(resolve(cfg$input$filled_table),check.names=FALSE)
    if (names(peak)[1] != "feature_id" || !setequal(names(peak)[-1],si$sample)) stop("filled_table 的样本列必须与纳入样本表完全一致")
  } else {
    if (is.null(cfg$input$mzml_manifest)) stop("需要 mzml_manifest")
    path <- resolve(cfg$input$mzml_manifest); files <- c(files,path)
    m <- read.delim(path,stringsAsFactors=FALSE)
    if (!all(c("file","sample","polarity") %in% names(m)) || !nrow(m) || anyNA(m) || anyDuplicated(m[c("sample","polarity")])) stop("mzML manifest 需要唯一 sample×polarity 及 file")
    if (!setequal(m$sample,si$sample) || !all(m$polarity %in% c("POS","NEG"))) stop("mzML manifest 样本或极性无效")
    for (pol in unique(m$polarity)) if (!setequal(m$sample[m$polarity==pol],si$sample)) stop("每个 polarity 必须覆盖全部纳入样本；不能将未采集极性填补为低丰度")
    mzfiles <- vapply(m$file,resolve,character(1)); if (any(!file.exists(mzfiles))) stop("缺少 mzML 文件")
    files <- c(files,mzfiles)
  }
  optfiles <- c(cfg$annotation$table,cfg$external_statistics_table,cfg$pathway$quantitative_table,cfg$pathway_enrichment_file)
  for (x in optfiles) if (!is.null(x) && nzchar(x)) files <- c(files,resolve(x))
  if (any(!file.exists(files))) stop("缺少输入: ",paste(files[!file.exists(files)],collapse=", "))
  result <- resolve(result_override %||% cfg$result_dir)
  list(config=cfg,files=files,result=result)
}
`%||%` <- function(x,y) if (is.null(x) || !length(x)) y else x
