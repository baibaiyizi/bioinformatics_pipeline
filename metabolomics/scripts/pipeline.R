#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i) || tail(i,1) == length(args)) default else args[[tail(i,1) + 1L]]
}
has_flag <- function(flag) flag %in% args
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "pipeline.R"
script_path <- gsub("~+~", " ", script_path, fixed=TRUE)
root_dir <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)

known <- c(
  "--config", "--module", "--result-dir", "--qc-role", "--qc-filter",
  "--qc-rsd-cutoff", "--drift-correction", "--normalization", "--run-order-adjustment",
  "--significance", "--company-comparison", "--detection-rate", "--fc-cutoff", "--p-cutoff",
  "--fdr-cutoff", "--force", "--dry-run", "--threads", "--check"
)
unknown <- args[grepl("^--", args) & !args %in% c(known, "--help")]
if (length(unknown)) stop("Unknown option: ", paste(unique(unknown), collapse = ", "), call. = FALSE)
if (has_flag("--help")) {
  cat("LC-MS 模板\n", "Usage: Rscript scripts/pipeline.R --config config/project.yml --module MODULE [options]\n", "Modules: audit inventory xcms preprocess plots annotate stats pathway report all\n", sep = "")
  quit(status = 0)
}

source(file.path(root_dir, "R", "template_contract.R"), local = TRUE)
source(file.path(root_dir, "R", "metabo_core.R"), local = TRUE)
opts <- list(
  root_dir = root_dir,
  config_file = normalizePath({ p <- value_after("--config", "config/project.yml"); if (grepl("^/", p)) p else file.path(root_dir, p) }, winslash = "/", mustWork = TRUE),
  module = value_after("--module", "all"),
  result_dir = value_after("--result-dir", NULL),
  qc_role = value_after("--qc-role", NULL),
  qc_filter = value_after("--qc-filter", NULL),
  qc_rsd_cutoff = value_after("--qc-rsd-cutoff", NULL),
  drift_correction = value_after("--drift-correction", NULL),
  normalization = value_after("--normalization", NULL),
  run_order_adjustment = value_after("--run-order-adjustment", NULL),
  significance = value_after("--significance", NULL),
  company_comparison = value_after("--company-comparison", NULL),
  detection_rate = value_after("--detection-rate", NULL),
  fc_cutoff = value_after("--fc-cutoff", NULL),
  p_cutoff = value_after("--p-cutoff", NULL),
  fdr_cutoff = value_after("--fdr-cutoff", NULL),
  threads = value_after("--threads", NULL),
  force = has_flag("--force"),
  dry_run = has_flag("--dry-run")
)
contract <- metabo_validate(opts$config_file,root_dir,opts$result_dir)
modules <- c("audit","inventory","xcms","preprocess","plots","annotate","stats","pathway","report","all","upstream","downstream")
if (!opts$module %in% modules) stop("未知模块: ",opts$module)
if (has_flag("--check") || opts$dry_run) { cat("PASS 输入配置检查；未创建输出。模块=",opts$module,"\n"); quit(status=0) }
if (!requireNamespace("digest",quietly=TRUE)) stop("需要 digest 包记录输入签名")
sig_opts <- opts[setdiff(names(opts),c("module","force","dry_run"))]
# 源码/配置/输入变化必须使用新的结果目录；不能仅凭文件存在复用。
small <- contract$files[!grepl("[.](mzML|mzXML)$",contract$files,ignore.case=TRUE)]
raw <- setdiff(contract$files,small)
sig <- digest::digest(list(opts=sig_opts,files=tools::md5sum(contract$files),code=tools::md5sum(c(script_path,list.files(file.path(root_dir,"R"),full.names=TRUE)))),algo="sha256")
sig_path <- file.path(contract$result,"input_signature.txt")
if (dir.exists(contract$result) && length(list.files(contract$result,all.files=FALSE))) {
  if (!file.exists(sig_path) || !identical(readLines(sig_path),sig)) stop("结果目录来源未验证或签名变化；请选择新的 result_dir")
}
dir.create(contract$result,recursive=TRUE,showWarnings=FALSE); writeLines(sig,sig_path)
metabo_run(opts)
