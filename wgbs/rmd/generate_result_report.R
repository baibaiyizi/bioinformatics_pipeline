#!/usr/bin/env Rscript

required_packages <- c("htmltools")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing_packages) > 0) stop("Missing report dependencies: ", paste(missing_packages, collapse = ", "))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if(length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else ""
default_root <- if(nzchar(script_path)) dirname(dirname(normalizePath(script_path, mustWork = FALSE))) else getwd()
root_dir <- normalizePath(Sys.getenv("RESULT_REPORT_ROOT", unset = default_root), winslash = "/", mustWork = TRUE)
result_dir <- file.path(root_dir, "result")
output_file <- file.path(result_dir, "index.html")
if(!dir.exists(result_dir)) stop("Result directory does not exist: ", result_dir)

path_parts <- strsplit(root_dir, "/", fixed = TRUE)[[1]]
assay <- "WGBS"
project_label <- basename(root_dir)

relative_path <- function(path) {
  root_norm <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root_norm, "/")
  if(!startsWith(path_norm, prefix)) stop("Report input is outside result/: ", path)
  substring(path_norm, nchar(prefix) + 1L)
}

file_url <- function(path) utils::URLencode(gsub("\\\\", "/", path), reserved = FALSE)
module_of <- function(path) strsplit(path, "/", fixed = TRUE)[[1]][[1]]
safe_id <- function(x) paste0("section-", gsub("[^A-Za-z0-9_-]+", "-", x))
pretty_name <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  stem <- gsub("_", " ", stem, fixed = TRUE)
  gsub("\\s+", " ", stem)
}

read_contrast <- function() {
  reference <- Sys.getenv("WGBS_REFERENCE_GROUP")
  comparison <- Sys.getenv("WGBS_COMPARISON_GROUP")
  if(!nzchar(reference) || !nzchar(comparison)) stop("报告需要通过入口加载显式分组配置")
  paste(comparison, reference, sep = " - ")
}

contrast <- read_contrast()
# 只有最新一次实际完成的模块才向报告提供图表；失败/中断留下的旧文件不纳入。
module_key <- function(x) {
  first <- vapply(strsplit(x, "/", fixed = TRUE), `[`, character(1), 1L)
  numbered <- grepl("^[0-9]+", first)
  first[numbered] <- as.character(as.integer(sub("^([0-9]+).*", "\\1", first[numbered])))
  first
}
execution_path <- file.path(result_dir, "execution.tsv")
execution <- if(file.exists(execution_path)) utils::read.delim(execution_path, colClasses = "character", check.names = FALSE) else data.frame()
if(nrow(execution) && !all(c("module", "status") %in% names(execution))) stop("execution.tsv 缺少 module/status 列")
latest <- if(nrow(execution)) {
  execution$key <- module_key(execution$module)
  execution[!duplicated(execution$key, fromLast = TRUE), , drop = FALSE]
} else data.frame(module = character(), status = character(), key = character())
latest <- latest[latest$key != "report", , drop = FALSE]
module_state <- setNames(latest$status, latest$key)
status_note <- htmltools::tags$aside(
  class = "report-note",
  htmltools::tags$h2("本次模块执行状态"),
  htmltools::tags$p("仅收录最新状态为 COMPLETED 的模块图表。FAILED、STARTED、SKIPPED 和未记录状态的旧产物不作为当前结果。输入资格、模块执行与真实数据复现分别记录。"),
  if(!nrow(latest)) htmltools::tags$p("UNVERIFIED: 尚无模块执行记录；不展示现存分析图表。"),
  htmltools::tags$ul(lapply(seq_len(nrow(latest)), function(i) htmltools::tags$li(paste(latest$module[i], latest$status[i], sep = ": ")))),
  if(file.exists(execution_path)) htmltools::tags$a(href = "execution.tsv", "查看完整执行历史")
)
all_files <- list.files(result_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
all_files <- all_files[normalizePath(all_files, winslash = "/", mustWork = FALSE) != normalizePath(output_file, winslash = "/", mustWork = FALSE)]
relative_files <- vapply(all_files, relative_path, character(1))
keys <- module_key(relative_files)
numbered <- grepl("^[0-9]+$", keys)
completed <- !is.na(module_state[keys]) & module_state[keys] == "COMPLETED"
# 顶层状态/日志是审计记录；所有编号模块产物都受执行状态门约束。
relative_files <- relative_files[completed | relative_files %in% c("execution.tsv", "optional_module_gates.tsv") | startsWith(relative_files, "logs/")]
png_paths <- sort(relative_files[grepl("\\.png$", relative_files, ignore.case = TRUE)])
pdf_paths <- sort(relative_files[grepl("\\.pdf$", relative_files, ignore.case = TRUE)])
table_paths <- sort(relative_files[grepl("\\.(csv|tsv|txt)$", relative_files, ignore.case = TRUE)])
other_html <- sort(relative_files[grepl("\\.html?$", relative_files, ignore.case = TRUE)])



modules <- sort(unique(c(
  vapply(png_paths, module_of, character(1)),
  if(length(table_paths) > 0) vapply(table_paths, module_of, character(1)) else character(0),
  if(length(other_html) > 0) vapply(other_html, module_of, character(1)) else character(0)
)))

module_titles <- c(
  "2_qc" = "质量控制 / Quality control",
  "3_dmr" = "差异甲基化区域 / Differentially methylated regions",
  "4_dmroverview" = "DMR 总览 / DMR overview",
  "5_enrichment" = "功能富集 / Functional enrichment",
  "6_genomecontext" = "基因组背景 / Genomic context",
  "7_regulatory" = "调控元件与 motif / Regulatory elements and motifs",
  "8_gene_methylation" = "基因甲基化汇总 / Gene methylation summary",
  "10_browser" = "基因组浏览轨迹 / Genome browser tracks",
  "11_sperm" = "精子特异分析 / Sperm-specific analysis",
  "cross_stage_analysis" = "跨组织比较 / Cross-context comparison",
  "02_project_inventory" = "项目与输入审计 / Project and input inventory",
  "03_atac_qc" = "质量控制 / Quality control",
  "04_peak_catalog_annotation" = "Peak 目录与注释 / Peak catalog and annotation",
  "05_differential_accessibility" = "差异可及性 / Differential accessibility",
  "06_genome_tracks" = "基因组轨迹 / Genome tracks",
  "07_peak_gene_enrichment" = "Peak-gene 富集 / Peak-gene enrichment",
  "08_motif_activity" = "Motif 活性 / Motif activity",
  "10_footprinting_tf_evidence" = "足迹与 TF 证据 / Footprinting and TF evidence",
  "12_tf_gene_network" = "TF-gene 网络 / TF-gene network",
  "13_sequence_conservation_gc" = "序列与保守性 / Sequence and conservation",
  "14_machine_learning" = "预测建模 / Predictive modeling",
  "15_keygenes_locus_reports" = "候选基因与位点报告 / Candidate genes and loci",
  "16_export_summary_manifest" = "导出审计 / Export and manifest audit",
  "IGV" = "IGV 文件 / IGV resources",
  "logs" = "运行日志 / Run logs"
)
module_title <- function(id) {
  value <- unname(module_titles[id])
  if(is.null(value) || is.na(value)) gsub("_", " ", id, fixed = TRUE) else value
}

scientific_notes <- if(identical(assay, "WGBS")) {
  c(
    "甲基化差异方向按 comparison - reference 定义；Hyper/Hypo 的数值与阈值以相邻结果表为准。",
    "DMR、区域注释和富集图分别回答区域差异、基因组位置及集合过度代表性问题，不能互相替代。",
    "基因关联和位点描述不构成调控因果证据。",
    "GO Chord/Cluster 最多显示 10 个术语和 20 个代表基因；treeplot 最多 8 个术语，cnetplot 默认 5 个术语。",
    "富集网络和 GOplot 标签仅为排版截断至 38 个字符；完整术语、统计量和成员关系以相邻结果表为准。"
  )
} else {
  c(
    "差异可及性方向按 comparison - reference 定义；效应量、p 值/FDR 和筛选口径以相邻结果表为准。",
    "Peak 注释、motif 富集、footprint 和 TF-gene 网络属于不同证据层，网络位置或中心性不等同于因果作用。",
    "ATAC 候选和 motif 位点关联不能替代功能实验。",
    "图内长术语可能为排版而缩写；完整名称、统计量和成员关系以 CSV/TSV 表为准。"
  )
}

make_link <- function(path, label = basename(path), class = NULL) {
  htmltools::tags$a(href = file_url(path), target = "_blank", rel = "noopener", class = class, label)
}

make_card <- function(path) {
  pdf_path <- sub("\\.png$", ".pdf", path, ignore.case = TRUE)
  if(!file.exists(file.path(result_dir, pdf_path))) {
    multipage_pdf <- sub("_page_[0-9]+\\.png$", ".pdf", path, ignore.case = TRUE)
    if(!identical(multipage_pdf, path) && file.exists(file.path(result_dir, multipage_pdf))) {
      pdf_path <- multipage_pdf
    }
  }
  has_pdf <- file.exists(file.path(result_dir, pdf_path))
  title <- pretty_name(path)
  htmltools::tags$article(
    class = "figure-card",
    `data-search` = tolower(paste(path, title, module_title(module_of(path)))),
    htmltools::tags$div(class = "card-head", htmltools::tags$h3(title), htmltools::tags$code(path)),
    htmltools::tags$a(
      href = file_url(path), target = "_blank", rel = "noopener",
      htmltools::tags$img(src = file_url(path), alt = title, loading = "lazy")
    ),
    htmltools::tags$div(
      class = "downloads",
      make_link(path, "PNG"),
      if(has_pdf) make_link(pdf_path, "PDF") else htmltools::tags$span(class = "missing", "PDF missing")
    ),
    htmltools::tags$p(class = "card-note", "图形用于描述当前模块结果；完整数值、筛选口径和未缩写字段以同目录表格为准。")
  )
}

make_module <- function(module_id) {
  figs <- png_paths[vapply(png_paths, module_of, character(1)) == module_id]
  tables <- table_paths[vapply(table_paths, module_of, character(1)) == module_id]
  reports <- other_html[vapply(other_html, module_of, character(1)) == module_id]
  htmltools::tags$section(
    id = safe_id(module_id), class = "module-section", `data-module` = module_id,
    htmltools::tags$div(
      class = "module-heading",
      htmltools::tags$h2(module_title(module_id)),
      htmltools::tags$span(class = "count", sprintf("%d figures · %d tables", length(figs), length(tables)))
    ),
    if(length(figs) > 0) htmltools::tags$div(class = "figure-grid", lapply(figs, make_card)),
    if(length(tables) + length(reports) > 0) htmltools::tags$details(
      class = "file-list",
      htmltools::tags$summary("表格与附加报告 / Tables and supplementary reports"),
      htmltools::tags$ul(c(
        lapply(tables, function(x) htmltools::tags$li(make_link(x, x))),
        lapply(reports, function(x) htmltools::tags$li(make_link(x, x)))
      ))
    )
  )
}

css <- htmltools::HTML("\n:root{--ink:#27323a;--muted:#64727c;--line:#dfe4e7;--panel:#f7f8f8;--blue:#4d7ea8;--red:#c96f67}\n*{box-sizing:border-box}body{margin:0;font-family:'Liberation Sans',Arial,sans-serif;color:var(--ink);background:white;line-height:1.5}\na{color:#315f84}.page{max-width:1500px;margin:auto;padding:28px}.hero{border-bottom:1px solid var(--line);padding-bottom:20px;margin-bottom:20px}\nh1{font-size:28px;margin:0 0 8px}.subtitle{color:var(--muted);margin:0}.summary{display:flex;gap:10px;flex-wrap:wrap;margin-top:14px}.badge{background:var(--panel);border:1px solid var(--line);border-radius:999px;padding:5px 10px;font-size:13px}\n.report-note{background:#fbfbf7;border-left:4px solid #d8c86c;padding:14px 18px;margin:20px 0}.report-note h2{font-size:17px;margin:0 0 8px}.report-note ul{margin:0;padding-left:20px}\n.toolbar{position:sticky;top:0;z-index:5;background:rgba(255,255,255,.96);border-bottom:1px solid var(--line);padding:12px 0;margin-bottom:18px}.toolbar input{width:min(620px,100%);padding:10px 12px;border:1px solid #b9c2c8;border-radius:6px;font-size:15px}\nnav{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0 20px}nav a{font-size:13px;text-decoration:none;border:1px solid var(--line);border-radius:5px;padding:5px 8px;background:white}\n.module-section{margin:30px 0 42px}.module-heading{display:flex;align-items:baseline;justify-content:space-between;border-bottom:2px solid var(--ink);gap:12px;margin-bottom:16px}.module-heading h2{font-size:21px;margin:0 0 7px}.count{color:var(--muted);font-size:13px;white-space:nowrap}\n.figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(390px,1fr));gap:18px}.figure-card{border:1px solid var(--line);border-radius:8px;background:white;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.05)}.card-head{padding:12px 14px 8px}.card-head h3{font-size:16px;margin:0 0 4px}.card-head code{display:block;color:var(--muted);font-size:11px;overflow-wrap:anywhere}.figure-card img{display:block;width:100%;height:330px;object-fit:contain;background:white;border-top:1px solid #eef0f1;border-bottom:1px solid #eef0f1}.downloads{display:flex;gap:8px;padding:10px 14px}.downloads a,.missing{border:1px solid var(--line);border-radius:4px;padding:3px 8px;text-decoration:none;font-size:12px}.missing{color:#9b4f49}.card-note{font-size:12px;color:var(--muted);margin:0;padding:0 14px 13px}.file-list{margin-top:14px;background:var(--panel);padding:10px 14px;border-radius:6px}.file-list ul{columns:2;padding-left:20px}.file-list li{break-inside:avoid;font-size:12px;margin:4px 0}.hidden{display:none!important}\nfooter{border-top:1px solid var(--line);color:var(--muted);font-size:12px;padding:18px 0}@media(max-width:700px){.page{padding:16px}.figure-grid{grid-template-columns:1fr}.figure-card img{height:auto}.file-list ul{columns:1}.module-heading{display:block}}\n")

js <- htmltools::HTML("\nconst search=document.getElementById('figure-search');\nsearch.addEventListener('input',()=>{const q=search.value.trim().toLowerCase();document.querySelectorAll('.figure-card').forEach(card=>card.classList.toggle('hidden',q && !card.dataset.search.includes(q)));document.querySelectorAll('.module-section').forEach(section=>{const cards=[...section.querySelectorAll('.figure-card')];section.classList.toggle('hidden',cards.length>0 && cards.every(x=>x.classList.contains('hidden')));});});\n")

page <- htmltools::tags$html(
  lang = "zh-CN",
  htmltools::tags$head(
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    htmltools::tags$title(paste(project_label, assay, "结果报告")),
    htmltools::tags$style(css)
  ),
  htmltools::tags$body(
    htmltools::tags$main(
      class = "page",
      htmltools::tags$header(
        class = "hero",
        htmltools::tags$h1(paste(project_label, assay, "结果报告")),
        htmltools::tags$p(class = "subtitle", "离线结果索引：按模块浏览图件，并下载同名 PDF、PNG 与结果表。"),
        htmltools::tags$div(
          class = "summary",
          htmltools::tags$span(class = "badge", paste("对比:", contrast)),
          htmltools::tags$span(class = "badge", paste("图件:", length(png_paths))),
          htmltools::tags$span(class = "badge", paste("PDF:", length(pdf_paths))),
          htmltools::tags$span(class = "badge", paste("表格:", length(table_paths))),
          htmltools::tags$span(class = "badge", "图件参数见各模块导出文件")
        )
      ),
      status_note,
      htmltools::tags$aside(
        class = "report-note",
        htmltools::tags$h2("解释边界 / Interpretation boundaries"),
        htmltools::tags$ul(lapply(scientific_notes, htmltools::tags$li))
      ),
      htmltools::tags$div(class = "toolbar", htmltools::tags$input(id = "figure-search", type = "search", placeholder = "按文件名、模块或关键词筛选图件…", `aria-label` = "筛选图件")),
      htmltools::tags$nav(lapply(modules, function(x) htmltools::tags$a(href = paste0("#", safe_id(x)), module_title(x)))),
      lapply(modules, make_module),
      htmltools::tags$footer(
        "本报告由固定脚本从 result/ 目录生成，不调用外部资源，不自动生成生物学结论。生成时间：",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      )
    ),
    htmltools::tags$script(js)
  )
)

htmltools::save_html(page, file = output_file, background = "white")
message("Result report written: ", output_file)
