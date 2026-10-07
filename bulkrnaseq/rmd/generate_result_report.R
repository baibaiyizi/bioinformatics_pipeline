#!/usr/bin/env Rscript

required_packages <- c("htmltools", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing report dependencies: ", paste(missing_packages, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else ""
script_path <- gsub("~+~", " ", script_path, fixed=TRUE)
default_root <- if (nzchar(script_path)) dirname(dirname(normalizePath(script_path, mustWork = FALSE))) else getwd()
root_dir <- normalizePath(Sys.getenv("RNASEQ_ROOT_DIR", unset = default_root), mustWork = TRUE)
result_dir <- file.path(root_dir, "result")
output_file <- file.path(result_dir, "index.html")

if (!dir.exists(result_dir)) {
  stop("Result directory does not exist: ", result_dir)
}

relative_to_result <- function(path) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root_norm, "/")
  if (!startsWith(path_norm, prefix)) {
    stop("Report input is outside result/: ", path)
  }
  substring(path_norm, nchar(prefix) + 1L)
}

natural_sort <- function(x) stringr::str_sort(x, numeric = TRUE)
file_url <- function(x) utils::URLencode(gsub("\\\\", "/", x), reserved = FALSE)
module_of <- function(path) strsplit(path, "/", fixed = TRUE)[[1]][[1]]
safe_id <- function(x) paste0("section-", gsub("[^A-Za-z0-9_-]+", "-", x))

all_files <- list.files(result_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
all_files <- all_files[normalizePath(all_files, winslash = "/", mustWork = FALSE) != normalizePath(output_file, winslash = "/", mustWork = FALSE)]
relative_files <- vapply(all_files, relative_to_result, character(1))

png_paths <- natural_sort(relative_files[grepl("\\.png$", relative_files, ignore.case = TRUE)])
table_paths <- natural_sort(relative_files[grepl("\\.(csv|tsv)$", relative_files, ignore.case = TRUE)])
other_report_paths <- natural_sort(relative_files[grepl("\\.html?$", relative_files, ignore.case = TRUE)])


read_context <- function() {
  path <- file.path(result_dir, "01_setup", "tables", "01_analysis_context.csv")
  if (!file.exists(path)) return(setNames(character(0), character(0)))
  tbl <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("metric", "value") %in% names(tbl))) return(setNames(character(0), character(0)))
  stats::setNames(as.character(tbl$value), as.character(tbl$metric))
}

context <- read_context()
context_value <- function(name, default = "Not available") {
  if (!name %in% names(context)) return(default)
  value <- unname(context[[name]])
  if (is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)) default else value
}

contrast_path <- file.path(result_dir, "05_de", "tables", "5_0_contrast_manifest.csv")
analysis_contrast <- context_value("analysis_contrast", "")
if (file.exists(contrast_path)) {
  contrast_tbl <- utils::read.csv(contrast_path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("analysis_contrast" %in% names(contrast_tbl) && nrow(contrast_tbl) > 0 && nzchar(contrast_tbl$analysis_contrast[[1]])) {
    analysis_contrast <- contrast_tbl$analysis_contrast[[1]]
  }
}
if (!nzchar(analysis_contrast)) {
  group1 <- Sys.getenv("RNASEQ_GROUP1"); group2 <- Sys.getenv("RNASEQ_GROUP2")
  if (!nzchar(group1) || !nzchar(group2)) stop("报告需要显式 RNASEQ_GROUP1/2")
  analysis_contrast <- paste(group2, group1, sep=" - ")
}

contrast_groups <- strsplit(analysis_contrast, " - ", fixed = TRUE)[[1]]
if (length(contrast_groups) != 2 || any(!nzchar(contrast_groups))) {
  stop("Invalid analysis contrast: ", analysis_contrast)
}
group2_label <- contrast_groups[[1]]
group1_label <- contrast_groups[[2]]
project_label <- basename(root_dir)
higher_group2 <- paste0(group2_label, " 较高")
higher_group1 <- paste0(group1_label, " 较高")

module_titles <- c(
  "00_audit" = "审计与质量记录 / Audit and quality records",
  "01_setup" = "分析环境 / Analysis context",
  "02_selection" = "样本选择 / Sample selection",
  "03_qc" = "质量控制 / Quality control",
  "04_pca" = "主成分分析 / Principal component analysis",
  "05_de" = "差异表达 / Differential expression",
  "06_enrichment" = "通路富集 / Pathway enrichment",
  "07_go" = "GO 富集 / Gene Ontology enrichment",
  "08_gsea" = "基因集富集 / Gene set enrichment analysis",
  "09_gsva" = "通路活性评分 / GSVA",
  "10_tf_regulator" = "转录因子活性 / TF regulator activity",
  "11_ppi" = "蛋白互作网络 / Protein-protein interaction network",
  "12_wgcna" = "共表达网络 / WGCNA",
  "13_apa" = "可变多聚腺苷酸化 / Alternative polyadenylation",
  "14_der" = "差异区域 / Differentially expressed regions",
  "15_splicing" = "可变剪接与转录本 / Alternative splicing and isoforms",
  "16_immune" = "免疫信号推断 / Immune signal inference",
  "17_pathway_activity" = "差异通路活性 / Differential pathway activity",
  "18_circle" = "富集圈图 / Enrichment circle plots"
)

module_title <- function(module_id) {
  title <- unname(module_titles[module_id])
  if (is.null(title) || length(title) == 0 || is.na(title)) module_id else title
}

directory_title <- function(relative_dir, module_id) {
  suffix <- sub(paste0("^", module_id, "/?"), "", relative_dir)
  suffix <- sub("(^|/)plots$", "", suffix)
  if (!nzchar(suffix)) return("图件 / Figures")
  labels <- c(
    "01_all" = "全部显著差异基因 / All significant DE genes",
    "02_up" = paste0(higher_group2, " / Higher in ", group2_label),
    "03_down" = paste0(higher_group1, " / Higher in ", group1_label),
    "human_ortholog" = "人源同源探索性结果 / Human ortholog (exploratory)",
    "01_event_level" = "事件层 / Event level",
    "02_gene_level" = "基因层 / Gene level",
    "03_visual_validation" = "可视验证 / Visual validation",
    "04_isoform_switch" = "转录本切换 / Isoform switching",
    "05_integration" = "多证据整合 / Multi-evidence integration",
    "06_sashimi_inputs" = "Sashimi 输入 / Sashimi inputs"
  )
  value <- unname(labels[suffix])
  if (is.null(value) || length(value) == 0 || is.na(value)) gsub("_", " ", suffix) else value
}

pretty_figure_name <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  stem <- sub("^[0-9]+_[0-9]+[A-Za-z]?_", "", stem)
  stem <- gsub("_", " ", stem, fixed = TRUE)
  replacements <- c(
    "barplot" = "bar plot",
    "dotplot" = "dot plot",
    "cnetplot" = "category-gene network",
    "emapplot" = "enrichment map",
    "treeplot" = "tree plot",
    "heatmap" = "heatmap",
    "boxplots" = "box plots"
  )
  for (from in names(replacements)) stem <- gsub(from, replacements[[from]], stem, ignore.case = TRUE)
  stem
}

note_rule <- function(id, pattern, method, reading, limitation) {
  list(id = id, pattern = pattern, method = method, reading = reading, limitation = limitation)
}

note_rules <- list(
  note_rule("sample-count", "selected_group_counts", "按最终纳入清单汇总每个分组的样本数。", "比较组间样本数是否平衡，并核对总数是否与样本清单一致。", "样本数平衡不等同于样本质量充分，仍需结合 QC 与样本级图件。"),
  note_rule("library-complexity", "library_size|library_complexity|top100_gene_fraction", "由原始计数计算文库总量、达到 CPM 阈值的基因数或 Top 100 基因占比。", "关注极端低复杂度、高集中度或与同组样本明显偏离的样本。", "这些指标用于技术 QC，不能单独证明或否定组间生物差异。"),
  note_rule("pca-scree", "pca_scree", "对归一化表达矩阵进行主成分分析，并汇总各主成分解释的方差比例。", "前几个主成分解释比例越高，说明主要样本变异可由较少维度概括。", "解释方差不提供组间显著性检验。"),
  note_rule("pca", "pca(\\.|_)|_pca_with_|PSI_sample_pca", "对相应的归一化表达或 PSI 矩阵中心化后进行主成分分析。", paste0("先看组内离散，再看 ", group2_label, " 与 ", group1_label, " 是否沿主要主成分分离，并检查是否由单个离群样本驱动。"), "PCA 是无监督描述，不等同于组间差异显著。"),
  note_rule("expression-boxplot", "expression_distribution_boxplot|rle_boxplot", "按样本汇总归一化表达分布；RLE 图以基因相对中位表达衡量样本偏移。", "各样本分布应具有可比中心与离散度，系统偏移提示归一化或样本质量问题。", "箱线图压缩了基因层信息，应与相关性、PCA 和原始 QC 联合判断。"),
  note_rule("sample-similarity", "sample_correlation_heatmap|sample_distance_heatmap|sample_dendrogram", "由样本表达向量计算 Pearson 相关、样本距离或层次聚类。", "关注组内相似性、组间结构和孤立样本，三种表示应相互印证。", "相似性结构受高变基因、变换方式和批次影响，不能替代设计矩阵中的统计建模。"),
  note_rule("fastp-qc", "fastp_", "从 fastp JSON 汇总过滤前后 reads、Q20/Q30 和 GC 等测序质量指标。", "检查各样本保留率和质量率是否一致，并识别异常低值。", "fastp 指标反映测序与过滤质量，不反映表达差异的统计显著性。"),
  note_rule("de-summary", "DE_summary", paste0("按既定 p-value 与 |log2FC| 阈值汇总 ", higher_group2, "、", higher_group1, "和不显著基因。"), "用于核对方向数量与阈值口径，不依据柱高推断机制。", "当前默认显示 nominal p-value；正式主张还应查看 padj/FDR。"),
  note_rule("volcano", "volcano", "横轴为效应量 log2FC，纵轴为当前图形显著性指标的 -log10 值，并按阈值标记方向。", paste0("同时关注效应量和显著性；右侧表示 ", higher_group2, "，左侧表示 ", higher_group1, "。"), "高显著性不等同于大效应；nominal p-value 图必须结合 padj/FDR 表解释。"),
  note_rule("ma-plot", "_MA\\.png$|ma_plot", "以平均表达量为横轴、log2FC 为纵轴展示基因效应。", "检查低表达区域的离散度、整体中心是否接近 0，以及显著基因的方向分布。", "低计数基因的效应更不稳定，不能只按 log2FC 极值排序。"),
  note_rule("de-concordance", "log2fc_concordance", "将主分析与全样本敏感性分析中匹配基因的 log2FC 对齐并比较。", "点越接近对角线表示效应方向和幅度越一致；同时观察相关系数和方向翻转。", "相关性高不保证每个基因均可复现，也不替代独立样本验证。"),
  note_rule("leading-edge-heatmap", "leading_edge_heatmap", "从显著 GSEA 条目的 leading-edge 成员中选取最多 8 条通路和全图最多 40 个代表基因；每条通路先保留高绝对效应基因，再按跨通路覆盖数和绝对 log2FC 补足。", "单元格表示成员基因的 log2FC，非成员位置为 0；0 使用中性色，通路名为便于显示而缩写。", "该图是展示性筛选，不改变 GSEA 结果；完整通路名、leading-edge 成员和统计量以相邻 CSV 表为准。"),
  note_rule("ranked-heatmap", "top_DE_heatmap|top_.*heatmap|PSI_heatmap|isoform_fraction_heatmap", "选取预先定义的高优先级特征，对样本级矩阵进行行标准化或相应比例变换后聚类展示。", paste0("看样本聚类、组内一致性以及特征在 ", group2_label, "/", group1_label, " 间的相对模式，不比较不同基因的绝对颜色。"), "Top 特征是选择后的展示，不能用于无偏估计整体分离能力。"),
  note_rule("ora-dot-bar", "ORA_.*(barplot|dotplot)|GO_.*(barplot|dotplot)|DAS_(GO|KEGG)_dotplot|directional_dotplot|combined_(barplot|dotplot)|node_enrichment_dotplot", "对预先筛选的基因集进行过度富集分析，并按图形显著性、GeneRatio/Rich factor 或基因数排序。", paste0("优先同时查看富集比例、基因数和显著性，并确认方向目录是全部、", higher_group2, "或", higher_group1, "。"), "ORA 依赖筛选阈值和背景基因集；正式解释应以对应库内 FDR 和源表为准。"),
  note_rule("gsea", "GSEA_.*(dotplot|ridgeplot|emapplot)|Mouse_M[25]_GSEA", "使用全基因排序运行 GSEA，并以 NES、集合大小及显著性展示富集方向；Reactome 方向 dot plot 最多显示 12 条结果。", "NES 为正和负代表排序列表两端的富集；图中数据库前缀和长标签仅为显示而缩写。", "GSEA 与 all-DE ORA 的 estimand 不同；完整名称、NES、leading edge 和 FDR 以相邻结果表为准。"),
  note_rule("chord", "Chord", "将显著通路与其成员基因连接，基因侧颜色编码 log2FC，通路侧表示富集条目；沿用分析脚本既定的通路和 Top 基因展示上限。", "用于识别多个通路共享的基因以及这些基因的效应方向；右侧图例名称最多显示 38 个字符。", "标签缩写和连线密度只影响展示，不改变成员关系；完整名称和成员以相邻表格为准，也不表示调控因果关系。"),
  note_rule("category-network", "cnetplot|emapplot|Cluster|Treeplot", "基于通路-基因归属或通路间基因重叠构建网络、聚类树或相似性结构；emapplot 和 treeplot 最多显示 10 条术语。", "关注共享基因形成的条目簇和冗余通路；标签最多显示 38 个字符，完整名称以相邻表格为准。", "网络边表示集合重叠或相似度，不表示分子相互作用或因果方向；显示上限不改变富集计算。"),
  note_rule("upset", "UpSet|upset", "把每个基因集编码为集合成员关系，并统计单集合及交集大小。", "顶部柱表示特定组合的交集数，左侧柱表示各集合总量，矩阵点定义组合。", "交集大小受输入阈值和集合定义影响；空白圆点不表示统计上的负相关。"),
  note_rule("gsva-activity", "GSVA|pathway_activity_(heatmap|barplot|boxplots)|activity_vs_GSEA", paste0("从样本表达矩阵计算通路级活性分数，再按 ", analysis_contrast, " 比较；一致性图将 activity logFC 与同来源 GSEA NES 对齐。"), "热图看样本级模式，条形/箱线图看组间方向，一致性散点看两种方法的方向是否相符。", "通路分数是基因集层代理量；人源同源结果仅作探索性补充，不能替代 mouse-native 主证据。"),
  note_rule("tf-activity", "TF_activity|TF_pathway_activity", paste0("基于 DoRothEA/decoupleR 调控网络从靶基因表达推断 TF 活性，并比较 ", group2_label, " 与 ", group1_label, "。"), "活性效应与 TF 自身表达不是同一量；优先查看方向、显著性和靶基因覆盖。", "推断依赖先验 regulon，不能直接证明 TF 的实验性激活或抑制。"),
  note_rule("interaction-network", "target_network|PPI_network|subnetwork_panels|adjacency_heatmap", "按已知 TF-target 或 STRING 互作边构建网络，并用节点属性编码方向、模块或中心性。", "看连接结构、边证据和候选节点是否在多个指标下稳定，而不是只看图中位置。", "网络拓扑来自数据库与阈值；中心节点不自动等同于生物学因果驱动因子。"),
  note_rule("ppi-topology", "PPI_(centrality|hub|edge_score|degree|bridge|subnetwork)", "在 STRING 子图上计算 degree、strength、betweenness、closeness、eigenvector 或 PageRank，并汇总网络稳健性。", "比较不同中心性指标的一致性、阈值敏感性和跨子网络连接比例。", "中心性依赖当前节点集合和 STRING 覆盖，不能外推到完整细胞网络。"),
  note_rule("event-summary", "event_summary|event_barplot|direction_barplot|overlap_barplot|summary_barplot", "按既定显著性和效应阈值汇总事件、方向或方法间重叠数量。", "先核对总数和方向定义，再看重叠比例是否由某一大集合主导。", "计数反映当前阈值下的发现数，不等同于效应强弱或可重复性。"),
  note_rule("splicing-effect", "deltaPSI|DSFA|JC_JCEC|DE_DAS|exon_usage|joint_effect_scatter|integrated_as_evidence", "整合 rMATS delta PSI、事件频率、exon usage 或基因表达效应，按相应显著性指标排序或比较。", "关注效应方向、方法一致性、事件类型分层以及是否由少数极端点驱动。", "不同剪接方法的统计单位不同；事件层与基因层证据不能直接相加解释为独立重复。"),
  note_rule("isoform-switch", "isoform_switch|isoform_consequence", "从转录本丰度和 isoform fraction 计算切换效应 dIF，并汇总显著性与预测功能后果。", "同时查看 dIF、q-value、基因内转录本方向及功能后果类别。", "预测后果依赖注释完整性，属于候选证据，需转录本或蛋白层验证。"),
  note_rule("immune", "immune_", "使用 ssGSEA 和适用时的 mMCPcounter 基因签名计算样本级免疫相关分数，并比较组间差异与方法一致性。", "把结果视为免疫相关表达信号；看组内分布、方法方向一致性和多重校正结果。", "早期胚胎 bulk RNA-seq 中这些分数不是细胞比例的直接测量，不能解释为真实免疫细胞浸润。"),
  note_rule("enrichment-circle", "enrichment_circle", "将富集条目的显著性、基因数、上下调构成和富集效应编码到同心轨道。", "由外向内依次核对条目标识、显著性、基因构成和效应方向，并与详细表对应。", "圈图是同一富集结果的二次展示，不构成独立证据。"),
  note_rule("generic-heatmap", "heatmap|correlation", "将相应矩阵按脚本指定的变换、排序或聚类绘制为热图。", "先阅读图例量纲，再判断行列模式、组内一致性和异常值。", "颜色是相对编码；不同热图间不能在未核对量纲时直接比较颜色深浅。"),
  note_rule("generic-box-violin", "boxplot|violin", "按组展示样本或事件分布，箱体表示四分位范围，violin 宽度表示核密度。", "比较中心、离散度、样本重叠和极端值，不只比较均值或中位数。", "样本量较小时分布形状不稳定，正式推断应以对应统计模型为准。"),
  note_rule("generic-scatter", "scatter|concordance", "将两个连续统计量对齐后绘制散点，并按需要标记方向、显著性或候选对象。", "看整体趋势、离群点、象限方向和是否由少量点驱动。", "相关或一致性不等同于因果，也不替代原始模型的显著性检验。"),
  note_rule("generic-ranking", "ranking|barplot|distribution", "按源表中的效应量、显著性、计数或网络指标排序并展示。", "先核对排序变量和单位，再比较方向与相对大小。", "Top 排名是展示性筛选，不应被解释为未经验证的优先级结论。")
)

classify_figure <- function(path) {
  for (rule in note_rules) {
    if (grepl(rule$pattern, path, ignore.case = TRUE, perl = TRUE)) return(rule)
  }
  NULL
}

figure_rules <- lapply(png_paths, classify_figure)
unclassified <- png_paths[vapply(figure_rules, is.null, logical(1))]
if (length(unclassified) > 0) {
  stop("Unclassified report figures:\n", paste(unclassified, collapse = "\n"))
}

figure_modules <- vapply(png_paths, module_of, character(1))
table_modules <- if (length(table_paths) > 0) vapply(table_paths, module_of, character(1)) else character(0)
report_modules <- if (length(other_report_paths) > 0) vapply(other_report_paths, module_of, character(1)) else character(0)
modules <- natural_sort(unique(c(figure_modules, table_modules, report_modules)))

make_link <- function(path, label = basename(path), class = NULL) {
  htmltools::tags$a(href = file_url(path), target = "_blank", rel = "noopener", class = class, label)
}

context_badges <- function(path) {
  badges <- character(0)
  if (grepl("/human_ortholog/", path, fixed = TRUE)) {
    badges <- c(badges, "Human ortholog · exploratory")
  } else if (module_of(path) %in% c("06_enrichment", "07_go", "08_gsea", "09_gsva", "17_pathway_activity", "18_circle")) {
    badges <- c(badges, "Mouse-native")
  }
  if (grepl("/01_all/", path, fixed = TRUE)) badges <- c(badges, "All significant DE genes")
  if (grepl("/02_up/", path, fixed = TRUE)) badges <- c(badges, paste("Higher in", group2_label))
  if (grepl("/03_down/", path, fixed = TRUE)) badges <- c(badges, paste("Higher in", group1_label))
  badges
}

make_figure_card <- function(path, rule) {
  pdf_path <- sub("\\.png$", ".pdf", path, ignore.case = TRUE)
  pdf_exists <- file.exists(file.path(result_dir, pdf_path))
  title <- pretty_figure_name(path)
  badges <- context_badges(path)
  search_text <- tolower(paste(path, title, module_title(module_of(path)), rule$id, paste(badges, collapse = " ")))

  htmltools::tags$article(
    class = "figure-card",
    `data-search` = search_text,
    `data-module` = module_of(path),
    htmltools::tags$div(
      class = "figure-card-head",
      htmltools::tags$h4(title),
      htmltools::tags$code(path)
    ),
    if (length(badges) > 0) htmltools::tags$div(
      class = "badges",
      lapply(badges, function(x) htmltools::tags$span(class = "badge", x))
    ),
    htmltools::tags$button(
      type = "button",
      class = "figure-open",
      `data-src` = file_url(path),
      `data-title` = title,
      `aria-label` = paste("查看大图", title),
      htmltools::tags$img(src = file_url(path), alt = title, loading = "lazy", decoding = "async")
    ),
    htmltools::tags$div(
      class = "file-actions",
      make_link(path, "PNG"),
      if (pdf_exists) make_link(pdf_path, "PDF")
    ),
    htmltools::tags$details(
      class = "figure-notes",
      htmltools::tags$summary("计算方法与通用读图"),
      htmltools::tags$dl(
        htmltools::tags$dt("计算方法"), htmltools::tags$dd(rule$method),
        htmltools::tags$dt("通用读图"), htmltools::tags$dd(rule$reading),
        htmltools::tags$dt("解释边界"), htmltools::tags$dd(rule$limitation)
      )
    )
  )
}

make_table_list <- function(paths) {
  if (length(paths) == 0) return(NULL)
  htmltools::tags$details(
    class = "module-files",
    htmltools::tags$summary(sprintf("结果表（%d）", length(paths))),
    htmltools::tags$ul(lapply(paths, function(path) htmltools::tags$li(make_link(path, path))))
  )
}

make_other_report_list <- function(paths) {
  if (length(paths) == 0) return(NULL)
  htmltools::tags$details(
    class = "module-files",
    htmltools::tags$summary(sprintf("其他 HTML 报告（%d）", length(paths))),
    htmltools::tags$ul(lapply(paths, function(path) htmltools::tags$li(make_link(path, path))))
  )
}

make_module_section <- function(module_id) {
  module_figures <- png_paths[figure_modules == module_id]
  module_tables <- table_paths[table_modules == module_id]
  module_reports <- other_report_paths[report_modules == module_id]
  figure_dirs <- natural_sort(unique(dirname(module_figures)))

  groups <- lapply(figure_dirs, function(one_dir) {
    these_paths <- module_figures[dirname(module_figures) == one_dir]
    these_rules <- figure_rules[match(these_paths, png_paths)]
    htmltools::tags$section(
      class = "figure-group",
      htmltools::tags$h3(directory_title(one_dir, module_id)),
      htmltools::tags$p(class = "directory-path", one_dir),
      htmltools::tags$div(
        class = "figure-grid",
        Map(make_figure_card, these_paths, these_rules)
      )
    )
  })

  htmltools::tags$section(
    id = safe_id(module_id),
    class = "module-section",
    htmltools::tags$div(
      class = "module-heading",
      htmltools::tags$div(
        htmltools::tags$h2(module_title(module_id)),
        htmltools::tags$code(module_id)
      ),
      htmltools::tags$span(class = "module-count", sprintf("%d figures", length(module_figures)))
    ),
    if (length(module_figures) == 0) htmltools::tags$p(class = "empty-note", "当前目录没有 PNG 图件。"),
    groups,
    make_table_list(module_tables),
    make_other_report_list(module_reports)
  )
}

css <- "
:root { --ink:#24323a; --muted:#66757f; --line:#d9e0e4; --soft:#f4f7f8; --accent:#315f70; --accent-soft:#e8f0f3; --paper:#ffffff; }
* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body { margin:0; color:var(--ink); background:var(--soft); font-family:'Liberation Sans','Noto Sans CJK SC','Microsoft YaHei',Arial,sans-serif; line-height:1.55; }
a { color:var(--accent); text-decoration:none; }
a:hover { text-decoration:underline; }
.layout { display:grid; grid-template-columns:270px minmax(0,1fr); min-height:100vh; }
.sidebar { position:sticky; top:0; height:100vh; overflow:auto; padding:24px 18px; background:#263840; color:#f7fafb; }
.sidebar h1 { margin:0 0 4px; font-size:20px; font-weight:600; }
.sidebar .subtitle { margin:0 0 20px; color:#c8d4d9; font-size:13px; }
.sidebar input { width:100%; border:1px solid #69808a; border-radius:4px; padding:9px 10px; background:#fff; color:#1f2b31; }
.sidebar nav { margin-top:20px; }
.sidebar nav a { display:block; padding:7px 8px; color:#e8eef1; border-left:2px solid transparent; font-size:13px; }
.sidebar nav a:hover { border-left-color:#9ec4d1; background:#314851; text-decoration:none; }
.sidebar .search-status { margin-top:8px; color:#c8d4d9; font-size:12px; }
main { width:min(1480px,100%); padding:32px 38px 64px; }
.report-header,.module-section { background:var(--paper); border:1px solid var(--line); border-radius:6px; }
.report-header { padding:28px 30px; margin-bottom:24px; }
.report-header h1 { margin:0 0 8px; font-size:26px; font-weight:600; }
.report-header .lede { margin:0; color:var(--muted); }
.summary-grid { display:grid; grid-template-columns:repeat(4,minmax(120px,1fr)); gap:12px; margin:22px 0; }
.summary-item { padding:14px; background:var(--soft); border-left:3px solid var(--accent); }
.summary-item strong { display:block; font-size:20px; font-weight:600; }
.summary-item span { color:var(--muted); font-size:12px; }
.report-note { margin:12px 0 0; padding:12px 14px; border:1px solid var(--line); background:#fbfcfc; }
.report-note ul { margin:6px 0 0; padding-left:20px; }
.module-section { padding:24px; margin-bottom:22px; scroll-margin-top:16px; }
.module-heading { display:flex; justify-content:space-between; gap:16px; align-items:flex-start; padding-bottom:14px; border-bottom:1px solid var(--line); }
.module-heading h2 { margin:0 0 3px; font-size:21px; font-weight:600; }
.module-count { color:var(--muted); white-space:nowrap; }
.figure-group h3 { margin:24px 0 0; font-size:17px; font-weight:600; }
.directory-path { margin:2px 0 12px; color:var(--muted); font-family:monospace; font-size:12px; }
.figure-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; }
.figure-card { min-width:0; border:1px solid var(--line); border-radius:5px; background:#fff; overflow:hidden; }
.figure-card[hidden] { display:none; }
.figure-card-head { padding:14px 14px 8px; }
.figure-card h4 { margin:0 0 3px; font-size:15px; font-weight:600; }
.figure-card code { display:block; overflow-wrap:anywhere; color:var(--muted); font-size:11px; }
.badges { display:flex; flex-wrap:wrap; gap:6px; padding:0 14px 10px; }
.badge { padding:2px 7px; border:1px solid #bdd0d7; border-radius:999px; background:var(--accent-soft); color:#31515d; font-size:11px; }
.figure-open { width:100%; min-height:230px; padding:8px; border:0; border-top:1px solid var(--line); border-bottom:1px solid var(--line); background:#fff; cursor:zoom-in; }
.figure-open img { width:100%; height:300px; object-fit:contain; display:block; }
.file-actions { display:flex; gap:8px; padding:10px 14px; }
.file-actions a { border:1px solid var(--line); border-radius:3px; padding:4px 9px; font-size:12px; }
.figure-notes { border-top:1px solid var(--line); padding:10px 14px 13px; background:#fbfcfc; }
.figure-notes summary,.module-files summary { cursor:pointer; color:var(--accent); font-weight:600; }
.figure-notes dl { display:grid; grid-template-columns:70px 1fr; gap:7px 12px; margin:12px 0 0; font-size:13px; }
.figure-notes dt { font-weight:600; }
.figure-notes dd { margin:0; color:#3d4c54; }
.module-files { margin-top:18px; padding-top:14px; border-top:1px solid var(--line); }
.module-files ul { columns:2; margin:10px 0 0; padding-left:20px; font-size:12px; }
.module-files li { break-inside:avoid; margin-bottom:4px; overflow-wrap:anywhere; }
.empty-note { color:var(--muted); }
dialog { width:min(94vw,1500px); max-height:94vh; padding:0; border:1px solid var(--line); border-radius:6px; }
dialog::backdrop { background:rgba(18,28,33,.72); }
.dialog-head { display:flex; justify-content:space-between; align-items:center; gap:12px; padding:10px 14px; border-bottom:1px solid var(--line); }
.dialog-head strong { font-size:14px; }
.dialog-head button { border:1px solid var(--line); background:#fff; border-radius:3px; padding:5px 10px; cursor:pointer; }
#dialog-image { display:block; max-width:100%; max-height:calc(94vh - 52px); margin:auto; object-fit:contain; }
@media (max-width:900px) { .layout{display:block}.sidebar{position:relative;height:auto}.sidebar nav{columns:2}main{padding:20px}.summary-grid{grid-template-columns:repeat(2,1fr)} }
@media (max-width:560px) { main{padding:12px}.report-header,.module-section{padding:16px}.figure-grid{grid-template-columns:1fr}.module-files ul{columns:1}.summary-grid{grid-template-columns:1fr}.figure-open img{height:240px} }
@media print { .sidebar,.file-actions,dialog{display:none!important}.layout{display:block}main{padding:0;width:100%}.module-section,.report-header{border:0;break-inside:avoid}.figure-card{break-inside:avoid}.figure-notes{display:block} }
"

js <- "
const cards = [...document.querySelectorAll('.figure-card')];
const search = document.getElementById('figure-search');
const status = document.getElementById('search-status');
function applySearch() {
  const query = search.value.trim().toLowerCase();
  let visible = 0;
  cards.forEach(card => {
    const show = !query || card.dataset.search.includes(query);
    card.hidden = !show;
    if (show) visible += 1;
  });
  status.textContent = `${visible} / ${cards.length} figures`;
}
search.addEventListener('input', applySearch);
applySearch();
const dialog = document.getElementById('figure-dialog');
const dialogImage = document.getElementById('dialog-image');
const dialogTitle = document.getElementById('dialog-title');
document.querySelectorAll('.figure-open').forEach(button => {
  button.addEventListener('click', () => {
    dialogImage.src = button.dataset.src;
    dialogImage.alt = button.dataset.title;
    dialogTitle.textContent = button.dataset.title;
    dialog.showModal();
  });
});
document.getElementById('dialog-close').addEventListener('click', () => dialog.close());
dialog.addEventListener('click', event => { if (event.target === dialog) dialog.close(); });
"

navigation <- htmltools::tags$nav(lapply(modules, function(module_id) {
  htmltools::tags$a(href = paste0("#", safe_id(module_id)), paste(module_id, sub(" /.*$", "", module_title(module_id))))
}))

page <- htmltools::tags$html(
  lang = "zh-CN",
  htmltools::tags$head(
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    htmltools::tags$title(paste(project_label, "RNA-seq 结果报告")),
    htmltools::tags$style(htmltools::HTML(css))
  ),
  htmltools::tags$body(
    htmltools::tags$div(
      class = "layout",
      htmltools::tags$aside(
        class = "sidebar",
        htmltools::tags$h1(paste(project_label, "RNA-seq")),
        htmltools::tags$p(class = "subtitle", "结果图与文件索引"),
        htmltools::tags$label(`for` = "figure-search", "搜索图名、路径或方法"),
        htmltools::tags$input(id = "figure-search", type = "search", placeholder = "例如 volcano、GO、splicing"),
        htmltools::tags$div(id = "search-status", class = "search-status"),
        navigation
      ),
      htmltools::tags$main(
        htmltools::tags$header(
          class = "report-header",
          htmltools::tags$h1(paste(project_label, "RNA-seq 分析结果")),
          htmltools::tags$p(class = "lede", "按 result/ 实际文件生成的离线图集与结果索引。页面不自动生成机制性或因果性结论。"),
          htmltools::tags$div(
            class = "summary-grid",
            htmltools::tags$div(class = "summary-item", htmltools::tags$strong(length(png_paths)), htmltools::tags$span("PNG figures")),
            htmltools::tags$div(class = "summary-item", htmltools::tags$strong(length(table_paths)), htmltools::tags$span("CSV/TSV tables")),
            htmltools::tags$div(class = "summary-item", htmltools::tags$strong(context_value("selected_samples_n")), htmltools::tags$span("Selected samples")),
            htmltools::tags$div(class = "summary-item", htmltools::tags$strong(analysis_contrast), htmltools::tags$span("Analysis contrast"))
          ),
          htmltools::tags$div(
            class = "report-note",
            htmltools::tags$strong("解释口径"),
            htmltools::tags$ul(
              htmltools::tags$li(paste0("log2FC = ", analysis_contrast, "；正值表示 ", higher_group2, "，负值表示 ", higher_group1, "。")),
              htmltools::tags$li(paste0("图形默认显著性列为 ", context_value("de_significance_column", "pvalue"), "；DE 阈值为 p < ", context_value("de_pvalue_cutoff", "0.05"), " 且 |log2FC| ≥ ", context_value("de_log2fc_cutoff", "1"), "。正式主张应同步检查 padj/FDR。")),
              htmltools::tags$li("Mouse-native 结果构成主证据；human_ortholog/ 仅为跨物种探索性补充。"),
              htmltools::tags$li("每张图的说明描述计算来源、通用读图方式和解释边界；当前数据的具体数值以图及配套表为准。")
            )
          ),
          htmltools::tags$p(class = "lede", paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
        ),
        lapply(modules, make_module_section)
      )
    ),
    htmltools::tags$dialog(
      id = "figure-dialog",
      htmltools::tags$div(class = "dialog-head", htmltools::tags$strong(id = "dialog-title"), htmltools::tags$button(id = "dialog-close", type = "button", "关闭")),
      htmltools::tags$img(id = "dialog-image", alt = "")
    ),
    htmltools::tags$script(htmltools::HTML(js))
  )
)

htmltools::save_html(page, file = output_file)
message("Result report written: ", output_file)
message("Report inventory: figures=", length(png_paths), ", tables=", length(table_paths), ", modules=", length(modules))
