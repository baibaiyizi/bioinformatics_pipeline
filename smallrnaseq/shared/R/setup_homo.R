# 人类设计适配；serum的匹配分装与过程QC队列明确隔离。
source(file.path(Sys.getenv("PANDORA_SHARED_ROOT"),"R","common.R"),local=TRUE)

knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 11, fig.height = 7, dpi = 300)

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(knitr)
})

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

hellinger_transform <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  if (any(totals <= 0)) stop("Hellinger 变换要求每个样本至少有一个非零特征。")
  sqrt(x / totals)
}

compute_pair_agreement_metrics <- function(
  abundance_df, sample_a_col, sample_b_col, feature_level,
  sample_a, sample_b, thresholds = c(0.5, 1, 5), pseudocount = 0.5
) {
  sample_a_cpm <- as.numeric(abundance_df[[sample_a_col]])
  sample_b_cpm <- as.numeric(abundance_df[[sample_b_col]])
  purrr::map_dfr(thresholds, function(threshold) {
    a_detected <- sample_a_cpm >= threshold
    b_detected <- sample_b_cpm >= threshold
    union <- a_detected | b_detected
    shared <- a_detected & b_detected
    shared_logfc <- log2((sample_b_cpm[shared] + pseudocount) / (sample_a_cpm[shared] + pseudocount))
    tibble(
      FeatureLevel = feature_level,
      CPMThreshold = threshold,
      SampleA = sample_a,
      SampleB = sample_b,
      SampleADetected = sum(a_detected),
      SampleBDetected = sum(b_detected),
      SampleAOnly = sum(a_detected & !b_detected),
      SharedDetected = sum(shared),
      SampleBOnly = sum(!a_detected & b_detected),
      UnionDetected = sum(union),
      Jaccard = if_else(sum(union) > 0, sum(shared) / sum(union), NA_real_),
      SharedSpearman = if (sum(shared) >= 3) safe_cor(log2(sample_a_cpm[shared] + pseudocount), log2(sample_b_cpm[shared] + pseudocount)) else NA_real_,
      UnionSpearman = if (sum(union) >= 3) safe_cor(log2(sample_a_cpm[union] + pseudocount), log2(sample_b_cpm[union] + pseudocount)) else NA_real_,
      SharedMedianLog2FC = if (sum(shared) > 0) stats::median(shared_logfc) else NA_real_,
      SharedMedianAbsLog2FC = if (sum(shared) > 0) stats::median(abs(shared_logfc)) else NA_real_,
      SharedWithinTwoFoldPct = if (sum(shared) > 0) 100 * mean(abs(shared_logfc) <= 1) else NA_real_
    )
  })
}

ROOT_DIR <- normalizePath(Sys.getenv("PANDORA_PROJECT_ROOT"), mustWork = TRUE)

RAW_DIR           <- file.path(ROOT_DIR, "01rawdata")

SPORTS_RESULT_DIR <- file.path(ROOT_DIR, "sportsresult")


METADATA_FILE     <- file.path(ROOT_DIR, "config", "samples.tsv")

RESULT_DIR        <- Sys.getenv("PANDORA_RESULT_DIR", file.path(ROOT_DIR, "result", "describe"))

INTERNAL_CACHE_DIR <- file.path(ROOT_DIR, ".cache", "shared_analysis")

PROJECT_SPECIES <- "Homo sapiens"

PROJECT_HOST_LABEL <- "human"

PROJECT_SAMPLE_TYPE <- Sys.getenv("PROJECT_SAMPLE_TYPE", Sys.getenv("PANDORA_PROJECT_NAME"))

PROJECT_MATRIX <- Sys.getenv("PROJECT_MATRIX", Sys.getenv("PANDORA_PROJECT_NAME"))

PROJECT_CONDITION <- Sys.getenv("PROJECT_CONDITION", "unspecified")

SPIKEIN_SEQ        <- "TCACCGGGTGTAAATCAGCTTG"

SPIKEIN_NAME       <- "cel-miR-39"

SPIKEIN_AMOUNT_FMOL <- 60

SPIKEIN_SAMPLE     <- "Normal3_1_spikein"

SPIKEIN_CONTROL_SAMPLE <- "Normal3_1_1"

METHOD_COMPARISON_BIOLOGICAL_SAMPLE <- "2"

SPIKEIN_BIOLOGICAL_SAMPLE <- "3"

ANALYSIS_SCOPE_LEVELS <- c("all", "method_comparison", "spikein_pair")

spikein_edit1_sequences <- function(target = SPIKEIN_SEQ) {
  chars <- strsplit(as.character(target), "", fixed = TRUE)[[1]]
  alphabet <- c(LETTERS, "N")
  substitutions <- unlist(lapply(seq_along(chars), function(i) {
    vapply(setdiff(alphabet, chars[[i]]), function(base) {
      left <- if (i > 1) chars[seq_len(i - 1)] else character()
      right <- if (i < length(chars)) chars[(i + 1):length(chars)] else character()
      paste0(c(left, base, right), collapse = "")
    }, character(1))
  }), use.names = FALSE)
  deletions <- vapply(seq_along(chars), function(i) {
    paste0(chars[-i], collapse = "")
  }, character(1))
  insertions <- unlist(lapply(seq_len(length(chars) + 1L), function(i) {
    vapply(alphabet, function(base) {
      left <- if (i > 1) chars[seq_len(i - 1)] else character()
      right <- if (i <= length(chars)) chars[i:length(chars)] else character()
      paste0(c(left, base, right), collapse = "")
    }, character(1))
  }), use.names = FALSE)
  unique(c(target, substitutions, deletions, insertions))
}

SPIKEIN_DIR        <- file.path(RESULT_DIR, "13_spikein_qc")

SPIKEIN_FIG        <- file.path(SPIKEIN_DIR, "plots")

SPIKEIN_TAB        <- file.path(SPIKEIN_DIR, "tables")

SETUP_DIR         <- file.path(RESULT_DIR, "1_setup")

EXPRESSION_MATRIX_DIR <- file.path(SETUP_DIR, "1_expression_matrix")

SINGLE_DIR        <- file.path(RESULT_DIR, "2_single")

SINGLE_FIG        <- file.path(SINGLE_DIR, "plots")

SINGLE_TAB        <- file.path(SINGLE_DIR, "tables")

PCA_DIR           <- file.path(RESULT_DIR, "3_pca")

PCA_FIG           <- file.path(PCA_DIR, "plots")

PCA_TAB           <- file.path(PCA_DIR, "tables")

SNC_COMPARE_DIR   <- file.path(RESULT_DIR, "4_sncRNA_compare")

SNC_COMPARE_FIG   <- file.path(SNC_COMPARE_DIR, "plots")

SNC_COMPARE_TAB   <- file.path(SNC_COMPARE_DIR, "tables")

TSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "5_tsrna_compare")

TSRNA_COMPARE_FIG <- file.path(TSRNA_COMPARE_DIR, "plots")

TSRNA_COMPARE_TAB <- file.path(TSRNA_COMPARE_DIR, "tables")

RSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "6_rsrna_compare")

RSRNA_COMPARE_FIG <- file.path(RSRNA_COMPARE_DIR, "plots")

RSRNA_COMPARE_TAB <- file.path(RSRNA_COMPARE_DIR, "tables")

MIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "7_mirna_compare")

MIRNA_COMPARE_FIG <- file.path(MIRNA_COMPARE_DIR, "plots")

MIRNA_COMPARE_TAB <- file.path(MIRNA_COMPARE_DIR, "tables")

PIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "8_pirna_compare")

PIRNA_COMPARE_FIG <- file.path(PIRNA_COMPARE_DIR, "plots")

PIRNA_COMPARE_TAB <- file.path(PIRNA_COMPARE_DIR, "tables")

CAPTURE_DIR        <- file.path(RESULT_DIR, "9_differential_capture")

CAPTURE_FIG        <- file.path(CAPTURE_DIR, "plots")

CAPTURE_TAB        <- file.path(CAPTURE_DIR, "tables")

CAPTURE_BIAS_DIR   <- file.path(RESULT_DIR, "10_capture_bias_diagnostics")

CAPTURE_BIAS_FIG   <- file.path(CAPTURE_BIAS_DIR, "plots")

CAPTURE_BIAS_TAB   <- file.path(CAPTURE_BIAS_DIR, "tables")

SNRNA_COMPARE_DIR  <- file.path(RESULT_DIR, "12_snRNA_compare")

SNRNA_COMPARE_FIG  <- file.path(SNRNA_COMPARE_DIR, "plots")

SNRNA_COMPARE_TAB  <- file.path(SNRNA_COMPARE_DIR, "tables")

RUN_FASTQC <- tolower(Sys.getenv("RUN_FASTQC", "auto"))

RUN_FASTQC <- if (RUN_FASTQC %in% c("auto", "1", "0", "true", "false", "yes", "no")) RUN_FASTQC else "auto"

DE_SIGNIFICANCE <- toupper(Sys.getenv("DE_SIGNIFICANCE", "FDR"))

DE_SIGNIFICANCE <- if (DE_SIGNIFICANCE %in% c("FDR", "PVALUE", "P_VALUE", "P")) DE_SIGNIFICANCE else "FDR"

REPORT_VIEW <- "sense"

ALL_VIEW_LEVELS <- c("merged", "sense", "antisense")

VIEW_LEVELS <- c(REPORT_VIEW)

OUTPUT_FIG_DIR <- NULL

OUTPUT_TAB_DIR <- NULL

capture_command_version <- function(command, args = c("--version")) {
  command_path <- Sys.which(command)
  if (!nzchar(command_path)) {
    return(tibble(tool = command, available = FALSE, path = "", version = "not available"))
  }

  version_text <- tryCatch(
    {
      timeout_path <- Sys.which("timeout")
      if (nzchar(timeout_path)) {
        out <- suppressWarnings(system2(timeout_path, args = c("15", command_path, args), stdout = TRUE, stderr = TRUE))
      } else {
        out <- suppressWarnings(system2(command_path, args = args, stdout = TRUE, stderr = TRUE))
      }
      status <- attr(out, "status")
      if (!is.null(status) && identical(as.integer(status), 124L)) "version check timed out"
      else paste(head(out, 3), collapse = " | ")
    },
    error = function(e) paste("version check failed:", conditionMessage(e))
  )

  tibble(tool = command, available = TRUE, path = unname(command_path), version = version_text)
}

build_runtime_config_table <- function() {
  tibble(
    key = c(
      "ROOT_DIR", "RAW_DIR", "SPORTS_RESULT_DIR", "RESULT_DIR",
      "PROJECT_SPECIES", "PROJECT_SAMPLE_TYPE", "PROJECT_MATRIX", "PROJECT_CONDITION",
      "RUN_FASTQC", "DE_SIGNIFICANCE", "REPORT_VIEW"
    ),
    value = c(
      ROOT_DIR, RAW_DIR, SPORTS_RESULT_DIR, RESULT_DIR,
      PROJECT_SPECIES, PROJECT_SAMPLE_TYPE, PROJECT_MATRIX, PROJECT_CONDITION,
      RUN_FASTQC, DE_SIGNIFICANCE, REPORT_VIEW
    )
  )
}

parse_tsrna_amino_acid <- function(sub_class) { out <- str_match(sub_class,"(?:mt_)?tRNA-([A-Za-z]+)-[A-Z]{3}")[,2]; if_else(is.na(out)|out=="","Unknown",out) }

parse_tsrna_anticodon <- function(sub_class) { out <- str_match(sub_class,"(?:mt_)?tRNA-[A-Za-z]+-([A-Z]{3})")[,2]; if_else(is.na(out)|out=="","Unknown",out) }

make_group_palette <- function(groups) {
  base <- c("#8ECFC9", "#FFBE7A", "#FF9F9B", "#A1C9F4", "#BEB8DC", "#FFED6F", "#FDCDAC", "#CBD5E8", "#F4CAE4")
  groups <- as.character(groups)
  if (length(groups) > length(base)) {
    extra <- grDevices::hcl.colors(length(groups) - length(base), palette = "Pastel 1")
    base <- c(base, extra)
  }
  setNames(base[seq_along(groups)], groups)
}

category_levels <- c("rsRNA","tsRNA","miRNA","piRNA","ysRNA","lncRNA","snRNA","snoRNA","scaRNA","misc_RNA","vault_RNA","ribozyme","other","unanno_MG","unanno_UMG")

sports_main_palette <- c(
  unanno_UMG="#8ECFC9", miRNA="#FFBE7A", piRNA="#FF9F9B", tsRNA="#A1C9F4", 
  rsRNA="#BEB8DC", other="#B3E2CD", unanno_MG="#FDCDAC", antisense="#CBD5E8", 
  ysRNA="#F4CAE4"
)

main_palette <- c(
  sports_main_palette, 
  lncRNA="#E6F5C9", snRNA="#FFF2AE", snoRNA="#F1E2CC", scaRNA="#CCCCCC", 
  misc_RNA="#B2DF8A", ribozyme="#FB8072", mt_tsRNA="#FDB462", mt_rsRNA="#FCCDE5", 
  DNA="#D9D9D9", mt_DNA="#BC80BD", scRNA="#CCEBC5", sRNA="#FFED6F", vault_RNA="#8DD3C7", 
  mt="#FF9F9B", cyto="#A1C9F4", Genomic="#A1C9F4", Mitochondrial="#FF9F9B"
)

category_palette <- main_palette[category_levels]

radar_levels  <- c("rsRNA","tsRNA","miRNA","piRNA","ysRNA","lncRNA","snRNA","snoRNA","other","unanno_MG","unanno_UMG")

radar_labels  <- c(other="other", unanno_MG="unanno MG", unanno_UMG="unanno UMG",rsRNA="rsRNA", tsRNA="tsRNA", miRNA="miRNA", piRNA="piRNA", ysRNA="ysRNA", lncRNA="lncRNA", snRNA="snRNA", snoRNA="snoRNA")

radar_sector_palette <- main_palette[radar_levels]

radar_label_palette  <- make_sector_label_palette(radar_sector_palette, amount = 0.28)

sports_length_categories <- c("miRNA","piRNA","ysRNA","tsRNA","mt_tsRNA","rsRNA","mt_rsRNA","lncRNA","snRNA","snoRNA","scaRNA","misc_RNA","vault_RNA","ribozyme","scRNA","sRNA","DNA","mt_DNA","other","unanno_MG","unanno_UMG")

sports_length_palette <- main_palette[sports_length_categories]

rrna_display_order <- c("5S","5.8S","12S","16S","18S","28S","45S","4.5S","other")

rrna_coverage_subtypes <- setdiff(rrna_display_order, "other")

rrna_palette <- c("12S"="#8ECFC9","16S"="#FFBE7A","5S"="#FF9F9B","5.8S"="#A1C9F4","18S"="#BEB8DC","28S"="#B3E2CD","45S"="#FDCDAC","4.5S"="#F4CAE4","other"="#CBD5E8")

rsrna_visual_levels  <- c("5S-rRNA","5.8S-rRNA","12S-rRNA","16S-rRNA","18S-rRNA","28S-rRNA","45S-rRNA","4.5S-rRNA","2S-rRNA","mt-rRNA","Ensembl-rRNA","Rfam-rRNA-like","other-rRNA")

rsrna_visual_palette <- c(
  "12S-rRNA" = "#8ECFC9",
  "16S-rRNA" = "#FFBE7A",
  "5S-rRNA" = "#FF9F9B",
  "5.8S-rRNA" = "#A1C9F4",
  "18S-rRNA" = "#BEB8DC",
  "28S-rRNA" = "#B3E2CD",
  "45S-rRNA" = "#FDCDAC",
  "other-rRNA" = "#CBD5E8",
  "4.5S-rRNA" = "#F4CAE4",
  "2S-rRNA" = "#E6F5C9",
  "mt-rRNA" = "#FFF2AE",
  "Ensembl-rRNA" = "#F1E2CC",
  "Rfam-rRNA-like" = "#CCCCCC"
)

rsrna_visual_labels  <- c("5S-rRNA"="5S","5.8S-rRNA"="5.8S","12S-rRNA"="12S","16S-rRNA"="16S","18S-rRNA"="18S","28S-rRNA"="28S","45S-rRNA"="45S","4.5S-rRNA"="4.5S","2S-rRNA"="2S","mt-rRNA"="mt-rRNA","Ensembl-rRNA"="Ensembl","Rfam-rRNA-like"="Rfam","other-rRNA"="other")

rsrna_radar_levels  <- c("12S-rRNA","16S-rRNA","5S-rRNA","5.8S-rRNA","18S-rRNA","28S-rRNA","45S-rRNA","other-rRNA")

rsrna_radar_palette <- rsrna_visual_palette[rsrna_radar_levels]

rsrna_radar_labels  <- rsrna_visual_labels[rsrna_radar_levels]

ts_origin_palette    <- c("Genomic tRNA"="#A1C9F4","Mitochondrial tRNA"="#FF9F9B")

ts_genomic_palette   <- c("Genomic tRNA 5'end"="#8ECFC9","Genomic tRNA 3'end"="#FFBE7A","Genomic tRNA CCA-end"="#FF9F9B","Genomic tRNA internal"="#A1C9F4")

ts_mito_palette      <- c("Mitochondrial tRNA 5'end"="#8ECFC9","Mitochondrial tRNA 3'end"="#FFBE7A","Mitochondrial tRNA CCA-end"="#FF9F9B","Mitochondrial tRNA internal"="#A1C9F4")

ts_cleavage_palette  <- c("5'end"="#8ECFC9","3'end"="#FFBE7A","CCA-end"="#FF9F9B","Internal/full"="#A1C9F4")

ts_maturity_palette  <- c("Genomic precursor tRNA"="#B3E2CD","Genomic mature tRNA"="#A1C9F4","Mito mature tRNA"="#FF9F9B")

mt_cyto_palette      <- c(mt="#FF9F9B", cyto="#A1C9F4")

rsrna_source_palette <- c("rRNAdb known subtypes"="#A1C9F4","rRNAdb other"="#CBD5E8","Ensembl rRNA"="#FFBE7A","Ensembl Mt rRNA"="#FF9F9B")

canonical_aa_order <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","SeC","iMet")

amino_acid_palette <- c(
  Ala = "#8ECFC9", Arg = "#FFBE7A", Asn = "#FF9F9B", Asp = "#A1C9F4", Cys = "#BEB8DC",
  Gln = "#B3E2CD", Glu = "#FDCDAC", Gly = "#CBD5E8", His = "#F4CAE4", Ile = "#E6F5C9",
  Leu = "#FFF2AE", Lys = "#F1E2CC", Met = "#CCCCCC", Phe = "#B2DF8A", Pro = "#FB8072",
  Ser = "#FDB462", Thr = "#FCCDE5", Trp = "#D9D9D9", Tyr = "#BC80BD", Val = "#CCEBC5",
  SeC = "#FFED6F", iMet = "#8DD3C7", Unknown = "#E5D8BD"
)

serum_blue <- "#A1C9F4"

serum_orange <- "#FFBE7A"

serum_green <- "#8ECFC9"

serum_red <- "#FF9F9B"

serum_purple <- "#BEB8DC"

serum_neutral_line <- "#C9CED6"

serum_neutral_text <- "#9AA1AC"

serum_axis_label_color <- "#6D7784"

serum_text_primary <- "#1F2329"

serum_text_muted <- "#47596E"

serum_transition_yellow <- "#FDFFB6"

serum_qc_palette <- c(CleanReads = serum_blue, MatchGenome = serum_green)

serum_diverging_palette <- c(low = serum_blue, mid = serum_transition_yellow, high = serum_red)

serum_direction_palette <- c(up = serum_red, down = serum_blue, ns = serum_neutral_line)

serum_direction_binary_palette <- c(up = serum_red, down = serum_blue)

serum_target_group_palette <- c(all = serum_purple, up = serum_red, down = serum_blue)

serum_enrichment_gradient <- c(serum_blue, serum_transition_yellow, serum_red)

serum_heatmap_gradient <- c(serum_blue, serum_transition_yellow, serum_red)

serum_theme_bw <- function(...) {
  theme_bw() + theme(panel.grid.minor = element_blank(), ...)
}

serum_theme_minimal <- function(...) {
  theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    ...
  )
}

load_metadata <- function() {
  src <- METADATA_FILE
  meta <- read_tsv_flex(src)
  if (nrow(meta)==0) stop("未检测到样本信息文件内容: ", src)
  names(meta) <- tolower(names(meta))
  required_cols <- c(
    "sample", "group", "method", "biological_sample", "replicate", "batch",
    "analysis_role", "pair_id", "spikein_sequence", "spikein_input_fmol"
  )
  missing_cols <- setdiff(required_cols, names(meta))
  if (length(missing_cols) > 0) {
    stop("样本元数据缺少必需列: ", paste(missing_cols, collapse = ", "), ". 文件: ", src)
  }
  if (!"raw_fastq_gz" %in% names(meta)) stop("样本信息文件必须包含 raw_fastq_gz 列: ", src)
  if (anyNA(meta$sample) || anyNA(meta$group) || any(!nzchar(trimws(as.character(meta$sample)))) || any(!nzchar(trimws(as.character(meta$group)))) || anyDuplicated(meta$sample)) stop("sample/group 不能为空，sample 必须唯一: ",src)
  out <- tibble(
    Group=trimws(as.character(meta$group)), Sample=trimws(as.character(meta$sample)),
    Method=trimws(as.character(meta$method)),
    Condition=if("condition"%in%names(meta)) trimws(as.character(meta$condition)) else PROJECT_CONDITION,
    SampleType=if("sample_type"%in%names(meta)) trimws(as.character(meta$sample_type)) else PROJECT_SAMPLE_TYPE,
    Species=if("species"%in%names(meta)) trimws(as.character(meta$species)) else PROJECT_SPECIES,
    Matrix=if("matrix"%in%names(meta)) trimws(as.character(meta$matrix)) else PROJECT_MATRIX,
    BiologicalSample=trimws(as.character(meta$biological_sample)),
    Replicate=trimws(as.character(meta$replicate)),
    Batch=trimws(as.character(meta$batch)),
    AnalysisRole=trimws(as.character(meta$analysis_role)),
    PairID=trimws(as.character(meta$pair_id)),
    SpikeinSequence=trimws(as.character(meta$spikein_sequence)),
    SpikeinInputFmol=clean_numeric(meta$spikein_input_fmol),
    RawFastqOriginal=trimws(as.character(meta$raw_fastq_gz))
  ) %>% filter(Group!="", Sample!="") %>%
    mutate(
      Condition = if_else(is.na(Condition) | Condition == "", PROJECT_CONDITION, Condition),
      SampleType = if_else(is.na(SampleType) | SampleType == "", PROJECT_SAMPLE_TYPE, SampleType),
      Species = if_else(is.na(Species) | Species == "", PROJECT_SPECIES, Species),
      Matrix = if_else(is.na(Matrix) | Matrix == "", PROJECT_MATRIX, Matrix),
      Batch = if_else(is.na(Batch), "", Batch),
      AnalysisRole = if_else(is.na(AnalysisRole), "", AnalysisRole),
      PairID = if_else(is.na(PairID), "", PairID),
      SpikeinSequence = if_else(is.na(SpikeinSequence), "", SpikeinSequence),
      SpikeinInputFmol = coalesce(SpikeinInputFmol, 0),
      Spikein = AnalysisRole == "spikein_pair_test",
      GroupOrder=suppressWarnings(as.integer(Group)),
      GroupOrder=if_else(is.na(GroupOrder),row_number(),GroupOrder),
      Sample=factor(Sample,levels=Sample),
      Group=factor(Group,levels=unique(Group[order(GroupOrder)]))
    ) %>%
    arrange(GroupOrder, Sample) %>%
    mutate(
      RawFastq = vapply(RawFastqOriginal, resolve_raw_fastq_path, character(1)),
      RawFastqExists = file.exists(RawFastq),
      SummaryFile=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_summary.txt")),
      LengthFile=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_length_distribution.txt")),
      OutputFile=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_output.txt")),
      PDF_sncRNA=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_sncRNA_distribution.pdf")),
      PDF_rsRNA=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_rsRNA_distribution.pdf")),
      PDF_rRNA=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_rRNA_mapping.pdf")),
      PDF_tRNA=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_tRNA_mature_mapping.pdf")),
      RRNAMatchFile=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_processed"),paste0(as.character(Sample),"_output_rRNA_match_genome")),
      RRNAUnmatchFile=file.path(SPORTS_RESULT_DIR,as.character(Sample),paste0(as.character(Sample),"_processed"),paste0(as.character(Sample),"_output_rRNA_unmatch_genome"))
    )
  missing <- out %>% mutate(Missing=!file.exists(SummaryFile)|!file.exists(LengthFile)|!file.exists(OutputFile)) %>% filter(Missing)
  if (nrow(missing)>0) stop("未检测到完整 SPORTS 结果，请通过所属分支 run.sh --stage upstream 生成输入。缺失样本: ", paste(as.character(missing$Sample),collapse=", "))
  if (Sys.getenv("PANDORA_DESIGN") == "method_capture") {
  allowed_roles <- c("method_comparison", "spikein_pair_control", "spikein_pair_test")
  if (any(!out$AnalysisRole %in% allowed_roles)) {
    stop("样本包含无效 AnalysisRole: ", paste(unique(out$AnalysisRole[!out$AnalysisRole %in% allowed_roles]), collapse = ", "))
  }
  method_rows <- out %>% filter(AnalysisRole == "method_comparison")
  if (nrow(method_rows) != 10L || any(as.character(method_rows$BiologicalSample) != METHOD_COMPARISON_BIOLOGICAL_SAMPLE)) {
    stop("method_comparison 队列必须恰好包含 2 号供体的 10 个库。")
  }
  method12_rows <- method_rows %>% filter(as.character(Method) %in% c("1", "2"))
  method12_pairs <- method12_rows %>%
    group_by(PairID) %>%
    summarise(
      SampleN = n(), MethodN = n_distinct(Method), ReplicateN = n_distinct(Replicate),
      BatchN = n_distinct(Batch), DonorN = n_distinct(BiologicalSample), .groups = "drop"
    )
  if (nrow(method12_rows) != 6L || nrow(method12_pairs) != 3L || any(method12_pairs$PairID == "") ||
      any(method12_pairs$SampleN != 2L) || any(method12_pairs$MethodN != 2L) ||
      any(method12_pairs$ReplicateN != 1L) || any(method12_pairs$BatchN != 1L) || any(method12_pairs$DonorN != 1L)) {
    stop("方法 1/2 必须形成三个供体、replicate 和 batch 均一致的一一匹配分装。")
  }
  if (sum(out$AnalysisRole == "spikein_pair_test") != 1L || sum(out$AnalysisRole == "spikein_pair_control") != 1L) {
    stop("spike-in 配对队列必须恰好包含一个 control 和一个 spike-in 样本。")
  }
  pair_rows <- out %>% filter(AnalysisRole %in% c("spikein_pair_control", "spikein_pair_test"))
  if (any(as.character(pair_rows$BiologicalSample) != SPIKEIN_BIOLOGICAL_SAMPLE) ||
      any(as.character(pair_rows$Method) != "1") ||
      any(as.character(pair_rows$Group) != "1") ||
      n_distinct(pair_rows$PairID) != 1L || any(pair_rows$PairID == "") ||
      n_distinct(pair_rows$Batch) != 1L) {
    stop("spike-in 配对样本必须来自 3 号供体、方法 1 和第 1 组。")
  }
  spike_row <- pair_rows %>% filter(AnalysisRole == "spikein_pair_test")
  if (nrow(spike_row) != 1L ||
      !identical(as.character(spike_row$SpikeinSequence[[1]]), SPIKEIN_SEQ) ||
      !isTRUE(as.numeric(spike_row$SpikeinInputFmol[[1]]) == SPIKEIN_AMOUNT_FMOL) ||
      !identical(SPIKEIN_AMOUNT_FMOL, 60) ||
      !identical(SPIKEIN_SEQ, "TCACCGGGTGTAAATCAGCTTG")) {
    stop("spike-in 元数据必须明确记录 cel-miR-39 序列和 60 fmol 投入量。")
  }
  }
  out
}

load_serum_context <- function(report_view = REPORT_VIEW, analysis_scope = "all") {
  if (!analysis_scope %in% ANALYSIS_SCOPE_LEVELS) {
    stop("Unsupported analysis_scope: ", analysis_scope, ". Choose one of: ", paste(ANALYSIS_SCOPE_LEVELS, collapse = ", "))
  }
  sample_info_all <- load_metadata()
  sample_info <- switch(
    analysis_scope,
    all = sample_info_all,
    method_comparison = sample_info_all %>% filter(AnalysisRole == "method_comparison"),
    spikein_pair = sample_info_all %>% filter(AnalysisRole %in% c("spikein_pair_control", "spikein_pair_test"))
  )
  if (nrow(sample_info) == 0) stop("analysis_scope 没有可用样本: ", analysis_scope)
  if (analysis_scope == "spikein_pair") {
    sample_info <- sample_info %>%
      arrange(match(as.character(Sample), c(SPIKEIN_CONTROL_SAMPLE, SPIKEIN_SAMPLE))) %>%
      mutate(Group = as.character(Sample))
  }
  sample_info <- sample_info %>%
    mutate(
      Sample = factor(as.character(Sample), levels = unique(as.character(Sample))),
      Group = factor(as.character(Group), levels = unique(as.character(Group)))
    )
  group_levels <- levels(sample_info$Group)
  group_palette <- make_group_palette(group_levels)

  # 原始 summary 和 length 表先统一加上方向、数据库来源、主类别等派生字段。
  summary_all_raw <- load_summary_tables(sample_info) %>%
    mutate(
      Direction = class_direction(Class, SubClass),
      GenomeStatus = class_genome_status(Class),
      SourceDB = class_source_db(Class),
      MainCategory = class_main_category(Class, SubClass)
    )
  length_all_raw <- load_length_tables(sample_info) %>%
    mutate(
      Direction = class_direction(Class),
      GenomeStatus = class_genome_status(Class),
      SourceDB = class_source_db(Class),
      MainCategory = class_main_category(Class),
      LengthCategory = length_category_from_class(Class)
    )
  output_all <- load_output_tables(sample_info)
  integrated <- integrate_endosirna(summary_all_raw, length_all_raw, output_all, sample_info)
  summary_all_raw <- integrated$summary
  length_all_raw <- integrated$length
  output_all <- integrated$output
  endo_evidence <- integrated$evidence
  endo_transfer_audit <- integrated$audit
  endo_reference_audit <- integrated$references

  # Clean_Reads / Match_Genome 是很多百分比与 RPM 计算的分母，因此单独提前抽出来。
  clean_reads_tbl <- summary_all_raw %>%
    filter(Class == "Clean_Reads") %>%
    transmute(Sample, CleanReads = Reads)
  match_reads_tbl <- summary_all_raw %>%
    filter(Class == "Match_Genome") %>%
    transmute(Sample, MatchGenome = Reads)

  # 根据 report_view 过滤方向。当前默认是 sense，但函数本身保留了 merged/antisense 扩展能力。
  summary_all <- filter_direction_view(summary_all_raw, report_view)
  summary_all_view <- summary_all
  length_all <- filter_direction_view(length_all_raw, report_view)
  length_all_view <- length_all

  # summary_top 仅保留真正的主类别汇总项，去掉 Clean_Reads / Match_Genome 等技术性条目。
  # tRNA 末端类是基础 tRNA 类的子集，也必须排除，否则类别组成会重复计数。
  summary_top <- summary_all %>%
    filter(
      SubClass == "-",
      !Class %in% c("Clean_Reads", "Match_Genome"),
      !is_tsrna_nested_class(Class)
    ) %>%
    mutate(MainCategory = factor(MainCategory, levels = category_levels))
  summary_top_view <- summary_top

  # 组级 clean/match 读数在大量组间比较图中复用，因此单独聚合一次。
  group_clean_reads_tbl <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group)) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    group_by(Group) %>%
    summarise(CleanReads = sum(CleanReads, na.rm = TRUE), .groups = "drop")
  group_match_tbl <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group)) %>%
    left_join(match_reads_tbl, by = "Sample") %>%
    group_by(Group) %>%
    summarise(MatchGenome = sum(MatchGenome, na.rm = TRUE), .groups = "drop")

  # 样本级百分比先计算；组间百分比用样本等权均值，避免重复数或测序深度改变比较尺度。
  sample_comp_total <- summary_top %>%
    group_by(Sample, Group, MainCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(Percent = if_else(CleanReads > 0, 100 * Reads / CleanReads, 0))
  sample_comp_view <- sample_comp_total
  composition_integrity_check <- sample_comp_total %>%
    group_by(Sample) %>%
    summarise(
      ExclusiveClassReads = sum(Reads, na.rm = TRUE),
      CleanReads = first(CleanReads),
      ReadDifference = ExclusiveClassReads - CleanReads,
      PercentSum = sum(Percent, na.rm = TRUE),
      PercentDifference = PercentSum - 100,
      Pass = abs(ReadDifference) <= 2,
      .groups = "drop"
    )
  if (any(!composition_integrity_check$Pass)) {
    bad_samples <- composition_integrity_check %>% filter(!Pass) %>% pull(Sample)
    stop("互斥 RNA 类别组成未能回收到 100%: ", paste(bad_samples, collapse = ", "))
  }
  sample_comp_direction <- summary_top %>%
    group_by(Sample, Group, Direction, MainCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(Percent = if_else(CleanReads > 0, 100 * Reads / CleanReads, 0))
  group_comp_total <- sample_comp_total %>%
    summarise_group_mean_percent("MainCategory") %>%
    left_join(group_clean_reads_tbl, by = "Group")
  group_comp_view <- group_comp_total
  group_comp_direction <- sample_comp_direction %>%
    summarise_group_mean_percent(c("Direction", "MainCategory")) %>%
    left_join(group_clean_reads_tbl, by = "Group")
  group_comp_direction_view <- group_comp_direction %>% select(-Direction)

  sample_qc_table <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), BiologicalSample = BiologicalSample, Replicate = Replicate) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    left_join(match_reads_tbl, by = "Sample") %>%
    left_join(
      sample_comp_view %>%
        arrange(Sample, desc(Percent)) %>%
        group_by(Sample) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        transmute(Sample, DominantCategory = as.character(MainCategory), DominantPercent = round(Percent, 2)),
      by = "Sample"
    ) %>%
    mutate(GenomeRate = round(if_else(CleanReads > 0, 100 * MatchGenome / CleanReads, 0), 2))

  group_qc_table <- group_clean_reads_tbl %>%
    left_join(group_match_tbl, by = "Group") %>%
    left_join(
      sample_qc_table %>%
        group_by(Group = as.character(Group)) %>%
        summarise(SampleN = n_distinct(Sample), GenomeRate = round(mean(GenomeRate, na.rm = TRUE), 2), .groups = "drop"),
      by = "Group"
    )

  # length_detail 是长度分布图的统一底表，RPM 也在这里一次性计算好。
  length_detail <- length_all %>%
    filter(
      !Class %in% c("Clean_Reads", "Match_Genome"),
      !is_tsrna_nested_class(Class),
      !is_nested_length_class(Class)
    ) %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0))
  length_detail_view <- length_detail
  length_focus_categories <- length_detail_view %>%
    group_by(MainCategory) %>%
    summarise(TotalRPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalRPM)) %>%
    pull(MainCategory) %>%
    as.character()
  length_focus_categories <- unique(c(c("rsRNA", "tsRNA", "miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA"), length_focus_categories)) %>%
    intersect(category_levels)
  sample_length_display <- length_detail_view %>%
    filter(MainCategory %in% length_focus_categories) %>%
    group_by(Sample, Group, Direction, MainCategory, Length) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop")
  sample_length_view <- sample_length_display %>% select(-Direction)
  group_length_display <- sample_length_display %>%
    summarise_group_mean_rpm(c("Direction", "MainCategory", "Length"))
  group_length_view <- group_length_display %>% select(-Direction)

  # rsRNA 先做 subtype 解析，再做样本级/组级占比统计。
  # visual 版本会把低频 subtype 合并，方便雷达图和 alluvial 展示。
  rsrna_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "rsRNA") %>%
    mutate(RsRNAClass = parse_rsrna_class(Class, SubClass))
  sample_rsrna <- rsrna_detail %>%
    group_by(Sample, Group, Direction, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_rsrna <- sample_rsrna %>%
    summarise_group_mean_percent(c("Direction", "RsRNAClass"))
  sample_rsrna_visual <- sample_rsrna %>%
    mutate(RsRNAClass = collapse_rsrna_visual_class(as.character(RsRNAClass))) %>%
    group_by(Sample, Group, Direction, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup() %>%
    mutate(RsRNAClass = factor(RsRNAClass, levels = rsrna_visual_levels))
  sample_rsrna_visual_by_view <- setNames(
    list(sample_rsrna_visual %>% select(-Direction)),
    report_view
  )
  sample_rsrna_visual_view <- sample_rsrna_visual_by_view[[report_view]]
  group_rsrna_visual <- sample_rsrna_visual %>%
    summarise_group_mean_percent(c("Direction", "RsRNAClass")) %>%
    mutate(RsRNAClass = factor(RsRNAClass, levels = rsrna_visual_levels))
  group_rsrna_visual_view <- group_rsrna_visual %>% select(-Direction)
  group_rsrna_visual_merged <- group_rsrna_visual_view %>%
    group_by(Group, RsRNAClass) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), Percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
    ungroup()
  group_rsrna_view <- group_rsrna %>% select(-Direction)

  # SPORTS 的 tRNA 末端行是基础 tRNA 行的子集。基础表用于总量/origin/AA/anticodon，
  # cleavage 表则把 Internal/full 重建为“基础总量 - 5' - 3' - CCA”，保证四类互斥。
  tsrna_base_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "tsRNA", !is_tsrna_nested_class(Class)) %>%
    mutate(
      BaseClass = as.character(Class),
      BaseSubClass = as.character(SubClass),
      Origin = parse_tsrna_origin(Class, SubClass),
      Maturity = parse_tsrna_maturity(Class, SubClass),
      AminoAcid = parse_tsrna_amino_acid(SubClass),
      Anticodon = parse_tsrna_anticodon(SubClass)
    )

  tsrna_end_detail <- summary_all %>%
    filter(SubClass != "-", MainCategory == "tsRNA", is_tsrna_nested_class(Class)) %>%
    mutate(
      BaseClass = strip_tsrna_end_from_class(Class),
      BaseSubClass = strip_tsrna_end_from_subclass(SubClass),
      Origin = parse_tsrna_origin(BaseClass, BaseSubClass),
      Maturity = parse_tsrna_maturity(BaseClass, BaseSubClass),
      Cleavage = parse_tsrna_cleavage(Class, SubClass),
      AminoAcid = parse_tsrna_amino_acid(BaseSubClass),
      Anticodon = parse_tsrna_anticodon(BaseSubClass)
    )

  tsrna_reconstruction_keys <- c(
    "Sample", "Group", "Direction", "BaseClass", "BaseSubClass",
    "Origin", "Maturity", "AminoAcid", "Anticodon"
  )
  orphan_tsrna_ends <- tsrna_end_detail %>%
    anti_join(tsrna_base_detail, by = tsrna_reconstruction_keys)
  if (nrow(orphan_tsrna_ends) > 0) {
    stop("SPORTS tRNA 末端子集缺少对应基础总量行，无法重建互斥 cleavage 计数。")
  }

  tsrna_end_totals <- tsrna_end_detail %>%
    group_by(across(all_of(tsrna_reconstruction_keys))) %>%
    summarise(EndReads = sum(Reads, na.rm = TRUE), .groups = "drop")

  tsrna_cleavage_reconstruction_check <- tsrna_base_detail %>%
    transmute(
      across(all_of(tsrna_reconstruction_keys)),
      BaseReads = Reads
    ) %>%
    left_join(tsrna_end_totals, by = tsrna_reconstruction_keys) %>%
    mutate(
      EndReads = coalesce(EndReads, 0),
      InternalReadsRaw = BaseReads - EndReads,
      ExcessEndReads = pmax(EndReads - BaseReads, 0),
      Pass = ExcessEndReads <= 1,
      InternalReads = pmax(InternalReadsRaw, 0),
      ReconstructedReads = EndReads + InternalReads,
      ResidualReads = ReconstructedReads - BaseReads
    )
  if (any(!tsrna_cleavage_reconstruction_check$Pass)) {
    bad <- tsrna_cleavage_reconstruction_check %>% filter(!Pass) %>% pull(Sample) %>% unique()
    stop("tRNA 末端子集和超过基础总量 1 read: ", paste(bad, collapse = ", "))
  }

  tsrna_internal_detail <- tsrna_cleavage_reconstruction_check %>%
    transmute(
      Sample, Group, Direction,
      Class = BaseClass,
      SubClass = BaseSubClass,
      Reads = InternalReads,
      MainCategory = "tsRNA",
      BaseClass, BaseSubClass, Origin, Maturity,
      Cleavage = "Internal/full", AminoAcid, Anticodon
    )
  tsrna_cleavage_detail <- bind_rows(
    tsrna_end_detail %>%
      transmute(
        Sample, Group, Direction, Class, SubClass, Reads,
        MainCategory = "tsRNA", BaseClass, BaseSubClass, Origin, Maturity,
        Cleavage, AminoAcid, Anticodon
      ),
    tsrna_internal_detail
  ) %>%
    filter(Reads > 0)

  sample_tsrna_origin <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, Origin) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_origin <- sample_tsrna_origin %>%
    summarise_group_mean_percent(c("Direction", "Origin"))
  group_tsrna_origin_view <- group_tsrna_origin %>% select(-Direction)
  sample_tsrna_cleavage <- tsrna_cleavage_detail %>%
    group_by(Sample, Group, Direction, Origin, Cleavage) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction, Origin) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_cleavage <- sample_tsrna_cleavage %>%
    summarise_group_mean_percent(c("Direction", "Origin", "Cleavage"))
  group_tsrna_cleavage_view <- group_tsrna_cleavage %>% select(-Direction)
  sample_tsrna_aa <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, AminoAcid) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_aa <- sample_tsrna_aa %>%
    summarise_group_mean_percent(c("Direction", "AminoAcid"))
  group_tsrna_aa_view <- group_tsrna_aa %>% select(-Direction)
  sample_tsrna_anticodon <- tsrna_base_detail %>%
    group_by(Sample, Group, Direction, Anticodon) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample, Group, Direction) %>%
    mutate(Percent = group_percent(Reads)) %>%
    ungroup()
  group_tsrna_anticodon <- sample_tsrna_anticodon %>%
    summarise_group_mean_percent(c("Direction", "Anticodon"))
  group_tsrna_anticodon_view <- group_tsrna_anticodon %>% select(-Direction)

  # output_all 保留逐条注释粒度，主要用于 top feature / 详细条目类分析。
  output_all <- output_all %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(
      Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", Annotation),
      PrimaryAnnotation = str_split_fixed(Annotation, ";", 2)[, 1],
      FeatureCategory = if_else(EndoCandidate, "endo-siRNA", annotation_main_category(PrimaryAnnotation)),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)
    )
  group_top_features <- bind_rows(
    summary_all %>% filter(SubClass != "-", !is_tsrna_nested_class(Class), !Class %in% c("Clean_Reads","Match_Genome")) %>%
      transmute(Sample,Group,PrimaryAnnotation=SubClass,FeatureCategory=as.character(MainCategory),Reads),
    output_all %>% filter(EndoCandidate) %>% transmute(Sample,Group,PrimaryAnnotation=paste0("endo-siRNA:",Sequence),FeatureCategory="endo-siRNA",Reads)
  ) %>%
    group_by(Sample, Group, PrimaryAnnotation, FeatureCategory) %>% summarise(Reads=sum(Reads),.groups="drop") %>%
    left_join(clean_reads_tbl,by="Sample") %>% mutate(RPM=1e6*Reads/CleanReads) %>%
    summarise_group_mean_rpm(c("PrimaryAnnotation","FeatureCategory")) %>%
    filter(RPM > 0 | Reads > 0) %>% arrange(Group,desc(RPM),desc(Reads),PrimaryAnnotation) %>%
    group_by(Group) %>% mutate(Rank=row_number()) %>% filter(Rank<=20) %>% ungroup()

  # 每样本 cel-miR-39 外参检出表（模块 13 使用）。
  # 基于 SPORTS output 矩阵做精确序列 + 编辑距离 1 的替换/插入/缺失近邻计数。
  # 注意：spike-in 序列无法比对到人源参考，SPORTS 会把它归到 NO_Annotation/other。
  compute_spikein_counts <- function(output_all_df, sample_info_df) {
    candidate_set <- spikein_edit1_sequences(SPIKEIN_SEQ)
    si_seqs <- output_all_df %>%
      filter(
        Sample %in% as.character(sample_info_df$Sample),
        Sequence %in% candidate_set
      )
    exact <- si_seqs %>%
      filter(Sequence == SPIKEIN_SEQ) %>%
      group_by(Sample) %>%
      summarise(SpikeinExactReads = sum(Reads, na.rm = TRUE), .groups = "drop")
    near <- si_seqs %>%
      filter(Sequence != SPIKEIN_SEQ) %>%
      group_by(Sample) %>%
      summarise(SpikeinNearReads = sum(Reads, na.rm = TRUE), .groups = "drop")
    sample_info_df %>%
      transmute(Sample = as.character(Sample), Group = as.character(Group), BiologicalSample = BiologicalSample, Spikein = Spikein) %>%
      left_join(clean_reads_tbl, by = "Sample") %>%
      left_join(exact, by = "Sample") %>%
      left_join(near, by = "Sample") %>%
      mutate(
        SpikeinExactReads = if_else(is.na(SpikeinExactReads), 0, SpikeinExactReads),
        SpikeinNearReads = if_else(is.na(SpikeinNearReads), 0, SpikeinNearReads),
        SpikeinReads = SpikeinExactReads + SpikeinNearReads,
        SpikeinRPM = if_else(CleanReads > 0, 1e6 * SpikeinExactReads / CleanReads, 0)
      )
  }
  spikein_counts_by_sample <- compute_spikein_counts(output_all, sample_info)

  # legacy_pdf_index 用来保留 SPORTS 原生 PDF 的索引，便于在结果表里回链查看。
  ratio_all_sample_df <- sample_comp_view %>%
    filter(Reads > 0) %>%
    mutate(MainCategory = factor(as.character(MainCategory), levels = category_levels))
  legacy_pdf_index <- sample_info %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      sncRNA_distribution = PDF_sncRNA,
      rsRNA_distribution = PDF_rsRNA,
      rRNA_mapping = PDF_rRNA,
      tRNA_mature_mapping = PDF_tRNA
    )

  # keep_names 决定哪些对象会暴露给后续模块环境。
  # 这样可以明确共享边界，避免把太多临时变量泄露出去。
  keep_names <- c(
    "sample_info", "analysis_scope", "group_levels", "group_palette", "endo_evidence", "endo_transfer_audit", "endo_reference_audit",
    "summary_all_raw", "summary_all", "summary_all_view", "summary_top", "summary_top_view",
    "length_all_raw", "length_all", "length_all_view",
    "clean_reads_tbl", "match_reads_tbl", "group_clean_reads_tbl", "group_match_tbl",
    "group_comp_total", "group_comp_view", "group_comp_direction_view",
    "sample_comp_view", "composition_integrity_check", "sample_qc_table", "group_qc_table",
    "length_detail", "length_detail_view", "length_focus_categories", "sample_length_display", "sample_length_view", "group_length_display", "group_length_view",
    "rsrna_detail", "sample_rsrna_visual_by_view", "sample_rsrna_visual_view", "group_rsrna_visual_view", "group_rsrna_visual_merged", "group_rsrna_view",
    "tsrna_base_detail", "tsrna_cleavage_detail", "tsrna_cleavage_reconstruction_check",
    "group_tsrna_origin", "group_tsrna_origin_view", "group_tsrna_cleavage", "group_tsrna_cleavage_view", "group_tsrna_aa_view", "group_tsrna_anticodon_view",
    "output_all", "group_top_features", "ratio_all_sample_df", "legacy_pdf_index",
    "spikein_counts_by_sample"
  )
  mget(keep_names, inherits = FALSE)
}

activate_serum_context <- function(env = parent.frame(), report_view = REPORT_VIEW, analysis_scope = "all") {
  if (PANDORA_DESIGN != "method_capture") analysis_scope <- "all"
  list2env(load_serum_context(report_view = report_view, analysis_scope = analysis_scope), envir = env)
  assign("PANDORA_ANALYSIS_SCOPE",analysis_scope,envir=env)
  invisible(env)
}

build_sector_radar_df <- function(input_df, radar_max, sector_levels=radar_levels, sector_labels=radar_labels, sector_palette=radar_sector_palette, bottom_anchor_category="rsRNA", angle_offset=0) {
  n <- length(sector_levels); sw <- 2*pi/n
  bai <- match(bottom_anchor_category,sector_levels); if(is.na(bai)) bai <- ceiling(n/2)
  label_palette <- make_sector_label_palette(sector_palette[sector_levels], amount = 0.28)
  input_df %>% mutate(PieCategory=factor(as.character(PieCategory),levels=sector_levels)) %>%
    complete(PieCategory=factor(sector_levels,levels=sector_levels),fill=list(Reads=0,Percent=0)) %>% arrange(PieCategory) %>%
    mutate(SectorId=row_number(), CenterAngle=-pi/2+(bai-SectorId)*sw+angle_offset, StartAngle=CenterAngle+sw/2, EndAngle=CenterAngle-sw/2,
           Radius=Percent/radar_max, DisplayLabel=unname(sector_labels[as.character(PieCategory)]),
           DisplayValue=sprintf("%.2f",Percent), ValueColor=serum_text_primary, LabelColor=unname(label_palette[as.character(PieCategory)]),
           UpperHalf=sin(CenterAngle)>=0,
           TextAngle=((if_else(UpperHalf,CenterAngle*180/pi-90,CenterAngle*180/pi+90)+180)%%360)-180)
}

plot_sector_radar_core <- function(plot_df, title_text, radar_max, radar_breaks, sector_levels=radar_levels, sector_palette=radar_sector_palette,
  axis_label="Relative Abundance (%)", label_radius=1.24, value_radius=NULL, coord_xlim=c(-1.92,1.72), coord_ylim=c(-1.54,1.54), label_size=3.35, value_size=3.55) {
  rir<-1.02; ror<-1.16; if(is.null(value_radius)) value_radius<-0.93
  sw <- 2*pi/length(sector_levels)
  cdf <- build_circle_path_df(radar_breaks[radar_breaks>0]/radar_max)
  wdf <- build_sector_polygon_df(plot_df,ir=0,or=NA_real_,use_data=TRUE)
  rdf <- build_sector_polygon_df(plot_df,ir=rir,or=ror)
  sdf <- tibble(x=0,y=0,xend=ror*cos(plot_df$StartAngle),yend=ror*sin(plot_df$StartAngle))
  ldf <- build_sector_arc_text_df(plot_df,"DisplayLabel",label_radius,sw,"LabelColor")
  vdf <- build_sector_text_df(plot_df,"DisplayValue",value_radius,"ValueColor")
  ggplot() +
    geom_path(data=cdf,aes(x=x,y=y,group=r),color=serum_neutral_line,linewidth=0.55) +
    geom_polygon(data=wdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.4,linejoin="round",alpha=0.86) +
    geom_polygon(data=rdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.45,linejoin="round",alpha=1) +
    geom_segment(data=sdf,aes(x=x,y=y,xend=xend,yend=yend),color=serum_neutral_line,linewidth=0.75,lineend="round") +
    geom_text(data=ldf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=label_size,fontface="bold",inherit.aes=FALSE) +
    geom_text(data=vdf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=value_size,fontface="bold",inherit.aes=FALSE) +
    annotate("text",x=-1.60,y=0,label=axis_label,angle=90,size=4.0,color=serum_axis_label_color,fontface="bold") +
    scale_fill_manual(values=sector_palette,drop=FALSE) + scale_color_identity() +
    coord_fixed(xlim=coord_xlim,ylim=coord_ylim,clip="off") + theme_void() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5,face="bold",size=15,color=serum_text_muted),
          plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
          plot.margin=margin(10,20,10,28)) + labs(title=title_text)
}

prepare_group_sector_radar <- function(input_df, radar_max, sector_levels = radar_levels, sector_labels = radar_labels,
  bottom_anchor_category = "rsRNA", sector_palette = radar_sector_palette, angle_offset = 0,
  value_digits = 2, min_value_label = 0) {
  n_sector <- length(sector_levels)
  sector_width <- 2 * pi / n_sector
  bai <- match(bottom_anchor_category, sector_levels)
  if (is.na(bai)) bai <- ceiling(n_sector / 2)
  label_palette <- make_sector_label_palette(sector_palette[sector_levels], amount = 0.28)

  sector_df <- tibble(PieCategory = factor(sector_levels, levels = sector_levels)) %>%
    mutate(
      SectorId = row_number(),
      SectorCenterAngle = -pi / 2 + (bai - SectorId) * sector_width + angle_offset,
      SectorStartAngle = SectorCenterAngle + sector_width / 2,
      SectorEndAngle = SectorCenterAngle - sector_width / 2,
      CenterAngle = SectorCenterAngle,
      StartAngle = SectorStartAngle,
      EndAngle = SectorEndAngle,
      DisplayLabel = unname(sector_labels[as.character(PieCategory)]),
      LabelColor = unname(label_palette[as.character(PieCategory)]),
      UpperHalf = sin(CenterAngle) >= 0,
      TextAngle = ((if_else(UpperHalf, CenterAngle * 180 / pi - 90, CenterAngle * 180 / pi + 90) + 180) %% 360) - 180
    )

  radar_df <- expand_grid(
    Group = factor(group_levels, levels = group_levels),
    PieCategory = factor(sector_levels, levels = sector_levels)
  ) %>%
    left_join(
      input_df %>%
        mutate(
          Group = factor(as.character(Group), levels = group_levels),
          PieCategory = factor(as.character(PieCategory), levels = sector_levels)
        ) %>%
        select(Group, PieCategory, Percent),
      by = c("Group", "PieCategory")
    ) %>%
    mutate(Percent = replace_na(Percent, 0)) %>%
    left_join(sector_df, by = "PieCategory") %>%
    mutate(
      GroupIndex = match(as.character(Group), group_levels),
      GroupSliceWidth = sector_width / length(group_levels),
      StartAngle = SectorStartAngle - (GroupIndex - 1) * GroupSliceWidth,
      EndAngle = SectorStartAngle - GroupIndex * GroupSliceWidth,
      CenterAngle = (StartAngle + EndAngle) / 2,
      Radius = Percent / radar_max,
      DisplayValue = if_else(Percent >= min_value_label, sprintf(paste0("%.", value_digits, "f"), Percent), ""),
      ValueColor = serum_text_primary,
      UpperHalf = sin(CenterAngle) >= 0,
      TextAngle = ((if_else(UpperHalf, CenterAngle * 180 / pi - 90, CenterAngle * 180 / pi + 90) + 180) %% 360) - 180
    )

  list(sectors = sector_df, slices = radar_df)
}

plot_group_sector_radar <- function(radar_data, title_text, radar_max, radar_breaks, axis_label = "Relative Abundance (%)",
  label_radius = 1.24, value_radius = 0.965, coord_xlim = c(-1.92, 1.72), coord_ylim = c(-1.54, 1.54),
  label_size = 3.35, value_size = 2.7, sector_palette = radar_sector_palette) {
  sector_df <- radar_data$sectors
  slice_df <- radar_data$slices
  ring_inner <- 1.02
  ring_outer <- 1.16

  cdf <- build_circle_path_df(radar_breaks[radar_breaks > 0] / radar_max)
  wedge_df <- build_sector_polygon_df(
    slice_df %>% mutate(SectorId = paste0(PieCategory, "_", Group), PieCategory = Group),
    ir = 0,
    or = NA_real_,
    use_data = TRUE
  ) %>%
    mutate(FillColor = unname(group_palette[as.character(PieCategory)]))
  ring_df <- build_sector_polygon_df(
    sector_df %>% mutate(Radius = ring_outer),
    ir = ring_inner,
    or = ring_outer
  )
  ring_df <- ring_df %>% mutate(FillColor = unname(sector_palette[as.character(PieCategory)]))
  sdf <- tibble(
    x = 0,
    y = 0,
    xend = ring_outer * cos(sector_df$SectorStartAngle),
    yend = ring_outer * sin(sector_df$SectorStartAngle)
  )
  ldf <- build_sector_arc_text_df(sector_df, "DisplayLabel", label_radius, 2 * pi / nrow(sector_df), "LabelColor")
  vdf <- build_sector_text_df(slice_df, "DisplayValue", value_radius, "ValueColor")

  ggplot() +
    geom_path(data = cdf, aes(x = x, y = y, group = r), color = serum_neutral_line, linewidth = 0.55) +
    geom_polygon(data = wedge_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.38, linejoin = "round", alpha = 0.86) +
    geom_polygon(data = ring_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.45, linejoin = "round", alpha = 1) +
    geom_segment(data = sdf, aes(x = x, y = y, xend = xend, yend = yend), color = serum_neutral_line, linewidth = 0.75, lineend = "round") +
    geom_text(data = ldf, aes(x = TextX, y = TextY, label = Label, angle = TextAngle, color = TextColor), size = label_size, fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = vdf %>% filter(Label != ""), aes(x = TextX, y = TextY, label = Label, angle = TextAngle), size = value_size, fontface = "bold", inherit.aes = FALSE) +
    annotate("text", x = -1.60, y = 0, label = axis_label, angle = 90, size = 4.0, color = serum_axis_label_color, fontface = "bold") +
    scale_fill_identity(
      guide = "legend",
      breaks = unname(group_palette[group_levels]),
      labels = group_levels,
      name = NULL
    ) +
    scale_color_identity() +
    coord_fixed(xlim = coord_xlim, ylim = coord_ylim, clip = "off") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, color = serum_text_muted),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 20, 10, 28)
    ) +
    labs(title = title_text, fill = NULL)
}

build_trend_bar_plot <- function(df, levels, value_col, title, y_label, group_levels, fill_palette = main_palette) {
  plot_df <- prepare_trend_plot_df(df, levels, group_levels)
  display_df <- plot_df %>% filter(.data[[value_col]] > 0)
  if (nrow(display_df) == 0) return(plot_note(title, "No trend data available."))
  p <- ggplot(display_df, aes(x = Group, y = .data[[value_col]], fill = TrendCategory)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = fill_palette[levels], drop = TRUE) +
    labs(x = "Group", y = y_label, title = title) +
    serum_theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
    guides(fill = guide_legend(ncol = if (any(nchar(levels) > 22, na.rm = TRUE)) 2 else 4, byrow = TRUE))
  if (identical(value_col, "Reads")) p <- p + scale_y_continuous(labels = scales::comma)
  p
}

build_two_group_summary_bar <- function(sample_df, value_col, title, y_label, group_levels, fill_palette = group_palette,
                                        sample_info = sample_info, comma_labels = FALSE, label_digits = 1,
                                        show_pvalue = FALSE) {
  plot_df <- complete_sample_metric(sample_df, value_col, sample_info = sample_info, group_levels = group_levels)
  if (nrow(plot_df) == 0) return(plot_note(title, "No sample-level data available."))

  summary_df <- plot_df %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = stats::sd(Value, na.rm = TRUE),
      N = sum(is.finite(Value)),
      SEM = if_else(N > 1, SD / sqrt(N), 0),
      .groups = "drop"
    ) %>%
    mutate(
      SD = replace_na(SD, 0),
      SEM = replace_na(SEM, 0),
      Group = factor(as.character(Group), levels = group_levels)
    )

  descriptive_singletons <- nrow(summary_df) > 0 && all(summary_df$N == 1)
  if (descriptive_singletons) {
    y_label <- sub("^Mean ", "Observed ", y_label)
  }

  max_y <- max(c(summary_df$Mean + summary_df$SEM, plot_df$Value), na.rm = TRUE)
  if (!is.finite(max_y)) max_y <- 0
  y_pad <- if (max_y > 0) max_y * 0.12 else 1
  bracket_y <- max_y + y_pad
  tip_y <- max_y + y_pad * 0.45
  label_y <- bracket_y + y_pad * 0.35
  upper_y <- label_y + y_pad * 0.45
  pvalue <- calculate_two_group_pvalue(plot_df, value_col = "Value", group_col = "Group", group_levels = group_levels)
  p_label <- format_pvalue_label(pvalue, show_pvalue = show_pvalue)

  bar_plot <- ggplot(summary_df, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.64, color = "white", linewidth = 0.25)
  if (!descriptive_singletons) {
    bar_plot <- bar_plot +
      geom_errorbar(aes(ymin = pmax(Mean - SEM, 0), ymax = Mean + SEM), width = 0.16, linewidth = 0.45, color = serum_text_primary)
  }
  bar_plot <- bar_plot +
    geom_point(
      data = plot_df,
      aes(x = Group, y = Value),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.065, height = 0, seed = 11),
      size = 2.2,
      shape = 21,
      stroke = 0.25,
      color = serum_text_primary,
      fill = "white",
      alpha = 0.9
    ) +
    geom_text(
      aes(label = sprintf(paste0("%.", label_digits, "f"), Mean), y = Mean + SEM + y_pad * 0.12),
      size = 3.2,
      color = serum_text_primary,
      vjust = 0
    ) +
    scale_fill_manual(values = fill_palette, drop = FALSE) +
    scale_y_continuous(
      labels = if (comma_labels) scales::comma else waiver(),
      limits = c(0, upper_y),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL, y = y_label, title = title,
      subtitle = if (descriptive_singletons) "One library per group; descriptive only" else NULL
    ) +
    serum_theme_bw() +
    theme(legend.position = "none")

  if (length(group_levels) == 2 && !is.na(p_label)) {
    bracket_df <- tibble(
      x = 1,
      xend = 2,
      y = bracket_y,
      yend = bracket_y,
      label = p_label,
      label_x = 1.5,
      label_y = label_y
    )
    tick_df <- tibble(x = c(1, 2), xend = c(1, 2), y = c(tip_y, tip_y), yend = c(bracket_y, bracket_y))
    bar_plot <- bar_plot +
      geom_segment(data = bracket_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.45) +
      geom_segment(data = tick_df, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.45) +
      geom_text(data = bracket_df, aes(x = label_x, y = label_y, label = label), inherit.aes = FALSE, size = 4.0, fontface = "bold")
  }

  bar_plot
}

build_trend_alluvial_plot <- function(df, levels, value_col, title, y_label, group_levels, fill_palette = main_palette) {
  plot_df <- prepare_trend_plot_df(df, levels, group_levels)
  display_df <- plot_df %>% filter(.data[[value_col]] > 0)
  if (!pkg_has("ggalluvial")) return(plot_note(title, "ggalluvial package is not available."))
  if (nrow(display_df) == 0) return(plot_note(title, "No trend data available."))
  # 这里显式压缩组轴间距，让“柱宽 = 柱间距 = 0.1”。
  # flow 宽度与 stratum 宽度保持一致，这样流带会从柱子的边界出发，到下一根柱子的边界结束。
  stratum_width <- 0.1
  flow_width <- stratum_width
  group_axis_spacing <- stratum_width * 2
  group_axis_positions <- seq(1, by = group_axis_spacing, length.out = length(group_levels))
  names(group_axis_positions) <- group_levels
  display_df <- display_df %>%
    mutate(GroupX = unname(group_axis_positions[as.character(Group)]))

  left_group <- group_levels[[1]]
  if (identical(value_col, "Reads")) {
    label_layer <- geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(scales::comma(count))),
      size = 3,
      color = serum_text_primary,
      check_overlap = TRUE
    )
    axis_label_layer <- NULL
    y_scale <- scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.03, 0.03)))
  } else {
    label_layer <- geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(ifelse(count > 3.5, sprintf("%.1f%%", count), ""))),
      size = 3,
      color = serum_text_primary,
      check_overlap = FALSE
    )
    # 左侧柱子上的刻度标注 0/25/50/75/100
    axis_label_df <- data.frame(
      x = group_axis_positions[[1]],
      y = c(0, 25, 50, 75, 100),
      label = c("0", "25", "50", "75", "100")
    )
    axis_label_layer <- geom_text(
      data = axis_label_df,
      aes(x = x - stratum_width / 2 - 0.02, y = y, label = label,
          vjust = ifelse(y == 0, 0, ifelse(y == 100, 1, 0.5))),
      size = 3, color = serum_text_primary, hjust = 1, inherit.aes = FALSE
    )
    y_scale <- scale_y_continuous(
      breaks = seq(0, 100, 25),
      labels = NULL,
      expand = expansion(mult = c(0.03, 0.03))
    )
  }

  p <- ggplot(display_df, aes(x = GroupX, y = .data[[value_col]], stratum = TrendCategory, alluvium = TrendCategory, fill = TrendCategory)) +
    ggalluvial::geom_flow(width = flow_width, alpha = 0.32) +
    ggalluvial::geom_stratum(width = stratum_width, color = "white", linewidth = 0.25) +
    label_layer +
    axis_label_layer +
    scale_fill_manual(values = fill_palette[levels], drop = TRUE) +
    scale_x_continuous(
      breaks = unname(group_axis_positions),
      labels = names(group_axis_positions),
      expand = expansion(add = stratum_width / 2)
    ) +
    y_scale +
    labs(x = "Group", y = y_label, title = title) +
    serum_theme_bw() + theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = if (any(nchar(levels) > 22, na.rm = TRUE)) 2 else 4, byrow = TRUE))
  p
}

build_capture_volcano_plot <- function(de_tbl, title_text, lfc_cutoff = log2(1.5), pvalue_cutoff = 0.05, label_n = 15) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No capture results available."))

  significance_values <- resolve_de_significance_values(de_tbl)
  significance_methods <- if ("SignificanceMethod" %in% names(de_tbl)) as.character(de_tbl$SignificanceMethod) else character(0)
  significance_methods <- significance_methods[!is.na(significance_methods) & nzchar(significance_methods)]
  significance_label <- if (length(significance_methods) > 0) {
    unique(significance_methods)[[1]]
  } else {
    de_significance_label()
  }

  plot_df <- de_tbl %>%
    mutate(
      SignificanceMetric = significance_values,
      direction = classify_de_direction(SignificanceMetric, logFC, pvalue_cutoff = pvalue_cutoff, lfc_cutoff = lfc_cutoff),
      label_score = abs(coalesce(logFC, 0)) * -log10(coalesce(SignificanceMetric, 1) + 1e-300),
      label = if_else(direction != "ns" & rank(-label_score, ties.method = "first") <= label_n, as.character(Feature), NA_character_)
    )

  p <- ggplot(plot_df, aes(x = logFC, y = -log10(SignificanceMetric + 1e-300), color = direction)) +
    geom_point(alpha = 0.75, size = 1.6) +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = serum_neutral_text) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = serum_neutral_text) +
    scale_color_manual(values = serum_direction_palette, drop = FALSE) +
    labs(x = "log2 Fold Change", y = paste0("-log10(", significance_label, ")"), title = title_text) +
    serum_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  if (pkg_has("ggrepel")) {
    p <- p + ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, na.rm = TRUE, show.legend = FALSE)
  } else {
    p <- p + geom_text(aes(label = label), size = 3, na.rm = TRUE, show.legend = FALSE, vjust = -0.4)
  }

  p
}

build_capture_ma_plot <- function(de_tbl, title_text) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No capture results available."))

  ggplot(de_tbl, aes(x = log10(baseMean + 1), y = logFC, color = direction)) +
    geom_point(alpha = 0.72, size = 1.5) +
    scale_color_manual(values = serum_direction_palette, drop = FALSE) +
    labs(x = "log10(baseMean + 1)", y = "log2 Fold Change", title = title_text) +
    serum_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")
}

write_setup_provenance <- function(ctx, output_dir = SETUP_DIR) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  tool_manifest <- bind_rows(
    tibble(tool = "R", available = TRUE, path = R.home(), version = paste(R.version.string, collapse = " ")),
    capture_command_version("Rscript", "--version"),
    capture_command_version("perl", "-v"),
    capture_command_version("python3", "--version"),
    capture_command_version("cutadapt", "--version"),
    capture_command_version("fastqc", "--version"),
    capture_command_version("multiqc", "--version"),
    capture_command_version("bowtie", "--version"),
    capture_command_version("bowtie2", "--version"),
    capture_command_version("samtools", "--version")
  )

  input_audit <- ctx$sample_info %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      BiologicalSample = BiologicalSample,
      Method = Method,
      Replicate = Replicate,
      Batch = Batch,
      PairID = PairID,
      AnalysisRole = AnalysisRole,
      Spikein = Spikein,
      SpikeinSequence = SpikeinSequence,
      SpikeinInputFmol = SpikeinInputFmol,
      Condition = Condition,
      SampleType = SampleType,
      Species = Species,
      Matrix = Matrix,
      RawFastqOriginal = RawFastqOriginal,
      RawFastq = RawFastq,
      RawFastqExists = RawFastqExists,
      SummaryFile = SummaryFile,
      SummaryFileExists = file.exists(SummaryFile),
      LengthFile = LengthFile,
      LengthFileExists = file.exists(LengthFile),
      OutputFile = OutputFile,
      OutputFileExists = file.exists(OutputFile)
    )

  sports_script <- file.path(Sys.getenv("PANDORA_SHARED_ROOT"), "upstream", paste0("sports_",tolower(Sys.getenv("PANDORA_SPECIES")),".sh"))
  sports_config <- if (file.exists(sports_script)) {
    script_lines <- readLines(sports_script, warn = FALSE)
    config_keys <- c(
      "DB_ROOT", "GENOME_BASE", "MIRNA_BASE", "RRNA_BASE", "TRNA_BASE",
      "PIRNA_BASE", "ENSEMBL_BASE", "RFAM_BASE",
      "ADAPTER3", "TRIM_ADAPTER", "MISMATCH", "RUN_FASTQC",
      "THREADS", "PARALLEL_JOBS", "MIN_LEN", "MAX_LEN"
    )
    purrr::map_dfr(config_keys, function(key) {
      hit <- script_lines[str_detect(script_lines, paste0("^", key, "="))][1]
      tibble(key = key, script_value = hit %||% "", environment_value = Sys.getenv(key, unset = ""))
    })
  } else {
    tibble(key = character(), script_value = character(), environment_value = character())
  }

  write_tsv_safe(build_runtime_config_table(), file.path(output_dir, "analysis_runtime_config.tsv"))
  write_tsv_safe(tool_manifest, file.path(output_dir, "tool_version_manifest.tsv"))
  write_tsv_safe(input_audit, file.path(output_dir, "input_path_audit.tsv"))
  write_tsv_safe(sports_config, file.path(output_dir, "sports_pipeline_config.tsv"))
  write_tsv_safe(
    ctx$sample_info %>%
      transmute(
        Sample = as.character(Sample),
        BiologicalSample = as.character(BiologicalSample),
        Method = as.character(Method),
        Group = as.character(Group),
        Replicate = as.character(Replicate),
        Batch = as.character(Batch),
        PairID = as.character(PairID),
        AnalysisRole,
        Spikein,
        SpikeinSequence,
        SpikeinInputFmol
      ),
    file.path(output_dir, "analysis_cohort_manifest.tsv")
  )
  utils::write.csv(
    ctx$sample_comp_view %>%
      transmute(
        Sample = as.character(Sample), Group = as.character(Group),
        MainCategory = as.character(MainCategory), Reads, CleanReads, Percent
      ),
    file.path(output_dir, "exclusive_class_composition.csv"), row.names = FALSE
  )
  utils::write.csv(
    ctx$composition_integrity_check,
    file.path(output_dir, "composition_integrity_check.csv"), row.names = FALSE
  )
  utils::write.csv(
    ctx$tsrna_cleavage_reconstruction_check,
    file.path(output_dir, "tsrna_cleavage_reconstruction_check.csv"), row.names = FALSE
  )
  invisible(output_dir)
}

write_serum_setup_expression_matrices <- function(report_view = REPORT_VIEW) {
  ctx <- load_serum_context(report_view = report_view)
  sample_order <- as.character(ctx$sample_info$Sample)
  output_dir <- EXPRESSION_MATRIX_DIR
  seq_dir <- file.path(output_dir, "seq_count")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)
  write_setup_provenance(ctx, output_dir = SETUP_DIR)

  spikein_tab_path <- file.path(SETUP_DIR, "spikein_counts_by_sample.tsv")
  dir.create(SETUP_DIR, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(ctx$spikein_counts_by_sample, spikein_tab_path, sep = "\t", row.names = FALSE, quote = FALSE)

  non_tsrna_annotation <- ctx$summary_all %>%
    filter(SubClass != "-", MainCategory != "tsRNA", !Class %in% c("Clean_Reads", "Match_Genome")) %>%
    transmute(
      Sample = as.character(Sample),
      Feature = if_else(is.na(SubClass) | SubClass == "", Class, SubClass),
      FeatureCategory = as.character(MainCategory),
      Reads = Reads
    )

  tsrna_annotation <- ctx$tsrna_cleavage_detail %>%
    transmute(
      Sample = as.character(Sample),
      Feature = make_company_tsrna_feature(SubClass, Cleavage),
      FeatureCategory = if_else(Origin == "Mitochondrial", "mt_tsRNA", "tsRNA"),
      Reads = Reads
    )

  annotation_input <- bind_rows(non_tsrna_annotation, tsrna_annotation)
  annotation_sets <- list(
    all = annotation_input,
    miRNA = annotation_input %>% filter(FeatureCategory == "miRNA"),
    tsRNA = annotation_input %>% filter(FeatureCategory == "tsRNA"),
    mt_tsRNA = annotation_input %>% filter(FeatureCategory == "mt_tsRNA")
  )
  walk(names(annotation_sets), ~write_annotation_expression_set(annotation_sets[[.x]], .x, output_dir, sample_order, ctx$clean_reads_tbl))

  if (pkg_has("data.table")) {
    sequence_input <- data.table::as.data.table(ctx$output_all)
    sequence_input <- sequence_input[!is.na(Sequence) & Sequence != ""]
    sequence_input[, FeatureCategory := as.character(FeatureCategory)]
    sequence_input[, PrimaryAnnotation := as.character(PrimaryAnnotation)]
    sequence_input[is.na(PrimaryAnnotation) | PrimaryAnnotation == "", PrimaryAnnotation := "NO_Annotation"]
    sequence_input[, SequenceSet := FeatureCategory]
    sequence_input[
      FeatureCategory == "tsRNA" & grepl("mt_tRNA|mitochond", PrimaryAnnotation, ignore.case = TRUE),
      SequenceSet := "mt_tsRNA"
    ]
    sequence_input <- sequence_input[, .(Sample, Sequence, Reads, RPM, PrimaryAnnotation, OriginalAnnotation, FeatureCategory, SequenceSet)]
    sequence_sets <- list(
      all = sequence_input,
      miRNA = sequence_input[SequenceSet == "miRNA"],
      tsRNA = sequence_input[SequenceSet == "tsRNA"],
      mt_tsRNA = sequence_input[SequenceSet == "mt_tsRNA"],
      piRNA = sequence_input[SequenceSet == "piRNA"],
      rsRNA = sequence_input[SequenceSet == "rsRNA"],
      `endo-siRNA` = sequence_input[SequenceSet == "endo-siRNA"]
    )
  } else {
    sequence_input <- ctx$output_all %>%
      filter(!is.na(Sequence), Sequence != "") %>%
      mutate(
        FeatureCategory = as.character(FeatureCategory),
        PrimaryAnnotation = if_else(is.na(PrimaryAnnotation) | PrimaryAnnotation == "", "NO_Annotation", PrimaryAnnotation),
        IsMtTsRNA = FeatureCategory == "tsRNA" & str_detect(PrimaryAnnotation, regex("mt_tRNA|mitochond", ignore_case = TRUE)),
        SequenceSet = case_when(
          IsMtTsRNA ~ "mt_tsRNA",
          TRUE ~ FeatureCategory
        )
      ) %>%
      select(Sample, Sequence, Reads, RPM, PrimaryAnnotation, OriginalAnnotation, FeatureCategory, SequenceSet)

    sequence_sets <- list(
      all = sequence_input,
      miRNA = sequence_input %>% filter(SequenceSet == "miRNA"),
      tsRNA = sequence_input %>% filter(SequenceSet == "tsRNA"),
      mt_tsRNA = sequence_input %>% filter(SequenceSet == "mt_tsRNA"),
      piRNA = sequence_input %>% filter(SequenceSet == "piRNA"),
      rsRNA = sequence_input %>% filter(SequenceSet == "rsRNA"),
      `endo-siRNA` = sequence_input %>% filter(SequenceSet == "endo-siRNA")
    )
  }
  if (!ENDO_ENABLED) sequence_sets[["endo-siRNA"]] <- NULL
  walk(names(sequence_sets), ~write_sequence_expression_set(sequence_sets[[.x]], .x, output_dir, seq_dir, sample_order, duplicate_root = .x == "all"))

  invisible(output_dir)
}

