# 小鼠设计适配；共享函数不在项目内复制。
source(file.path(Sys.getenv("PANDORA_SHARED_ROOT"),"R","common.R"),local=TRUE)

knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.width = 11, fig.height = 7, dpi = 300)

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(knitr)
})

ROOT_DIR <- normalizePath(Sys.getenv("PANDORA_PROJECT_ROOT"), mustWork = TRUE)

RAW_DIR           <- file.path(ROOT_DIR, "01rawdata")

SPORTS_RESULT_DIR <- file.path(ROOT_DIR, "sportsresult")


SPORTS_METADATA_DIR <- file.path(SPORTS_RESULT_DIR, "02_metadata")

SPORTS_QC_DIR <- file.path(SPORTS_RESULT_DIR, "03_qc")

SPORTS_ANNOTATION_DIR <- file.path(SPORTS_RESULT_DIR, "04_sports_annotation")

SPORTS_SAMPLE_ROOT <- file.path(SPORTS_ANNOTATION_DIR, "samples")

SPORTS_ENDOSIRNA_DIR <- file.path(SPORTS_RESULT_DIR, "05_endosirna")

SPORTS_ENDOSIRNA_SUMMARY_DIR <- file.path(SPORTS_ENDOSIRNA_DIR, "summary")

SPORTS_LOG_DIR <- file.path(SPORTS_RESULT_DIR, "06_logs")

METADATA_FILE     <- file.path(ROOT_DIR, "config", "samples.tsv")

RESULT_DIR        <- Sys.getenv("PANDORA_RESULT_DIR", file.path(ROOT_DIR, "result", "describe"))

INTERNAL_CACHE_DIR <- file.path(ROOT_DIR, ".cache", "shared_analysis")

PROJECT_SPECIES <- "Mus musculus"

PROJECT_HOST_LABEL <- "mouse"

PROJECT_STAGE <- Sys.getenv("PROJECT_STAGE")

PROJECT_CONTRAST_BASELINE <- Sys.getenv("CONTRAST_BASELINE")

PROJECT_CONTRAST_CASE <- Sys.getenv("CONTRAST_CASE")

PROJECT_GROUP_LABELS <- c(MYMII = "mouse-young-MII", MOMII = "mouse-old-MII")

MIRNA_SPECIES <- "Mus musculus"

MIRNA_PREFIX <- "mmu"

MULTIMIR_ORG <- "mmu"

KEGG_ORGANISM <- "mmu"

MIRNA_TARGET_ORGDB_PACKAGE <- "org.Mm.eg.db"

MIRNA_TARGET_ORGDB_OBJECT <- "org.Mm.eg.db"

SETUP_DIR         <- file.path(RESULT_DIR, "1_setup")

EXPRESSION_MATRIX_DIR <- file.path(SETUP_DIR, "1_expression_matrix")

SINGLE_DIR        <- file.path(RESULT_DIR, "2_single_sample_qc")

SINGLE_FIG        <- file.path(SINGLE_DIR, "plots")

SINGLE_TAB        <- file.path(SINGLE_DIR, "tables")

PCA_DIR           <- file.path(RESULT_DIR, "3_cohort_qc_pca")

PCA_FIG           <- file.path(PCA_DIR, "plots")

PCA_TAB           <- file.path(PCA_DIR, "tables")

SNC_COMPARE_DIR   <- file.path(RESULT_DIR, "4_global_sncRNA")

SNC_COMPARE_FIG   <- file.path(SNC_COMPARE_DIR, "plots")

SNC_COMPARE_TAB   <- file.path(SNC_COMPARE_DIR, "tables")

TSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "5_tsrna")

TSRNA_COMPARE_FIG <- file.path(TSRNA_COMPARE_DIR, "plots")

TSRNA_COMPARE_TAB <- file.path(TSRNA_COMPARE_DIR, "tables")

RSRNA_COMPARE_DIR <- file.path(RESULT_DIR, "6_rsrna")

RSRNA_COMPARE_FIG <- file.path(RSRNA_COMPARE_DIR, "plots")

RSRNA_COMPARE_TAB <- file.path(RSRNA_COMPARE_DIR, "tables")

MIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "7_mirna")

MIRNA_COMPARE_FIG <- file.path(MIRNA_COMPARE_DIR, "plots")

MIRNA_COMPARE_TAB <- file.path(MIRNA_COMPARE_DIR, "tables")

PIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "8_pirna")

PIRNA_COMPARE_FIG <- file.path(PIRNA_COMPARE_DIR, "plots")

PIRNA_COMPARE_TAB <- file.path(PIRNA_COMPARE_DIR, "tables")

ENDOSIRNA_COMPARE_DIR <- file.path(RESULT_DIR, "9_endosirna")

ENDOSIRNA_COMPARE_FIG <- file.path(ENDOSIRNA_COMPARE_DIR, "plots")

ENDOSIRNA_COMPARE_TAB <- file.path(ENDOSIRNA_COMPARE_DIR, "tables")

RUN_FASTQC <- tolower(Sys.getenv("RUN_FASTQC", "auto"))

RUN_FASTQC <- if (RUN_FASTQC %in% c("auto", "1", "0", "true", "false", "yes", "no")) RUN_FASTQC else "auto"

RUN_EXOGENOUS <- tolower(Sys.getenv("RUN_EXOGENOUS", "0"))

RUN_EXOGENOUS <- if (RUN_EXOGENOUS %in% c("1", "true", "yes")) "1" else "0"

EXOGENOUS_DB_ROOT <- Sys.getenv("EXOGENOUS_DB_ROOT", "")

MIRNA_TARGET_SOURCE <- tolower(Sys.getenv("MIRNA_TARGET_SOURCE", "validated"))

MIRNA_TARGET_SOURCE <- if (MIRNA_TARGET_SOURCE %in% c("validated", "predicted", "msigdb_fallback")) MIRNA_TARGET_SOURCE else "validated"

DE_SIGNIFICANCE <- toupper(Sys.getenv("DE_SIGNIFICANCE", "FDR"))

DE_SIGNIFICANCE <- if (DE_SIGNIFICANCE %in% c("FDR", "PVALUE", "P_VALUE", "P")) DE_SIGNIFICANCE else "FDR"

RAW_AUDIT_MAX_READS <- suppressWarnings(as.integer(Sys.getenv("RAW_AUDIT_MAX_READS", "100000")))

if (!is.finite(RAW_AUDIT_MAX_READS) || RAW_AUDIT_MAX_READS <= 0) RAW_AUDIT_MAX_READS <- 100000L

RAW_AUDIT_KMER <- suppressWarnings(as.integer(Sys.getenv("RAW_AUDIT_KMER", "12")))

if (!is.finite(RAW_AUDIT_KMER) || RAW_AUDIT_KMER < 4) RAW_AUDIT_KMER <- 12L

RAW_AUDIT_REFRESH <- tolower(Sys.getenv("RAW_AUDIT_REFRESH", "0"))

RAW_AUDIT_REFRESH <- if (RAW_AUDIT_REFRESH %in% c("1", "true", "yes")) "1" else "0"

REPORT_VIEW <- "sense"

ALL_VIEW_LEVELS <- c("merged", "sense", "antisense")

VIEW_LEVELS <- c(REPORT_VIEW)

OUTPUT_FIG_DIR <- NULL

OUTPUT_TAB_DIR <- NULL

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

radar_levels  <- c("unanno_UMG","miRNA","piRNA","tsRNA","rsRNA","other","unanno_MG")

radar_labels  <- c(unanno_UMG="unanno UMG", miRNA="miRNA", piRNA="piRNA", tsRNA="tsRNA", rsRNA="rsRNA", other="other", unanno_MG="unanno MG")

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

mii_blue <- "#A1C9F4"

mii_orange <- "#FFBE7A"

mii_green <- "#8ECFC9"

mii_red <- "#FF9F9B"

mii_purple <- "#BEB8DC"

mii_neutral_line <- "#C9CED6"

mii_neutral_text <- "#9AA1AC"

mii_axis_label_color <- "#6D7784"

mii_text_primary <- "#1F2329"

mii_text_muted <- "#47596E"

mii_transition_yellow <- "#FDFFB6"

mii_qc_palette <- c(CleanReads = mii_blue, MatchGenome = mii_green)

mii_diverging_palette <- c(low = mii_blue, mid = mii_transition_yellow, high = mii_red)

mii_direction_palette <- c(up = mii_red, down = mii_blue, ns = mii_neutral_line)

mii_direction_binary_palette <- c(up = mii_red, down = mii_blue)

mii_target_group_palette <- c(all = mii_purple, up = mii_red, down = mii_blue)

mii_enrichment_gradient <- c(mii_blue, mii_transition_yellow, mii_red)

mii_heatmap_gradient <- c(mii_blue, mii_transition_yellow, mii_red)

load_metadata <- function() {
  src <- METADATA_FILE
  meta <- read_tsv_flex(src)
  if (nrow(meta)==0) stop("未检测到样本信息文件内容: ", src)
  names(meta) <- tolower(names(meta))
  if (!all(c("group","sample") %in% names(meta))) stop("样本信息文件必须至少包含 group 和 sample 列: ", src)
  if (!"raw_fastq_gz" %in% names(meta)) stop("样本信息文件必须包含 raw_fastq_gz 列: ", src)
  if (anyNA(meta$sample) || anyNA(meta$group) || any(!nzchar(trimws(as.character(meta$sample)))) || any(!nzchar(trimws(as.character(meta$group)))) || anyDuplicated(meta$sample)) stop("sample/group 不能为空，sample 必须唯一: ",src)
  out <- tibble(
    Group=trimws(as.character(meta$group)), Sample=trimws(as.character(meta$sample)),
    Condition=if("condition"%in%names(meta)) trimws(as.character(meta$condition)) else NA_character_,
    AgeGroup=if("age_group"%in%names(meta)) trimws(as.character(meta$age_group)) else rep("unspecified", nrow(meta)),
    Species=if("species"%in%names(meta)) trimws(as.character(meta$species)) else PROJECT_SPECIES,
    Stage=if("stage"%in%names(meta)) trimws(as.character(meta$stage)) else PROJECT_STAGE,
    Replicate=if("replicate"%in%names(meta)) trimws(as.character(meta$replicate)) else NA_character_,
    RawFastqOriginal=trimws(as.character(meta$raw_fastq_gz))
  ) %>% filter(Group!="", Sample!="") %>%
    mutate(
      AgeGroup = if_else(is.na(AgeGroup) | AgeGroup == "", "unspecified", AgeGroup),
      Species = if_else(is.na(Species) | Species == "", PROJECT_SPECIES, Species),
      Stage = if_else(is.na(Stage) | Stage == "", PROJECT_STAGE, Stage),
      GroupOrder = match(Group, unique(Group)),
      Sample=factor(Sample,levels=Sample),
      Group=factor(Group,levels=unique(Group[order(GroupOrder)]))
    ) %>%
    arrange(GroupOrder, Sample) %>%
    mutate(
      RawFastq = vapply(RawFastqOriginal, resolve_raw_fastq_path, character(1)),
      RawFastqExists = file.exists(RawFastq),
      SampleSportsDir = vapply(as.character(Sample), sports_sample_dir, character(1)),
      SummaryFile=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_summary.txt")),
      LengthFile=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_length_distribution.txt")),
      OutputFile=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_output.txt")),
      PDF_sncRNA=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_sncRNA_distribution.pdf")),
      PDF_rsRNA=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_rsRNA_distribution.pdf")),
      PDF_rRNA=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_rRNA_mapping.pdf")),
      PDF_tRNA=sports_sample_file(as.character(Sample), paste0(as.character(Sample), "_tRNA_mature_mapping.pdf")),
      RRNAMatchFile=file.path(sports_sample_processed_dir(as.character(Sample)), paste0(as.character(Sample), "_output_rRNA_match_genome")),
      RRNAUnmatchFile=file.path(sports_sample_processed_dir(as.character(Sample)), paste0(as.character(Sample), "_output_rRNA_unmatch_genome"))
    )
  missing <- out %>% mutate(Missing=!file.exists(SummaryFile)|!file.exists(LengthFile)|!file.exists(OutputFile)) %>% filter(Missing)
  if (nrow(missing)>0) stop("未检测到完整 SPORTS 结果，请通过所属分支 run.sh --stage upstream 生成输入。缺失样本: ", paste(as.character(missing$Sample),collapse=", "))
  out
}

load_mii_context <- function(report_view = REPORT_VIEW) {
  sample_info <- load_metadata()
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

  # summary_top 仅保留互斥主类别；tRNA 末端类是基础 tRNA 总量的子集，不能重复相加。
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
    stop("Mutually exclusive RNA classes did not recover 100% for: ", paste(bad_samples, collapse = ", "))
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
    transmute(Sample = as.character(Sample), Group = as.character(Group), Replicate = Replicate) %>%
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

  # 长度表同时包含父级总量与 tRNA/rRNA 子集；全局长度构成只保留父级总量。
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
  length_composition_integrity_check <- sample_length_display %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    group_by(Sample) %>%
    summarise(
      ExclusiveLengthReads = sum(Reads, na.rm = TRUE),
      CleanReads = first(CleanReads),
      ReadDifference = ExclusiveLengthReads - CleanReads,
      PercentSum = if_else(CleanReads > 0, 100 * ExclusiveLengthReads / CleanReads, NA_real_),
      Pass = abs(ReadDifference) <= 2,
      .groups = "drop"
    )
  if (any(!length_composition_integrity_check$Pass)) {
    bad_samples <- length_composition_integrity_check %>% filter(!Pass) %>% pull(Sample)
    stop("Mutually exclusive length classes did not recover clean reads for: ", paste(bad_samples, collapse = ", "))
  }
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

  # SPORTS 的 tRNA 末端行是基础总量的子集。基础表用于总量/origin/AA/anticodon；
  # cleavage 表以“基础总量 - 5' - 3' - CCA”重建 Internal/full，保证四类互斥。
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
    stop("SPORTS tRNA end subsets lack matching base-total rows.")
  }

  tsrna_end_totals <- tsrna_end_detail %>%
    group_by(across(all_of(tsrna_reconstruction_keys))) %>%
    summarise(EndReads = sum(Reads, na.rm = TRUE), .groups = "drop")

  tsrna_cleavage_reconstruction_check <- tsrna_base_detail %>%
    transmute(across(all_of(tsrna_reconstruction_keys)), BaseReads = Reads) %>%
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
    bad_samples <- tsrna_cleavage_reconstruction_check %>% filter(!Pass) %>% pull(Sample) %>% unique()
    stop("tRNA end subsets exceed base totals by more than 1 read for: ", paste(bad_samples, collapse = ", "))
  }

  tsrna_internal_detail <- tsrna_cleavage_reconstruction_check %>%
    transmute(
      Sample, Group, Direction,
      Class = BaseClass, SubClass = BaseSubClass, Reads = InternalReads,
      MainCategory = "tsRNA", BaseClass, BaseSubClass, Origin, Maturity,
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

  # output_all 保留逐条序列粒度。PrimaryAnnotation 只用于展示和逐序列分析，
  # tRNA 总量及 isoacceptor 统计使用上面的 SPORTS 分数分配结果，不由首个注释决定。
  output_all <- output_all %>%
    left_join(clean_reads_tbl, by = "Sample") %>%
    mutate(
      Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", Annotation),
      PrimaryAnnotation = str_split_fixed(Annotation, ";", 2)[, 1],
      FeatureCategory = if_else(EndoCandidate, "endo-siRNA", annotation_main_category(PrimaryAnnotation)),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)
    )

  sequence_complexity_tables <- build_sequence_complexity_tables(output_all, sample_info, clean_reads_tbl)
  sequence_complexity_by_sample <- sequence_complexity_tables$sequence_complexity_by_sample
  sequence_complexity_by_group <- sequence_complexity_tables$sequence_complexity_by_group
  sequence_rank_by_sample <- sequence_complexity_tables$sequence_rank_by_sample
  sequence_length_complexity <- sequence_complexity_tables$sequence_length_complexity
  raw_fastq_audit_tables <- load_or_build_raw_fastq_audit_tables(sample_info)
  raw_fastq_audit_by_sample <- raw_fastq_audit_tables$raw_fastq_audit_by_sample
  raw_fastq_audit_by_group <- raw_fastq_audit_tables$raw_fastq_audit_by_group
  raw_fastq_length_distribution <- raw_fastq_audit_tables$raw_fastq_length_distribution
  raw_fastq_end_base <- raw_fastq_audit_tables$raw_fastq_end_base
  raw_fastq_kmer <- raw_fastq_audit_tables$raw_fastq_kmer

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
    "sample_info", "group_levels", "group_palette", "endo_evidence", "endo_transfer_audit", "endo_reference_audit",
    "summary_all_raw", "summary_all", "summary_all_view", "summary_top", "summary_top_view",
    "length_all_raw", "length_all", "length_all_view",
    "clean_reads_tbl", "match_reads_tbl", "group_clean_reads_tbl", "group_match_tbl",
    "group_comp_total", "group_comp_view", "group_comp_direction_view",
    "sample_comp_view", "composition_integrity_check", "sample_qc_table", "group_qc_table",
    "length_detail", "length_detail_view", "length_focus_categories", "sample_length_display", "sample_length_view", "group_length_display", "group_length_view", "length_composition_integrity_check",
    "rsrna_detail", "sample_rsrna_visual_by_view", "sample_rsrna_visual_view", "group_rsrna_visual_view", "group_rsrna_visual_merged", "group_rsrna_view",
    "tsrna_base_detail", "tsrna_cleavage_detail", "tsrna_cleavage_reconstruction_check",
    "group_tsrna_origin", "group_tsrna_origin_view", "group_tsrna_cleavage", "group_tsrna_cleavage_view", "group_tsrna_aa_view", "group_tsrna_anticodon_view",
    "output_all", "group_top_features", "ratio_all_sample_df", "legacy_pdf_index",
    "sequence_complexity_by_sample", "sequence_complexity_by_group", "sequence_rank_by_sample", "sequence_length_complexity",
    "raw_fastq_audit_by_sample", "raw_fastq_audit_by_group", "raw_fastq_length_distribution", "raw_fastq_end_base", "raw_fastq_kmer"
  )
  mget(keep_names, inherits = FALSE)
}

activate_mii_context <- function(env = parent.frame(), report_view = REPORT_VIEW) {
  list2env(load_mii_context(report_view = report_view), envir = env)
  invisible(env)
}

write_mii_setup_expression_matrices <- function(report_view = REPORT_VIEW) {
  ctx <- load_mii_context(report_view = report_view)
  sample_order <- as.character(ctx$sample_info$Sample)
  output_dir <- EXPRESSION_MATRIX_DIR
  seq_dir <- file.path(output_dir, "seq_count")
  dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)
  write_setup_provenance(ctx, output_dir = SETUP_DIR)

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

