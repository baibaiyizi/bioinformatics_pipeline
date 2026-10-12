# SPORTS公共读取、分配、聚合与绘图函数。项目路径和设计在物种adapter中定义。
pkg_has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

hellinger_transform <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  if (any(totals <= 0)) stop("Hellinger transform requires at least one non-zero feature per sample.")
  sqrt(x / totals)
}

read_tsv_flex <- function(file, select = NULL) {
  if (!file.exists(file)) return(tibble())
  if (pkg_has("data.table")) {
    data.table::fread(file, sep="\t", select=select, data.table=FALSE, showProgress=FALSE)
  } else {
    out <- read.delim(file, sep="\t", header=TRUE, check.names=FALSE, quote="", comment.char="", stringsAsFactors=FALSE)
    if (!is.null(select)) out <- out[, select, drop=FALSE]
    out
  } |> as_tibble()
}

read_alignment_tsv_flex <- function(file, col_names) {
  if (!file.exists(file)) return(tibble())
  if (pkg_has("data.table")) {
    out <- data.table::fread(file, sep="\t", header=FALSE, data.table=FALSE, showProgress=FALSE, quote="", fill=TRUE)
  } else {
    out <- read.delim(file, sep="\t", header=FALSE, check.names=FALSE, quote="", comment.char="", fill=TRUE, stringsAsFactors=FALSE)
  }
  if (nrow(out) == 0) return(tibble())
  if (ncol(out) < length(col_names)) for (i in seq.int(ncol(out)+1L, length(col_names))) out[[i]] <- NA_character_
  out <- out[, seq_along(col_names), drop=FALSE]; names(out) <- col_names; as_tibble(out)
}

write_csv_safe <- function(df, file_name, tab_dir = OUTPUT_TAB_DIR) {
  if (is.null(tab_dir) || !nzchar(tab_dir)) stop("tab_dir is not set for write_csv_safe().")
  target_file <- file.path(tab_dir, file_name)
  dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, target_file, row.names = FALSE)
}

write_tsv_safe <- function(df, file_path) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(df, file_path, sep = "\t", quote = FALSE, na = "")
}

normalize_views <- function(view) {
  view <- unique(as.character(view)); invalid <- setdiff(view, VIEW_LEVELS)
  if (length(invalid) > 0) stop("Unsupported view: ", paste(invalid, collapse=", "))
  view
}

safe_dir_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

sncrna_length_axis <- function() {
  ggplot2::scale_x_continuous(
    breaks=function(limits) unique(round(pretty(limits,n=4))),
    limits=function(limits) if(diff(limits)==0) limits+c(-1,1) else limits,
    expand=ggplot2::expansion(mult=0.02))
}

prepare_module_outputs <- function(fig_dir, tab_dir = NULL, sections = character(), stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  if (is.null(fig_dir) || !nzchar(fig_dir)) stop("fig_dir is not set for prepare_module_outputs().")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(tab_dir) && nzchar(tab_dir)) dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(sections) > 0) cleanup_plot_sections(fig_dir, unique(sections), stale_views = unique(c(ALL_VIEW_LEVELS, stale_views)))
}

build_nested_dir <- function(base_dir, section, sample = NULL, view = "merged") {
  parts <- c(base_dir)
  if (length(VIEW_LEVELS) > 1) parts <- c(parts, normalize_views(view)[[1]])
  if (!is.null(section)) parts <- c(parts, section)
  if (!is.null(sample)) parts <- c(parts, safe_dir_name(sample))
  out <- do.call(file.path, as.list(parts)); dir.create(out, recursive=TRUE, showWarnings=FALSE); out
}

save_plot <- function(plot_obj, file_name, width=11, height=7, view="merged") {
  if (length(VIEW_LEVELS) == 1) {
    dir.create(OUTPUT_FIG_DIR, recursive=TRUE, showWarnings=FALSE)
    save_ggplot_file(plot_obj, file.path(OUTPUT_FIG_DIR, file_name), width=width, height=height)
  } else {
    for (v in normalize_views(view)) {
      d <- file.path(OUTPUT_FIG_DIR, v); dir.create(d, recursive=TRUE, showWarnings=FALSE)
      save_ggplot_file(plot_obj, file.path(d, file_name), width=width, height=height)
    }
  }
}

save_plot_nested <- function(plot_obj, section, file_name, sample=NULL, width=11, height=7, view="merged") {
  for (v in normalize_views(view)) {
    d <- build_nested_dir(OUTPUT_FIG_DIR, section, sample=sample, view=v)
    save_ggplot_file(plot_obj, file.path(d, file_name), width=width, height=height)
  }
}

save_pdf <- function(plot_obj, filepath, width=11, height=7) {
  save_ggplot_file(plot_obj, filepath, width=width, height=height)
}

filter_direction_view <- function(df, view, direction_col="Direction") {
  if (view == "merged" || !direction_col %in% names(df)) return(df)
  df %>% filter(as.character(.data[[direction_col]]) == view)
}

view_title_suffix <- function(view) {
  switch(as.character(view), merged="", sense="", antisense=" (antisense only)", "")
}

plot_note <- function(title, note) {
  ggplot() + annotate("text", x=0, y=0, label=paste(title, note, sep="\n"), size=5) +
    xlim(-1,1) + ylim(-1,1) + theme_void()
}

clean_numeric <- function(x) suppressWarnings(as.numeric(x))

group_percent <- function(x) { t <- sum(x, na.rm=TRUE); if(t>0) 100*x/t else rep(0, length(x)) }

calculate_gini <- function(x) {
  x <- sort(as.numeric(x[is.finite(x) & x >= 0]))
  n <- length(x)
  total <- sum(x, na.rm = TRUE)
  if (n == 0 || total <= 0) return(NA_real_)
  (2 * sum(seq_len(n) * x, na.rm = TRUE) / (n * total)) - ((n + 1) / n)
}

sequence_entropy <- function(x) {
  x <- as.numeric(x)
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  p <- x[x > 0 & is.finite(x)] / total
  -sum(p * log(p), na.rm = TRUE)
}

top_read_sum <- function(x, n = 1L) {
  x <- sort(as.numeric(x), decreasing = TRUE)
  sum(head(x[is.finite(x) & x > 0], n), na.rm = TRUE)
}

get_col <- function(df, col_name) {
  if (col_name %in% colnames(df)) df[[col_name]] else rep(0, nrow(df))
}

summarise_group_mean_rpm <- function(sample_level_df, key_cols, sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$RPM <- numeric()
    return(out)
  }
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!"Reads" %in% names(sample_level_df)) sample_level_df$Reads <- 0
  if (!"RPM" %in% names(sample_level_df)) stop("sample_level_df must contain an RPM column.")

  sample_summary <- sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      RPM = sum(RPM, na.rm = TRUE),
      .groups = "drop"
    )
  key_tbl <- sample_level_df %>%
    select(all_of(key_cols)) %>%
    distinct()
  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")
  group_summary <- sample_summary %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      RPM = sum(RPM, na.rm = TRUE),
      .groups = "drop"
    )

  tidyr::expand_grid(group_grid, key_tbl) %>%
    left_join(group_summary, by = c("Group", key_cols)) %>%
    mutate(
      Reads = replace_na(Reads, 0),
      RPM = if_else(sample_n > 0, replace_na(RPM, 0) / sample_n, 0)
    ) %>%
    select(-sample_n)
}

summarise_sparse_group_mean_rpm <- function(sample_level_df, key_cols, sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$RPM <- numeric()
    return(out)
  }
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!"Reads" %in% names(sample_level_df)) sample_level_df$Reads <- 0
  if (!"RPM" %in% names(sample_level_df)) stop("sample_level_df must contain an RPM column.")

  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")

  sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    left_join(group_grid, by = "Group") %>%
    mutate(RPM = if_else(sample_n > 0, RPM / sample_n, 0)) %>%
    select(-sample_n)
}

summarise_group_mean_percent <- function(sample_level_df, key_cols, percent_col = "Percent", reads_col = "Reads", sample_info_df = NULL) {
  if (is.null(sample_info_df)) sample_info_df <- get("sample_info", envir = parent.frame())
  if (nrow(sample_level_df) == 0) {
    out <- tibble(Group = character())
    for (key_col in key_cols) out[[key_col]] <- (sample_level_df[[key_col]] %||% character())[0]
    out$Reads <- numeric()
    out$Percent <- numeric()
    return(out)
  }
  if (!percent_col %in% names(sample_level_df)) stop("sample_level_df must contain a percent column.")
  sample_level_df <- sample_level_df %>%
    ungroup() %>%
    mutate(Sample = as.character(Sample), Group = as.character(Group))
  if (!reads_col %in% names(sample_level_df)) sample_level_df[[reads_col]] <- 0

  sample_summary <- sample_level_df %>%
    group_by(Sample, Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(.data[[reads_col]], na.rm = TRUE),
      Percent = sum(.data[[percent_col]], na.rm = TRUE),
      .groups = "drop"
    )
  key_tbl <- sample_level_df %>%
    select(all_of(key_cols)) %>%
    distinct()
  group_grid <- sample_info_df %>%
    transmute(Group = as.character(Group)) %>%
    count(Group, name = "sample_n")
  group_summary <- sample_summary %>%
    group_by(Group, across(all_of(key_cols))) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      Percent = sum(Percent, na.rm = TRUE),
      .groups = "drop"
    )

  tidyr::expand_grid(group_grid, key_tbl) %>%
    left_join(group_summary, by = c("Group", key_cols)) %>%
    mutate(
      Reads = replace_na(Reads, 0),
      Percent = if_else(sample_n > 0, replace_na(Percent, 0) / sample_n, 0)
    ) %>%
    select(-sample_n)
}

resolve_raw_fastq_path <- function(path_value) {
  if (is.na(path_value) || !nzchar(trimws(path_value))) return(NA_character_)
  if (grepl("^/", path_value)) path_value else file.path(ROOT_DIR, path_value)
}

capture_command_version <- function(command, args = c("--version")) {
  command_path <- Sys.which(command)
  if (!nzchar(command_path)) {
    return(tibble(tool = command, available = FALSE, path = "", version = "not available"))
  }

  version_text <- tryCatch(
    {
      out <- suppressWarnings(system2(command_path, args = args, stdout = TRUE, stderr = TRUE))
      paste(head(out, 3), collapse = " | ")
    },
    error = function(e) paste("version check failed:", conditionMessage(e))
  )

  tibble(tool = command, available = TRUE, path = unname(command_path), version = version_text)
}

build_runtime_config_table <- function() {
  tibble(
    key = c(
      "ROOT_DIR", "RAW_DIR", "SPORTS_RESULT_DIR", "SPORTS_METADATA_DIR", "SPORTS_QC_DIR",
      "SPORTS_SAMPLE_ROOT", "SPORTS_ENDOSIRNA_DIR", "SPORTS_LOG_DIR", "RESULT_DIR",
      "PROJECT_SPECIES", "PROJECT_STAGE", "PROJECT_CONTRAST_BASELINE", "PROJECT_CONTRAST_CASE",
      "RUN_FASTQC", "RUN_EXOGENOUS", "EXOGENOUS_DB_ROOT",
      "MIRNA_TARGET_SOURCE", "MIRNA_SPECIES", "MIRNA_PREFIX", "MULTIMIR_ORG",
      "MIRNA_TARGET_ORGDB_PACKAGE", "KEGG_ORGANISM", "DE_SIGNIFICANCE",
      "RAW_AUDIT_MAX_READS", "RAW_AUDIT_KMER", "RAW_AUDIT_REFRESH",
      "REPORT_VIEW"
    ),
    value = c(
      ROOT_DIR, RAW_DIR, SPORTS_RESULT_DIR, SPORTS_METADATA_DIR, SPORTS_QC_DIR,
      SPORTS_SAMPLE_ROOT, SPORTS_ENDOSIRNA_DIR, SPORTS_LOG_DIR, RESULT_DIR,
      PROJECT_SPECIES, PROJECT_STAGE, PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE,
      RUN_FASTQC, RUN_EXOGENOUS, EXOGENOUS_DB_ROOT,
      MIRNA_TARGET_SOURCE, MIRNA_SPECIES, MIRNA_PREFIX, MULTIMIR_ORG,
      MIRNA_TARGET_ORGDB_PACKAGE, KEGG_ORGANISM, DE_SIGNIFICANCE,
      RAW_AUDIT_MAX_READS, RAW_AUDIT_KMER, RAW_AUDIT_REFRESH,
      REPORT_VIEW
    )
  )
}

is_mito_label <- function(x) {
  str_detect(replace_na(as.character(x),""), regex("(^|[-_])(mt|mito)|mitochond|chrM|12S-rRNA|16S-rRNA", ignore_case=TRUE))
}

remove_path_if_exists <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  } else if (file.exists(path)) {
    unlink(path, force = TRUE)
  }
}

cleanup_stale_view_dirs <- function(base_dir, stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  walk(stale_views, ~remove_path_if_exists(file.path(base_dir, .x)))
}

cleanup_plot_sections <- function(base_dir, sections, views = VIEW_LEVELS, stale_views = setdiff(ALL_VIEW_LEVELS, VIEW_LEVELS)) {
  cleanup_stale_view_dirs(base_dir, stale_views = stale_views)
  for (section_name in sections) {
    remove_path_if_exists(file.path(base_dir, section_name))
  }
  for (view_name in unique(c(views, ALL_VIEW_LEVELS))) {
    for (section_name in sections) {
      remove_path_if_exists(file.path(base_dir, view_name, section_name))
    }
  }
}

build_plot_collage <- function(plot_list, title_text = NULL, cols = 4) {
  valid_plots <- compact(plot_list)
  if (length(valid_plots) == 0) {
    return(plot_note(title_text %||% "Plot collage", "No plots available."))
  }
  if (!pkg_has("patchwork")) {
    return(plot_note(title_text %||% "Plot collage", "patchwork package is not available."))
  }

  patchwork::wrap_plots(valid_plots, ncol = cols) +
    patchwork::plot_annotation(
      title = title_text,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
      )
    )
}

collage_dimensions <- function(n_plots, cols = 4, base_width = 4, base_height = 3.1) {
  rows <- max(1, ceiling(n_plots / cols))
  list(
    width = max(10, min(18, cols * base_width)),
    height = max(6, rows * base_height + 0.9)
  )
}

class_direction <- function(class_name, sub_class="") {
  if_else(str_detect(paste(class_name, sub_class), "antisense"), "antisense", "sense")
}

class_genome_status <- function(class_name) {
  case_when(str_detect(class_name,"_Match_Genome$")~"Match_Genome", str_detect(class_name,"_Unmatch_Genome$")~"Unmatch_Genome", TRUE~"Unknown")
}

class_source_db <- function(class_name) {
  case_when(
    str_starts(class_name,"miRBase-")~"miRBase", str_starts(class_name,"rRNAdb-")~"rRNAdb",
    str_starts(class_name,"YRNAdb-")~"YRNAdb", str_starts(class_name,"GtRNAdb-")~"GtRNAdb",
    str_starts(class_name,"mitotRNAdb-")~"mitotRNAdb", str_starts(class_name,"piRNAdb-")~"piRNA",
    str_starts(class_name,"ensembl-")~"Ensembl", str_starts(class_name,"Rfam-")~"Rfam",
    str_starts(class_name,"Unannotated_")~"Unannotated", TRUE~"Other")
}

class_main_category <- function(class_name, sub_class="") {
  token <- paste(class_name, sub_class)
  case_when(
    str_detect(token,"^endo-siRNA")~"endo-siRNA",
    str_detect(token,"^Unannotated_Match_Genome")~"unanno_MG", str_detect(token,"^Unannotated_Unmatch_Genome")~"unanno_UMG",
    str_detect(token,"YRNA")~"ysRNA", str_detect(token,"miRNA")~"miRNA", str_detect(token,"piRNA")~"piRNA",
    str_detect(token,"scaRNA")~"scaRNA", str_detect(token,"snoRNA")~"snoRNA", str_detect(token,"snRNA")~"snRNA",
    str_detect(token,"lncRNA|lincRNA")~"lncRNA", str_detect(token,"vault_RNA")~"vault_RNA",
    str_detect(token,"misc_RNA")~"misc_RNA", str_detect(token,"ribozyme")~"ribozyme",
    str_detect(token,"tRNA")~"tsRNA", str_detect(token,"rRNA")~"rsRNA", TRUE~"other")
}

annotation_main_category <- function(annotation) {
  token <- if_else(is.na(annotation)|annotation=="", "NO_Annotation", annotation)
  case_when(
    str_detect(token,"^endo-siRNA")~"endo-siRNA",
    str_detect(token,"YRNA")~"ysRNA", str_detect(token,regex("mir|let-",ignore_case=TRUE))~"miRNA",
    str_detect(token,regex("piRNA",ignore_case=TRUE))~"piRNA", str_detect(token,"scaRNA")~"scaRNA",
    str_detect(token,"snoRNA")~"snoRNA", str_detect(token,"snRNA")~"snRNA",
    str_detect(token,"lncRNA|lincRNA")~"lncRNA", str_detect(token,"vault_RNA")~"vault_RNA",
    str_detect(token,"misc_RNA")~"misc_RNA", str_detect(token,"ribozyme")~"ribozyme",
    str_detect(token,"tRNA")~"tsRNA", str_detect(token,"rRNA")~"rsRNA",
    str_detect(token,"NO_Annotation")~"other", TRUE~"other")
}

length_category_from_class <- function(class_name) {
  case_when(
    str_detect(class_name,"^endo-siRNA")~"endo-siRNA",
    str_detect(class_name,"miRNA")~"miRNA", str_detect(class_name,"piRNA")~"piRNA",
    str_detect(class_name,"YRNA")~"ysRNA", str_detect(class_name,"^mitotRNAdb")~"mt_tsRNA",
    str_detect(class_name,"Mt_rRNA")~"mt_rsRNA", str_detect(class_name,"tRNA")~"tsRNA",
    str_detect(class_name,"rRNA")~"rsRNA", str_detect(class_name,"lncRNA")~"lncRNA",
    str_detect(class_name,"snRNA")~"snRNA", str_detect(class_name,"snoRNA")~"snoRNA",
    str_detect(class_name,"scaRNA")~"scaRNA", str_detect(class_name,"misc_RNA")~"misc_RNA",
    str_detect(class_name,"vault_RNA")~"vault_RNA", str_detect(class_name,"ribozyme")~"ribozyme",
    str_detect(class_name,"scRNA")~"scRNA", str_detect(class_name,"(^|-)sRNA_")~"sRNA",
    str_detect(class_name,"mt_DNA")~"mt_DNA", str_detect(class_name,"DNA")~"DNA", TRUE~NA_character_)
}

parse_rsrna_class <- function(class_name, sub_class) {
  token <- if_else(sub_class=="-"|is.na(sub_class), class_name, sub_class) |> str_remove("-antisense$")
  case_when(
    str_detect(token,"45S-rRNA")~"45S-rRNA", str_detect(token,"28S-rRNA")~"28S-rRNA",
    str_detect(token,"26S-rRNA")~"26S-rRNA", str_detect(token,"25S-rRNA")~"25S-rRNA",
    str_detect(token,"23S-rRNA")~"23S-rRNA", str_detect(token,"18S-rRNA")~"18S-rRNA",
    str_detect(token,"17S-rRNA")~"17S-rRNA", str_detect(token,"16S-rRNA")~"16S-rRNA",
    str_detect(token,"12S-rRNA")~"12S-rRNA", str_detect(token,"5\\.8S-rRNA")~"5.8S-rRNA",
    str_detect(token,"5\\.3S-rRNA")~"5.3S-rRNA", str_detect(token,"5S-rRNA")~"5S-rRNA",
    str_detect(token,"4\\.5S-rRNA")~"4.5S-rRNA", str_detect(token,"2S-rRNA")~"2S-rRNA",
    str_detect(token,"mt_rRNA|Mt_rRNA")~"mt-rRNA", str_detect(token,"other-rRNA")~"other-rRNA",
    str_starts(class_name,"ensembl-rRNA")~"Ensembl-rRNA", str_starts(class_name,"Rfam-")~"Rfam-rRNA-like",
    TRUE~"other-rRNA")
}

parse_tsrna_origin <- function(class_name, sub_class) if_else(str_detect(paste(class_name,sub_class),"mitotRNAdb|mt_tRNA"), "Mitochondrial", "Genomic")

parse_tsrna_maturity <- function(class_name, sub_class) if_else(str_detect(paste(class_name,sub_class),"pre"), "Pre-tRNA", "Mature-tRNA")

parse_tsrna_cleavage <- function(class_name, sub_class) {
  token <- paste(class_name, sub_class)
  case_when(str_detect(token,"5_end")~"5'end", str_detect(token,"3_end")~"3'end", str_detect(token,"CCA_end")~"CCA-end", TRUE~"Internal/full")
}

parse_tsrna_amino_acid <- function(sub_class) {
  out <- str_match(sub_class, "(?:mt_)?tRNA-([A-Za-z]+[0-9]*)-[A-Z]{3}")[, 2]
  if_else(is.na(out) | out == "", "Unknown", out)
}

parse_tsrna_anticodon <- function(sub_class) {
  out <- str_match(sub_class, "(?:mt_)?tRNA-[A-Za-z]+[0-9]*-([A-Z]{3})")[, 2]
  if_else(is.na(out) | out == "", "Unknown", out)
}

is_tsrna_nested_class <- function(class_name) {
  str_detect(as.character(class_name), "(?:tRNA|mt_tRNA)_(?:5_end|3_end|CCA_end)_")
}

is_nested_length_class <- function(class_name) {
  token <- as.character(class_name)
  (str_starts(token, "rRNAdb-") & !str_starts(token, "rRNAdb-rRNA_")) |
    (str_starts(token, "YRNAdb-") & !str_starts(token, "YRNAdb-YRNA_"))
}

strip_tsrna_end_from_class <- function(class_name) {
  str_replace(as.character(class_name), "_(?:5_end|3_end|CCA_end)(?=_)", "")
}

strip_tsrna_end_from_subclass <- function(sub_class) {
  str_replace(as.character(sub_class), "_(?:5_end|3_end|CCA_end)$", "")
}

collapse_rsrna_visual_class <- function(x) case_when(x %in% c("17S-rRNA","23S-rRNA","25S-rRNA","26S-rRNA","5.3S-rRNA")~"other-rRNA", TRUE~x)

darken_color <- function(cols, amount = 0.22) {
  if (length(cols) == 0) return(character())
  rgb_mat <- grDevices::col2rgb(cols) / 255
  darker <- rgb_mat * (1 - amount)
  darker[darker < 0] <- 0
  grDevices::rgb(darker[1, ], darker[2, ], darker[3, ])
}

order_by_reference <- function(levels, reference_levels) {
  levels <- unique(as.character(levels))
  c(reference_levels[reference_levels %in% levels], setdiff(levels, reference_levels))
}

make_sector_label_palette <- function(sector_palette, amount = 0.28) {
  sector_palette <- sector_palette[!is.na(sector_palette)]
  setNames(darken_color(unname(sector_palette), amount = amount), names(sector_palette))
}

make_group_palette <- function(groups) {
  base <- c("#8ECFC9", "#FFBE7A", "#FF9F9B", "#A1C9F4", "#BEB8DC", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4")
  groups <- as.character(groups)
  if (length(groups) > length(base)) {
    extra <- grDevices::hcl.colors(length(groups) - length(base), palette = "Pastel 1")
    base <- c(base, extra)
  }
  setNames(base[seq_along(groups)], groups)
}

collapse_rsrna_radar_class <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% setdiff(rsrna_radar_levels, "other-rRNA") ~ x,
    TRUE ~ "other-rRNA"
  )
}

mii_theme_bw <- function(...) {
  theme_bw() + theme(panel.grid.minor = element_blank(), ...)
}

mii_theme_minimal <- function(...) {
  theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    ...
  )
}

sports_sample_dir <- function(sample_id) {
  file.path(SPORTS_SAMPLE_ROOT, as.character(sample_id))
}

sports_sample_file <- function(sample_id, ...) {
  file.path(sports_sample_dir(sample_id), ...)
}

sports_sample_processed_dir <- function(sample_id) {
  sports_sample_file(sample_id, paste0(as.character(sample_id), "_processed"))
}

sports_endosirna_summary_file <- function() {
  file.path(SPORTS_ENDOSIRNA_SUMMARY_DIR, "endo_siRNA_summary.tsv")
}

load_summary_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$SummaryFile[i]); if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]), Replicate=as.character(si$Replicate[i]),
           Class=as.character(df$Class), SubClass=as.character(df$Sub_Class), Reads=clean_numeric(df$Reads))
  })
}

load_length_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$LengthFile[i]); if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]),
           Class=as.character(df$Class), Length=clean_numeric(df$Length), Reads=clean_numeric(df$Reads))
  })
}

load_output_tables <- function(si) {
  map_dfr(seq_len(nrow(si)), function(i) {
    df <- read_tsv_flex(si$OutputFile[i], select=c("ID","Sequence","Length","Reads","Match_Genome","Annotation"))
    if(nrow(df)==0) return(tibble())
    tibble(Sample=as.character(si$Sample[i]), Group=as.character(si$Group[i]), ID=as.character(df$ID),
           Sequence=as.character(df$Sequence), Length=clean_numeric(df$Length), Reads=clean_numeric(df$Reads),
           MatchGenome=as.character(df$Match_Genome), Annotation=as.character(df$Annotation))
  })
}

normalize_filter_sequence <- function(x) {
  toupper(gsub("U", "T", trimws(as.character(x))))
}

build_sequence_complexity_tables <- function(output_all_tbl, sample_info_tbl, clean_reads_tbl, top_rank_limit = 50L) {
  sample_base <- sample_info_tbl %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group))

  empty_result <- list(
    sequence_complexity_by_sample = sample_base %>%
      left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
      mutate(
        CleanReads = replace_na(CleanReads, 0),
        SequenceReads = 0,
        UniqueSequences = 0,
        Top1Reads = 0,
        Top5Reads = 0,
        Top10Reads = 0,
        Top1CleanFraction = 0,
        Top5CleanFraction = 0,
        Top10CleanFraction = 0,
        Top1OutputFraction = 0,
        Top5OutputFraction = 0,
        Top10OutputFraction = 0,
        ShannonEntropy = NA_real_,
        EffectiveSequences = NA_real_,
        GiniIndex = NA_real_,
        DominantSequence = NA_character_,
        DominantAnnotation = NA_character_,
        DominantCategory = NA_character_,
        DominantLength = NA_real_,
        DominantReads = 0,
        DominantRPM = 0
      ),
    sequence_complexity_by_group = tibble(),
    sequence_rank_by_sample = tibble(),
    sequence_length_complexity = tibble()
  )

  if (is.null(output_all_tbl) || nrow(output_all_tbl) == 0) return(empty_result)

  if (pkg_has("data.table")) {
    dt <- data.table::as.data.table(output_all_tbl)
    dt <- data.table::copy(dt)
    if (!"PrimaryAnnotation" %in% names(dt)) {
      dt[, Annotation := ifelse(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation))]
      dt[, PrimaryAnnotation := sub(";.*$", "", Annotation)]
    }
    if (!"FeatureCategory" %in% names(dt)) {
      dt[, FeatureCategory := annotation_main_category(PrimaryAnnotation)]
    }

    dt[, Sample := as.character(Sample)]
    dt[, Group := as.character(Group)]
    dt[, Sequence := normalize_filter_sequence(Sequence)]
    dt[, Length := clean_numeric(Length)]
    dt[, Reads := clean_numeric(Reads)]
    dt[, Annotation := ifelse(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation))]
    dt[, PrimaryAnnotation := as.character(PrimaryAnnotation)]
    dt[, FeatureCategory := as.character(FeatureCategory)]
    dt <- dt[!is.na(Sequence) & Sequence != "" & is.finite(Reads) & Reads > 0]
    if (nrow(dt) == 0) return(empty_result)
    dt[is.na(Length) | Length <= 0, Length := nchar(Sequence)]

    sequence_counts_dt <- dt[, .(Reads = sum(Reads, na.rm = TRUE)), by = .(Sample, Group, Sequence, Length, PrimaryAnnotation, FeatureCategory)]
    data.table::setorder(sequence_counts_dt, Sample, Group, Sequence, Length, -Reads, PrimaryAnnotation)
    sequence_counts_dt <- sequence_counts_dt[, .(
      Reads = sum(Reads, na.rm = TRUE),
      PrimaryAnnotation = PrimaryAnnotation[1L],
      FeatureCategory = FeatureCategory[1L]
    ), by = .(Sample, Group, Sequence, Length)]

    data.table::setorder(sequence_counts_dt, Sample, -Reads, Sequence)
    sequence_counts_dt[, Rank := seq_len(.N), by = Sample]
    sequence_counts_dt[, SequenceReads := sum(Reads, na.rm = TRUE), by = Sample]
    sequence_counts_dt[, CumulativeReads := cumsum(Reads), by = Sample]
    sequence_counts_dt[, FractionOfOutput := ifelse(SequenceReads > 0, Reads / SequenceReads, 0)]
    sequence_counts_dt[, CumulativeFractionOfOutput := ifelse(SequenceReads > 0, CumulativeReads / SequenceReads, 0)]

    clean_dt <- data.table::as.data.table(clean_reads_tbl %>% mutate(Sample = as.character(Sample)))
    sequence_rank_dt <- merge(
      sequence_counts_dt[Rank <= top_rank_limit],
      clean_dt,
      by = "Sample",
      all.x = TRUE,
      sort = FALSE
    )
    sequence_rank_dt[is.na(CleanReads), CleanReads := 0]
    sequence_rank_dt[, FractionOfClean := ifelse(CleanReads > 0, Reads / CleanReads, 0)]
    sequence_rank_dt[, CumulativeFractionOfClean := ifelse(CleanReads > 0, CumulativeReads / CleanReads, 0)]
    sequence_rank_dt[, RPM := ifelse(CleanReads > 0, 1e6 * Reads / CleanReads, 0)]
    data.table::setorder(sequence_rank_dt, Sample, Rank)

    dominant_dt <- sequence_rank_dt[Rank == 1, .(
      Sample,
      DominantSequence = Sequence,
      DominantAnnotation = PrimaryAnnotation,
      DominantCategory = FeatureCategory,
      DominantLength = Length,
      DominantReads = Reads,
      DominantRPM = RPM
    )]

    sequence_complexity_dt <- sequence_counts_dt[, .(
      SequenceReads = sum(Reads, na.rm = TRUE),
      UniqueSequences = .N,
      Top1Reads = top_read_sum(Reads, 1),
      Top5Reads = top_read_sum(Reads, 5),
      Top10Reads = top_read_sum(Reads, 10),
      ShannonEntropy = sequence_entropy(Reads),
      EffectiveSequences = exp(sequence_entropy(Reads)),
      GiniIndex = calculate_gini(Reads)
    ), by = .(Sample, Group)]

    sample_base_dt <- data.table::as.data.table(sample_base)
    sequence_complexity_dt <- merge(sample_base_dt, sequence_complexity_dt, by = c("Sample", "Group"), all.x = TRUE, sort = FALSE)
    sequence_complexity_dt <- merge(sequence_complexity_dt, clean_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    sequence_complexity_dt <- merge(sequence_complexity_dt, dominant_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    numeric_zero_cols <- c("CleanReads", "SequenceReads", "UniqueSequences", "Top1Reads", "Top5Reads", "Top10Reads", "DominantReads", "DominantRPM")
    for (col_name in numeric_zero_cols) {
      if (col_name %in% names(sequence_complexity_dt)) {
        data.table::set(sequence_complexity_dt, which(is.na(sequence_complexity_dt[[col_name]])), col_name, 0)
      }
    }
    sequence_complexity_dt[, Top1CleanFraction := ifelse(CleanReads > 0, Top1Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top5CleanFraction := ifelse(CleanReads > 0, Top5Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top10CleanFraction := ifelse(CleanReads > 0, Top10Reads / CleanReads, 0)]
    sequence_complexity_dt[, Top1OutputFraction := ifelse(SequenceReads > 0, Top1Reads / SequenceReads, 0)]
    sequence_complexity_dt[, Top5OutputFraction := ifelse(SequenceReads > 0, Top5Reads / SequenceReads, 0)]
    sequence_complexity_dt[, Top10OutputFraction := ifelse(SequenceReads > 0, Top10Reads / SequenceReads, 0)]
    data.table::setorder(sequence_complexity_dt, Group, Sample)

    complexity_metric_cols <- c(
      "CleanReads", "SequenceReads", "UniqueSequences",
      "Top1CleanFraction", "Top5CleanFraction", "Top10CleanFraction",
      "Top1OutputFraction", "Top5OutputFraction", "Top10OutputFraction",
      "ShannonEntropy", "EffectiveSequences", "GiniIndex", "DominantRPM"
    )
    for (col_name in complexity_metric_cols) {
      if (col_name %in% names(sequence_complexity_dt)) {
        sequence_complexity_dt[, (col_name) := as.numeric(.SD[[1L]]), .SDcols = col_name]
      }
    }
    sequence_group_dt <- sequence_complexity_dt[, c(
      list(SampleN = .N),
      unlist(lapply(.SD, function(x) {
        list(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
      }), recursive = FALSE)
    ), by = Group, .SDcols = complexity_metric_cols]
    expected_group_names <- c(
      "Group", "SampleN",
      unlist(lapply(complexity_metric_cols, function(x) paste0(x, c("_mean", "_sd", "_median"))))
    )
    if (length(expected_group_names) == ncol(sequence_group_dt)) names(sequence_group_dt) <- expected_group_names

    length_complexity_dt <- sequence_counts_dt[, .(Reads = sum(Reads, na.rm = TRUE), UniqueSequences = .N), by = .(Sample, Group, Length)]
    length_complexity_dt <- merge(length_complexity_dt, clean_dt, by = "Sample", all.x = TRUE, sort = FALSE)
    length_complexity_dt[is.na(CleanReads), CleanReads := 0]
    length_complexity_dt[, RPM := ifelse(CleanReads > 0, 1e6 * Reads / CleanReads, 0)]
    length_complexity_dt[, FractionOfClean := ifelse(CleanReads > 0, Reads / CleanReads, 0)]
    data.table::setorder(length_complexity_dt, Sample, Length)

    return(list(
      sequence_complexity_by_sample = as_tibble(sequence_complexity_dt),
      sequence_complexity_by_group = as_tibble(sequence_group_dt),
      sequence_rank_by_sample = as_tibble(sequence_rank_dt),
      sequence_length_complexity = as_tibble(length_complexity_dt)
    ))
  }

  input_tbl <- output_all_tbl
  if (!"PrimaryAnnotation" %in% names(input_tbl)) {
    input_tbl <- input_tbl %>%
      mutate(
        Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation)),
        PrimaryAnnotation = str_split_fixed(Annotation, ";", 2)[, 1]
      )
  }
  if (!"FeatureCategory" %in% names(input_tbl)) {
    input_tbl <- input_tbl %>%
      mutate(FeatureCategory = annotation_main_category(PrimaryAnnotation))
  }

  output_clean <- input_tbl %>%
    mutate(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Sequence = normalize_filter_sequence(Sequence),
      Length = clean_numeric(Length),
      Reads = clean_numeric(Reads),
      Annotation = if_else(is.na(Annotation) | Annotation == "", "NO_Annotation", as.character(Annotation)),
      PrimaryAnnotation = as.character(PrimaryAnnotation),
      FeatureCategory = as.character(FeatureCategory)
    ) %>%
    filter(!is.na(Sequence), Sequence != "", Reads > 0) %>%
    mutate(Length = if_else(is.na(Length) | Length <= 0, as.numeric(nchar(Sequence)), Length))

  if (nrow(output_clean) == 0) return(empty_result)

  sequence_counts <- output_clean %>%
    group_by(Sample, Group, Sequence, Length, PrimaryAnnotation, FeatureCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    arrange(Sample, Group, Sequence, Length, desc(Reads), PrimaryAnnotation) %>%
    group_by(Sample, Group, Sequence, Length) %>%
    summarise(
      PrimaryAnnotation = first(PrimaryAnnotation),
      FeatureCategory = first(FeatureCategory),
      Reads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      PrimaryAnnotation = replace_na(PrimaryAnnotation, "NO_Annotation"),
      FeatureCategory = replace_na(FeatureCategory, "unanno_UMG")
    )

  sequence_rank_by_sample <- sequence_counts %>%
    group_by(Sample) %>%
    arrange(desc(Reads), Sequence, .by_group = TRUE) %>%
    mutate(
      Rank = row_number(),
      SequenceReads = sum(Reads, na.rm = TRUE),
      CumulativeReads = cumsum(Reads),
      FractionOfOutput = if_else(SequenceReads > 0, Reads / SequenceReads, 0),
      CumulativeFractionOfOutput = if_else(SequenceReads > 0, CumulativeReads / SequenceReads, 0)
    ) %>%
    ungroup() %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      FractionOfClean = if_else(CleanReads > 0, Reads / CleanReads, 0),
      CumulativeFractionOfClean = if_else(CleanReads > 0, CumulativeReads / CleanReads, 0),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)
    ) %>%
    filter(Rank <= top_rank_limit) %>%
    arrange(Sample, Rank)

  dominant_by_sample <- sequence_rank_by_sample %>%
    filter(Rank == 1) %>%
    transmute(
      Sample,
      DominantSequence = Sequence,
      DominantAnnotation = PrimaryAnnotation,
      DominantCategory = FeatureCategory,
      DominantLength = Length,
      DominantReads = Reads,
      DominantRPM = RPM
    )

  sequence_complexity_by_sample <- sequence_counts %>%
    group_by(Sample, Group) %>%
    summarise(
      SequenceReads = sum(Reads, na.rm = TRUE),
      UniqueSequences = n(),
      Top1Reads = top_read_sum(Reads, 1),
      Top5Reads = top_read_sum(Reads, 5),
      Top10Reads = top_read_sum(Reads, 10),
      ShannonEntropy = sequence_entropy(Reads),
      EffectiveSequences = exp(ShannonEntropy),
      GiniIndex = calculate_gini(Reads),
      .groups = "drop"
    ) %>%
    right_join(sample_base, by = c("Sample", "Group")) %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    left_join(dominant_by_sample, by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      SequenceReads = replace_na(SequenceReads, 0),
      UniqueSequences = replace_na(UniqueSequences, 0),
      Top1Reads = replace_na(Top1Reads, 0),
      Top5Reads = replace_na(Top5Reads, 0),
      Top10Reads = replace_na(Top10Reads, 0),
      Top1CleanFraction = if_else(CleanReads > 0, Top1Reads / CleanReads, 0),
      Top5CleanFraction = if_else(CleanReads > 0, Top5Reads / CleanReads, 0),
      Top10CleanFraction = if_else(CleanReads > 0, Top10Reads / CleanReads, 0),
      Top1OutputFraction = if_else(SequenceReads > 0, Top1Reads / SequenceReads, 0),
      Top5OutputFraction = if_else(SequenceReads > 0, Top5Reads / SequenceReads, 0),
      Top10OutputFraction = if_else(SequenceReads > 0, Top10Reads / SequenceReads, 0),
      DominantReads = replace_na(DominantReads, 0),
      DominantRPM = replace_na(DominantRPM, 0)
    ) %>%
    arrange(Group, Sample)

  complexity_metric_cols <- c(
    "CleanReads", "SequenceReads", "UniqueSequences",
    "Top1CleanFraction", "Top5CleanFraction", "Top10CleanFraction",
    "Top1OutputFraction", "Top5OutputFraction", "Top10OutputFraction",
    "ShannonEntropy", "EffectiveSequences", "GiniIndex", "DominantRPM"
  )

  sequence_complexity_by_group <- sequence_complexity_by_sample %>%
    group_by(Group) %>%
    summarise(
      SampleN = n(),
      across(
        all_of(complexity_metric_cols),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~stats::sd(.x, na.rm = TRUE),
          median = ~stats::median(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  sequence_length_complexity <- sequence_counts %>%
    group_by(Sample, Group, Length) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), UniqueSequences = n(), .groups = "drop") %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(
      CleanReads = replace_na(CleanReads, 0),
      RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0),
      FractionOfClean = if_else(CleanReads > 0, Reads / CleanReads, 0)
    ) %>%
    arrange(Sample, Length)

  list(
    sequence_complexity_by_sample = sequence_complexity_by_sample,
    sequence_complexity_by_group = sequence_complexity_by_group,
    sequence_rank_by_sample = sequence_rank_by_sample,
    sequence_length_complexity = sequence_length_complexity
  )
}

empty_raw_fastq_audit_tables <- function(sample_info_tbl = tibble()) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) {
    sample_base <- tibble(Sample = character(), Group = character(), RawFastq = character(), RawFastqExists = logical())
  } else {
    sample_info_tbl <- as_tibble(sample_info_tbl)
    if (!"RawFastq" %in% names(sample_info_tbl)) sample_info_tbl$RawFastq <- ""
    if (!"Group" %in% names(sample_info_tbl)) sample_info_tbl$Group <- ""
    if (!"Sample" %in% names(sample_info_tbl)) sample_info_tbl$Sample <- ""
  sample_base <- sample_info_tbl %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
        RawFastq = as.character(RawFastq),
      RawFastqExists = file.exists(RawFastq)
    )
  }

  list(
    raw_fastq_audit_by_sample = sample_base %>%
      mutate(
        Status = if_else(RawFastqExists, "not_run", "missing_fastq"),
        MaxReadsRequested = RAW_AUDIT_MAX_READS,
        ReadsExamined = 0L,
        KmerK = RAW_AUDIT_KMER,
        MinLength = NA_real_,
        Q1Length = NA_real_,
        MedianLength = NA_real_,
        MeanLength = NA_real_,
        Q3Length = NA_real_,
        MaxLength = NA_real_,
        FractionLt18 = NA_real_,
        Fraction18to32 = NA_real_,
        FractionGt32 = NA_real_,
        AdapterSeedDetectedFraction = NA_real_,
        PolyA3pFraction = NA_real_,
        PolyT3pFraction = NA_real_,
        PolyG3pFraction = NA_real_,
        PolyC3pFraction = NA_real_,
        MostCommon5pKmer = NA_character_,
        MostCommon5pKmerFraction = NA_real_,
        MostCommon3pKmer = NA_character_,
        MostCommon3pKmerFraction = NA_real_
      ),
    raw_fastq_audit_by_group = tibble(
      Group = character(), SampleN = integer(), ReadsExamined_mean = numeric(),
      MeanLength_mean = numeric(), MeanLength_sd = numeric(), MedianLength_mean = numeric(),
      FractionLt18_mean = numeric(), Fraction18to32_mean = numeric(), FractionGt32_mean = numeric(),
      AdapterSeedDetectedFraction_mean = numeric(), PolyA3pFraction_mean = numeric(),
      PolyT3pFraction_mean = numeric(), PolyG3pFraction_mean = numeric(), PolyC3pFraction_mean = numeric()
    ),
    raw_fastq_length_distribution = tibble(
      Sample = character(), Group = character(), Length = integer(), ReadN = integer(), Fraction = numeric()
    ),
    raw_fastq_end_base = tibble(
      Sample = character(), Group = character(), End = character(), Base = character(), ReadN = integer(), Fraction = numeric()
    ),
    raw_fastq_kmer = tibble(
      Sample = character(), Group = character(), KmerSide = character(), Kmer = character(), ReadN = integer(), Fraction = numeric()
    )
  )
}

read_fastq_sequences_limited <- function(file_path, max_reads = RAW_AUDIT_MAX_READS, chunk_reads = 10000L) {
  if (!file.exists(file_path) || max_reads <= 0) return(character())
  con <- if (str_ends(file_path, fixed(".gz"))) gzfile(file_path, open = "rt") else file(file_path, open = "rt")
  on.exit(close(con), add = TRUE)

  out <- character()
  while (length(out) < max_reads) {
    read_target <- min(chunk_reads, max_reads - length(out))
    lines <- readLines(con, n = read_target * 4L, warn = FALSE)
    if (length(lines) == 0) break
    seq_idx <- seq(2L, length(lines), by = 4L)
    seq_lines <- lines[seq_idx]
    if (length(seq_lines) > 0) out <- c(out, seq_lines)
    if (length(lines) < read_target * 4L) break
  }

  normalize_filter_sequence(out)
}

top_kmer_table <- function(seqs, sample_id, group_id, side = c("5prime", "3prime"), k = RAW_AUDIT_KMER, top_n = 25L) {
  side <- match.arg(side)
  seq_len_vec <- nchar(seqs)
  seqs <- seqs[seq_len_vec >= k]
  if (length(seqs) == 0) {
    return(tibble(Sample = character(), Group = character(), KmerSide = character(), Kmer = character(), ReadN = integer(), Fraction = numeric()))
  }
  kmers <- if (side == "5prime") substr(seqs, 1L, k) else substr(seqs, nchar(seqs) - k + 1L, nchar(seqs))
  kmers <- kmers[str_detect(kmers, "^[ACGTN]+$")]
  if (length(kmers) == 0) {
    return(tibble(Sample = character(), Group = character(), KmerSide = character(), Kmer = character(), ReadN = integer(), Fraction = numeric()))
  }

  kmer_tbl <- tibble(Kmer = kmers) %>%
    count(Kmer, name = "ReadN") %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      KmerSide = side,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    arrange(desc(ReadN), Kmer) %>%
    slice_head(n = top_n) %>%
    select(Sample, Group, KmerSide, Kmer, ReadN, Fraction)

  kmer_tbl
}

audit_one_raw_fastq <- function(sample_id, group_id, raw_fastq, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  raw_fastq <- as.character(raw_fastq %||% "")
  exists_flag <- nzchar(raw_fastq) && file.exists(raw_fastq)
  adapter_seed <- substr(Sys.getenv("ADAPTER3", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"), 1L, max(8L, min(k, 12L)))

  if (!exists_flag) {
    return(list(
      summary = empty_raw_fastq_audit_tables(tibble(Sample = sample_id, Group = group_id, RawFastq = raw_fastq))$raw_fastq_audit_by_sample %>%
        mutate(Status = "missing_fastq"),
      length_distribution = tibble(),
      end_base = tibble(),
      kmer = tibble()
    ))
  }

  seqs <- tryCatch(
    read_fastq_sequences_limited(raw_fastq, max_reads = max_reads),
    error = function(e) character()
  )
  seqs <- seqs[nzchar(seqs)]
  read_n <- length(seqs)
  if (read_n == 0) {
    return(list(
      summary = tibble(
        Sample = sample_id,
        Group = group_id,
        RawFastq = raw_fastq,
        RawFastqExists = TRUE,
        Status = "empty_or_unreadable",
        MaxReadsRequested = max_reads,
        ReadsExamined = 0L,
        KmerK = k,
        MinLength = NA_real_,
        Q1Length = NA_real_,
        MedianLength = NA_real_,
        MeanLength = NA_real_,
        Q3Length = NA_real_,
        MaxLength = NA_real_,
        FractionLt18 = NA_real_,
        Fraction18to32 = NA_real_,
        FractionGt32 = NA_real_,
        AdapterSeedDetectedFraction = NA_real_,
        PolyA3pFraction = NA_real_,
        PolyT3pFraction = NA_real_,
        PolyG3pFraction = NA_real_,
        PolyC3pFraction = NA_real_,
        MostCommon5pKmer = NA_character_,
        MostCommon5pKmerFraction = NA_real_,
        MostCommon3pKmer = NA_character_,
        MostCommon3pKmerFraction = NA_real_
      ),
      length_distribution = tibble(),
      end_base = tibble(),
      kmer = tibble()
    ))
  }

  len <- nchar(seqs)
  first_base <- substr(seqs, 1L, 1L)
  last_base <- substr(seqs, len, len)
  kmer_tbl <- bind_rows(
    top_kmer_table(seqs, sample_id, group_id, side = "5prime", k = k),
    top_kmer_table(seqs, sample_id, group_id, side = "3prime", k = k)
  )
  top5 <- kmer_tbl %>% filter(KmerSide == "5prime") %>% arrange(desc(ReadN), Kmer) %>% slice(1)
  top3 <- kmer_tbl %>% filter(KmerSide == "3prime") %>% arrange(desc(ReadN), Kmer) %>% slice(1)

  summary_tbl <- tibble(
    Sample = sample_id,
    Group = group_id,
    RawFastq = raw_fastq,
    RawFastqExists = TRUE,
    Status = "completed",
    MaxReadsRequested = max_reads,
    ReadsExamined = read_n,
    KmerK = k,
    MinLength = min(len, na.rm = TRUE),
    Q1Length = as.numeric(stats::quantile(len, 0.25, na.rm = TRUE, names = FALSE)),
    MedianLength = stats::median(len, na.rm = TRUE),
    MeanLength = mean(len, na.rm = TRUE),
    Q3Length = as.numeric(stats::quantile(len, 0.75, na.rm = TRUE, names = FALSE)),
    MaxLength = max(len, na.rm = TRUE),
    FractionLt18 = mean(len < 18, na.rm = TRUE),
    Fraction18to32 = mean(len >= 18 & len <= 32, na.rm = TRUE),
    FractionGt32 = mean(len > 32, na.rm = TRUE),
    AdapterSeedDetectedFraction = if (nzchar(adapter_seed)) mean(str_detect(seqs, fixed(adapter_seed)), na.rm = TRUE) else NA_real_,
    PolyA3pFraction = mean(str_detect(seqs, "A{5,}$"), na.rm = TRUE),
    PolyT3pFraction = mean(str_detect(seqs, "T{5,}$"), na.rm = TRUE),
    PolyG3pFraction = mean(str_detect(seqs, "G{5,}$"), na.rm = TRUE),
    PolyC3pFraction = mean(str_detect(seqs, "C{5,}$"), na.rm = TRUE),
    MostCommon5pKmer = if (nrow(top5) > 0) top5$Kmer[[1]] else NA_character_,
    MostCommon5pKmerFraction = if (nrow(top5) > 0) top5$Fraction[[1]] else NA_real_,
    MostCommon3pKmer = if (nrow(top3) > 0) top3$Kmer[[1]] else NA_character_,
    MostCommon3pKmerFraction = if (nrow(top3) > 0) top3$Fraction[[1]] else NA_real_
  )

  length_distribution <- tibble(Length = len) %>%
    count(Length, name = "ReadN") %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    select(Sample, Group, Length, ReadN, Fraction) %>%
    arrange(Sample, Length)

  end_base <- bind_rows(
    tibble(End = "5prime", Base = first_base),
    tibble(End = "3prime", Base = last_base)
  ) %>%
    filter(str_detect(Base, "^[ACGTN]$")) %>%
    count(End, Base, name = "ReadN") %>%
    group_by(End) %>%
    mutate(
      Sample = sample_id,
      Group = group_id,
      Fraction = ReadN / sum(ReadN, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Sample, Group, End, Base, ReadN, Fraction)

  list(summary = summary_tbl, length_distribution = length_distribution, end_base = end_base, kmer = kmer_tbl)
}

build_raw_fastq_audit_tables <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) return(empty_raw_fastq_audit_tables())
  audits <- purrr::map(seq_len(nrow(sample_info_tbl)), function(i) {
    audit_one_raw_fastq(
      sample_id = as.character(sample_info_tbl$Sample[[i]]),
      group_id = as.character(sample_info_tbl$Group[[i]]),
      raw_fastq = as.character(sample_info_tbl$RawFastq[[i]]),
      max_reads = max_reads,
      k = k
    )
  })

  empty_tables <- empty_raw_fastq_audit_tables()
  raw_fastq_audit_by_sample <- bind_rows(empty_tables$raw_fastq_audit_by_sample, map(audits, "summary"))
  raw_fastq_length_distribution <- bind_rows(empty_tables$raw_fastq_length_distribution, map(audits, "length_distribution"))
  raw_fastq_end_base <- bind_rows(empty_tables$raw_fastq_end_base, map(audits, "end_base"))
  raw_fastq_kmer <- bind_rows(empty_tables$raw_fastq_kmer, map(audits, "kmer"))

  raw_fastq_audit_by_group <- raw_fastq_audit_by_sample %>%
    filter(Status == "completed") %>%
    group_by(Group) %>%
    summarise(
      SampleN = n(),
      ReadsExamined_mean = mean(ReadsExamined, na.rm = TRUE),
      MeanLength_mean = mean(MeanLength, na.rm = TRUE),
      MeanLength_sd = stats::sd(MeanLength, na.rm = TRUE),
      MedianLength_mean = mean(MedianLength, na.rm = TRUE),
      FractionLt18_mean = mean(FractionLt18, na.rm = TRUE),
      Fraction18to32_mean = mean(Fraction18to32, na.rm = TRUE),
      FractionGt32_mean = mean(FractionGt32, na.rm = TRUE),
      AdapterSeedDetectedFraction_mean = mean(AdapterSeedDetectedFraction, na.rm = TRUE),
      PolyA3pFraction_mean = mean(PolyA3pFraction, na.rm = TRUE),
      PolyT3pFraction_mean = mean(PolyT3pFraction, na.rm = TRUE),
      PolyG3pFraction_mean = mean(PolyG3pFraction, na.rm = TRUE),
      PolyC3pFraction_mean = mean(PolyC3pFraction, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    raw_fastq_audit_by_sample = raw_fastq_audit_by_sample,
    raw_fastq_audit_by_group = raw_fastq_audit_by_group,
    raw_fastq_length_distribution = raw_fastq_length_distribution,
    raw_fastq_end_base = raw_fastq_end_base,
    raw_fastq_kmer = raw_fastq_kmer
  )
}

build_raw_fastq_audit_manifest <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  if (is.null(sample_info_tbl) || nrow(sample_info_tbl) == 0) {
    return(tibble(Sample = character(), Group = character(), RawFastq = character(), RawFastqExists = logical()))
  }

  sample_info_tbl %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      RawFastq = as.character(RawFastq),
      RawFastqExists = file.exists(RawFastq),
      RawFastqSize = map_dbl(RawFastq, ~if (file.exists(.x)) file.info(.x)$size else NA_real_),
      RawFastqMTime = map_chr(RawFastq, ~if (file.exists(.x)) as.character(file.info(.x)$mtime) else ""),
      MaxReadsRequested = max_reads,
      KmerK = k
    ) %>%
    arrange(Group, Sample)
}

raw_fastq_audit_manifest_matches <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  a <- as.data.frame(a, stringsAsFactors = FALSE)
  b <- as.data.frame(b, stringsAsFactors = FALSE)
  rownames(a) <- NULL
  rownames(b) <- NULL
  identical(a, b)
}

load_or_build_raw_fastq_audit_tables <- function(sample_info_tbl, max_reads = RAW_AUDIT_MAX_READS, k = RAW_AUDIT_KMER) {
  manifest <- build_raw_fastq_audit_manifest(sample_info_tbl, max_reads = max_reads, k = k)
  cache_dir <- file.path(INTERNAL_CACHE_DIR, "raw_fastq_audit")
  cache_file <- file.path(cache_dir, "raw_fastq_audit_cache.rds")

  if (RAW_AUDIT_REFRESH != "1" && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(cached) && raw_fastq_audit_manifest_matches(cached$manifest, manifest) && is.list(cached$tables)) {
      return(cached$tables)
    }
  }

  tables <- build_raw_fastq_audit_tables(sample_info_tbl, max_reads = max_reads, k = k)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    saveRDS(list(manifest = manifest, tables = tables), cache_file),
    error = function(e) NULL
  )
  tables
}

complete_length_series <- function(df, len, category, clean_total) {
  if (nrow(df)==0) out <- tibble(Length=len,Reads=0) else
    out <- df %>% group_by(Length) %>% summarise(Reads=sum(Reads,na.rm=TRUE),.groups="drop") %>%
      right_join(tibble(Length=len),by="Length") %>% mutate(Reads=replace_na(Reads,0)) %>% arrange(Length)
  out %>% mutate(Category=category, RPM=if(clean_total>0) 1e6*Reads/clean_total else rep(0,n())) %>% select(Category,Length,Reads,RPM)
}

compute_pie_tables <- function(summary_top_df) {
  summary_top_df %>% filter(!is.na(MainCategory)) %>%
    mutate(PieCategory=if_else(as.character(MainCategory)%in%radar_levels,as.character(MainCategory),"other")) %>%
    group_by(Sample,Group,PieCategory) %>% summarise(Reads=sum(Reads),.groups="drop") %>%
    group_by(Sample,Group) %>% complete(PieCategory=unique(c(radar_levels,"other")),fill=list(Reads=0)) %>%
    mutate(Percent=100*Reads/sum(Reads),PieCategory=factor(PieCategory,levels=unique(c(radar_levels,"other")))) %>% ungroup()
}

build_sector_radar_df <- function(input_df, radar_max, sector_levels=radar_levels, sector_labels=radar_labels, sector_palette=radar_sector_palette, bottom_anchor_category="rsRNA", angle_offset=0) {
  n <- length(sector_levels); sw <- 2*pi/n
  bai <- match(bottom_anchor_category,sector_levels); if(is.na(bai)) bai <- ceiling(n/2)
  label_palette <- make_sector_label_palette(sector_palette[sector_levels], amount = 0.28)
  input_df %>% mutate(PieCategory=factor(as.character(PieCategory),levels=sector_levels)) %>%
    complete(PieCategory=factor(sector_levels,levels=sector_levels),fill=list(Reads=0,Percent=0)) %>% arrange(PieCategory) %>%
    mutate(SectorId=row_number(), CenterAngle=-pi/2+(bai-SectorId)*sw+angle_offset, StartAngle=CenterAngle+sw/2, EndAngle=CenterAngle-sw/2,
           Radius=Percent/radar_max, DisplayLabel=unname(sector_labels[as.character(PieCategory)]),
           DisplayValue=sprintf("%.2f",Percent), ValueColor=mii_text_primary, LabelColor=unname(label_palette[as.character(PieCategory)]),
           UpperHalf=sin(CenterAngle)>=0,
           TextAngle=((if_else(UpperHalf,CenterAngle*180/pi-90,CenterAngle*180/pi+90)+180)%%360)-180)
}

build_sector_text_df <- function(pdf,tcol,rv,ccol) tibble(Label=pdf[[tcol]],TextX=rv*cos(pdf$CenterAngle),TextY=rv*sin(pdf$CenterAngle),TextAngle=pdf$TextAngle,TextColor=pdf[[ccol]])

build_sector_arc_text_df <- function(pdf,tcol,rv,sw,ccol) {
  pmap_dfr(list(pdf$CenterAngle,pdf[[tcol]],pdf[[ccol]],pdf$UpperHalf), function(ca,tv,tc,uh) {
    ch <- strsplit(as.character(tv),"",fixed=TRUE)[[1]]; ch <- ch[nzchar(ch)]; if(length(ch)==0) return(tibble())
    ts <- min(sw*0.72, sw*0.068*max(length(ch)-1,1))
    if(uh) { ca2<-seq(ca+ts/2,ca-ts/2,length.out=length(ch)); ta<-ca2*180/pi-90 }
    else   { ca2<-seq(ca-ts/2,ca+ts/2,length.out=length(ch)); ta<-ca2*180/pi+90 }
    tibble(Label=ch,TextX=rv*cos(ca2),TextY=rv*sin(ca2),TextAngle=((ta+180)%%360)-180,TextColor=tc)
  })
}

build_sector_polygon_df <- function(pdf,ir,or,use_data=FALSE,pc=160) {
  pmap_dfr(list(pdf$SectorId,pdf$PieCategory,pdf$StartAngle,pdf$EndAngle,pdf$Radius), function(sid,pcat,sa,ea,rv) {
    cor <- if(use_data) rv else or; if(is.na(cor)||cor<=ir+1e-8) return(tibble(x=double(),y=double(),Group=character(),PieCategory=character()))
    to <- seq(sa,ea,length.out=pc); od <- tibble(x=cor*cos(to),y=cor*sin(to))
    pd <- if(ir<=0) bind_rows(tibble(x=0,y=0),od,tibble(x=0,y=0)) else {
      ti <- seq(ea,sa,length.out=pc); id <- tibble(x=ir*cos(ti),y=ir*sin(ti)); bind_rows(od,id) }
    pd %>% mutate(Group=paste0("sector_",sid,"_",if(use_data)"wedge"else"ring"), PieCategory=pcat)
  })
}

build_circle_path_df <- function(rv,pc=240) { th<-seq(0,2*pi,length.out=pc); map_dfr(rv,function(r) tibble(r=r,x=r*cos(th),y=r*sin(th))) }

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
    geom_path(data=cdf,aes(x=x,y=y,group=r),color=mii_neutral_line,linewidth=0.55) +
    geom_polygon(data=wdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.4,linejoin="round",alpha=0.86) +
    geom_polygon(data=rdf,aes(x=x,y=y,group=Group,fill=PieCategory),color="white",linewidth=0.45,linejoin="round",alpha=1) +
    geom_segment(data=sdf,aes(x=x,y=y,xend=xend,yend=yend),color=mii_neutral_line,linewidth=0.75,lineend="round") +
    geom_text(data=ldf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=label_size,fontface="bold",inherit.aes=FALSE) +
    geom_text(data=vdf,aes(x=TextX,y=TextY,label=Label,angle=TextAngle,color=TextColor),size=value_size,fontface="bold",inherit.aes=FALSE) +
    annotate("text",x=-1.60,y=0,label=axis_label,angle=90,size=4.0,color=mii_axis_label_color,fontface="bold") +
    scale_fill_manual(values=sector_palette,drop=FALSE) + scale_color_identity() +
    coord_fixed(xlim=coord_xlim,ylim=coord_ylim,clip="off") + theme_void() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5,face="bold",size=15,color=mii_text_muted),
          plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
          plot.margin=margin(10,20,10,28)) + labs(title=title_text)
}

plot_sector_radar_sample <- function(sample_id, pie_df, rmax, rbrk, view="merged") {
  pdf <- pie_df %>% filter(as.character(Sample)==sample_id) %>% build_sector_radar_df(radar_max=rmax)
  plot_sector_radar_core(pdf, paste0("sncRNA Proportion of ",sample_id,view_title_suffix(view)), rmax, rbrk)
}

prepare_group_sector_radar <- function(input_df, radar_max, sector_levels = radar_levels, sector_labels = radar_labels,
  bottom_anchor_category = "rsRNA", sector_palette = radar_sector_palette, angle_offset = 0) {
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
      DisplayValue = if_else(Percent > 0, sprintf("%.2f", Percent), ""),
      ValueColor = mii_text_primary,
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
    geom_path(data = cdf, aes(x = x, y = y, group = r), color = mii_neutral_line, linewidth = 0.55) +
    geom_polygon(data = wedge_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.38, linejoin = "round", alpha = 0.86) +
    geom_polygon(data = ring_df, aes(x = x, y = y, group = Group, fill = FillColor), color = "white", linewidth = 0.45, linejoin = "round", alpha = 1) +
    geom_segment(data = sdf, aes(x = x, y = y, xend = xend, yend = yend), color = mii_neutral_line, linewidth = 0.75, lineend = "round") +
    geom_text(data = ldf, aes(x = TextX, y = TextY, label = Label, angle = TextAngle, color = TextColor), size = label_size, fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = vdf %>% filter(Label != ""), aes(x = TextX, y = TextY, label = Label, angle = TextAngle), size = value_size, fontface = "bold", inherit.aes = FALSE) +
    annotate("text", x = -1.60, y = 0, label = axis_label, angle = 90, size = 4.0, color = mii_axis_label_color, fontface = "bold") +
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
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, color = mii_text_muted),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 20, 10, 28)
    ) +
    labs(title = title_text, fill = NULL)
}

prepare_trend_plot_df <- function(df, levels, group_levels) {
  if (nrow(df) == 0) {
    return(tibble(Group = factor(character(), levels = group_levels), TrendCategory = factor(character(), levels = levels), Reads = numeric(), Percent = numeric()))
  }
  df %>% mutate(Group = factor(as.character(Group), levels = group_levels), TrendCategory = factor(as.character(TrendCategory), levels = levels)) %>%
    filter(!is.na(TrendCategory)) %>%
    complete(Group = factor(group_levels, levels = group_levels), TrendCategory = factor(levels, levels = levels), fill = list(Reads = 0, Percent = 0)) %>%
    arrange(Group, TrendCategory)
}

build_trend_bar_plot <- function(df, levels, value_col, title, y_label, group_levels, fill_palette = main_palette) {
  plot_df <- prepare_trend_plot_df(df, levels, group_levels)
  display_df <- plot_df %>% filter(.data[[value_col]] > 0)
  if (nrow(display_df) == 0) return(plot_note(title, "No trend data available."))
  p <- ggplot(display_df, aes(x = Group, y = .data[[value_col]], fill = TrendCategory)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = fill_palette[levels], drop = TRUE) +
    labs(x = "Group", y = y_label, title = title) +
    mii_theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal") +
    guides(fill = guide_legend(ncol = if (any(nchar(levels) > 22, na.rm = TRUE)) 2 else 4, byrow = TRUE))
  if (identical(value_col, "Reads")) p <- p + scale_y_continuous(labels = scales::comma)
  p
}

pvalue_to_stars <- function(pvalue) {
  case_when(
    is.na(pvalue) ~ NA_character_,
    pvalue < 0.0001 ~ "****",
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

format_pvalue_label <- function(pvalue, show_pvalue = FALSE) {
  stars <- pvalue_to_stars(pvalue)
  if (is.na(stars)) return(NA_character_)
  if (!show_pvalue || is.na(pvalue)) return(stars)
  paste0(stars, "\n", format.pval(pvalue, digits = 2, eps = 0.001))
}

calculate_two_group_pvalue <- function(df, value_col = "Value", group_col = "Group", group_levels = NULL, min_replicates = 2) {
  if (get0("PANDORA_ANALYSIS_STAGE",inherits=TRUE,ifnotfound="describe") != "differential" || Sys.getenv("PANDORA_DESIGN") != "independent") return(NA_real_)
  if (nrow(df) == 0 || !all(c(value_col, group_col) %in% names(df))) return(NA_real_)
  stat_df <- df %>%
    mutate(
      .group = as.character(.data[[group_col]]),
      .value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    filter(!is.na(.group), is.finite(.value))

  groups <- if (is.null(group_levels)) unique(stat_df$.group) else as.character(group_levels)
  groups <- groups[groups %in% unique(stat_df$.group)]
  if (length(groups) != 2) return(NA_real_)

  x <- stat_df$.value[stat_df$.group == groups[[1]]]
  y <- stat_df$.value[stat_df$.group == groups[[2]]]
  if (length(x) < min_replicates || length(y) < min_replicates) return(NA_real_)
  if (length(unique(c(x, y))) < 2) return(1)

  tryCatch(
    suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

complete_sample_metric <- function(df, value_col, sample_info, group_levels) {
  base_df <- sample_info %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group))

  if (nrow(df) == 0 || !value_col %in% names(df)) {
    return(base_df %>% mutate(Value = 0, Group = factor(Group, levels = group_levels)))
  }

  metric_df <- df %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Value = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    group_by(Sample, Group) %>%
    summarise(Value = sum(Value, na.rm = TRUE), .groups = "drop")

  base_df %>%
    left_join(metric_df, by = c("Sample", "Group")) %>%
    mutate(
      Value = replace_na(Value, 0),
      Group = factor(Group, levels = group_levels)
    )
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
    geom_col(width = 0.64, color = "white", linewidth = 0.25) +
    geom_errorbar(aes(ymin = pmax(Mean - SEM, 0), ymax = Mean + SEM), width = 0.16, linewidth = 0.45, color = mii_text_primary) +
    geom_point(
      data = plot_df,
      aes(x = Group, y = Value),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.065, height = 0, seed = 11),
      size = 2.2,
      shape = 21,
      stroke = 0.25,
      color = mii_text_primary,
      fill = "white",
      alpha = 0.9
    ) +
    geom_text(
      aes(label = sprintf(paste0("%.", label_digits, "f"), Mean), y = Mean + SEM + y_pad * 0.12),
      size = 3.2,
      color = mii_text_primary,
      vjust = 0
    ) +
    scale_fill_manual(values = fill_palette, drop = FALSE) +
    scale_y_continuous(
      labels = if (comma_labels) scales::comma else waiver(),
      limits = c(0, upper_y),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = y_label, title = title) +
    mii_theme_bw() +
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
      color = mii_text_primary,
      check_overlap = TRUE
    )
    axis_label_layer <- NULL
    y_scale <- scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.03, 0.03)))
  } else {
    label_layer <- geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(ifelse(count > 3.5, sprintf("%.1f%%", count), ""))),
      size = 3,
      color = mii_text_primary,
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
      size = 3, color = mii_text_primary, hjust = 1, inherit.aes = FALSE
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
    mii_theme_bw() + theme(
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

plot_percent_pie <- function(df, category_col, palette, title_text, empty_note) {
  df <- df %>% filter(.data[[category_col]] != "Unknown", Reads > 0 | Percent > 0)
  if (nrow(df) == 0) return(plot_note(title_text, empty_note))
  df <- df %>% mutate(!!category_col := factor(.data[[category_col]], levels = names(palette)))
  ggplot(df, aes(x = "", y = Percent, fill = .data[[category_col]])) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = palette, drop = FALSE) +
    geom_text(aes(label = ifelse(Percent >= 1, sprintf("%.1f%%", Percent), "")), position = position_stack(vjust = 0.5), size = 3.2) +
    labs(x = NULL, y = NULL, title = title_text, fill = NULL) +
    theme_void() + theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
}

safe_scale_rows <- function(x) {
  x <- as.matrix(x)
  if (length(x) == 0 || nrow(x) == 0 || ncol(x) == 0) return(x)
  scaled <- tryCatch(t(scale(t(x))), error = function(e) matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

build_feature_count_matrix <- function(df, feature_col, sample_col = "Sample", count_col = "Reads", sample_order = NULL, feature_order = NULL, min_total = 0) {
  if (is.null(df) || nrow(df) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = character(), feature_order = character()))
  }

  required_cols <- c(feature_col, sample_col, count_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns for build_feature_count_matrix(): ", paste(missing_cols, collapse = ", "))
  }

  count_tbl <- df %>%
    transmute(
      Sample = as.character(.data[[sample_col]]),
      Feature = as.character(.data[[feature_col]]),
      Reads = clean_numeric(.data[[count_col]])
    ) %>%
    filter(!is.na(Sample), Sample != "", !is.na(Feature), Feature != "", !is.na(Reads)) %>%
    group_by(Sample, Feature) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop")

  if (nrow(count_tbl) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = character(), feature_order = character()))
  }

  if (is.null(sample_order)) {
    sample_order <- count_tbl %>% distinct(Sample) %>% pull(Sample)
  } else {
    sample_order <- unique(as.character(sample_order))
  }

  if (is.null(feature_order)) {
    feature_order <- count_tbl %>%
      group_by(Feature) %>%
      summarise(TotalReads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
      filter(TotalReads > min_total) %>%
      arrange(desc(TotalReads), Feature) %>%
      pull(Feature)
  } else {
    feature_order <- unique(as.character(feature_order))
  }

  if (length(sample_order) == 0 || length(feature_order) == 0) {
    empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(list(matrix = empty_matrix, table = tibble(), sample_order = sample_order, feature_order = feature_order))
  }

  full_tbl <- tidyr::expand_grid(Sample = sample_order, Feature = feature_order) %>%
    left_join(count_tbl, by = c("Sample", "Feature")) %>%
    mutate(Reads = replace_na(Reads, 0), Sample = factor(Sample, levels = sample_order), Feature = factor(Feature, levels = feature_order))

  count_matrix <- xtabs(Reads ~ Feature + Sample, data = full_tbl)
  count_matrix <- as.matrix(count_matrix)
  count_matrix <- count_matrix[feature_order, sample_order, drop = FALSE]

  list(
    matrix = count_matrix,
    table = full_tbl %>% arrange(Feature, Sample),
    sample_order = sample_order,
    feature_order = feature_order
  )
}

resolve_pairwise_contrasts <- function(group_levels, contrast_specs = NULL) {
  group_levels <- unique(as.character(group_levels))
  if (length(group_levels) < 2) return(tibble(group1 = character(), group2 = character(), contrast = character()))

  if (is.null(contrast_specs)) {
    pairs <- utils::combn(group_levels, 2, simplify = FALSE)
    return(map_dfr(pairs, function(pair) tibble(group1 = pair[[1]], group2 = pair[[2]], contrast = paste0(pair[[2]], "_vs_", pair[[1]]))))
  }

  if (is.data.frame(contrast_specs)) {
    if (!all(c("group1", "group2") %in% names(contrast_specs))) {
      stop("contrast_specs data frame must contain group1 and group2 columns.")
    }
    return(as_tibble(contrast_specs) %>% transmute(group1 = as.character(group1), group2 = as.character(group2), contrast = if ("contrast" %in% names(.)) as.character(contrast) else paste0(group2, "_vs_", group1)))
  }

  if (is.list(contrast_specs)) {
    return(map_dfr(contrast_specs, function(pair) {
      pair <- as.character(pair)
      if (length(pair) < 2) return(tibble())
      tibble(group1 = pair[[1]], group2 = pair[[2]], contrast = if (length(pair) >= 3 && nzchar(pair[[3]])) pair[[3]] else paste0(pair[[2]], "_vs_", pair[[1]]))
    }))
  }

  if (is.character(contrast_specs)) {
    if (length(contrast_specs) == 2) {
      return(tibble(group1 = contrast_specs[[1]], group2 = contrast_specs[[2]], contrast = paste0(contrast_specs[[2]], "_vs_", contrast_specs[[1]])))
    }
  }

  stop("Unsupported contrast_specs format.")
}

classify_de_direction <- function(pvalue, log_fc, pvalue_cutoff = 0.05, lfc_cutoff = log2(1.5)) {
  case_when(
    is.na(pvalue) | is.na(log_fc) ~ "ns",
    pvalue < pvalue_cutoff & log_fc >= lfc_cutoff ~ "up",
    pvalue < pvalue_cutoff & log_fc <= -lfc_cutoff ~ "down",
    TRUE ~ "ns"
  )
}

de_significance_label <- function() {
  if (DE_SIGNIFICANCE %in% c("PVALUE", "P_VALUE", "P")) "PValue" else "FDR"
}

resolve_de_significance_values <- function(de_tbl) {
  if ("SignificanceMetric" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$SignificanceMetric))
  }
  preferred_col <- de_significance_label()
  if (preferred_col %in% names(de_tbl)) {
    return(as.numeric(de_tbl[[preferred_col]]))
  }
  if ("FDR" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$FDR))
  }
  if ("PValue" %in% names(de_tbl)) {
    return(as.numeric(de_tbl$PValue))
  }
  rep(NA_real_, nrow(de_tbl))
}

run_edgeR_pairwise_de <- function(count_matrix, sample_groups, group1, group2, feature_ids = NULL, min_replicates = 2) {
  if (!pkg_has("edgeR")) {
    return(list(
      skipped = TRUE,
      message = "edgeR package is not available.",
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  sample_groups <- as.character(sample_groups)
  keep_samples <- sample_groups %in% c(group1, group2)
  if (sum(keep_samples) < 4) {
    return(list(
      skipped = TRUE,
      message = paste0("Not enough samples for contrast ", group2, " vs ", group1, "."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  count_matrix <- as.matrix(count_matrix[, keep_samples, drop = FALSE])
  sample_groups <- sample_groups[keep_samples]
  if (min(table(sample_groups)) < min_replicates) {
    return(list(
      skipped = TRUE,
      message = paste0("Contrast ", group2, " vs ", group1, " needs at least ", min_replicates, " replicates per group."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  if (nrow(count_matrix)==0 || any(colSums(count_matrix)<=0)) {
    return(list(skipped=TRUE,message="insufficient_class_counts: zero-count library; biological replicates retained",result=tibble(),sample_groups=sample_groups,feature_ids=feature_ids))
  }
  group_factor <- factor(sample_groups, levels = c(group1, group2))
  base_mean_values <- rowMeans(count_matrix, na.rm = TRUE)
  base_mean_names <- rownames(count_matrix)
  if (is.null(base_mean_names) || anyNA(base_mean_names) || !all(nzchar(base_mean_names))) {
    if (!is.null(feature_ids) && length(feature_ids) == nrow(count_matrix)) {
      base_mean_names <- as.character(feature_ids)
    } else {
      base_mean_names <- as.character(seq_len(nrow(count_matrix)))
    }
  }
  names(base_mean_values) <- base_mean_names
  dge <- edgeR::DGEList(counts = count_matrix, group = group_factor)
  keep <- as.vector(edgeR::filterByExpr(dge, group = group_factor))
  if (!any(keep) || any(colSums(count_matrix[keep,,drop=FALSE])<=0)) {
    return(list(skipped=TRUE,message="insufficient_class_counts after expression filtering; biological replicates retained",result=tibble(),sample_groups=sample_groups,feature_ids=feature_ids))
  }
  dge <- edgeR::DGEList(counts = count_matrix[keep, , drop = FALSE], group = group_factor)

  if (nrow(dge) == 0) {
    return(list(
      skipped = TRUE,
      message = paste0("No features passed expression filtering for contrast ", group2, " vs ", group1, "."),
      result = tibble(),
      sample_groups = sample_groups,
      feature_ids = feature_ids
    ))
  }

  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~ 0 + group_factor)
  safe_group_names <- setNames(make.names(levels(group_factor)), levels(group_factor))
  colnames(design) <- unname(safe_group_names)
  de_method <- "edgeR_robust_glmQLFTest"
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  contrast <- limma::makeContrasts(
    contrasts = paste0(safe_group_names[[as.character(group2)]], "-", safe_group_names[[as.character(group1)]]),
    levels = design
  )
  test <- edgeR::glmQLFTest(fit, contrast = contrast)
  de_table <- edgeR::topTags(test, n = Inf)$table
  result <- de_table
  result <- tibble::rownames_to_column(as.data.frame(result), var = "Feature")

  if (!is.null(feature_ids) && length(feature_ids) == nrow(count_matrix)) {
    result <- result %>% left_join(tibble(Feature = as.character(feature_ids), FeatureID = as.character(feature_ids)), by = "Feature")
  }

  result <- result %>%
    mutate(baseMean = unname(base_mean_values[Feature])) %>%
    mutate(
      SignificanceMetric = if (de_significance_label() == "FDR") FDR else PValue,
      SignificanceMethod = de_significance_label(),
      DEMethod = de_method,
      direction = classify_de_direction(SignificanceMetric, logFC),
      rank_score = sign(coalesce(logFC, 0)) * -log10(coalesce(SignificanceMetric, 1) + 1e-300),
      contrast = paste0(group2, "_vs_", group1)
    ) %>%
    arrange(SignificanceMetric, PValue, desc(abs(logFC)))

  list(
    skipped = FALSE,
    message = NULL,
    result = result,
    sample_groups = sample_groups,
    feature_ids = feature_ids
  )
}

build_de_volcano_plot <- function(de_tbl, title_text, lfc_cutoff = log2(1.5), pvalue_cutoff = 0.05, label_n = 15) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No DE results available."))

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
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = mii_neutral_text) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = mii_neutral_text) +
    scale_color_manual(values = mii_direction_palette, drop = FALSE) +
    labs(x = "log2 Fold Change", y = paste0("-log10(", significance_label, ")"), title = title_text) +
    mii_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  if (pkg_has("ggrepel")) {
    p <- p + ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, na.rm = TRUE, show.legend = FALSE)
  } else {
    p <- p + geom_text(aes(label = label), size = 3, na.rm = TRUE, show.legend = FALSE, vjust = -0.4)
  }

  p
}

build_de_ma_plot <- function(de_tbl, title_text) {
  if (nrow(de_tbl) == 0) return(plot_note(title_text, "No DE results available."))

  ggplot(de_tbl, aes(x = log10(baseMean + 1), y = logFC, color = direction)) +
    geom_point(alpha = 0.72, size = 1.5) +
    scale_color_manual(values = mii_direction_palette, drop = FALSE) +
    labs(x = "log10(baseMean + 1)", y = "log2 Fold Change", title = title_text) +
    mii_theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")
}

build_de_robustness_table <- function(result_tbl, count_matrix, sample_order, sample_groups, group1, group2) {
  if (nrow(result_tbl) == 0 || nrow(count_matrix) == 0 || length(sample_order) == 0) return(tibble())
  sample_order <- as.character(sample_order)
  sample_groups <- as.character(sample_groups)
  if (!all(sample_order %in% colnames(count_matrix))) return(tibble())

  top_features <- result_tbl %>%
    arrange(SignificanceMetric, PValue, desc(abs(logFC))) %>%
    pull(Feature) %>%
    as.character() %>%
    unique()
  top_features <- top_features[top_features %in% rownames(count_matrix)]
  if (length(top_features) == 0) return(tibble())

  baseline_samples <- sample_order[sample_groups == group1]
  case_samples <- sample_order[sample_groups == group2]
  if (length(baseline_samples) == 0 || length(case_samples) == 0) return(tibble())

  purrr::map_dfr(top_features, function(feature_name) {
    values <- as.numeric(count_matrix[feature_name, sample_order, drop = TRUE])
    names(values) <- sample_order
    baseline_values <- values[baseline_samples]
    case_values <- values[case_samples]
    baseline_total <- sum(baseline_values, na.rm = TRUE)
    case_total <- sum(case_values, na.rm = TRUE)
    all_total <- baseline_total + case_total
    max_sample <- names(values)[which.max(values)]
    max_reads <- max(values, na.rm = TRUE)
    max_group <- if (max_sample %in% baseline_samples) group1 else group2
    max_baseline_fraction <- if (baseline_total > 0) max(baseline_values, na.rm = TRUE) / baseline_total else NA_real_
    max_case_fraction <- if (case_total > 0) max(case_values, na.rm = TRUE) / case_total else NA_real_
    max_group_fraction <- max(c(max_baseline_fraction, max_case_fraction), na.rm = TRUE)
    if (!is.finite(max_group_fraction)) max_group_fraction <- NA_real_

    de_row <- result_tbl %>% filter(as.character(Feature) == feature_name) %>% slice(1)
    direction <- as.character(de_row$direction[[1]])
    baseline_median <- stats::median(baseline_values, na.rm = TRUE)
    case_median <- stats::median(case_values, na.rm = TRUE)
    case_support_fraction <- case_when(
      direction == "up" ~ mean(case_values > baseline_median, na.rm = TRUE),
      direction == "down" ~ mean(case_values < baseline_median, na.rm = TRUE),
      TRUE ~ NA_real_
    )
    baseline_support_fraction <- case_when(
      direction == "up" ~ mean(baseline_values < case_median, na.rm = TRUE),
      direction == "down" ~ mean(baseline_values > case_median, na.rm = TRUE),
      TRUE ~ NA_real_
    )
    single_sample_dominated <- is.finite(max_group_fraction) & max_group_fraction >= 0.80
    weak_replicate_support <- direction %in% c("up", "down") & is.finite(case_support_fraction) & case_support_fraction < 0.50
    robustness_flag <- case_when(
      direction == "ns" ~ "not_significant",
      single_sample_dominated ~ "single_sample_dominated",
      weak_replicate_support ~ "weak_replicate_support",
      TRUE ~ "replicate_supported"
    )

    tibble(
      Feature = feature_name,
      direction = direction,
      logFC = de_row$logFC[[1]],
      SignificanceMetric = de_row$SignificanceMetric[[1]],
      PValue = de_row$PValue[[1]],
      FDR = if ("FDR" %in% names(de_row)) de_row$FDR[[1]] else NA_real_,
      BaselineGroup = group1,
      CaseGroup = group2,
      BaselineTotalReads = baseline_total,
      CaseTotalReads = case_total,
      TotalReads = all_total,
      BaselineDetectedSamples = sum(baseline_values > 0, na.rm = TRUE),
      CaseDetectedSamples = sum(case_values > 0, na.rm = TRUE),
      BaselineSampleN = length(baseline_samples),
      CaseSampleN = length(case_samples),
      BaselineMedianReads = baseline_median,
      CaseMedianReads = case_median,
      CaseDirectionSupportFraction = case_support_fraction,
      BaselineDirectionSupportFraction = baseline_support_fraction,
      MaxSample = max_sample,
      MaxSampleGroup = max_group,
      MaxSampleReads = max_reads,
      MaxSampleFractionTotal = if_else(all_total > 0, max_reads / all_total, NA_real_),
      MaxBaselineSampleFraction = max_baseline_fraction,
      MaxCaseSampleFraction = max_case_fraction,
      MaxGroupSampleFraction = max_group_fraction,
      SingleSampleDominated = single_sample_dominated,
      RobustnessFlag = robustness_flag,
      BaselineSampleReads = paste(paste0(baseline_samples, ":", baseline_values), collapse = ";"),
      CaseSampleReads = paste(paste0(case_samples, ":", case_values), collapse = ";")
    )
  })
}

build_de_robustness_plot <- function(robustness_tbl, title_text, label_n = 30) {
  if (nrow(robustness_tbl) == 0) return(plot_note(title_text, "No robustness table available."))
  plot_df <- robustness_tbl %>%
    arrange(desc(MaxGroupSampleFraction), SignificanceMetric) %>%
    slice_head(n = label_n) %>%
    mutate(
      Feature = factor(Feature, levels = rev(unique(Feature))),
      RobustnessFlag = factor(
        RobustnessFlag,
        levels = c("replicate_supported", "weak_replicate_support", "single_sample_dominated", "not_significant")
      )
    )

  ggplot(plot_df, aes(x = Feature, y = MaxGroupSampleFraction, fill = RobustnessFlag)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.15) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(
      values = c(
        replicate_supported = mii_green,
        weak_replicate_support = mii_orange,
        single_sample_dominated = mii_red,
        not_significant = mii_neutral_line
      ),
      drop = FALSE
    ) +
    labs(title = title_text, x = NULL, y = "Largest single-sample fraction within a group", fill = "Flag") +
    mii_theme_bw() +
    theme(legend.position = "top")
}

resolve_mii_contrast_specs <- function(sample_info_tbl = sample_info, group_levels_vec = group_levels) {
  contrast_file <- file.path(ROOT_DIR, "metadata", "contrast_pairs.tsv")
  age_group_design <- sample_info_tbl %>%
    transmute(Group = as.character(Group), AgeGroup = tolower(trimws(as.character(AgeGroup)))) %>%
    distinct() %>%
    filter(AgeGroup %in% c("young", "old"))
  young_groups <- age_group_design %>% filter(AgeGroup == "young") %>% pull(Group) %>% unique()
  old_groups <- age_group_design %>% filter(AgeGroup == "old") %>% pull(Group) %>% unique()

  if (file.exists(contrast_file)) {
    contrast_specs_raw <- read_tsv_flex(contrast_file)
  } else if (length(young_groups) == 1 && length(old_groups) == 1) {
    contrast_specs_raw <- tibble(
      group1 = young_groups[[1]],
      group2 = old_groups[[1]],
      contrast = paste0(old_groups[[1]], "_vs_", young_groups[[1]])
    )
  } else if (all(c(PROJECT_CONTRAST_BASELINE, PROJECT_CONTRAST_CASE) %in% as.character(group_levels_vec))) {
    contrast_specs_raw <- tibble(
      group1 = PROJECT_CONTRAST_BASELINE,
      group2 = PROJECT_CONTRAST_CASE,
      contrast = paste0(PROJECT_CONTRAST_CASE, "_vs_", PROJECT_CONTRAST_BASELINE)
    )
  } else if (all(c("1", "2") %in% as.character(group_levels_vec))) {
    contrast_specs_raw <- tibble(group1 = "1", group2 = "2")
  } else {
    contrast_specs_raw <- resolve_pairwise_contrasts(group_levels_vec)
  }

  resolve_pairwise_contrasts(group_levels_vec, contrast_specs_raw)
}

write_de_run_config <- function(contrast_specs, section = "00_de_config") {
  write_csv_safe(contrast_specs, file.path(section, "contrast_manifest.csv"))
  write_csv_safe(
    tibble(
      significance_method = de_significance_label(),
      significance_cutoff = 0.05,
      log2_fold_change_cutoff = log2(1.5),
      rule = paste0("direction != ns requires ", de_significance_label(), " < 0.05 and |logFC| >= log2(1.5)")
    ),
    file.path(section, "de_threshold_manifest.csv")
  )
  invisible(contrast_specs)
}

sports_annotation_counts <- function(summary_df = summary_all, clean_df = clean_reads_tbl) {
  # 使用上游已分配的 annotation counts，不把一条序列全额归给首个标签。
  summary_df %>%
    filter(
      !is.na(SubClass), SubClass != "-", SubClass != "",
      !Class %in% c("Clean_Reads", "Match_Genome"),
      !is_tsrna_nested_class(Class)
    ) %>%
    transmute(
      Sample = as.character(Sample), Group = as.character(Group),
      PrimaryAnnotation = as.character(SubClass),
      FeatureCategory = as.character(MainCategory), Reads
    ) %>%
    group_by(Sample, Group, PrimaryAnnotation, FeatureCategory) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_df %>% transmute(Sample = as.character(Sample), CleanReads), by = "Sample") %>%
    mutate(RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0)) %>%
    select(-CleanReads)
}

build_global_de_inputs <- function() {
  general_annotation_categories <- c("miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "ribozyme", "other")
  general_annotation_categories <- general_annotation_categories[general_annotation_categories %in% category_levels]
  if (length(general_annotation_categories) == 0) {
    general_annotation_categories <- c("miRNA", "piRNA", "ysRNA", "lncRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "ribozyme", "other")
  }

  category_input <- summary_top %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), Feature = as.character(MainCategory), Reads = Reads)

  subtype_input <- summary_all %>%
    filter(SubClass != "-", !MainCategory %in% c("tsRNA", "rsRNA"), !Class %in% c("Clean_Reads", "Match_Genome")) %>%
    transmute(Sample = as.character(Sample), Group = as.character(Group), Feature = paste0(as.character(MainCategory), "::", as.character(SubClass)), Reads = Reads)

  annotation_input <- build_annotation_de_input(general_annotation_categories)

  list(
    category_input = category_input,
    subtype_input = subtype_input,
    annotation_input = annotation_input,
    category_meta = tibble(Feature = as.character(category_levels), FeatureCategory = as.character(category_levels)),
    subtype_meta = subtype_input %>% distinct(Feature) %>% mutate(FeatureCategory = sub("::.*$", "", Feature)),
    annotation_meta = annotation_input %>% distinct(Feature, FeatureCategory)
  )
}

build_annotation_de_input <- function(category_name) {
  # 注释丰度沿用 SPORTS 的多命中分数分配；首个序列注释仅可作标签。
  sports_annotation_counts() %>%
    filter(FeatureCategory %in% category_name) %>%
    rename(Feature = PrimaryAnnotation)
}

build_tsrna_annotation_de_input <- function() {
  tsrna_cleavage_detail %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = make_company_tsrna_feature(SubClass, Cleavage),
      FeatureCategory = if_else(Origin == "Mitochondrial", "mt_tsRNA", "tsRNA"),
      Reads = Reads
    )
}

build_rsrna_subtype_de_input <- function() {
  rsrna_detail %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature = as.character(RsRNAClass),
      FeatureCategory = "rsRNA",
      Reads = Reads
    )
}

build_sequence_de_base <- function() {
  output_all %>%
    filter(!is.na(Sequence), Sequence != "", !is.na(Reads), Reads > 0) %>%
    mutate(
      Feature = as.character(Sequence),
      AnnotationLabel = if_else(
        is.na(PrimaryAnnotation) | PrimaryAnnotation == "" | PrimaryAnnotation == "NO_Annotation",
        "unannotated",
        as.character(PrimaryAnnotation)
      ),
      SequenceCategory = as.character(FeatureCategory),
      IsMito = is_mito_label(paste(Annotation, PrimaryAnnotation, sep = " "))
    ) %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(Group),
      Feature,
      Reads,
      AnnotationLabel,
      SequenceCategory,
      IsMito
    )
}

sequence_de_input_for <- function(sequence_de_base, category_name) {
  switch(
    category_name,
    miRNA = sequence_de_base %>% filter(SequenceCategory == "miRNA"),
    tsRNA = sequence_de_base %>% filter(SequenceCategory == "tsRNA", !IsMito),
    mt_tsRNA = sequence_de_base %>% filter(SequenceCategory == "tsRNA", IsMito),
    piRNA = sequence_de_base %>% filter(SequenceCategory == "piRNA"),
    rsRNA = sequence_de_base %>% filter(SequenceCategory == "rsRNA"),
    tibble()
  )
}

run_sequence_de_categories <- function(categories, layer_prefix = "sequence_de", sequence_de_base = build_sequence_de_base()) {
  categories <- unique(as.character(categories))
  category_labels <- c(miRNA = "miRNA", tsRNA = "tsRNA", mt_tsRNA = "mt-tsRNA", piRNA = "piRNA", rsRNA = "rsRNA")
  runs <- purrr::map(categories, function(category) {
    seq_input_raw <- sequence_de_input_for(sequence_de_base, category)
    seq_input <- seq_input_raw %>%
      group_by(Sample, Group, Feature) %>%
      summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop")

    seq_meta <- seq_input_raw %>%
      group_by(Feature) %>%
      summarise(
        FeatureCategory = category,
        TopAnnotation = {
          ann_counts <- sort(table(AnnotationLabel), decreasing = TRUE)
          if (length(ann_counts) == 0) "unannotated" else names(ann_counts)[[1]]
        },
        TotalReads = sum(Reads, na.rm = TRUE),
        .groups = "drop"
      )

    run_layer_de(
      layer_name = file.path(layer_prefix, paste0(category, "_seq")),
      layer_title = paste0(category_labels[[category]] %||% category, " sequence-level DE"),
      layer_input = seq_input,
      feature_order = NULL,
      feature_meta = seq_meta,
      min_total = 5
    )
  })
  names(runs) <- categories

  list(
    runs = runs,
    results = bind_rows(map(runs, "results")),
    summary = bind_rows(map(runs, "summary")),
    manifest = bind_rows(map(runs, "manifest")),
    robustness_summary = bind_rows(map(runs, "robustness_summary"))
  )
}

summarise_de_robustness_flags <- function(robustness_tbl) {
  if (is.null(robustness_tbl) || nrow(robustness_tbl) == 0 || !"RobustnessFlag" %in% names(robustness_tbl)) {
    return(tibble())
  }

  robustness_tbl %>%
    mutate(
      RobustnessFlag = as.character(RobustnessFlag),
      direction = as.character(direction),
      is_significant = direction != "ns"
    ) %>%
    group_by(layer, layer_title, contrast, RobustnessFlag) %>%
    summarise(
      FeatureN = n(),
      SignificantFeatureN = sum(is_significant, na.rm = TRUE),
      UpN = sum(direction == "up", na.rm = TRUE),
      DownN = sum(direction == "down", na.rm = TRUE),
      MedianMaxGroupSampleFraction = median(MaxGroupSampleFraction, na.rm = TRUE),
      MaxMaxGroupSampleFraction = max(MaxGroupSampleFraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(layer, contrast, desc(SignificantFeatureN), RobustnessFlag)
}

run_layer_de <- function(layer_name, layer_title, layer_input, feature_order = NULL, feature_meta = NULL, min_total = 0) {
  layer_results <- list()
  layer_summary <- list()
  layer_manifest <- list()
  layer_robustness <- list()
  layer_robustness_summary <- list()

  for (i in seq_len(nrow(contrast_specs))) {
    group1 <- as.character(contrast_specs$group1[[i]])
    group2 <- as.character(contrast_specs$group2[[i]])
    contrast_label <- as.character(contrast_specs$contrast[[i]])
    pair_samples <- sample_info %>%
      filter(as.character(Group) %in% c(group1, group2)) %>%
      mutate(Group = factor(as.character(Group), levels = c(group1, group2))) %>%
      arrange(Group, Sample)

    pair_sample_order <- as.character(pair_samples$Sample)
    pair_sample_groups <- as.character(pair_samples$Group)
    pair_section <- file.path(layer_name, safe_dir_name(contrast_label))

    layer_subset <- layer_input %>%
      filter(Sample %in% pair_sample_order) %>%
      select(Sample, Group, Feature, Reads)

    matrix_info <- build_feature_count_matrix(
      layer_subset,
      feature_col = "Feature",
      sample_order = pair_sample_order,
      feature_order = feature_order,
      min_total = min_total
    )

    if (nrow(matrix_info$matrix) == 0 || ncol(matrix_info$matrix) == 0) {
      skip_message <- paste0("No features available for ", layer_title, " in contrast ", contrast_label, ".")
      layer_manifest[[length(layer_manifest) + 1]] <- tibble(
        layer = layer_name,
        contrast = contrast_label,
        group1 = group1,
        group2 = group2,
        status = "skipped",
        message = skip_message,
        sample_n = length(pair_sample_order),
        feature_n = 0,
        significant_n = 0,
        up_n = 0,
        down_n = 0
      )
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), skip_message),
        pair_section,
        paste0("skipped_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
      next
    }

    feature_lookup <- tibble(
      FeatureInternal = sprintf("feature_%06d", seq_len(nrow(matrix_info$matrix))),
      Feature = rownames(matrix_info$matrix)
    )
    de_matrix <- matrix_info$matrix
    rownames(de_matrix) <- feature_lookup$FeatureInternal

    de_run <- run_edgeR_pairwise_de(
      count_matrix = de_matrix,
      sample_groups = pair_sample_groups,
      group1 = group1,
      group2 = group2,
      feature_ids = feature_lookup$FeatureInternal,
      min_replicates = 2
    )

    if (isTRUE(de_run$skipped)) {
      layer_manifest[[length(layer_manifest) + 1]] <- tibble(
        layer = layer_name,
        contrast = contrast_label,
        group1 = group1,
        group2 = group2,
        status = "skipped",
        message = de_run$message %||% "Skipped.",
        sample_n = length(pair_sample_order),
        feature_n = nrow(matrix_info$matrix),
        significant_n = 0,
        up_n = 0,
        down_n = 0
      )
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), de_run$message %||% "Skipped."),
        pair_section,
        paste0("skipped_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
      next
    }

    result_tbl <- de_run$result %>%
      mutate(FeatureInternal = as.character(Feature)) %>%
      select(-any_of(c("FeatureID", "Feature"))) %>%
      left_join(feature_lookup, by = "FeatureInternal") %>%
      select(Feature, everything(), -FeatureInternal) %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        group1 = group1,
        group2 = group2,
        contrast = contrast_label
      )

    if (!is.null(feature_meta) && nrow(feature_meta) > 0) {
      result_tbl <- result_tbl %>% left_join(feature_meta, by = c("Feature" = "Feature"))
    }

    summary_tbl <- result_tbl %>%
      count(direction, name = "feature_count") %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        group1 = group1,
        group2 = group2,
        contrast = contrast_label
      )

    sig_tbl <- result_tbl %>% filter(direction != "ns")
    up_tbl <- result_tbl %>% filter(direction == "up") %>% arrange(SignificanceMetric, PValue, desc(logFC))
    down_tbl <- result_tbl %>% filter(direction == "down") %>% arrange(SignificanceMetric, PValue, logFC)
    robustness_tbl <- build_de_robustness_table(
      result_tbl = result_tbl,
      count_matrix = matrix_info$matrix,
      sample_order = pair_sample_order,
      sample_groups = pair_sample_groups,
      group1 = group1,
      group2 = group2
    ) %>%
      mutate(
        layer = layer_name,
        layer_title = layer_title,
        contrast = contrast_label
      )
    robustness_summary_tbl <- summarise_de_robustness_flags(robustness_tbl)

    layer_results[[length(layer_results) + 1]] <- result_tbl
    layer_summary[[length(layer_summary) + 1]] <- summary_tbl
    layer_robustness[[length(layer_robustness) + 1]] <- robustness_tbl
    layer_robustness_summary[[length(layer_robustness_summary) + 1]] <- robustness_summary_tbl
    layer_manifest[[length(layer_manifest) + 1]] <- tibble(
      layer = layer_name,
      contrast = contrast_label,
      group1 = group1,
      group2 = group2,
      status = "completed",
      message = "",
      sample_n = length(pair_sample_order),
      feature_n = nrow(matrix_info$matrix),
      significant_n = nrow(sig_tbl),
      up_n = sum(result_tbl$direction == "up", na.rm = TRUE),
      down_n = sum(result_tbl$direction == "down", na.rm = TRUE)
    )

    write_csv_safe(result_tbl, file.path(pair_section, paste0("de_results_", safe_dir_name(contrast_label), ".csv")))
    write_csv_safe(summary_tbl, file.path(pair_section, paste0("de_summary_", safe_dir_name(contrast_label), ".csv")))
    if (nrow(up_tbl) > 0) {
      write_csv_safe(up_tbl, file.path(pair_section, paste0("de_up_", safe_dir_name(contrast_label), ".csv")))
    }
    if (nrow(down_tbl) > 0) {
      write_csv_safe(down_tbl, file.path(pair_section, paste0("de_down_", safe_dir_name(contrast_label), ".csv")))
    }
    if (nrow(robustness_tbl) > 0) {
      write_csv_safe(robustness_tbl, file.path(pair_section, paste0("de_robustness_", safe_dir_name(contrast_label), ".csv")))
      write_csv_safe(robustness_summary_tbl, file.path(pair_section, paste0("de_robustness_summary_", safe_dir_name(contrast_label), ".csv")))
      save_plot_nested(
        build_de_robustness_plot(
          robustness_tbl,
          paste0(layer_title, " robustness | ", contrast_label)
        ),
        pair_section,
        paste0("robustness_", safe_dir_name(contrast_label), ".pdf"),
        width = 8.8,
        height = 6.4,
        view = REPORT_VIEW
      )
    }

    volcano_plot <- build_de_volcano_plot(
      result_tbl,
      paste0(layer_title, " | ", contrast_label),
      label_n = 15
    )
    save_plot_nested(
      volcano_plot,
      pair_section,
      paste0("volcano_", safe_dir_name(contrast_label), ".pdf"),
      width = 8.6,
      height = 6.2,
      view = REPORT_VIEW
    )

    ma_plot <- build_de_ma_plot(
      result_tbl,
      paste0(layer_title, " MA | ", contrast_label)
    )
    save_plot_nested(
      ma_plot,
      pair_section,
      paste0("ma_", safe_dir_name(contrast_label), ".pdf"),
      width = 8.6,
      height = 6.2,
      view = REPORT_VIEW
    )

    top_heat_features <- result_tbl %>%
      arrange(SignificanceMetric, PValue, desc(abs(logFC))) %>%
      slice_head(n = 30) %>%
      pull(Feature) %>%
      unique()
    top_heat_features <- top_heat_features[top_heat_features %in% rownames(matrix_info$matrix)]

    if (length(top_heat_features) >= 2 && pkg_has("pheatmap")) {
      heat_dge <- edgeR::DGEList(counts = matrix_info$matrix[, pair_sample_order, drop = FALSE])
      heat_dge <- edgeR::calcNormFactors(heat_dge)
      heat_mat <- edgeR::cpm(heat_dge, log = TRUE, prior.count = 1)[top_heat_features, , drop = FALSE]
      heat_mat <- safe_scale_rows(heat_mat)
      heat_mat <- pmax(pmin(heat_mat, 2.5), -2.5)
      heat_anno <- data.frame(Group = factor(pair_sample_groups, levels = c(group1, group2)))
      rownames(heat_anno) <- pair_sample_order

      heat_plot_path <- file.path(OUTPUT_FIG_DIR, pair_section, paste0("top_heatmap_", safe_dir_name(contrast_label), ".pdf"))
      save_pheatmap_file(
        heat_mat,
        annotation_col = heat_anno,
        annotation_colors = list(Group = group_palette[c(group1, group2)]),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        color = grDevices::colorRampPalette(mii_heatmap_gradient)(101),
        breaks = seq(-2.5, 2.5, length.out = 102),
        legend_breaks = seq(-2.5, 2.5, by = 1.25),
        fontsize_row = 7,
        fontsize_col = 9,
        main = paste0(layer_title, " top features | ", contrast_label),
        filename = heat_plot_path,
        width = 9,
        height = max(5.5, 0.25 * nrow(heat_mat) + 2.5)
      )
    } else {
      save_plot_nested(
        plot_note(paste0(layer_title, " | ", contrast_label), "No heatmap available."),
        pair_section,
        paste0("heatmap_note_", safe_dir_name(contrast_label), ".pdf"),
        width = 8,
        height = 5,
        view = REPORT_VIEW
      )
    }
  }

  combined_results <- bind_rows(layer_results)
  combined_summary <- bind_rows(layer_summary)
  combined_manifest <- bind_rows(layer_manifest)
  combined_robustness <- bind_rows(layer_robustness)
  combined_robustness_summary <- bind_rows(layer_robustness_summary)

  write_csv_safe(combined_results, paste0(layer_name, "_de_results.csv"))
  write_csv_safe(combined_summary, paste0(layer_name, "_de_summary.csv"))
  write_csv_safe(combined_manifest, paste0(layer_name, "_de_manifest.csv"))
  if (nrow(combined_robustness) > 0) {
    write_csv_safe(combined_robustness, paste0(layer_name, "_de_robustness.csv"))
    write_csv_safe(combined_robustness_summary, paste0(layer_name, "_de_robustness_summary.csv"))
  }

  list(
    results = combined_results,
    summary = combined_summary,
    manifest = combined_manifest,
    robustness = combined_robustness,
    robustness_summary = combined_robustness_summary
  )
}

normalize_mirna_name <- function(x) {
  x <- as.character(x)
  x <- tolower(trimws(x))
  x <- gsub("^[a-z]{2,3}[-_]+", "", x)
  x <- gsub("^mirna", "mir", x)
  x <- gsub("[-_. ]+(3p|5p)$", "", x)
  x <- gsub("[-_. ]+", "", x)
  x <- sub("^mirlet", "let", x)
  toupper(x)
}

extract_mirna_tokens_from_gs_name <- function(gs_name) {
  tokens <- unlist(strsplit(as.character(gs_name), "_", fixed = TRUE))
  tokens <- tokens[grepl("^(MIR|LET)", toupper(tokens))]
  unique(normalize_mirna_name(tokens))
}

empty_mirna_target_sets <- function(query_tokens = character(), source = "not_available", status = "skipped", message = "") {
  list(
    mapping = tibble(
      mirna_token = character(),
      target_gene = character(),
      source_database = character(),
      evidence_source = character(),
      gs_name = character()
    ),
    summary = tibble(
      mirna_token = character(),
      target_gene_count = integer(),
      target_set_count = integer(),
      total_links = integer(),
      target_genes = character(),
      target_sets = character()
    ),
    target_genes = tibble(target_gene = character()),
    query_tokens = query_tokens,
    source_summary = tibble(
      target_source = source,
      status = status,
      message = message,
      query_miRNA_n = length(query_tokens),
      target_gene_n = 0L,
      target_link_n = 0L
    )
  )
}

summarise_mirna_target_mapping <- function(mapping, query_tokens, source, status = "completed", message = "") {
  if (nrow(mapping) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = source, status = "empty", message = message))
  }

  mapping <- mapping %>%
    mutate(
      mirna_token = normalize_mirna_name(mirna_token),
      target_gene = trimws(as.character(target_gene)),
      gs_name = if ("gs_name" %in% names(.)) as.character(gs_name) else paste(source_database, mirna_token, sep = ":"),
      source_database = if ("source_database" %in% names(.)) as.character(source_database) else source,
      evidence_source = if ("evidence_source" %in% names(.)) as.character(evidence_source) else source
    ) %>%
    filter(mirna_token %in% query_tokens, !is.na(target_gene), nzchar(target_gene)) %>%
    group_by(mirna_token, target_gene, source_database, evidence_source) %>%
    summarise(
      target_set_n = n_distinct(gs_name),
      gs_name = paste(head(sort(unique(gs_name)), 20), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(mirna_token, target_gene, source_database)

  if (nrow(mapping) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = source, status = "empty", message = message))
  }

  summary_tbl <- mapping %>%
    group_by(mirna_token) %>%
    summarise(
      target_gene_count = n_distinct(target_gene),
      target_set_count = sum(target_set_n, na.rm = TRUE),
      total_links = n(),
      target_genes = paste(sort(unique(target_gene)), collapse = ";"),
      target_sets = paste(sort(unique(gs_name)), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(target_gene_count), mirna_token)

  target_genes <- mapping %>%
    distinct(target_gene) %>%
    arrange(target_gene)

  list(
    mapping = mapping,
    summary = summary_tbl,
    target_genes = target_genes,
    query_tokens = query_tokens,
    source_summary = tibble(
      target_source = source,
      status = status,
      message = message,
      query_miRNA_n = length(query_tokens),
      target_gene_n = nrow(target_genes),
      target_link_n = nrow(mapping)
    )
  )
}

format_multimir_query <- function(mirna_names) {
  x <- trimws(as.character(mirna_names))
  x <- x[!is.na(x) & nzchar(x)]
  prefix <- paste0(MIRNA_PREFIX, "-")
  x <- str_replace(x, regex("^[a-z]{2,3}-mir", ignore_case = TRUE), paste0(prefix, "miR"))
  x <- str_replace(x, regex("^mir", ignore_case = TRUE), paste0(prefix, "miR"))
  x <- str_replace(x, regex("^[a-z]{2,3}-let", ignore_case = TRUE), paste0(prefix, "let"))
  x <- str_replace(x, regex("^let", ignore_case = TRUE), paste0(prefix, "let"))
  unique(x)
}

collect_multimir_target_sets <- function(mirna_names, query_tokens, target_source) {
  if (!pkg_has("multiMiR")) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "skipped", message = "multiMiR package is not available."))
  }

  target_tables <- if (identical(target_source, "predicted")) c("validated", "predicted") else "validated"
  query_mirnas <- format_multimir_query(mirna_names)
  if (length(query_mirnas) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "empty", message = "No miRNA names were available for multiMiR."))
  }

  mapping <- purrr::map_dfr(target_tables, function(target_table) {
    mm_obj <- tryCatch(
      multiMiR::get_multimir(org = MULTIMIR_ORG, mirna = query_mirnas, table = target_table, legacy.out = FALSE),
      error = function(e) e
    )
    if (inherits(mm_obj, "error")) {
      return(tibble(.message = conditionMessage(mm_obj), .table = target_table))
    }

    mm_df <- tryCatch(
      as_tibble(mm_obj@data),
      error = function(e) tryCatch(as_tibble(as.data.frame(mm_obj)), error = function(e2) tibble())
    )
    if (nrow(mm_df) == 0) return(tibble())

    mirna_col <- intersect(c("mature_mirna_id", "mirna", "mirna_id", "mature_mirna_acc"), names(mm_df))[1]
    target_col <- intersect(c("target_symbol", "target_gene", "symbol", "target"), names(mm_df))[1]
    db_col <- intersect(c("database", "source", "support_type"), names(mm_df))[1]
    if (is.na(mirna_col) || is.na(target_col)) return(tibble())

    tibble(
      mirna_token = normalize_mirna_name(mm_df[[mirna_col]]),
      target_gene = as.character(mm_df[[target_col]]),
      source_database = if (!is.na(db_col)) as.character(mm_df[[db_col]]) else "multiMiR",
      evidence_source = target_table,
      gs_name = paste0("multiMiR:", target_table, ":", if (!is.na(db_col)) as.character(mm_df[[db_col]]) else "database")
    )
  })

  if (".message" %in% names(mapping) && all(!all(c("mirna_token", "target_gene") %in% names(mapping)))) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "skipped", message = paste(unique(mapping$.message), collapse = "; ")))
  }
  mapping <- mapping %>% select(any_of(c("mirna_token", "target_gene", "source_database", "evidence_source", "gs_name")))

  summarise_mirna_target_mapping(mapping, query_tokens, source = paste0("multiMiR_", target_source), status = "completed")
}

collect_msigdb_mir_target_sets <- function(mirna_names, query_tokens, species = MIRNA_SPECIES, subcollections = c("MIR:MIRDB", "MIR:MIR_LEGACY"), fallback_message = "") {
  if (!pkg_has("msigdbr")) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "skipped", message = "msigdbr package is not available."))
  }

  target_db <- tryCatch(
    purrr::map_dfr(subcollections, function(subcollection) {
      msigdbr::msigdbr(species = species, collection = "C3", subcollection = subcollection) %>%
        transmute(
          gs_collection,
          gs_subcollection,
          gs_name,
          gs_description,
          target_gene = gene_symbol,
          mirna_token = list(extract_mirna_tokens_from_gs_name(gs_name)),
          source_database = "MSigDB_C3_MIR",
          evidence_source = "msigdb_fallback"
        ) %>%
        filter(purrr::map_lgl(mirna_token, ~any(.x %in% query_tokens)))
    }),
    error = function(e) tibble(.message = conditionMessage(e))
  )

  if (".message" %in% names(target_db)) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "skipped", message = unique(target_db$.message)[[1]]))
  }
  if (nrow(target_db) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = "msigdb_fallback", status = "empty", message = "No MSigDB C3 MIR records were returned."))
  }

  mapping <- target_db %>%
    tidyr::unnest_longer(mirna_token, values_to = "mirna_token") %>%
    filter(!is.na(mirna_token), mirna_token != "", mirna_token %in% query_tokens) %>%
    transmute(
      mirna_token,
      target_gene,
      source_database,
      evidence_source,
      gs_name
    )

  status <- if (nzchar(fallback_message)) "fallback" else "completed"
  summarise_mirna_target_mapping(mapping, query_tokens, source = "msigdb_fallback", status = status, message = fallback_message)
}

collect_mirna_target_sets <- function(mirna_names, species = MIRNA_SPECIES, subcollections = c("MIR:MIRDB", "MIR:MIR_LEGACY"), target_source = MIRNA_TARGET_SOURCE) {
  query_tokens <- unique(normalize_mirna_name(mirna_names))
  query_tokens <- query_tokens[nzchar(query_tokens)]
  if (length(query_tokens) == 0) {
    return(empty_mirna_target_sets(query_tokens, source = target_source, status = "empty", message = "No query miRNA tokens were available."))
  }

  if (!identical(target_source, "msigdb_fallback")) {
    multimir_sets <- collect_multimir_target_sets(mirna_names, query_tokens, target_source = target_source)
    if (nrow(multimir_sets$mapping) > 0) {
      return(multimir_sets)
    }
    multimir_sets$source_summary <- multimir_sets$source_summary %>%
      mutate(
        status = "skipped",
        message = paste0(
          "Primary source unavailable or empty: ",
          message %||% "multiMiR returned no targets.",
          " Exploratory MSigDB C3 MIR fallback is disabled unless MIRNA_TARGET_SOURCE=msigdb_fallback."
        )
      )
    return(multimir_sets)
  }

  collect_msigdb_mir_target_sets(mirna_names, query_tokens, species = species, subcollections = subcollections)
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
      "RESULT_DIR", "METADATA_DIR", "QC_DIR", "SPORTS_ANNOTATION_DIR", "SAMPLE_ROOT",
      "ENDOSIRNA_DIR", "ENDOSIRNA_SAMPLE_ROOT", "ENDOSIRNA_SUMMARY_DIR",
      "LOG_DIR", "SAMPLE_LOG_DIR", "CUTADAPT_LOG_DIR", "ENDOSIRNA_LOG_DIR", "WORK_DIR", "TMP_DIR",
      "DB_ROOT", "GENOME_BASE", "MIRNA_BASE", "RRNA_BASE", "TRNA_BASE",
      "TRNA_MITO_BASE", "PIRNA_BASE", "ENSEMBL_BASE", "RFAM_BASE",
      "ADAPTER3", "PRETRIM_FASTQ", "TRIM_ADAPTER", "CUTADAPT_NEXTSEQ_TRIM", "CUTADAPT_QUALITY_CUTOFF", "MISMATCH", "RUN_FASTQC",
      "RUN_ENDOSIRNA", "ENDOSIRNA_REFS", "ENDOSIRNA_MISMATCH", "ENDOSIRNA_MAX_ALIGNMENTS",
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
  utils::write.csv(
    ctx$sample_comp_view %>%
      transmute(
        Sample = as.character(Sample), Group = as.character(Group),
        MainCategory = as.character(MainCategory), Reads, CleanReads, Percent
      ),
    file.path(output_dir, "exclusive_class_composition.csv"), row.names = FALSE
  )
  utils::write.csv(ctx$composition_integrity_check, file.path(output_dir, "composition_integrity_check.csv"), row.names = FALSE)
  utils::write.csv(ctx$length_composition_integrity_check, file.path(output_dir, "length_composition_integrity_check.csv"), row.names = FALSE)
  utils::write.csv(ctx$tsrna_cleavage_reconstruction_check, file.path(output_dir, "tsrna_cleavage_reconstruction_check.csv"), row.names = FALSE)
  write_tsv_safe(if (is.null(ctx$raw_fastq_audit_by_sample)) empty_raw_fastq_audit_tables()$raw_fastq_audit_by_sample else ctx$raw_fastq_audit_by_sample, file.path(output_dir, "raw_fastq_audit_by_sample.tsv"))
  write_tsv_safe(if (is.null(ctx$raw_fastq_audit_by_group)) empty_raw_fastq_audit_tables()$raw_fastq_audit_by_group else ctx$raw_fastq_audit_by_group, file.path(output_dir, "raw_fastq_audit_by_group.tsv"))
  invisible(output_dir)
}

collapse_unique_text <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else paste(sort(x), collapse = ";")
}

make_company_tsrna_feature <- function(sub_class, cleavage) {
  base_name <- as.character(sub_class) %>%
    str_remove("^mature-") %>%
    str_remove("^pre-") %>%
    str_replace("_(5_end|3_end|CCA_end)$", "")
  suffix <- case_when(
    cleavage == "5'end" ~ "_5_end",
    cleavage == "3'end" ~ "_3_end",
    cleavage == "CCA-end" ~ "_3_CCA",
    TRUE ~ "_inner"
  )
  paste0(base_name, suffix)
}

empty_expression_matrix <- function(output_col, sample_order) {
  out <- tibble(.feature = character())
  names(out)[1] <- output_col
  for (sample_name in sample_order) out[[sample_name]] <- numeric()
  out
}

build_expression_wide_matrix <- function(df, feature_col, value_col, output_col, sample_order) {
  if (nrow(df) == 0) return(empty_expression_matrix(output_col, sample_order))

  long_df <- df %>%
    transmute(
      .feature = as.character(.data[[feature_col]]),
      Sample = as.character(Sample),
      .value = clean_numeric(.data[[value_col]])
    ) %>%
    filter(!is.na(.feature), .feature != "", Sample %in% sample_order) %>%
    group_by(.feature, Sample) %>%
    summarise(.value = sum(.value, na.rm = TRUE), .groups = "drop")

  if (nrow(long_df) == 0) return(empty_expression_matrix(output_col, sample_order))

  wide_df <- long_df %>%
    mutate(Sample = factor(Sample, levels = sample_order)) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = .value, values_fill = 0)

  for (sample_name in sample_order) {
    if (!sample_name %in% names(wide_df)) wide_df[[sample_name]] <- 0
  }

  wide_df <- wide_df %>%
    mutate(.total = rowSums(across(all_of(sample_order)), na.rm = TRUE)) %>%
    arrange(desc(.total), .feature) %>%
    select(.feature, all_of(sample_order))
  names(wide_df)[1] <- output_col
  wide_df
}

prepare_annotation_expression_values <- function(annotation_df, clean_reads_tbl) {
  annotation_df %>%
    group_by(Feature, Sample) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    left_join(clean_reads_tbl %>% mutate(Sample = as.character(Sample)), by = "Sample") %>%
    mutate(RPM = if_else(CleanReads > 0, 1e6 * Reads / CleanReads, 0))
}

write_annotation_expression_set <- function(annotation_df, set_name, output_dir, sample_order, clean_reads_tbl) {
  values_df <- prepare_annotation_expression_values(annotation_df, clean_reads_tbl)
  count_mat <- build_expression_wide_matrix(values_df, "Feature", "Reads", "Sub_Class", sample_order)
  rpm_mat <- build_expression_wide_matrix(values_df, "Feature", "RPM", "Sub_Class", sample_order)
  write_tsv_safe(count_mat, file.path(output_dir, paste0("count_", set_name, ".txt")))
  write_tsv_safe(rpm_mat, file.path(output_dir, paste0("rpm_", set_name, ".txt")))
}

write_sequence_expression_set <- function(sequence_df, set_name, output_dir, seq_dir, sample_order, duplicate_root = FALSE) {
  if (pkg_has("data.table")) {
    sequence_dt <- data.table::as.data.table(sequence_df)
    if (nrow(sequence_dt) == 0) {
      count_mat <- empty_expression_matrix("Sequence", sample_order)
      rpm_mat <- empty_expression_matrix("Sequence", sample_order)
      type_tbl <- tibble(Sequence = character(), Annotation = character(), FeatureCategory = character())
    } else {
      sequence_dt <- sequence_dt[!is.na(Sequence) & Sequence != ""]
      sequence_dt[, Sample := as.character(Sample)]
      sequence_dt <- sequence_dt[Sample %in% sample_order]
      if (nrow(sequence_dt) == 0) {
        count_mat <- empty_expression_matrix("Sequence", sample_order)
        rpm_mat <- empty_expression_matrix("Sequence", sample_order)
        type_tbl <- tibble(Sequence = character(), Annotation = character(), FeatureCategory = character())
      } else {
        sequence_dt[, Reads := clean_numeric(Reads)]
        sequence_dt[, RPM := clean_numeric(RPM)]

        values_dt <- sequence_dt[, .(
          Reads = sum(Reads, na.rm = TRUE),
          RPM = sum(RPM, na.rm = TRUE)
        ), by = .(Sequence, Sample)]

        count_dt <- data.table::dcast(values_dt, Sequence ~ Sample, value.var = "Reads", fill = 0)
        rpm_dt <- data.table::dcast(values_dt, Sequence ~ Sample, value.var = "RPM", fill = 0)
        for (sample_name in sample_order) {
          if (!sample_name %in% names(count_dt)) count_dt[, (sample_name) := 0]
          if (!sample_name %in% names(rpm_dt)) rpm_dt[, (sample_name) := 0]
        }

        if (nrow(count_dt) > 0) {
          count_dt[, .total := rowSums(.SD, na.rm = TRUE), .SDcols = sample_order]
          data.table::setorderv(count_dt, c(".total", "Sequence"), c(-1L, 1L))
          sequence_order <- count_dt$Sequence
          rpm_dt <- rpm_dt[match(sequence_order, Sequence)]
          count_mat <- as_tibble(count_dt[, c("Sequence", sample_order), with = FALSE])
          rpm_mat <- as_tibble(rpm_dt[, c("Sequence", sample_order), with = FALSE])
        } else {
          sequence_order <- character()
          count_mat <- empty_expression_matrix("Sequence", sample_order)
          rpm_mat <- empty_expression_matrix("Sequence", sample_order)
        }

        collapse_dt_text <- function(x) {
          x <- unique(as.character(x))
          x <- x[!is.na(x) & nzchar(x)]
          if (length(x) == 0) "" else paste(sort(x), collapse = ";")
        }
        type_unique <- unique(sequence_dt[, .(Sequence, Annotation=OriginalAnnotation, FeatureCategory)])
        multiple <- type_unique$Sequence[duplicated(type_unique$Sequence)]
        type_dt <- type_unique[!Sequence %in% multiple]
        if (length(multiple)) type_dt <- data.table::rbindlist(list(type_dt,
          type_unique[Sequence %in% multiple, .(
            Annotation=collapse_dt_text(Annotation),FeatureCategory=collapse_dt_text(FeatureCategory)
          ),by=Sequence]))
        type_dt[, .order := match(Sequence, sequence_order)]
        data.table::setorderv(type_dt, c(".order", "Sequence"), c(1L, 1L), na.last = TRUE)
        type_dt[, .order := NULL]
        type_tbl <- as_tibble(type_dt)
      }
    }
  } else {
    values_df <- sequence_df %>%
      group_by(Sequence, Sample) %>%
      summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop")

    count_mat <- build_expression_wide_matrix(values_df, "Sequence", "Reads", "Sequence", sample_order)
    rpm_mat <- build_expression_wide_matrix(values_df, "Sequence", "RPM", "Sequence", sample_order)
    sequence_order <- count_mat$Sequence

    type_tbl <- sequence_df %>%
      group_by(Sequence) %>%
      summarise(
        Annotation = collapse_unique_text(OriginalAnnotation),
        FeatureCategory = collapse_unique_text(FeatureCategory),
        .groups = "drop"
      ) %>%
      mutate(.order = match(Sequence, sequence_order)) %>%
      arrange(.order, Sequence) %>%
      select(-.order)
  }

  write_tsv_safe(count_mat, file.path(seq_dir, paste0("seq_count_", set_name, ".txt")))
  write_tsv_safe(rpm_mat, file.path(seq_dir, paste0("seq_rpm_", set_name, ".txt")))
  write_tsv_safe(type_tbl, file.path(seq_dir, paste0("seq_type_", set_name, ".txt")))

  if (isTRUE(duplicate_root)) {
    write_tsv_safe(count_mat, file.path(output_dir, paste0("seq_count_", set_name, ".txt")))
    write_tsv_safe(rpm_mat, file.path(output_dir, paste0("seq_rpm_", set_name, ".txt")))
    write_tsv_safe(type_tbl, file.path(output_dir, paste0("seq_type_", set_name, ".txt")))
  }
}


compute_sports_snc_for_sample <- function(sample_id, length_df, sample_info_tbl) {
  sample_meta <- sample_info_tbl %>% filter(as.character(Sample) == sample_id) %>% slice(1)
  sample_df <- length_df %>% filter(as.character(Sample) == sample_id)
  if (nrow(sample_df) == 0 || nrow(sample_meta) == 0) {
    return(list(length = tibble(), rsrna_status = tibble(), tsrna_cleavage = tibble()))
  }

  len <- sample_df %>% filter(Class == "Clean_Reads") %>% pull(Length) %>% unique() %>% sort()
  if (length(len) == 0) len <- seq(min(sample_df$Length, na.rm = TRUE), max(sample_df$Length, na.rm = TRUE))
  clean_total <- (sample_df %>% filter(Class == "Clean_Reads") %>% summarise(Total = sum(Reads, na.rm = TRUE), .groups = "drop") %>% pull(Total)) %||% 0
  class_text <- sample_df$Class
  is_tsrna_total <- str_detect(class_text, regex("tRNA_Match_Genome|tRNA_Unmatch_Genome", ignore_case = TRUE))
  is_tsrna_5 <- str_detect(class_text, regex("tRNA_5_end_Match_Genome|tRNA_5_end_Unmatch_Genome", ignore_case = TRUE))
  is_tsrna_3 <- str_detect(class_text, regex("tRNA_3_end_Match_Genome|tRNA_3_end_Unmatch_Genome", ignore_case = TRUE))
  is_tsrna_cca <- str_detect(class_text, regex("tRNA_CCA_end_Match_Genome|tRNA_CCA_end_Unmatch_Genome", ignore_case = TRUE))
  is_rsrna <- str_detect(class_text, regex("-rRNA|_rRNA|-YRNA", ignore_case = TRUE)) &
    !str_detect(class_text, regex("antisense|S-rRNA|other-rRNA|YRNA", ignore_case = TRUE))
  is_rsrna_match <- is_rsrna & str_detect(class_text, regex("_Match_Genome", ignore_case = TRUE))
  is_rsrna_unmatch <- is_rsrna & str_detect(class_text, regex("_Unmatch_Genome", ignore_case = TRUE))
  is_unanno_mg <- str_detect(class_text, regex("Unannotated_Match_Genome", ignore_case = TRUE))
  is_unanno_umg <- str_detect(class_text, regex("Unannotated_Unmatch_Genome", ignore_case = TRUE))
  exclusive_sample_df <- sample_df %>%
    filter(!is_tsrna_nested_class(Class), !is_nested_length_class(Class))

  classified <- map_dfr(
    setdiff(sports_length_categories, c("other", "unanno_MG", "unanno_UMG")),
    function(category_name) complete_length_series(exclusive_sample_df %>% filter(LengthCategory == category_name), len, category_name, clean_total)
  )
  classified_total <- classified %>%
    group_by(Length) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    right_join(tibble(Length = len), by = "Length") %>%
    mutate(Reads = replace_na(Reads, 0)) %>%
    arrange(Length)
  umg_series <- complete_length_series(sample_df[is_unanno_umg, ], len, "unanno_UMG", clean_total)
  mg_series <- complete_length_series(sample_df[is_unanno_mg, ], len, "unanno_MG", clean_total)
  clean_series <- complete_length_series(sample_df %>% filter(Class == "Clean_Reads"), len, "Clean_Reads", clean_total)
  other_reads <- pmax(clean_series$Reads - classified_total$Reads - umg_series$Reads - mg_series$Reads, 0)
  other_series <- tibble(
    Category = "other",
    Length = len,
    Reads = other_reads,
    RPM = if (clean_total > 0) 1e6 * other_reads / clean_total else rep(0, length(other_reads))
  )
  length_out <- bind_rows(classified, other_series, mg_series, umg_series) %>%
    mutate(
      Sample = sample_id,
      Group = as.character(sample_meta$Group[[1]]),
      Category = factor(Category, levels = sports_length_categories)
    ) %>%
    select(Sample, Group, Category, Length, Reads, RPM)

  rsrna_status <- bind_rows(
    complete_length_series(sample_df[is_rsrna_match, ], len, "Match_Genome", clean_total),
    complete_length_series(sample_df[is_rsrna_unmatch, ], len, "Unmatch_Genome", clean_total)
  ) %>%
    group_by(Category) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Sample = sample_id,
      Group = as.character(sample_meta$Group[[1]]),
      GenomeStatus = factor(Category, levels = c("Match_Genome", "Unmatch_Genome")),
      Percent = group_percent(Reads)
    ) %>%
    select(Sample, Group, GenomeStatus, Reads, RPM, Percent)

  tsrna_5 <- complete_length_series(sample_df[is_tsrna_5, ], len, "5'end", clean_total)
  tsrna_3 <- complete_length_series(sample_df[is_tsrna_3, ], len, "3'end", clean_total)
  tsrna_cca <- complete_length_series(sample_df[is_tsrna_cca, ], len, "CCA-end", clean_total)
  tsrna_total <- complete_length_series(sample_df[is_tsrna_total, ], len, "tsRNA", clean_total)
  tsrna_end_sum <- tsrna_5$Reads + tsrna_3$Reads + tsrna_cca$Reads
  if (any(tsrna_end_sum > tsrna_total$Reads + 1, na.rm = TRUE)) {
    stop("tRNA end subsets exceed base totals by more than 1 read in sample ", sample_id)
  }
  tsrna_other_reads <- pmax(tsrna_total$Reads - tsrna_end_sum, 0)
  tsrna_other <- tibble(
    Category = "Internal-other",
    Length = len,
    Reads = tsrna_other_reads,
    RPM = if (clean_total > 0) 1e6 * tsrna_other_reads / clean_total else rep(0, length(tsrna_other_reads))
  )
  tsrna_cleavage <- bind_rows(tsrna_5, tsrna_3, tsrna_cca, tsrna_other) %>%
    group_by(Category) %>%
    summarise(Reads = sum(Reads, na.rm = TRUE), RPM = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Sample = sample_id,
      Group = as.character(sample_meta$Group[[1]]),
      Cleavage = factor(Category, levels = c("5'end", "3'end", "CCA-end", "Internal-other")),
      Percent = group_percent(Reads)
    ) %>%
    select(Sample, Group, Cleavage, Reads, RPM, Percent)

  list(length = length_out, rsrna_status = rsrna_status, tsrna_cleavage = tsrna_cleavage)
}
