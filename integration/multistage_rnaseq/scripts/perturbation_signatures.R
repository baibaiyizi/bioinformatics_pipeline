#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_packages <- c("data.table", "fgsea", "ggplot2", "jsonlite", "scales", "digest")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  library(data.table)
  library(fgsea)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

parse_cli <- function(args) {
  allowed <- c("ranked", "vst", "samples", "obox-signatures", "pabpn-effects", "output", "config")
  result <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% allowed || i == length(args)) stop("Unknown or incomplete argument: ", args[[i]])
    result[[gsub("-", "_", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- gsub("-", "_", setdiff(allowed, "config"))
  missing <- required[!vapply(required, function(x) nzchar(result[[x]] %||% ""), logical(1))]
  if (length(missing)) stop("Missing arguments: ", paste(paste0("--", gsub("_", "-", missing)), collapse = ", "))
  result
}

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input does not exist: ", path)
  data.table::fread(path, sep = "\t", na.strings = c("", "NA", "NaN"), data.table = TRUE, check.names = FALSE)
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(data), path, sep = "\t", quote = FALSE, na = "NA")
}

star_label <- function(value) {
  ifelse(!is.finite(value), "NT", ifelse(value < 0.001, "***", ifelse(value < 0.01, "**", ifelse(value < 0.05, "*", "ns"))))
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
md5_file <- function(path) digest::digest(file = path, algo = "md5", serialize = FALSE)

publication_theme <- function() {
  ggplot2::theme_classic(base_family = "Liberation Sans", base_size = 7) +
    ggplot2::theme(
      text = element_text(family = "Liberation Sans", color = "#202124"),
      axis.title = element_text(size = 8), axis.text = element_text(size = 7, color = "#202124"),
      legend.title = element_text(size = 7.5), legend.text = element_text(size = 7),
      plot.title = element_text(size = 9, face = "plain", hjust = 0),
      plot.subtitle = element_text(size = 7, color = "#4D5156", hjust = 0),
      axis.line = element_line(linewidth = 0.35, color = "#30343B"),
      axis.ticks = element_line(linewidth = 0.35, color = "#30343B"),
      legend.key = element_blank(), legend.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(4, 5, 4, 4, unit = "mm")
    )
}

save_figure <- function(plot, stem, width_mm, height_mm) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  ggplot2::ggsave(pdf_path, plot = plot, width = width_mm, height = height_mm, units = "mm",
                  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(png_path, plot = plot, width = width_mm, height = height_mm, units = "mm",
                  dpi = 300, bg = "white", limitsize = FALSE)
  list(pdf = pdf_path, png = png_path, width_mm = width_mm, height_mm = height_mm)
}

relative_to <- function(path, root) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- paste0(normalizePath(root, winslash = "/", mustWork = FALSE), "/")
  if (startsWith(p, r)) substring(p, nchar(r) + 1L) else p
}

analysis_script_path <- function() {
  configured <- Sys.getenv("PERTURBATION_SCRIPT_PATH", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, winslash = "/", mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) stop("Cannot determine the executing perturbation-signature script")
  normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
}

make_figure_manifest_row <- function(asset_id, figure_id, media, output_dir, source_table,
                                     question, effect_metric, comparison, test, caption,
                                     experimental_unit, n_definition, layout_rationale, script_path) {
  data.table(
    asset_id = asset_id, asset_class = "standalone", figure_id = figure_id,
    panel_label = "", panel_ids = "", coordinates_pdf = 1L, coordinates_png = 1L,
    pdf = relative_to(media$pdf, output_dir), png = relative_to(media$png, output_dir),
    pdf_md5 = md5_file(media$pdf), png_md5 = md5_file(media$png),
    width_mm = media$width_mm, height_mm = media$height_mm,
    layout_rationale = layout_rationale, minimum_font_pt = 7,
    rng_unchanged = "true", encoding = "Direction is encoded by sign and text in addition to color; significance is encoded by stars.",
    experimental_unit = experimental_unit, n_definition = n_definition,
    effect_metric = effect_metric, uncertainty = "No uncertainty interval is inferred beyond the supplied model or descriptive summary.",
    test = test, comparison = comparison,
    source_tables = relative_to(source_table, file.path(output_dir, "tables")),
    source_md5 = md5_file(source_table), data_file = relative_to(source_table, output_dir),
    data_md5 = md5_file(source_table), script = normalizePath(script_path, winslash = "/"),
    code_sha256 = sha256_file(script_path), caption = caption,
    scientific_question = question, source_level = "Project model plus public perturbation reference",
    contract_version = "sci-source-data-composition-2026-10-12-v6"
  )
}

collapse_ranked_symbols <- function(stage_data) {
  data <- copy(stage_data)
  data[, gene_symbol := trimws(as.character(gene_symbol))]
  eligible <- data[!is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(stat)]
  eligible[, abs_stat_for_mapping := abs(stat)]
  setorder(eligible, gene_symbol, -abs_stat_for_mapping, gene_id)
  duplicate_symbols <- eligible[duplicated(gene_symbol), uniqueN(gene_symbol)]
  retained <- eligible[!duplicated(gene_symbol)]
  retained[, abs_stat_for_mapping := NULL]
  setorder(retained, -stat, gene_symbol)
  list(
    retained = retained,
    audit = data.table(
      input_rows = nrow(data), missing_symbol_or_stat = nrow(data) - nrow(eligible),
      duplicate_symbols = duplicate_symbols, retained_symbols = nrow(retained),
      tied_statistic_rows = sum(duplicated(retained$stat))
    )
  )
}

run_signature_projection <- function(ranked, signatures, stages, cfg) {
  signature_meta <- unique(signatures[, .(signature_id, source_stage, direction)])
  signature_order <- signature_meta$signature_id
  all_results <- list()
  all_mapping <- list()
  all_leading <- list()
  rank_audit <- list()
  collapsed <- list()
  for (stage_index in seq_along(stages)) {
    current_stage <- stages[[stage_index]]
    stage_rank <- collapse_ranked_symbols(ranked[ranked$stage == current_stage])
    stage_rank$audit[, stage := current_stage]
    rank_audit[[current_stage]] <- stage_rank$audit
    collapsed[[current_stage]] <- stage_rank$retained
    stat_vector <- stage_rank$retained$stat
    names(stat_vector) <- stage_rank$retained$gene_symbol
    pathway_list <- list()
    result_shell <- signature_meta[, .(signature_id, source_stage, direction)]
    result_shell[, `:=`(comparison_stage = current_stage, source_size = 0L, mapped_size = 0L,
                        test_status = "not_tested", ES = NA_real_, NES = NA_real_,
                        pvalue = NA_real_, n_more_extreme = NA_integer_, leading_edge = NA_character_)]
    for (signature_index in seq_len(nrow(signature_meta))) {
      current_signature <- signature_meta$signature_id[[signature_index]]
      members <- unique(signatures[signatures$signature_id == current_signature, gene_symbol])
      mapped <- intersect(members, names(stat_vector))
      missing <- setdiff(members, names(stat_vector))
      result_shell[result_shell$signature_id == current_signature, `:=`(source_size = length(members), mapped_size = length(mapped))]
      map_rows <- data.table(
        signature_id = current_signature, source_stage = signature_meta$source_stage[[signature_index]],
        direction = signature_meta$direction[[signature_index]], comparison_stage = current_stage,
        gene_symbol = members, mapped = members %in% mapped
      )
      map_rows[, project_gene_id := stage_rank$retained$gene_id[match(gene_symbol, stage_rank$retained$gene_symbol)]]
      map_rows[, project_stat := stage_rank$retained$stat[match(gene_symbol, stage_rank$retained$gene_symbol)]]
      map_rows[, exclusion_reason := ifelse(mapped, "", "no_exact_project_gene_symbol_match")]
      all_mapping[[paste(current_stage, current_signature, sep = "::")]] <- map_rows
      if (length(mapped) < cfg$min_size) {
        result_shell[result_shell$signature_id == current_signature, test_status := "not_tested_below_min_size"]
      } else if (length(mapped) > cfg$max_size) {
        result_shell[result_shell$signature_id == current_signature, test_status := "not_tested_above_max_size"]
      } else {
        pathway_list[[current_signature]] <- mapped
        result_shell[result_shell$signature_id == current_signature, test_status := "tested"]
      }
      if (length(missing) == 0L && length(members) == 0L) stop("Empty signature: ", current_signature)
    }
    if (length(pathway_list)) {
      set.seed(as.integer(cfg$seed) + stage_index)
      fg <- suppressWarnings(fgsea::fgseaSimple(
        pathways = pathway_list, stats = stat_vector, minSize = cfg$min_size,
        maxSize = cfg$max_size, nperm = cfg$nperm_simple, nproc = 1
      ))
      fg <- as.data.table(fg)
      fg[, leading_edge := vapply(leadingEdge, paste, character(1), collapse = ";")]
      fg[, leadingEdge := NULL]
      for (current_signature in fg$pathway) {
        record <- fg[fg$pathway == current_signature][1L]
        result_shell[result_shell$signature_id == current_signature, `:=`(
          ES = record$ES, NES = record$NES, pvalue = record$pval,
          n_more_extreme = record$nMoreExtreme, leading_edge = record$leading_edge
        )]
        if (nzchar(record$leading_edge %||% "")) {
          genes <- strsplit(record$leading_edge, ";", fixed = TRUE)[[1L]]
          all_leading[[paste(current_stage, current_signature, sep = "::")]] <- data.table(
            comparison_stage = current_stage, signature_id = current_signature, gene_symbol = genes,
            leading_edge_rank = seq_along(genes), project_stat = unname(stat_vector[genes])
          )
        }
      }
    }
    all_results[[current_stage]] <- result_shell
  }
  result <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  result[test_status == "tested" & (!is.finite(pvalue) | !is.finite(NES)), test_status := "not_estimable"]
  expected <- as.integer(cfg$expected_family_size %||% nrow(result))
  if (expected < nrow(result)) stop("expected_family_size is smaller than the planned signature-by-stage family")
  result[, padj_family := stats::p.adjust(pvalue, method = "BH", n = expected)]
  result[, display_stars := ifelse(test_status == "not_tested_below_min_size", "NT",
                                   ifelse(test_status == "not_tested_above_max_size", "NT",
                                          ifelse(test_status == "not_estimable", "NE", star_label(padj_family))))]
  result[, concordance_NES := ifelse(direction == "down_in_mzKO", -NES, NES)]
  result[, matching_stage := comparison_stage == source_stage]
  result[, family_size := expected]
  list(
    results = result,
    mapping = rbindlist(all_mapping, use.names = TRUE, fill = TRUE),
    leading = rbindlist(all_leading, use.names = TRUE, fill = TRUE),
    ranking_audit = rbindlist(rank_audit, use.names = TRUE, fill = TRUE),
    collapsed = collapsed,
    signature_order = signature_order
  )
}

select_matching_leading_edges <- function(projection, ranked, top_n) {
  if (!nrow(projection$leading)) return(data.table())
  meta <- projection$results[matching_stage == TRUE & test_status == "tested" & is.finite(padj_family) & padj_family < 0.05,
                             .(comparison_stage, signature_id, source_stage, direction, NES,
                               concordance_NES, padj_family, signature_stars = display_stars)]
  values <- merge(projection$leading, meta, by = c("comparison_stage", "signature_id"), all = FALSE)
  values <- values[order(signature_id, -abs(project_stat), gene_symbol), head(.SD, top_n), by = signature_id]
  effect_cols <- ranked[, .(comparison_stage = stage, gene_id, gene_symbol, log2FC, padj)]
  values <- merge(values, effect_cols, by = c("comparison_stage", "gene_symbol"), all.x = TRUE)
  values[, gene_stars := star_label(padj)]
  values[, display_gene := paste0(gene_symbol, " (", comparison_stage, ")")]
  effect_span <- diff(range(values$log2FC, finite = TRUE))
  values[, star_x := max(log2FC, na.rm = TRUE) + max(0.5, 0.08 * effect_span)]
  values[]
}

make_expression_table <- function(selected, vst, samples) {
  if (!nrow(selected)) return(data.table())
  selected_ids <- unique(selected[!is.na(gene_id), .(gene_id, gene_symbol, comparison_stage, signature_id, direction, display_gene)])
  value_rows <- vst[gene_id %in% selected_ids$gene_id]
  sample_ids <- intersect(samples$sample_id, setdiff(names(value_rows), "gene_id"))
  if (!length(sample_ids)) stop("No sample IDs in VST matrix match sample metadata")
  long <- melt(value_rows, id.vars = "gene_id", measure.vars = sample_ids,
               variable.name = "sample_id", value.name = "vst_expression")
  long <- merge(long, samples[, .(sample_id, stage, age_group)], by = "sample_id", all.x = TRUE)
  long <- merge(long, selected_ids, by.x = c("gene_id", "stage"),
                by.y = c("gene_id", "comparison_stage"), all = FALSE)
  setnames(long, "stage", "comparison_stage")
  long[]
}

prepare_pabpn_concordance <- function(pabpn, ranked) {
  eligible <- pabpn[as.logical(pass_expression) & is.finite(effect_eps1) & is.finite(effect_eps0_1)]
  project <- ranked[, .(stage, gene_id, gene_symbol, project_log2fc = log2FC, project_padj = padj)]
  merged <- merge(eligible, project, by = c("stage", "gene_symbol"), all = FALSE)
  merged <- merged[is.finite(project_log2fc)]
  merged[, concordant_sign_eps1 := sign(effect_eps1) == sign(project_log2fc) & sign(effect_eps1) != 0]
  merged[, concordant_sign_eps0_1 := sign(effect_eps0_1) == sign(project_log2fc) & sign(effect_eps0_1) != 0]
  summary <- merged[, .(
    mapped_genes = .N,
    spearman_rho_eps1 = suppressWarnings(cor(effect_eps1, project_log2fc, method = "spearman")),
    spearman_rho_eps0_1 = suppressWarnings(cor(effect_eps0_1, project_log2fc, method = "spearman")),
    source_effect_sensitivity_rho = suppressWarnings(cor(effect_eps1, effect_eps0_1, method = "spearman")),
    concordant_fraction_eps1 = mean(concordant_sign_eps1),
    concordant_fraction_eps0_1 = mean(concordant_sign_eps0_1),
    source_sign_stability = mean(sign(effect_eps1) == sign(effect_eps0_1))
  ), by = stage]
  list(full = merged, summary = summary)
}

plot_obox_matrix <- function(data, stages, signature_meta) {
  labels <- signature_meta[, .(
    signature_id,
    signature_label = paste0(source_stage, " mzKO ", ifelse(direction == "up_in_mzKO", "up", "down"))
  )]
  plot_data <- merge(copy(data), labels, by = "signature_id", all.x = TRUE)
  plot_data[, comparison_stage := factor(comparison_stage, levels = stages)]
  plot_data[, signature_label := factor(signature_label, levels = rev(labels$signature_label))]
  ggplot(plot_data, aes(comparison_stage, signature_label, fill = concordance_NES)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = display_stars), size = 2.5, family = "Liberation Sans", color = "#202124") +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                         na.value = "#D9DDE1", name = "Concordance\nNES") +
    labs(
      title = "OBOX mzKO signatures across maternal-age contrasts",
      subtitle = "Positive concordance indicates that MO−MY effects resemble the signed mzKO response",
      x = "Embryonic stage", y = NULL
    ) + publication_theme() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")
}

plot_leading_effects <- function(data) {
  plot_data <- copy(data)
  plot_data[, display_gene := factor(display_gene, levels = rev(unique(display_gene)))]
  ggplot(plot_data, aes(log2FC, display_gene, shape = direction, color = direction)) +
    geom_vline(xintercept = 0, color = "#B5BBC2", linewidth = 0.35) +
    geom_point(size = 2.2, stroke = 0.5) +
    geom_text(aes(x = star_x, label = gene_stars), hjust = 0, size = 2.5,
              family = "Liberation Sans", color = "#202124", show.legend = FALSE) +
    scale_color_manual(values = c(up_in_mzKO = "#D55E00", down_in_mzKO = "#0072B2"),
                       labels = c(up_in_mzKO = "Up in mzKO", down_in_mzKO = "Down in mzKO"), name = NULL) +
    scale_shape_manual(values = c(up_in_mzKO = 17, down_in_mzKO = 16),
                       labels = c(up_in_mzKO = "Up in mzKO", down_in_mzKO = "Down in mzKO"), name = NULL) +
    labs(title = "Matching-stage leading-edge genes",
         subtitle = "Stars show the project stage-specific DE FDR; effects are MO−MY",
         x = "MO−MY log2 fold change", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    publication_theme() + theme(legend.position = "top")
}

plot_leading_expression <- function(data) {
  plot_data <- copy(data)
  gene_levels <- rev(unique(plot_data$display_gene))
  plot_data[, gene_index := match(display_gene, gene_levels)]
  plot_data[, within_group_rank := frank(sample_id, ties.method = "first"), by = .(display_gene, age_group)]
  plot_data[, y_position := gene_index + ifelse(age_group == "MY", -0.12, 0.12) +
              (within_group_rank - mean(within_group_rank)) * 0.035,
            by = .(display_gene, age_group)]
  medians <- plot_data[, .(vst_expression = median(vst_expression),
                           y_position = unique(gene_index) + ifelse(unique(age_group) == "MY", -0.12, 0.12)),
                       by = .(display_gene, age_group)]
  ggplot(plot_data, aes(vst_expression, y_position, color = age_group, shape = age_group)) +
    geom_point(size = 1.6, alpha = 0.82, stroke = 0.45) +
    geom_point(data = medians, aes(vst_expression, y_position, color = age_group),
               inherit.aes = FALSE, shape = 95, size = 5.5, show.legend = FALSE) +
    scale_color_manual(values = c(MY = "#0072B2", MO = "#D55E00"), name = "Maternal age") +
    scale_shape_manual(values = c(MY = 16, MO = 17), name = "Maternal age") +
    scale_y_continuous(breaks = seq_along(gene_levels), labels = gene_levels,
                       limits = c(0.5, length(gene_levels) + 0.5), expand = expansion(mult = 0)) +
    labs(title = "Observed expression of matching-stage leading-edge genes",
         subtitle = "Each point is one Smart-seq library; horizontal marks denote medians",
         x = "VST expression", y = NULL) +
    publication_theme() + theme(legend.position = "top")
}

plot_pabpn_concordance <- function(data, summary) {
  annotation <- summary[, paste0(stage, ": rho = ", sprintf("%.2f", spearman_rho_eps1),
                                  ", mapped genes = ", mapped_genes)]
  stages <- unique(as.character(data$stage))
  color_values <- c(L2C = "#009E73", `8C` = "#E69F00")
  missing_colors <- setdiff(stages, names(color_values))
  if (length(missing_colors)) {
    extras <- c("#0072B2", "#CC79A7", "#56B4E9", "#D55E00")
    color_values <- c(color_values, setNames(extras[seq_along(missing_colors)], missing_colors))
  }
  shape_values <- setNames(rep(c(16, 17, 15, 18), length.out = length(stages)), stages)
  ggplot(data, aes(effect_eps1, project_log2fc, color = stage, shape = stage)) +
    geom_hline(yintercept = 0, color = "#C6CBD1", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#C6CBD1", linewidth = 0.3) +
    geom_point(size = 0.65, alpha = 0.22, stroke = 0.25) +
    annotate("text", x = -Inf, y = Inf, label = paste(annotation, collapse = "\n"),
             hjust = -0.05, vjust = 1.15, size = 2.5, family = "Liberation Sans", color = "#202124") +
    scale_color_manual(values = color_values[stages], name = "Stage") +
    scale_shape_manual(values = shape_values, name = "Stage") +
    labs(title = "PABPN1-depletion and maternal-age effect concordance",
         subtitle = "Descriptive gene-level comparison; no external differential P or FDR was inferred",
         x = "PABPN1 depletion log2 effect (pseudocount = 1)",
         y = "Project MO−MY log2 fold change") +
    publication_theme() + theme(legend.position = "top")
}

run_perturbation_signatures <- function(args) {
  script_path <- analysis_script_path()
  cfg <- list(
    min_size = 10L, max_size = 5000L, nperm_simple = 10000L, seed = 20261015L,
    top_leading_edge_per_signature = 5L, expected_family_size = NULL,
    stage_order = NULL,
    core_contract_version = "sci-core-2026-10-15-v2",
    core_contract_sha256 = "e443ec5bc51799b498f507ea6284acfcbb56b2d6690a9c693f8f1652d6395f9a",
    analysis_contract_version = "sci-analysis-routing-2026-10-13-v9",
    analysis_contract_sha256 = "830d7b2e97febd9c33f263a338a3fef452ec3171a07144443129eea65c836070",
    figure_contract_version = "sci-figure-task-loading-2026-10-13-v5",
    figure_contract_sha256 = "980e66aadb4d8f3e30c92b786ef4c3c464f443da69cf3c4ff8da406d9284fff1"
  )
  if (nzchar(args$config %||% "")) cfg <- utils::modifyList(cfg, jsonlite::fromJSON(args$config, simplifyVector = TRUE))
  cfg$min_size <- as.integer(cfg$min_size)
  cfg$max_size <- as.integer(cfg$max_size)
  cfg$nperm_simple <- as.integer(cfg$nperm_simple)
  cfg$seed <- as.integer(cfg$seed)
  cfg$top_leading_edge_per_signature <- as.integer(cfg$top_leading_edge_per_signature)
  output_dir <- normalizePath(args$output, winslash = "/", mustWork = FALSE)
  dir.create(file.path(output_dir, "tables", "plot_data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

  ranked <- read_tsv(args$ranked)
  vst <- read_tsv(args$vst)
  samples <- read_tsv(args$samples)
  signatures <- read_tsv(args$obox_signatures)
  pabpn <- read_tsv(args$pabpn_effects)
  assert_columns(ranked, c("stage", "gene_id", "gene_symbol", "stat", "log2FC", "padj"), "ranked table")
  assert_columns(vst, "gene_id", "VST matrix")
  assert_columns(samples, c("sample_id", "stage", "age_group"), "sample metadata")
  assert_columns(signatures, c("signature_id", "source_stage", "direction", "gene_symbol"), "OBOX signatures")
  assert_columns(pabpn, c("stage", "gene_symbol", "effect_eps1", "effect_eps0_1", "mean_KD", "mean_WT", "pass_expression"), "PABPN1 effects")
  ranked[, stage := as.character(stage)]
  signatures[, gene_symbol := trimws(as.character(gene_symbol))]
  stages <- cfg$stage_order %||% unique(ranked$stage)
  stages <- as.character(stages[stages %in% unique(ranked$stage)])
  if (!length(stages)) stop("No configured stages occur in the ranked table")
  planned <- length(stages) * uniqueN(signatures$signature_id)
  if (is.null(cfg$expected_family_size)) cfg$expected_family_size <- planned
  if (as.integer(cfg$expected_family_size) != planned) stop("expected_family_size must equal stages × signatures (", planned, ")")

  projection <- run_signature_projection(ranked, signatures, stages, cfg)
  selected <- select_matching_leading_edges(projection, ranked, cfg$top_leading_edge_per_signature)
  expression <- make_expression_table(selected, vst, samples)
  pabpn_concordance <- prepare_pabpn_concordance(pabpn, ranked)

  tables_dir <- file.path(output_dir, "tables")
  write_tsv(projection$results, file.path(tables_dir, "obox_fgsea_results.tsv"))
  write_tsv(projection$mapping, file.path(tables_dir, "obox_signature_mapping.tsv"))
  write_tsv(projection$leading, file.path(tables_dir, "obox_leading_edges.tsv"))
  write_tsv(projection$ranking_audit, file.path(tables_dir, "ranking_symbol_audit.tsv"))
  write_tsv(selected, file.path(tables_dir, "obox_matching_leading_edge_display.tsv"))
  write_tsv(expression, file.path(tables_dir, "obox_matching_leading_edge_expression.tsv"))
  write_tsv(pabpn_concordance$full, file.path(tables_dir, "pabpn1_project_effect_concordance.tsv"))
  write_tsv(pabpn_concordance$summary, file.path(tables_dir, "pabpn1_concordance_summary.tsv"))

  signature_meta <- unique(signatures[, .(signature_id, source_stage, direction)])
  source_a <- file.path(tables_dir, "plot_data", "14a_obox_signature_concordance.tsv")
  source_b <- file.path(tables_dir, "plot_data", "14b_matching_leading_edge_effects.tsv")
  source_c <- file.path(tables_dir, "plot_data", "14c_matching_leading_edge_expression.tsv")
  source_d <- file.path(tables_dir, "plot_data", "14d_pabpn1_effect_concordance.tsv")
  write_tsv(projection$results, source_a)
  write_tsv(selected, source_b)
  write_tsv(expression, source_c)
  write_tsv(pabpn_concordance$full, source_d)

  figures <- list()
  media_a <- save_figure(plot_obox_matrix(projection$results, stages, signature_meta),
                         file.path(output_dir, "figures", "14a_obox_signature_concordance"), 183, 105)
  figures[[1L]] <- make_figure_manifest_row(
    "14a_obox_signature_concordance", "Figure PS1", media_a, output_dir, source_a,
    "At which stages do maternal-age transcriptome effects resemble signed OBOX mzKO responses?",
    "Direction-oriented fgsea NES", "MO−MY within each stage", "fgseaSimple with 10,000 preranked permutations; BH over all 24 planned cells",
    "OBOX mzKO signatures across maternal-age contrasts. Positive concordance means the MO−MY ranking follows the signed mzKO response. Stars denote the 24-test family FDR; NT denotes a set below the prespecified minimum size.",
    "Smart-seq library", "Three pooled libraries per age group at each of six stages; public signature genes are features, not replicates.",
    "A 4×6 matrix retains all planned signature-stage cells and exposes non-testable cells in grey.", script_path
  )
  if (nrow(selected)) {
    media_b <- save_figure(plot_leading_effects(selected), file.path(output_dir, "figures", "14b_matching_leading_edge_effects"), 160, 110)
    figures[[2L]] <- make_figure_manifest_row(
      "14b_matching_leading_edge_effects", "Figure PS2", media_b, output_dir, source_b,
      "Which genes drive matching-stage OBOX signature enrichment?", "Stage-specific DE log2 fold change",
      "MO−MY at the OBOX source-matched stage", "Project DESeq2 stage contrast; displayed stars use project DE FDR",
      "Matching-stage leading-edge genes from family-FDR-supported signatures. Points show MO−MY log2 fold changes; stars denote project stage-specific differential-expression FDR and NT denotes unavailable FDR.",
      "Smart-seq library", "Three pooled libraries per age group at the matched stage; genes are displayed features.",
      "A horizontal effect plot preserves direction, gene identity, source stage and project DE support in one coordinate.",
      script_path
    )
    media_c <- save_figure(plot_leading_expression(expression), file.path(output_dir, "figures", "14c_matching_leading_edge_expression"), 183, 115)
    figures[[3L]] <- make_figure_manifest_row(
      "14c_matching_leading_edge_expression", "Figure PS3", media_c, output_dir, source_c,
      "Are leading-edge effects visible across the actual project libraries?", "VST expression",
      "MY and MO libraries at the OBOX source-matched stage", "No additional test; dots expose the model input libraries",
      "Observed expression of matching-stage leading-edge genes. Each point is one Smart-seq library and horizontal marks denote medians.",
      "Smart-seq library", "Three pooled libraries per age group at the matched stage.",
      "A horizontal library-level plot keeps all observations visible without treating genes as replicates.",
      script_path
    )
  }
  media_d <- save_figure(plot_pabpn_concordance(pabpn_concordance$full, pabpn_concordance$summary),
                         file.path(output_dir, "figures", "14d_pabpn1_effect_concordance"), 150, 105)
  figures[[length(figures) + 1L]] <- make_figure_manifest_row(
    "14d_pabpn1_effect_concordance", "Figure PS4", media_d, output_dir, source_d,
    "Do maternal-age effects show descriptive concordance with PABPN1 depletion?",
    "Public PABPN1-depletion log2 effect versus project MO−MY log2 fold change",
    "PABPN1 KD−WT compared with project MO−MY at L2C and 8C", "Descriptive Spearman rho; no external differential P or FDR inferred",
    "PABPN1-depletion and maternal-age effect concordance. Gene-level effects are descriptive references; rho summarizes rank concordance without treating genes as biological replicates.",
    "Smart-seq library and public perturbation library", "Project: three pooled libraries per age group per stage; public source: two KD and two WT libraries per stage.",
    "One shared coordinate compares both stages using color plus shape while retaining every expression-eligible mapped gene.",
    script_path
  )
  figure_manifest <- rbindlist(figures, use.names = TRUE, fill = TRUE)
  write_tsv(figure_manifest, file.path(tables_dir, "figure_manifest.tsv"))

  input_paths <- c(ranked = args$ranked, vst = args$vst, samples = args$samples,
                   obox_signatures = args$obox_signatures, pabpn_effects = args$pabpn_effects)
  input_manifest <- data.table(
    input = names(input_paths), path = normalizePath(input_paths, winslash = "/"),
    sha256 = vapply(input_paths, sha256_file, character(1)), bytes = file.info(input_paths)$size
  )
  write_tsv(input_manifest, file.path(tables_dir, "input_manifest.tsv"))
  claim_map <- data.table(
    claim_id = c("OBOX_STAGE_CONCORDANCE", "OBOX_MATCHING_LEADING_EDGE", "PABPN1_EFFECT_REFERENCE"),
    claim = c(
      "Signed OBOX mzKO gene sets can be projected across the six maternal-age contrasts.",
      "Matching-stage leading-edge genes locate the members driving the descriptive projection.",
      "PABPN1 depletion provides an effect-size reference without reconstructed external significance."
    ),
    input = c(args$obox_signatures, args$obox_signatures, args$pabpn_effects),
    script = script_path,
    table = c("tables/obox_fgsea_results.tsv", "tables/obox_matching_leading_edge_display.tsv", "tables/pabpn1_concordance_summary.tsv"),
    figure = c("figures/14a_obox_signature_concordance.pdf", "figures/14b_matching_leading_edge_effects.pdf", "figures/14d_pabpn1_effect_concordance.pdf"),
    test = c("fgseaSimple with 10,000 permutations and BH over 24 planned cells", "Existing project DESeq2 contrast", "Descriptive Spearman rho only"),
    scope = c("Exploratory perturbation similarity", "Exploratory driver localization", "Descriptive effect concordance"),
    status = "supported"
  )
  write_tsv(claim_map, file.path(tables_dir, "claim_result_map.tsv"))
  writeLines(c(
    "# Method plan",
    "",
    "## Knowledge reuse record",
    "",
    "knowledge_id | candidate_idea | decision | reason | adaptation | validation",
    "--- | --- | --- | --- | --- | ---",
    "doi-10.1038_s41586-023-06428-3 | Signed OBOX mzKO signatures | adapt | Direct early-embryo perturbation reference | Project full MO-MY Wald ranks; exact symbols; 24-cell BH family | Mapping audit, leading edges and actual-library display",
    "doi-10.1093_nar_gkab1213 | PABPN1 depletion effects | adapt | RNA-processing perturbation at L2C/8C | ERCC-FPKM effects only; no reconstructed external significance | Pseudocount sensitivity and sign/rank concordance"
  ), file.path(output_dir, "logs", "method_plan.md"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "logs", "sessionInfo.txt"))
  provenance <- data.table(
    provenance_key = c("generated", "command", "seed", "min_size", "max_size", "nperm_simple", "permutation_p_resolution", "family_size",
            "core_contract", "analysis_contract", "figure_contract", "inference_status", "comparison_direction"),
    value = c(
      format(Sys.time(), tz = "UTC", usetz = TRUE), paste(commandArgs(), collapse = " "), cfg$seed,
      cfg$min_size, cfg$max_size, cfg$nperm_simple,
      paste0("1/(n_same_sign+1), where n_same_sign <= ", cfg$nperm_simple,
             "; theoretical lower bound 1/(", cfg$nperm_simple, "+1) before family adjustment"),
      cfg$expected_family_size,
      paste(cfg$core_contract_version, cfg$core_contract_sha256, sep = " | "),
      paste(cfg$analysis_contract_version, cfg$analysis_contract_sha256, sep = " | "),
      paste(cfg$figure_contract_version, cfg$figure_contract_sha256, sep = " | "),
      "exploratory", "Positive project statistics and log2FC denote MO higher than MY"
    )
  )
  write_tsv(provenance, file.path(output_dir, "logs", "provenance.tsv"))
  write_tsv(data.table(status = "ready", message = "Perturbation signature projection completed; upstream differential model was not refit."),
            file.path(output_dir, "status.tsv"))
  invisible(list(projection = projection$results, pabpn = pabpn_concordance$summary, figures = figure_manifest))
}

if (sys.nframe() == 0L) {
  run_perturbation_signatures(parse_cli(commandArgs(trailingOnly = TRUE)))
}
