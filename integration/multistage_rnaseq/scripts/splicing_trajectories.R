#!/usr/bin/env Rscript

# Cross-stage rMATS event alignment and publication-ready PSI/Sashimi output.
# The script is intentionally self-contained so a project can freeze a local copy.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(as.character(x[[1L]]))) y else x
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)

parse_args <- function(x) {
  out <- list(skip_sashimi = FALSE)
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (key == "--skip-sashimi") {
      out$skip_sashimi <- TRUE
      i <- i + 1L
    } else {
      if (i == length(x)) stop("Missing value for ", key)
      value <- x[[i + 1L]]
      key <- sub("^--", "", key)
      key <- gsub("-", "_", key, fixed = TRUE)
      out[[key]] <- value
      i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("config_dir", "rnaseq_root", "output")
missing_args <- required_args[!vapply(required_args, function(x) nzchar(args[[x]] %||% ""), logical(1))]
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

normalize_existing <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
normalize_output <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)
config_dir <- normalize_existing(args$config_dir)
rnaseq_root <- normalize_existing(args$rnaseq_root)
output_root <- normalize_output(args$output)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(output_root, "tables")
plots_dir <- file.path(output_root, "plots")
logs_dir <- file.path(output_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

read_key_value <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(x), c("key", "value")) || anyDuplicated(x$key)) {
    stop("analysis.tsv must contain unique key/value rows")
  }
  stats::setNames(as.character(x$value), x$key)
}

split_csv <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
as_bool <- function(x) tolower(x) %in% c("true", "1", "yes")
analysis <- read_key_value(file.path(config_dir, "analysis.tsv"))
events <- read.delim(file.path(config_dir, "events.tsv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))

required_settings <- c("stages", "reference_group", "comparison_group", "count_modes", "plot_count_mode", "sashimi_required", "container_image")
if (!all(required_settings %in% names(analysis))) stop("analysis.tsv is missing required settings")
required_event_cols <- c("event_id", "event_type", "gene", "chrom", "strand", "sashimi_stage")
if (!all(required_event_cols %in% names(events)) || anyDuplicated(events$event_id)) stop("events.tsv has invalid or duplicate event definitions")

stages <- split_csv(analysis[["stages"]])
reference_group <- analysis[["reference_group"]]
comparison_group <- analysis[["comparison_group"]]
count_modes <- split_csv(analysis[["count_modes"]])
plot_count_mode <- analysis[["plot_count_mode"]]
if (!plot_count_mode %in% count_modes) stop("plot_count_mode must be one of count_modes")
if (reference_group == comparison_group) stop("reference and comparison groups must differ")
sashimi_required <- as_bool(analysis[["sashimi_required"]])

event_coordinate_columns <- list(
  SE = c("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
  A5SS = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"),
  A3SS = c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"),
  MXE = c("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
  RI = c("riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
)
if (!all(events$event_type %in% names(event_coordinate_columns))) stop("Unsupported event type in events.tsv")
event_types <- split_csv(analysis[["event_types"]] %||% paste(names(event_coordinate_columns), collapse = ","))
if (!all(event_types %in% names(event_coordinate_columns)) || !all(events$event_type %in% event_types)) {
  stop("event_types must be supported and include every configured target event type")
}
for (i in seq_len(nrow(events))) {
  needed_coords <- event_coordinate_columns[[events$event_type[[i]]]]
  if (!all(needed_coords %in% names(events)) || anyNA(events[i, needed_coords, drop = TRUE])) {
    stop("events.tsv lacks complete coordinates for ", events$event_id[[i]])
  }
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(out)) return(NA_character_)
  strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

contract_record <- function(component, filename) {
  root <- args$contract_root %||% ""
  if (!nzchar(root)) {
    return(data.frame(component = component, version = "not_provided", sha256 = NA_character_, stringsAsFactors = FALSE))
  }
  path <- file.path(normalize_existing(root), filename)
  if (!file.exists(path)) stop("Contract file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  version_lines <- grep("^contract_version:[[:space:]]*", lines, value = TRUE)
  if (length(version_lines) != 1L) stop("Expected one contract_version line in ", path)
  version <- trimws(sub("^contract_version:[[:space:]]*", "", version_lines[[1L]]))
  data.frame(component = component, version = version, sha256 = sha256_file(path), stringsAsFactors = FALSE)
}

parse_run_metadata <- function(path) {
  lines <- readLines(path, warn = FALSE)
  keys <- sub("=.*$", "", lines)
  values <- sub("^[^=]*=", "", lines)
  if (anyDuplicated(keys)) stop("Duplicate rMATS metadata field: ", path)
  stats::setNames(values, keys)
}

read_stage_design <- function(stage) {
  rmats_dir <- file.path(rnaseq_root, stage, "5.post_transcriptional_regulation", "splicing", "rmats")
  samples_file <- file.path(rnaseq_root, stage, "4.DE_analysis", "samples.txt")
  metadata_file <- file.path(rmats_dir, "run_metadata.tsv")
  if (!file.exists(metadata_file) || !file.exists(samples_file)) stop("Missing rMATS metadata or samples.txt for ", stage)
  md <- parse_run_metadata(metadata_file)
  needed <- c("group1", "group2", "group1_samples", "group2_samples", "sample_key")
  if (!all(needed %in% names(md)) || !setequal(md[c("group1", "group2")], c(reference_group, comparison_group))) {
    stop("Unrecoverable b1/b2 group definition for ", stage)
  }
  b1_samples <- split_csv(md[["group1_samples"]])
  b2_samples <- split_csv(md[["group2_samples"]])
  samples <- read.delim(samples_file, header = FALSE, col.names = c("group", "sample"), stringsAsFactors = FALSE)
  if (ncol(samples) != 2L || anyDuplicated(samples$sample) || !setequal(unique(samples$group), c(reference_group, comparison_group))) {
    stop("Invalid samples.txt for ", stage)
  }
  if (!setequal(b1_samples, samples$sample[samples$group == md[["group1"]]]) ||
      !setequal(b2_samples, samples$sample[samples$group == md[["group2"]]])) {
    stop("rMATS metadata and samples.txt disagree for ", stage)
  }
  expected_key <- paste(sort(c(b1_samples, b2_samples)), collapse = ",")
  if (!identical(expected_key, md[["sample_key"]])) stop("sample_key mismatch for ", stage)
  list(
    stage = stage, rmats_dir = rmats_dir, metadata_file = metadata_file, samples_file = samples_file,
    b1_group = md[["group1"]], b2_group = md[["group2"]], b1_samples = b1_samples, b2_samples = b2_samples,
    reference_samples = if (md[["group1"]] == reference_group) b1_samples else b2_samples,
    comparison_samples = if (md[["group1"]] == comparison_group) b1_samples else b2_samples,
    multiplier = if (md[["group1"]] == comparison_group) 1 else -1
  )
}

stage_designs <- stats::setNames(lapply(stages, read_stage_design), stages)
all_coordinate_columns <- unique(unlist(event_coordinate_columns, use.names = FALSE))
normalise_chrom <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
numeric_or_na <- function(x) suppressWarnings(as.numeric(as.character(x)))

make_event_key <- function(x, event_type) {
  coords <- event_coordinate_columns[[event_type]]
  paste(c(event_type, normalise_chrom(x[["chr"]]), as.character(x[["strand"]]), as.character(x[coords])), collapse = "|")
}

make_event_keys <- function(x, event_type, coord_limit = NULL) {
  coords <- event_coordinate_columns[[event_type]]
  if (!is.null(coord_limit)) coords <- head(coords, coord_limit)
  parts <- c(
    list(rep(event_type, nrow(x)), normalise_chrom(x[["chr"]]), as.character(x[["strand"]])),
    lapply(coords, function(column) as.character(x[[column]])),
    list(sep = "|")
  )
  do.call(paste, parts)
}

make_legacy_partial_key <- function(x, event_type) {
  coords <- head(event_coordinate_columns[[event_type]], 4L)
  paste(c(event_type, normalise_chrom(x[["chr"]]), as.character(x[["strand"]]), as.character(x[coords])), collapse = "|")
}

split_values <- function(x) suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))
stars_from_fdr <- function(x) {
  out <- rep("ns", length(x))
  out[is.na(x)] <- "NA"
  out[!is.na(x) & x < 0.05] <- "*"
  out[!is.na(x) & x < 0.01] <- "**"
  out[!is.na(x) & x < 0.001] <- "***"
  out
}

summary_rows <- list()
target_stats <- list()
target_psi <- list()
target_support <- list()
availability <- list()
raw_target_rows <- list()
input_manifest <- list()
validation <- list()
row_id <- 0L

add_validation <- function(check, status, detail) {
  validation[[length(validation) + 1L]] <<- data.frame(check = check, status = status, detail = detail, stringsAsFactors = FALSE)
}

for (stage in stages) {
  design <- stage_designs[[stage]]
  for (path in c(design$metadata_file, design$samples_file)) {
    info <- file.info(path)
    input_manifest[[length(input_manifest) + 1L]] <- data.frame(
      input_type = if (grepl("run_metadata", path)) "rmats_metadata" else "sample_table",
      stage = stage, path = normalize_existing(path), size_bytes = info$size, mtime = as.character(info$mtime),
      sha256 = sha256_file(path), hash_scope = "full_file", stringsAsFactors = FALSE
    )
  }
  for (event_type in event_types) {
    for (count_mode in count_modes) {
      path <- file.path(design$rmats_dir, sprintf("%s.MATS.%s.txt", event_type, count_mode))
      if (!file.exists(path)) stop("Missing rMATS table: ", path)
      info <- file.info(path)
      input_manifest[[length(input_manifest) + 1L]] <- data.frame(
        input_type = "rmats_event_table", stage = stage, path = normalize_existing(path), size_bytes = info$size,
        mtime = as.character(info$mtime), sha256 = sha256_file(path), hash_scope = "full_file", stringsAsFactors = FALSE
      )
      x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = c("NA", ""))
      coords <- event_coordinate_columns[[event_type]]
      needed <- c("GeneID", "geneSymbol", "chr", "strand", coords, "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference",
                  "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
      if (!all(needed %in% names(x))) stop("Malformed rMATS table: ", path)
      event_keys <- make_event_keys(x, event_type)
      legacy_keys <- make_event_keys(x, event_type, coord_limit = 4L)
      if (anyDuplicated(event_keys)) stop("Duplicate full-coordinate event key in ", path)
      base <- data.frame(
        stage = stage, event_type = event_type, count_mode = count_mode,
        event_key = event_keys, legacy_partial_key = legacy_keys,
        gene_id = as.character(x$GeneID), gene_symbol = as.character(x$geneSymbol),
        chrom = normalise_chrom(x$chr), source_chrom = as.character(x$chr), strand = as.character(x$strand),
        pvalue = numeric_or_na(x$PValue), fdr = numeric_or_na(x$FDR),
        delta_psi_raw_b1_minus_b2 = numeric_or_na(x$IncLevelDifference),
        delta_psi_comparison_minus_reference = numeric_or_na(x$IncLevelDifference) * design$multiplier,
        b1_group = design$b1_group, b2_group = design$b2_group,
        source_file = normalize_existing(path), stringsAsFactors = FALSE
      )
      for (column in all_coordinate_columns) base[[column]] <- if (column %in% names(x)) numeric_or_na(x[[column]]) else NA_real_
      summary_rows[[length(summary_rows) + 1L]] <- base

      selected_events <- events[events$event_type == event_type, , drop = FALSE]
      for (j in seq_len(nrow(selected_events))) {
        event <- selected_events[j, , drop = FALSE]
        keep <- normalise_chrom(x$chr) == normalise_chrom(event$chrom) & as.character(x$strand) == as.character(event$strand)
        for (column in coords) keep <- keep & numeric_or_na(x[[column]]) == numeric_or_na(event[[column]])
        hits <- which(keep)
        availability[[length(availability) + 1L]] <- data.frame(
          event_id = event$event_id, stage = stage, count_mode = count_mode, matched_rows = length(hits),
          status = if (length(hits) == 1L) "detected" else if (!length(hits)) "not_detected" else "coordinate_collision",
          source_file = normalize_existing(path), stringsAsFactors = FALSE
        )
        if (length(hits) > 1L) stop("Target coordinates matched multiple rows: ", event$event_id, " / ", stage, " / ", count_mode)
        if (!length(hits)) next
        i <- hits[[1L]]
        raw_target_rows[[paste(event$event_id, stage, count_mode, sep = "|")]] <- list(row = x[i, , drop = FALSE], path = path)
        stat <- base[i, , drop = FALSE]
        stat$event_id <- event$event_id
        stat$target_gene <- event$gene
        stat$significance <- stars_from_fdr(stat$fdr)
        target_stats[[length(target_stats) + 1L]] <- stat

        psi_b1 <- split_values(x$IncLevel1[[i]])
        psi_b2 <- split_values(x$IncLevel2[[i]])
        ijc_b1 <- split_values(x$IJC_SAMPLE_1[[i]])
        sjc_b1 <- split_values(x$SJC_SAMPLE_1[[i]])
        ijc_b2 <- split_values(x$IJC_SAMPLE_2[[i]])
        sjc_b2 <- split_values(x$SJC_SAMPLE_2[[i]])
        if (length(psi_b1) != length(design$b1_samples) || length(psi_b2) != length(design$b2_samples) ||
            length(ijc_b1) != length(design$b1_samples) || length(ijc_b2) != length(design$b2_samples)) {
          stop("Per-library rMATS values do not match metadata: ", event$event_id, " / ", stage)
        }
        append_group <- function(group, samples, psi, ijc, sjc) {
          target_psi[[length(target_psi) + 1L]] <<- data.frame(
            event_id = event$event_id, event_type = event_type, gene = event$gene, stage = stage,
            count_mode = count_mode, group = group, sample = samples, psi = psi,
            psi_status = ifelse(is.na(psi), "not_estimable", "observed"), stringsAsFactors = FALSE
          )
          target_support[[length(target_support) + 1L]] <<- data.frame(
            event_id = event$event_id, event_type = event_type, gene = event$gene, stage = stage,
            count_mode = count_mode, group = group, sample = samples,
            inclusion_junction_count = ijc, skipping_junction_count = sjc,
            read_support_status = ifelse(is.na(ijc) | is.na(sjc), "not_available", "observed"), stringsAsFactors = FALSE
          )
        }
        append_group(design$b1_group, design$b1_samples, psi_b1, ijc_b1, sjc_b1)
        append_group(design$b2_group, design$b2_samples, psi_b2, ijc_b2, sjc_b2)
        finite1 <- psi_b1[is.finite(psi_b1)]
        finite2 <- psi_b2[is.finite(psi_b2)]
        if (length(finite1) && length(finite2)) {
          reconstructed <- mean(finite1) - mean(finite2)
          if (abs(reconstructed - numeric_or_na(x$IncLevelDifference[[i]])) > 0.0011) {
            stop("IncLevelDifference is inconsistent with per-library PSI: ", event$event_id, " / ", stage)
          }
        }
      }
    }
  }
}

event_stage_summary <- do.call(rbind, summary_rows)
target_event_stats <- if (length(target_stats)) do.call(rbind, target_stats) else data.frame()
target_event_psi <- if (length(target_psi)) do.call(rbind, target_psi) else data.frame()
target_event_support <- if (length(target_support)) do.call(rbind, target_support) else data.frame()
event_availability <- do.call(rbind, availability)
input_manifest_df <- do.call(rbind, input_manifest)

write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_tsv(event_stage_summary, file.path(tables_dir, "event_stage_summary.tsv"))
write_tsv(target_event_stats, file.path(tables_dir, "target_event_stage_statistics.tsv"))
write_tsv(target_event_psi, file.path(tables_dir, "target_event_psi_long.tsv"))
write_tsv(target_event_support, file.path(tables_dir, "target_event_read_support.tsv"))
write_tsv(event_availability, file.path(tables_dir, "target_event_availability.tsv"))
write_tsv(input_manifest_df, file.path(tables_dir, "input_manifest.tsv"))

if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required for cross-stage event aggregation")
event_dt <- data.table::as.data.table(event_stage_summary)
meta_cols <- c("event_type", "count_mode", "event_key", "gene_id", "gene_symbol", "chrom", "source_chrom", "strand", all_coordinate_columns)
recurrence_meta <- unique(event_dt[, ..meta_cols])
if (recurrence_meta[, .N, by = .(event_key, count_mode)][, any(N > 1L)]) {
  stop("One full-coordinate event key maps to conflicting gene or coordinate metadata")
}
recurrence_counts <- event_dt[, .(
  detected_stage_n = data.table::uniqueN(stage),
  detected_stages = paste(stages[stages %in% unique(stage)], collapse = ";")
), by = .(event_key, count_mode)]
recurrence <- merge(recurrence_meta, recurrence_counts, by = c("event_key", "count_mode"), sort = FALSE)
recurrence <- as.data.frame(recurrence)
write_tsv(recurrence, file.path(tables_dir, "cross_stage_event_recurrence.tsv"))

legacy_map <- unique(event_dt[, .(event_type, count_mode, legacy_partial_key, event_key)])
legacy_collision_n <- legacy_map[, .(full_key_n = data.table::uniqueN(event_key)), by = .(event_type, count_mode, legacy_partial_key)][full_key_n > 1L, .N]
add_validation("full_coordinate_event_key_unique", "PASS", sprintf("%d standardized event-stage rows; no duplicate full keys within a stage/mode", nrow(event_stage_summary)))
add_validation("legacy_partial_coordinate_collisions", "PASS", sprintf("%d partial keys map to multiple full-coordinate keys and remain separated", legacy_collision_n))
add_validation("psi_na_preserved", "PASS", sprintf("%d target PSI values retained as NA", sum(is.na(target_event_psi$psi))))
add_validation("contrast_direction", "PASS", sprintf("All effects oriented as %s - %s from run_metadata.tsv", comparison_group, reference_group))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
suppressPackageStartupMessages(library(ggplot2))
group_colors <- stats::setNames(
  c(analysis[["reference_color"]] %||% "#0072B2", analysis[["comparison_color"]] %||% "#D55E00"),
  c(reference_group, comparison_group)
)
group_labels <- stats::setNames(
  c(analysis[["reference_label"]] %||% reference_group, analysis[["comparison_label"]] %||% comparison_group),
  c(reference_group, comparison_group)
)
group_shapes <- stats::setNames(c(16, 17), c(reference_group, comparison_group))

plot_rows <- target_event_psi[target_event_psi$count_mode == plot_count_mode, , drop = FALSE]
stat_rows <- target_event_stats[target_event_stats$count_mode == plot_count_mode, , drop = FALSE]
availability_plot <- event_availability[event_availability$count_mode == plot_count_mode, , drop = FALSE]
figure_manifest <- list()

for (event_id in events$event_id) {
  event <- events[events$event_id == event_id, , drop = FALSE]
  z <- plot_rows[plot_rows$event_id == event_id, , drop = FALSE]
  s <- stat_rows[stat_rows$event_id == event_id, , drop = FALSE]
  a <- availability_plot[availability_plot$event_id == event_id, , drop = FALSE]
  z$stage <- factor(z$stage, levels = stages)
  z$group <- factor(z$group, levels = c(reference_group, comparison_group))
  z$stage_index <- match(as.character(z$stage), stages)
  z$group_offset <- ifelse(z$group == reference_group, -0.13, 0.13)
  z$replicate_index <- ave(seq_len(nrow(z)), interaction(z$stage, z$group, drop = TRUE), FUN = seq_along)
  z$replicate_n <- ave(seq_len(nrow(z)), interaction(z$stage, z$group, drop = TRUE), FUN = length)
  z$jitter <- ifelse(z$replicate_n <= 1L, 0, (z$replicate_index - (z$replicate_n + 1) / 2) * 0.045)
  z$x_plot <- z$stage_index + z$group_offset + z$jitter
  means <- aggregate(psi ~ stage + group, data = z[is.finite(z$psi), , drop = FALSE], FUN = mean)
  means$stage_index <- match(as.character(means$stage), stages)
  means$group_offset <- ifelse(means$group == reference_group, -0.13, 0.13)
  means$x_plot <- means$stage_index + means$group_offset
  s$x_plot <- match(s$stage, stages)
  missing_stage <- stages[!stages %in% a$stage[a$status == "detected"]]
  missing_df <- data.frame(x_plot = match(missing_stage, stages), psi = rep(-0.055, length(missing_stage)))
  na_df <- z[is.na(z$psi), , drop = FALSE]
  na_df$psi <- rep(-0.055, nrow(na_df))
  title <- sprintf("%s %s inclusion", event$gene, event$event_type)
  p <- ggplot() +
    geom_line(data = means, aes(x = x_plot, y = psi, color = group, linetype = group, group = group), linewidth = 0.55, na.rm = TRUE) +
    geom_point(data = means, aes(x = x_plot, y = psi, color = group, shape = group), size = 2.1, stroke = 0.35, na.rm = TRUE) +
    geom_point(data = z[is.finite(z$psi), , drop = FALSE], aes(x = x_plot, y = psi, color = group, shape = group), size = 1.45, alpha = 0.58, stroke = 0.25) +
    geom_point(data = na_df, aes(x = x_plot, y = psi), shape = 4, color = "#777777", size = 1.8, stroke = 0.55) +
    geom_point(data = missing_df, aes(x = x_plot, y = psi), shape = 95, color = "#A0A0A0", size = 4.2) +
    geom_text(data = s, aes(x = x_plot, y = 1.055, label = significance), family = "Liberation Sans", size = 2.65, vjust = 0) +
    scale_x_continuous(breaks = seq_along(stages), labels = stages, limits = c(0.55, length(stages) + 0.45), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(-0.09, 1.12), expand = c(0, 0)) +
    scale_color_manual(values = group_colors, labels = group_labels, drop = FALSE) +
    scale_shape_manual(values = group_shapes, labels = group_labels, drop = FALSE) +
    scale_linetype_manual(values = c(2, 1), labels = group_labels, drop = FALSE) +
    labs(x = "Developmental stage", y = "Percent spliced in (PSI)", title = title, color = NULL, shape = NULL, linetype = NULL) +
    theme_classic(base_family = "Liberation Sans", base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "plain", hjust = 0), axis.title = element_text(size = 8),
      axis.text = element_text(size = 7, color = "black"), legend.text = element_text(size = 7),
      legend.position = "top", legend.justification = "left", legend.box.spacing = grid::unit(0.5, "mm"),
      axis.line = element_line(linewidth = 0.35, color = "black"), axis.ticks = element_line(linewidth = 0.35, color = "black"),
      plot.margin = margin(4, 4, 3, 4, unit = "mm")
    )
  stem <- sprintf("15_%s_PSI_trajectory", event_id)
  pdf_path <- file.path(plots_dir, paste0(stem, ".pdf"))
  png_path <- file.path(plots_dir, paste0(stem, ".png"))
  ggsave(pdf_path, p, width = 110, height = 78, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(png_path, p, width = 110, height = 78, units = "mm", dpi = 300, bg = "white")
  figure_manifest[[length(figure_manifest) + 1L]] <- data.frame(
    asset_id = paste0("splicing_", event_id, "_psi"), asset_class = "standalone", figure_id = stem, panel_label = "",
    pdf = normalizePath(pdf_path, winslash = "/"), png = normalizePath(png_path, winslash = "/"),
    source_table = normalizePath(file.path(tables_dir, "target_event_psi_long.tsv"), winslash = "/"),
    statistics_table = normalizePath(file.path(tables_dir, "target_event_stage_statistics.tsv"), winslash = "/"),
    question = "How does event inclusion differ between groups across observed developmental stages?",
    experimental_unit = "Smart-seq library", comparison = paste(comparison_group, reference_group, sep = " - "),
    significance_metric = "rMATS FDR within each stage", annotation_rule = "* FDR<0.05; ** FDR<0.01; *** FDR<0.001; ns otherwise",
    inference_scope = "Cross-sectional per-stage event inclusion; no stage-by-group interaction test",
    layout_rationale = "Single 110 x 78 mm coordinate retains all stages, individual libraries, missing events and NA PSI.",
    width_mm = 110, height_mm = 78, dpi = 300, coordinates_pdf = 1, coordinates_png = 1,
    stringsAsFactors = FALSE
  )
}

render_sashimi <- function(event) {
  event_id <- event$event_id
  stage <- as.character(event$sashimi_stage)
  status <- data.frame(
    event_id = event_id, stage = stage, status = "not_requested", message = "",
    b1_group = NA_character_, b2_group = NA_character_, b1_label = NA_character_, b2_label = NA_character_,
    b1_samples = NA_character_, b2_samples = NA_character_, b1_inc_levels = NA_character_, b2_inc_levels = NA_character_,
    group_order_check = NA_character_,
    track_aggregation = NA_character_, inc_level_aggregation = NA_character_, junction_arc_aggregation = NA_character_,
    pdf = NA_character_, png = NA_character_, stringsAsFactors = FALSE
  )
  if (is.na(stage) || !nzchar(stage)) return(status)
  key <- paste(event_id, stage, plot_count_mode, sep = "|")
  target <- raw_target_rows[[key]]
  if (is.null(target)) {
    status$status <- "event_not_detected"
    status$message <- "No matching event row at requested Sashimi stage"
    return(status)
  }
  design <- stage_designs[[stage]]
  status$b1_group <- design$b1_group
  status$b2_group <- design$b2_group
  status$b1_label <- group_labels[[design$b1_group]]
  status$b2_label <- group_labels[[design$b2_group]]
  status$b1_samples <- paste(design$b1_samples, collapse = ",")
  status$b2_samples <- paste(design$b2_samples, collapse = ",")
  status$b1_inc_levels <- as.character(target$row$IncLevel1)
  status$b2_inc_levels <- as.character(target$row$IncLevel2)
  status$track_aggregation <- "Pooled base coverage normalized by total mapped reads across BAMs in each rMATS group"
  status$inc_level_aggregation <- "Arithmetic mean of per-library rMATS IncLevel within each displayed group"
  status$junction_arc_aggregation <- "Rounded arithmetic mean of raw junction-supporting read counts across BAMs in each displayed group"
  bam_dir <- file.path(rnaseq_root, stage, analysis[["bam_subdir"]] %||% "1.Mapping")
  b1_bams <- file.path(bam_dir, paste0(design$b1_samples, ".bam"))
  b2_bams <- file.path(bam_dir, paste0(design$b2_samples, ".bam"))
  bams <- c(b1_bams, b2_bams)
  bais <- paste0(bams, ".bai")
  missing <- c(bams[!file.exists(bams)], bais[!file.exists(bais)])
  if (length(missing)) {
    status$status <- "missing_bam"
    status$message <- paste(missing, collapse = ";")
    return(status)
  }
  for (path in c(bams, bais)) {
    info <- file.info(path)
    input_manifest[[length(input_manifest) + 1L]] <<- data.frame(
      input_type = if (grepl("[.]bai$", path)) "bam_index" else "bam",
      stage = stage, path = normalize_existing(path), size_bytes = info$size, mtime = as.character(info$mtime),
      sha256 = if (grepl("[.]bai$", path)) sha256_file(path) else NA_character_,
      hash_scope = if (grepl("[.]bai$", path)) "full_file" else "size_mtime_only_large_bam", stringsAsFactors = FALSE
    )
  }
  samtools <- Sys.which("samtools")
  if (!nzchar(samtools)) {
    status$status <- "missing_samtools"
    status$message <- "samtools not found"
    return(status)
  }
  target_contig <- normalise_chrom(event$chrom)
  for (bam in bams) {
    contigs <- system2(samtools, c("idxstats", bam), stdout = TRUE, stderr = TRUE)
    present <- sub("\t.*$", "", contigs)
    if (!target_contig %in% present) {
      status$status <- "contig_not_found"
      status$message <- paste(target_contig, "not in", bam)
      return(status)
    }
  }
  singularity <- Sys.which("singularity")
  if (!nzchar(singularity)) {
    status$status <- "missing_executor"
    status$message <- "singularity not found"
    return(status)
  }
  event_input_dir <- file.path(tables_dir, "sashimi_events")
  raw_dir <- file.path(logs_dir, "sashimi_raw", event_id)
  dir.create(event_input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  event_file <- file.path(event_input_dir, sprintf("%s_%s_%s.MATS.%s.txt", event_id, stage, event$event_type, plot_count_mode))
  write.table(target$row, event_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  group_file <- file.path(event_input_dir, sprintf("%s_%s_groups.txt", event_id, stage))
  writeLines(c(sprintf("%s: 1-%d", group_labels[[design$b1_group]], length(b1_bams)),
               sprintf("%s: %d-%d", group_labels[[design$b2_group]], length(b1_bams) + 1L, length(bams))), group_file)
  mpl_dir <- tempfile(pattern = "mpl_", tmpdir = logs_dir)
  dir.create(mpl_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(mpl_dir, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines(c(
    "font.family: Liberation Sans", "font.sans-serif: Liberation Sans", "font.size: 9",
    "pdf.fonttype: 42", "ps.fonttype: 42", "figure.facecolor: white", "savefig.facecolor: white"
  ), file.path(mpl_dir, "matplotlibrc"))
  font_dir <- "/usr/share/fonts/truetype/liberation"
  bind_args <- c("--bind", paste0(rnaseq_root, ":", rnaseq_root), "--bind", paste0(output_root, ":", output_root))
  if (dir.exists(font_dir)) bind_args <- c(bind_args, "--bind", paste0(font_dir, ":", font_dir))
  command_args <- c(
    "exec", bind_args, paste0("docker://", analysis[["container_image"]]), "rmats2sashimiplot",
    "--b1", paste(b1_bams, collapse = ","), "--b2", paste(b2_bams, collapse = ","),
    "--event-type", event$event_type, "-e", event_file,
    "--l1", group_labels[[design$b1_group]], "--l2", group_labels[[design$b2_group]],
    "--group-info", group_file, "--exon_s", analysis[["exon_scale"]] %||% "1",
    "--intron_s", analysis[["intron_scale"]] %||% "20", "--min-counts", analysis[["min_counts"]] %||% "1",
    "--fig-width", "7.2047", "--fig-height", "3.9370", "--font-size", "9",
    "--color", shQuote(paste(group_colors[c(design$b1_group, design$b2_group)], collapse = ",")),
    "--remove-event-chr-prefix", "-o", raw_dir
  )
  log_file <- file.path(logs_dir, paste0("sashimi_", event_id, "_rmats2.log"))
  status_code <- system2(
    singularity, command_args, stdout = log_file, stderr = log_file,
    env = c(paste0("MPLCONFIGDIR=", mpl_dir), paste0("SINGULARITYENV_MPLCONFIGDIR=", mpl_dir))
  )
  if (!identical(status_code, 0L)) {
    status$status <- "render_failed"
    status$message <- sprintf("rmats2sashimiplot exit=%s; log=%s", status_code, log_file)
    return(status)
  }
  index_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
  index_dirs <- index_dirs[grepl("^Sashimi_index_", basename(index_dirs))]
  settings_files <- if (length(index_dirs) == 1L) file.path(index_dirs, "sashimi_plot_settings.txt") else character()
  pickle_files <- if (length(index_dirs) == 1L) list.files(index_dirs, pattern = "[.]pickle$", recursive = TRUE, full.names = TRUE) else character()
  if (length(pickle_files) != 1L || length(index_dirs) != 1L || !file.exists(settings_files)) {
    status$status <- "render_failed"
    status$message <- sprintf("Could not recover one native MISO event/index/settings set; log=%s", log_file)
    return(status)
  }
  event_name <- sub("[.]pickle$", "", basename(pickle_files[[1L]]))
  setting_values <- function(lines, key) {
    line <- grep(paste0("^", key, " = "), lines, value = TRUE)
    if (length(line) != 1L) return(character())
    hits <- regmatches(line, gregexpr('"[^"]*"', line))[[1L]]
    gsub('^"|"$', "", hits)
  }
  settings_lines <- readLines(settings_files, warn = FALSE)
  setting_bams <- setting_values(settings_lines, "bam_files")
  setting_labels <- setting_values(settings_lines, "sample_labels")
  setting_colors <- setting_values(settings_lines, "colors")
  source_inc1 <- split_values(target$row$IncLevel1)
  source_inc2 <- split_values(target$row$IncLevel2)
  setting_inc <- suppressWarnings(as.numeric(sub("^.*IncLevel:[[:space:]]*", "", setting_labels)))
  group_lines <- readLines(group_file, warn = FALSE)
  expected_group_lines <- c(
    sprintf("%s: 1-%d", group_labels[[design$b1_group]], length(b1_bams)),
    sprintf("%s: %d-%d", group_labels[[design$b2_group]], length(b1_bams) + 1L, length(bams))
  )
  b1_index <- seq_along(b1_bams)
  b2_index <- length(b1_bams) + seq_along(b2_bams)
  group_order_ok <-
    identical(setting_bams, bams) &&
    identical(group_lines, expected_group_lines) &&
    identical(setting_colors, unname(group_colors[c(design$b1_group, design$b2_group)])) &&
    length(setting_labels) == length(bams) &&
    all(grepl(paste0(" ", group_labels[[design$b1_group]], "-"), setting_labels[b1_index], fixed = TRUE)) &&
    all(grepl(paste0(" ", group_labels[[design$b2_group]], "-"), setting_labels[b2_index], fixed = TRUE)) &&
    length(source_inc1) == length(b1_bams) && length(source_inc2) == length(b2_bams) &&
    isTRUE(all.equal(setting_inc[b1_index], round(source_inc1, 2), check.attributes = FALSE)) &&
    isTRUE(all.equal(setting_inc[b2_index], round(source_inc2, 2), check.attributes = FALSE))
  status$group_order_check <- if (group_order_ok) "PASS" else "FAIL"
  if (!group_order_ok) {
    status$status <- "group_order_mismatch"
    status$message <- "Generated MISO settings do not preserve source rMATS b1/b2 BAM, label, color and IncLevel order"
    return(status)
  }
  final_dir <- file.path(raw_dir, "Sashimi_final")
  if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  native_log <- file.path(logs_dir, paste0("sashimi_", event_id, "_native.log"))
  native_script <- "/usr/local/lib/python3.10/site-packages/MISO/misopy/sashimi_plot/sashimi_plot.py"
  short_title <- sprintf("%s %s | %s", event$gene, event$event_type, stage)
  native_args <- c(
    "exec", bind_args, paste0("docker://", analysis[["container_image"]]), "python", native_script,
    "--plot-event", event_name, index_dirs[[1L]], settings_files[[1L]], "--no-posteriors",
    "--output-dir", final_dir, "--group-info", group_file,
    "--plot-title", shQuote(short_title)
  )
  native_status <- system2(
    singularity, native_args, stdout = native_log, stderr = native_log,
    env = c(paste0("MPLCONFIGDIR=", mpl_dir), paste0("SINGULARITYENV_MPLCONFIGDIR=", mpl_dir))
  )
  pdfs <- list.files(final_dir, pattern = "[.]pdf$", recursive = FALSE, full.names = TRUE)
  if (!identical(native_status, 0L) || length(pdfs) != 1L) {
    status$status <- "render_failed"
    status$message <- sprintf("native MISO exit=%s; pdf_count=%d; log=%s", native_status, length(pdfs), native_log)
    return(status)
  }
  stem <- sprintf("15_%s_%s_sashimi", event_id, stage)
  pdf_path <- file.path(plots_dir, paste0(stem, ".pdf"))
  png_path <- file.path(plots_dir, paste0(stem, ".png"))
  if (!file.copy(pdfs[[1L]], pdf_path, overwrite = TRUE)) stop("Could not copy Sashimi PDF")
  python <- Sys.which("python3")
  if (!nzchar(python)) stop("python3 with PyMuPDF is required to render Sashimi PNG")
  convert_code <- paste(
    "import fitz,sys",
    "doc=fitz.open(sys.argv[1])",
    "pix=doc[0].get_pixmap(dpi=300,alpha=False)",
    "pix.save(sys.argv[2])",
    sep = ";"
  )
  convert_status <- system2(python, c("-c", shQuote(convert_code), pdf_path, png_path))
  if (!identical(convert_status, 0L) || !file.exists(png_path)) stop("PyMuPDF failed to render Sashimi PNG")
  status$status <- "completed"
  status$message <- "Rendered from indexed BAMs with rMATS b1/b2 order preserved and native MISO title"
  status$pdf <- normalizePath(pdf_path, winslash = "/")
  status$png <- normalizePath(png_path, winslash = "/")
  figure_manifest[[length(figure_manifest) + 1L]] <<- data.frame(
    asset_id = paste0("splicing_", event_id, "_sashimi"), asset_class = "standalone", figure_id = stem, panel_label = "",
    pdf = status$pdf, png = status$png, source_table = normalizePath(event_file, winslash = "/"),
    statistics_table = normalizePath(file.path(tables_dir, "target_event_read_support.tsv"), winslash = "/"),
    question = "Do indexed reads support the selected splice junction structure in each age group?",
    experimental_unit = "Smart-seq library; junction curves grouped for display", comparison = paste(comparison_group, reference_group, sep = " - "),
    significance_metric = "none", annotation_rule = "Junction read counts shown by rmats2sashimiplot; no added significance label",
    inference_scope = "Representative locus-level read support at the configured stage",
    layout_rationale = "Single 183 x 100 mm locus coordinate with two rMATS-ordered grouped junction tracks and exon structure.",
    width_mm = 183, height_mm = 100, dpi = 300, coordinates_pdf = 1, coordinates_png = 1,
    stringsAsFactors = FALSE
  )
  status
}

sashimi_status <- if (isTRUE(args$skip_sashimi)) {
  data.frame(
    event_id = events$event_id, stage = events$sashimi_stage, status = "skipped_by_cli", message = "--skip-sashimi",
    b1_group = NA_character_, b2_group = NA_character_, b1_label = NA_character_, b2_label = NA_character_,
    b1_samples = NA_character_, b2_samples = NA_character_, b1_inc_levels = NA_character_, b2_inc_levels = NA_character_,
    group_order_check = NA_character_,
    track_aggregation = NA_character_, inc_level_aggregation = NA_character_, junction_arc_aggregation = NA_character_,
    pdf = NA_character_, png = NA_character_, stringsAsFactors = FALSE
  )
} else do.call(rbind, lapply(seq_len(nrow(events)), function(i) render_sashimi(events[i, , drop = FALSE])))
write_tsv(sashimi_status, file.path(tables_dir, "sashimi_status.tsv"))

if (sashimi_required && any(sashimi_status$status != "completed")) {
  write_tsv(do.call(rbind, validation), file.path(tables_dir, "validation.tsv"))
  write_tsv(do.call(rbind, input_manifest), file.path(tables_dir, "input_manifest.tsv"))
  stop("Required Sashimi rendering did not complete: ", paste(sashimi_status$event_id[sashimi_status$status != "completed"], collapse = ", "))
}
add_validation("sashimi_execution", if (all(sashimi_status$status %in% c("completed", "not_requested", "missing_bam", "skipped_by_cli"))) "PASS" else "FAIL",
               paste(paste(sashimi_status$event_id, sashimi_status$status, sep = "="), collapse = "; "))
rendered_sashimi <- sashimi_status$status == "completed"
add_validation(
  "sashimi_group_order",
  if (all(sashimi_status$group_order_check[rendered_sashimi] == "PASS")) "PASS" else "FAIL",
  paste(sprintf("%s:%s/%s=%s", sashimi_status$event_id, sashimi_status$b1_group, sashimi_status$b2_group, sashimi_status$group_order_check), collapse = "; ")
)

figure_manifest_df <- do.call(rbind, figure_manifest)
figure_manifest_df$pdf_sha256 <- vapply(figure_manifest_df$pdf, sha256_file, character(1))
figure_manifest_df$png_sha256 <- vapply(figure_manifest_df$png, sha256_file, character(1))
figure_manifest_df$source_sha256 <- vapply(figure_manifest_df$source_table, sha256_file, character(1))
write_tsv(figure_manifest_df, file.path(tables_dir, "figure_manifest.tsv"))
write_tsv(do.call(rbind, input_manifest), file.path(tables_dir, "input_manifest.tsv"))
write_tsv(do.call(rbind, validation), file.path(tables_dir, "validation.tsv"))

software <- rbind(
  data.frame(
    component = c("R", "ggplot2", "data.table", "rmats2sashimiplot_container"),
    version = c(as.character(getRversion()), as.character(utils::packageVersion("ggplot2")),
                as.character(utils::packageVersion("data.table")), analysis[["container_image"]]),
    sha256 = NA_character_, stringsAsFactors = FALSE
  ),
  contract_record("sci_core_contract", "sci_core_contract.md"),
  contract_record("sci_analysis_contract", "sci_analysis_contract.md"),
  contract_record("sci_figure_contract", "sci_figure_contract.md")
)
write_tsv(software, file.path(tables_dir, "software_versions.tsv"))

status <- data.frame(
  status = "completed", standardized_event_stage_rows = nrow(event_stage_summary),
  recurrent_event_rows = nrow(recurrence), target_psi_rows = nrow(target_event_psi),
  psi_figures = sum(grepl("_psi$", figure_manifest_df$asset_id)),
  sashimi_figures = sum(grepl("_sashimi$", figure_manifest_df$asset_id)),
  completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), stringsAsFactors = FALSE
)
write_tsv(status, file.path(output_root, "status.tsv"))

provenance <- c(
  "# Splicing trajectory provenance",
  "",
  sprintf("- Contrast: `%s - %s`.", comparison_group, reference_group),
  sprintf("- Stages: `%s`.", paste(stages, collapse = ", ")),
  sprintf("- rMATS count modes: `%s`; PSI figures use `%s`.", paste(count_modes, collapse = ", "), plot_count_mode),
  "- Direction is reconstructed from each stage's `run_metadata.tsv`; original rMATS columns are not rewritten.",
  "- Full event keys contain event type, normalized chromosome, strand and every event-specific coordinate.",
  "- Exact P values and FDR values remain in `target_event_stage_statistics.tsv`; figures show only `* / ** / *** / ns`.",
  "- Sashimi plots use real indexed BAMs. BAM provenance records size and modification time; BAI files and tabular inputs receive full SHA-256 hashes.",
  "- Inference is per-stage and cross-sectional. No stage-by-group interaction significance is created here.",
  "",
  "## Re-run",
  "",
  paste("Rscript", script_path,
        "--config-dir", config_dir, "--rnaseq-root", rnaseq_root, "--output", output_root)
)
writeLines(provenance, file.path(output_root, "provenance.md"))

claim_map <- data.frame(
  claim_id = c("splice_event_trajectory", "splice_locus_support"),
  claim = c("Configured events have stage-specific per-library PSI values with age contrasts oriented as comparison minus reference.",
            "Configured representative loci have read-level junction support in indexed BAM files where Sashimi status is completed."),
  input = c("rMATS MATS tables and run_metadata.tsv", "BAM/BAI and one-row rMATS event tables"),
  script = script_path,
  table = c("tables/target_event_psi_long.tsv;tables/target_event_stage_statistics.tsv", "tables/target_event_read_support.tsv;tables/sashimi_status.tsv"),
  figure = c("plots/*_PSI_trajectory.pdf", "plots/*_sashimi.pdf"),
  test = c("tables/validation.tsv", "tables/sashimi_status.tsv"),
  scope = c("exploratory per-stage event inclusion", "representative locus visualization"),
  status = c("supported", if (all(sashimi_status$status[sashimi_status$status != "not_requested"] == "completed")) "supported" else "inconclusive"),
  stringsAsFactors = FALSE
)
write_tsv(claim_map, file.path(output_root, "claim_result_map.tsv"))
