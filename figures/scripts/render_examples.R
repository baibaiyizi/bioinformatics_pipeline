#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(root, "scripts", "recipes.R"))
read <- function(name) read.delim(file.path(root, "examples", name), check.names = FALSE)
out <- function(name) file.path(root, "previews", paste0("r_", name), "figure")
plot_ora(read("marker_ora.tsv"), out("ora"))
plot_marker(read("marker_expression.tsv"), out("marker"))
plot_program_matrix(read("program_matrix.tsv"), read("composition.tsv"), out("program_matrix"))
plot_violin(read("gene_effect_violin.tsv"), out("violin"), pfos_palette("de")[c("up","down")] |> setNames(c("Hyper","Hypo")), "RNA log2FC")
receipts <- data.frame(recipe = c("r_ora", "r_marker", "r_program_matrix", "r_violin"), data_kind = "real_PFOS_frozen_source", statistics_recomputed = FALSE)
write.table(receipts, file.path(root, "previews", "r_examples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (recipe in receipts$recipe) {
  inputs <- switch(recipe, r_violin = "gene_effect_violin.tsv", r_ora = "marker_ora.tsv", r_marker = "marker_expression.tsv", r_program_matrix = c("program_matrix.tsv", "composition.tsv"))
  paths <- file.path(root, "examples", inputs)
  sha <- function(path) digest::digest(file = path, algo = "sha256")
  files <- file.path(root, "previews", recipe, c("figure.pdf", "figure.png"))
  record <- list(status = "COMPLETED", recipe = recipe, data_kind = "real_PFOS_frozen_source", statistics_recomputed = FALSE,
    input_sha256 = setNames(lapply(paths, sha), inputs), script_sha256 = sha(file.path(root, "scripts", "recipes.R")),
    runner_sha256 = sha(script), style_sha256 = sha(file.path(root, "scripts", "figure_style.R")),
    profile_sha256 = sha(file.path(root, "config", "figure_style.json")), output_sha256 = setNames(lapply(files, sha), basename(files)))
  jsonlite::write_json(record, file.path(root, "previews", recipe, "figure.json"), auto_unbox = TRUE, pretty = TRUE)
}
cat("COMPLETED 4 native R figure pairs\n")
