options(stringsAsFactors = FALSE)
options(timeout = 600)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(data.table)
})

first_nonempty <- function(...) {
  for (value in list(...)) {
    if (is.null(value) || length(value) == 0) next
    value <- trimws(as.character(value[[1]]))
    if (!is.na(value) && nzchar(value)) return(value)
  }
  ""
}

resolve_project_root <- function() {
  root_candidate <- first_nonempty(
    Sys.getenv("SPATIAL_ROOT_DIR", unset = ""),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
  )
  normalizePath(root_candidate, winslash = "/", mustWork = FALSE)
}

root_dir <- resolve_project_root()
config_file <- file.path(root_dir, "config", "spatial_config.yaml")
config <- if (file.exists(config_file)) yaml::read_yaml(config_file) else list()

result_dir <- file.path(root_dir, "result")
cache_dir <- file.path(root_dir, ".cache", "rmd")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

pandoraseq_palette <- c(
  "#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B",
  "#B39FDB", "#FDFFB6", "#D4D9DE", "#5A666F"
)

theme_spatial <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

module_dir <- function(module) {
  path <- file.path(result_dir, module)
  dir.create(file.path(path, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "plots"), recursive = TRUE, showWarnings = FALSE)
  path
}

table_path <- function(module, filename) {
  file.path(module_dir(module), "tables", filename)
}

plot_path <- function(module, filename) {
  file.path(module_dir(module), "plots", filename)
}

read_optional_tsv <- function(path) {
  if (!file.exists(path)) return(tibble())
  data.table::fread(path, sep = "\t", data.table = FALSE) |> as_tibble()
}

save_plot_versions <- function(plot, module, name, width = 7, height = 5) {
  pdf_file <- plot_path(module, paste0(name, ".pdf"))
  png_file <- plot_path(module, paste0(name, ".png"))
  ggsave(pdf_file, plot, width = width, height = height)
  ggsave(png_file, plot, width = width, height = height, dpi = 300)
  invisible(list(pdf = pdf_file, png = png_file))
}

write_review_table <- function(x, module, filename) {
  out <- table_path(module, filename)
  readr::write_csv(x, out)
  invisible(out)
}
