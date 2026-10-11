#!/usr/bin/env Rscript
run <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if(length(args) != 4L) stop("需要 step_id、root、rmd_dir、module_file")
  step_id <- as.integer(args[[1]])
  root_dir <- normalizePath(args[[2]], mustWork = TRUE)
  rmd_dir <- normalizePath(args[[3]], mustWork = TRUE)
  input_file <- file.path(rmd_dir, args[[4]])
  Sys.setenv(ATAC_ROOT_DIR = root_dir, ATAC_RMD_DIR = rmd_dir)
  setwd(root_dir)
  source(file.path(rmd_dir, "01_main_setup.R"), local = globalenv())
  if(step_id == 1L) return(invisible(NULL))
  load_atac_context(step_id, globalenv())
  knitr::opts_knit$set(root.dir = root_dir)
  output <- tempfile(sprintf("atac_step_%02d_", step_id), fileext = ".md")
  on.exit(unlink(output), add = TRUE)
  invisible(knitr::knit(input_file, output = output, envir = globalenv(), quiet = TRUE))
}
run()
