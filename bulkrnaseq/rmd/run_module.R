args <- commandArgs(trailingOnly=TRUE)
id <- as.integer(args[[1]])
root_dir <- Sys.getenv('RNASEQ_ROOT_DIR'); rmd_dir <- Sys.getenv('RNASEQ_RMD_DIR')
setwd(root_dir)
source(file.path(rmd_dir,'01_main_setup.R'),local=globalenv())
if (id != 1L) {
 if (id %in% c(9L,10L,13L,14L,15L,16L,17L)) load_post_context(globalenv())
 if (!requireNamespace('knitr',quietly=TRUE)) stop('MISSING R package knitr')
 knitr::opts_knit$set(root.dir=root_dir); knitr::opts_chunk$set(error=FALSE)
 target <- tempfile(sprintf('module_%02d_',id),fileext='.md')
 invisible(knitr::knit(file.path(rmd_dir,args[[2]]),output=target,envir=globalenv(),quiet=TRUE))
 unlink(target)
}
