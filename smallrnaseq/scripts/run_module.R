root <- Sys.getenv("SMALLRNA_ROOT_DIR")
setwd(root)
setup_file <- file.path(root,"rmd","01_setup.R")
module <- as.integer(Sys.getenv("SMALLRNA_MODULE"))
if (module == 1L) {
  source(setup_file, local=globalenv())
} else {
  input_file <- list.files(file.path(root,"rmd"),pattern=sprintf("^%02d_.*[.]Rmd$",module),full.names=TRUE)
  if(length(input_file)!=1L) stop("Missing or ambiguous Rmd module")
  source(setup_file, local=globalenv())
  knitr::opts_chunk$set(error=FALSE)
  knitr::knit(input_file,output=file.path(Sys.getenv("SMALLRNA_RESULT_DIR"),sprintf("module_%02d.md",module)),envir=globalenv(),quiet=TRUE)
}
