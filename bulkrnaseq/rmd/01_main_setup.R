install_and_load <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

try_install_and_load <- function(pkgs, bioc = FALSE) {
  loaded <- character(0)
  missing <- character(0)
  messages <- character(0)
  for (pkg in pkgs) {
    ok <- tryCatch({
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (bioc) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
          }
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
          install.packages(pkg, repos = "https://cloud.r-project.org")
        }
      }
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      TRUE
    }, error = function(e) {
      messages <<- c(messages, paste0(pkg, ": ", conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) {
      loaded <- c(loaded, pkg)
    } else {
      missing <- c(missing, pkg)
    }
  }
  list(ok = length(missing) == 0, loaded = loaded, missing = missing, messages = messages)
}

first_nonempty <- function(...) {
  for (value in list(...)) {
    if (is.null(value) || length(value) == 0) {
      next
    }
    value_chr <- trimws(as.character(value[[1]]))
    if (!is.na(value_chr) && nzchar(value_chr)) {
      return(value_chr)
    }
  }
  ""
}

normalize_rnaseq_species <- function(value) {
  species_key <- tolower(first_nonempty(value))
  if (species_key %in% c("", "mouse", "mus musculus", "mmu")) {
    return("mouse")
  }
  if (species_key %in% c("human", "homo sapiens", "hsa")) {
    return("human")
  }
  stop("不支持的 RNASEQ_SPECIES: ", value, "。当前仅支持 mouse 或 human。")
}

resolve_project_root <- function() {
  root_candidate <- first_nonempty(
    get0("root_dir", ifnotfound = "", inherits = TRUE),
    Sys.getenv("RNASEQ_ROOT_DIR", unset = ""),
    get0("working_dir", ifnotfound = "", inherits = TRUE),
    getwd()
  )
  normalizePath(root_candidate, winslash = "/", mustWork = FALSE)
}

resolve_orgdb_object <- function(object_name, package_name) {
  if (!exists(object_name, inherits = TRUE)) {
    stop("未找到注释库对象: ", object_name, "。请确认已安装并加载包 ", package_name, "。")
  }
  get(object_name, inherits = TRUE)
}

working_dir <- resolve_project_root()
if (!dir.exists(working_dir)) {
  stop("工作目录不存在: ", working_dir)
}
setwd(working_dir)
cat("工作目录设置为:", getwd(), "\n")

rnaseq_species <- normalize_rnaseq_species(first_nonempty(
  get0("rnaseq_species", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_SPECIES", unset = ""),
  "mouse"
))

species_defaults <- switch(
  rnaseq_species,
  mouse = list(
    orgdb_package = "org.Mm.eg.db",
    orgdb_object_name = "org.Mm.eg.db",
    kegg_organism_code = "mmu",
    reactome_organism = "mouse",
    msigdb_species = "Mus musculus"
  ),
  human = list(
    orgdb_package = "org.Hs.eg.db",
    orgdb_object_name = "org.Hs.eg.db",
    kegg_organism_code = "hsa",
    reactome_organism = "human",
    msigdb_species = "Homo sapiens"
  )
)

orgdb_package <- first_nonempty(
  get0("orgdb_package", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_ORGDB_PACKAGE", unset = ""),
  species_defaults$orgdb_package
)
orgdb_object_name <- first_nonempty(
  get0("orgdb_object_name", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_ORGDB_OBJECT", unset = ""),
  species_defaults$orgdb_object_name
)
kegg_organism_code <- first_nonempty(
  get0("kegg_organism_code", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_KEGG_ORGANISM", unset = ""),
  species_defaults$kegg_organism_code
)
reactome_organism <- first_nonempty(
  get0("reactome_organism", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_REACTOME_ORGANISM", unset = ""),
  species_defaults$reactome_organism
)
msigdb_species <- first_nonempty(
  get0("msigdb_species", ifnotfound = "", inherits = TRUE),
  Sys.getenv("RNASEQ_MSIGDB_SPECIES", unset = ""),
  species_defaults$msigdb_species
)

cran_pkgs <- c(
  "dplyr", "tibble", "tidyr", "readr", "stringr", "purrr",
  "ggplot2", "ggrepel", "ggdendro", "reshape2", "scales",
  "RColorBrewer", "patchwork", "pheatmap", "WGCNA",
  "GOplot", "UpSetR", "Cairo", "circlize", "msigdbr", "remotes",
  "jsonlite"
)
bioc_pkgs <- unique(c(
  "clusterProfiler", "enrichplot", orgdb_package,
  "AnnotationDbi", "edgeR", "ComplexHeatmap", "ReactomePA", "fgsea",
  "GSVA", "GSEABase", "DESeq2", "limma"
))

install_and_load(cran_pkgs, bioc = FALSE)
install_and_load(bioc_pkgs, bioc = TRUE)

ensure_github_package <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github(repo, upgrade = "never", dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

annotation_orgdb <- resolve_orgdb_object(orgdb_object_name, orgdb_package)

result_root_dir <- file.path(working_dir, "result")
result_01_dir <- file.path(result_root_dir, "01_setup")
result_02_dir <- file.path(result_root_dir, "02_selection")
result_03_dir <- file.path(result_root_dir, "03_qc")
result_04_dir <- file.path(result_root_dir, "04_pca")
result_05_dir <- file.path(result_root_dir, "05_de")
result_06_dir <- file.path(result_root_dir, "06_enrichment")
result_07_dir <- file.path(result_root_dir, "07_go")
result_08_dir <- file.path(result_root_dir, "08_gsea")
result_09_dir <- file.path(result_root_dir, "09_gsva")
result_10_dir <- file.path(result_root_dir, "10_tf_regulator")
result_11_dir <- file.path(result_root_dir, "11_ppi")
result_12_dir <- file.path(result_root_dir, "12_wgcna")
result_13_dir <- file.path(result_root_dir, "13_apa")
result_14_dir <- file.path(result_root_dir, "14_der")
result_15_dir <- file.path(result_root_dir, "15_splicing")
result_16_dir <- file.path(result_root_dir, "16_immune")
result_17_dir <- file.path(result_root_dir, "17_pathway_activity")
result_18_dir <- file.path(result_root_dir, "18_circle")

cache_root_dir <- file.path(working_dir, ".cache", "rmd")
cache_06_tables_path <- file.path(cache_root_dir, "06_main_enrichment_tables.rds")
cache_06_objects_path <- file.path(cache_root_dir, "06_main_enrichment_objects.rds")
cache_09_gsva_state_path <- file.path(cache_root_dir, "09_gsva_state.rds")
cache_10_tf_state_path <- file.path(cache_root_dir, "10_tf_regulator_state.rds")
cache_11_ppi_state_path <- file.path(cache_root_dir, "11_ppi_state.rds")
cache_12_wgcna_path <- file.path(cache_root_dir, "12_wgcna_state.rds")
cache_15_splicing_state_path <- file.path(cache_root_dir, "15_splicing_state.rds")
cache_16_immune_state_path <- file.path(cache_root_dir, "16_immune_state.rds")

# 兼容旧 step 变量命名，统一映射到新目录接口。
step1_input_dir <- result_02_dir
step1_overview_dir <- result_02_dir
step2_qc_dir <- result_03_dir
step3_pca_dir <- result_04_dir
step4_de_dir <- result_05_dir
step5_function_dir <- result_06_dir
step6_go_dir <- result_07_dir
step7_gsea_dir <- result_08_dir
step8_wgcna_dir <- result_12_dir
step9_circle_dir <- result_18_dir
step10_root_dir <- result_15_dir
step11_dir <- result_13_dir
step12_dir <- result_16_dir
step13_dir <- result_14_dir
step14_dir <- result_17_dir
step15_dir <- result_10_dir
step16_dir <- result_09_dir
step17_dir <- result_11_dir

ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  dirs <- unique(dirs[!is.na(dirs) & nzchar(dirs)])
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  if (is.atomic(x) || is.list(x)) {
    all_missing <- suppressWarnings(tryCatch(all(is.na(x)), error = function(e) FALSE))
    if (isTRUE(all_missing)) {
      return(y)
    }
  }
  x
}

DIRECTION_LEVELS <- c("all", "up", "down")
DIRECTION_LABELS <- c(all = "全部变化", up = "仅上调", down = "仅下调")
PLOT_FILE_EXTENSIONS <- c("pdf", "png", "jpg", "jpeg", "svg", "tif", "tiff", "bmp")

RESULT_LAYOUT_REGISTRY <- list(
  "01_setup" = list(type = "standard"),
  "02_selection" = list(type = "standard"),
  "03_qc" = list(type = "standard"),
  "04_pca" = list(type = "standard"),
  "05_de" = list(type = "standard"),
  "06_enrichment" = list(type = "direction", keys = c(all = "01_all", up = "02_up", down = "03_down"), default_key = "all"),
  "07_go" = list(type = "direction", keys = c(all = "01_all", up = "02_up", down = "03_down"), default_key = "all"),
  "08_gsea" = list(type = "direction", keys = c(all = "01_all", up = "02_up", down = "03_down"), default_key = "all"),
  "09_gsva" = list(type = "standard"),
  "10_tf_regulator" = list(type = "standard"),
  "11_ppi" = list(type = "standard"),
  "12_wgcna" = list(type = "standard"),
  "13_apa" = list(type = "standard"),
  "14_der" = list(type = "standard"),
  "15_splicing" = list(
    type = "module",
    keys = c(
      event_level = "01_event_level",
      gene_level = "02_gene_level",
      visual_validation = "03_visual_validation",
      isoform_switch = "04_isoform_switch",
      integration = "05_integration",
      sashimi_inputs = "06_sashimi_inputs"
    ),
    default_key = "event_level"
  ),
  "16_immune" = list(type = "standard"),
  "17_pathway_activity" = list(type = "standard"),
  "18_circle" = list(type = "standard")
)

MODULE_CACHE_REGISTRY <- list(
  "06_enrichment" = c(cache_06_tables_path, cache_06_objects_path),
  "09_gsva" = cache_09_gsva_state_path,
  "10_tf_regulator" = cache_10_tf_state_path,
  "11_ppi" = cache_11_ppi_state_path,
  "12_wgcna" = cache_12_wgcna_path,
  "15_splicing" = cache_15_splicing_state_path,
  "16_immune" = cache_16_immune_state_path
)

compose_path <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) {
    return("")
  }
  do.call(file.path, as.list(parts))
}

layout_entry_for_step <- function(step_dir) {
  RESULT_LAYOUT_REGISTRY[[step_dir]] %||% list(type = "standard", keys = character(0), default_key = NULL)
}

match_layout_key <- function(step_dir, key = NULL) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(NULL)
  }

  key_chr <- trimws(first_nonempty(key, entry$default_key, names(entry$keys)[1]))
  if (!nzchar(key_chr)) {
    return(NULL)
  }
  if (key_chr %in% names(entry$keys)) {
    return(unname(entry$keys[[key_chr]]))
  }
  if (key_chr %in% unname(entry$keys)) {
    return(key_chr)
  }
  NULL
}

layout_key_dir <- function(step_dir, key = NULL) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(NULL)
  }
  key_dir <- match_layout_key(step_dir, key)
  if (is.null(key_dir) || !nzchar(key_dir)) {
    stop("未识别的布局 key: ", key, " (step=", step_dir, ")")
  }
  key_dir
}

layout_key_name <- function(step_dir, key = NULL) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(NULL)
  }
  key_dir <- layout_key_dir(step_dir, key)
  key_name <- names(entry$keys)[match(key_dir, unname(entry$keys))]
  if (length(key_name) == 0 || is.na(key_name[[1]]) || !nzchar(key_name[[1]])) {
    key_dir
  } else {
    key_name[[1]]
  }
}

step_result_dir <- function(step_dir, key = NULL) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(file.path(result_root_dir, step_dir))
  }
  file.path(result_root_dir, step_dir, layout_key_dir(step_dir, key))
}

step_plot_dir <- function(step_dir, key = NULL) {
  file.path(step_result_dir(step_dir, key), "plots")
}

step_table_dir <- function(step_dir, key = NULL) {
  file.path(step_result_dir(step_dir, key), "tables")
}

step_result_path <- function(step_dir, kind = c("table", "plot"), key = NULL, ...) {
  kind <- match.arg(kind)
  sub_dir <- if (identical(kind, "plot")) "plots" else "tables"
  compose_path(step_result_dir(step_dir, key), sub_dir, ...)
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

is_plot_file_path <- function(path) {
  ext <- tolower(tools::file_ext(basename(as.character(path[[1]]))))
  ext %in% PLOT_FILE_EXTENSIONS
}

move_path_safe <- function(src, dst) {
  src_norm <- normalizePath(src, winslash = "/", mustWork = FALSE)
  dst_norm <- normalizePath(dst, winslash = "/", mustWork = FALSE)
  if (identical(src_norm, dst_norm)) {
    return(invisible(TRUE))
  }

  if (base::dir.exists(src)) {
    ensure_dirs(dst)
    children <- list.files(src, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    for (child in children) {
      move_path_safe(child, file.path(dst, basename(child)))
    }
    if (base::dir.exists(src) && length(list.files(src, all.files = TRUE, no.. = TRUE)) == 0) {
      unlink(src, recursive = TRUE, force = TRUE)
    }
    return(invisible(TRUE))
  }

  if (!base::file.exists(src) || base::file.exists(dst)) {
    return(invisible(FALSE))
  }

  ensure_parent_dir(dst)
  moved <- suppressWarnings(base::file.rename(src, dst))
  if (!isTRUE(moved)) {
    moved <- base::file.copy(src, dst, overwrite = FALSE)
    if (isTRUE(moved)) {
      unlink(src, force = TRUE)
    }
  }
  invisible(isTRUE(moved))
}

prune_empty_dirs <- function(path) {
  if (!base::dir.exists(path)) {
    return(invisible(NULL))
  }
  children <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  for (child in children) {
    if (base::dir.exists(child)) {
      prune_empty_dirs(child)
    }
  }
  if (base::dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) == 0) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  invisible(NULL)
}

truthy_env <- function(name, default = TRUE) {
  raw_value <- Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false")
  tolower(trimws(raw_value)) %in% c("1", "true", "yes", "y", "on")
}

cleaned_module_result_steps <- new.env(parent = emptyenv())

module_result_root_dir <- function(step_dir) {
  if (!step_dir %in% names(RESULT_LAYOUT_REGISTRY)) {
    stop("未知模块结果目录: ", step_dir)
  }
  file.path(result_root_dir, step_dir)
}

clean_active_module_result_once <- function(step_dir) {
  active_step <- Sys.getenv("RNASEQ_ACTIVE_RESULT_STEP", unset = "")
  if (!nzchar(active_step) || !identical(active_step, step_dir)) {
    return(invisible(FALSE))
  }
  if (isTRUE(cleaned_module_result_steps[[step_dir]])) {
    return(invisible(FALSE))
  }
  cleaned_module_result_steps[[step_dir]] <- TRUE

  if (truthy_env("RNASEQ_CLEAN_MODULE_RESULT", default = TRUE)) {
    target_dir <- module_result_root_dir(step_dir)
    result_root_norm <- normalizePath(result_root_dir, winslash = "/", mustWork = FALSE)
    target_norm <- normalizePath(target_dir, winslash = "/", mustWork = FALSE)
    if (!startsWith(paste0(target_norm, "/"), paste0(result_root_norm, "/"))) {
      stop("拒绝清理 result 目录之外的路径: ", target_dir)
    }
    if (base::dir.exists(target_dir)) {
      message("[cleanup] 删除旧结果目录: ", target_dir)
      unlink(target_dir, recursive = TRUE, force = TRUE)
    }
    ensure_dirs(target_dir)
  }

  if (truthy_env("RNASEQ_CLEAN_MODULE_CACHE", default = TRUE)) {
    cache_paths <- MODULE_CACHE_REGISTRY[[step_dir]] %||% character(0)
    cache_paths <- cache_paths[!is.na(cache_paths) & nzchar(cache_paths)]
    cache_hits <- cache_paths[base::file.exists(cache_paths)]
    if (length(cache_hits) > 0) {
      message("[cleanup] 删除旧模块缓存: ", paste(cache_hits, collapse = ", "))
      unlink(cache_hits, force = TRUE)
    }
  }

  invisible(TRUE)
}

migrate_layout_key_dir <- function(step_dir, key_token) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(invisible(NULL))
  }

  legacy_dir <- file.path(result_root_dir, step_dir, as.character(key_token))
  target_key <- layout_key_dir(step_dir, key_token)
  target_dir <- step_result_dir(step_dir, target_key)
  if (identical(
    normalizePath(legacy_dir, winslash = "/", mustWork = FALSE),
    normalizePath(target_dir, winslash = "/", mustWork = FALSE)
  )) {
    return(invisible(NULL))
  }
  if (!base::dir.exists(legacy_dir)) {
    return(invisible(NULL))
  }

  children <- list.files(legacy_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  for (child in children) {
    child_name <- basename(child)
    if (base::dir.exists(child) && child_name %in% c("plots", "tables")) {
      target_root <- if (identical(child_name, "plots")) step_plot_dir(step_dir, target_key) else step_table_dir(step_dir, target_key)
      move_path_safe(child, target_root)
    } else {
      target_root <- if (is_plot_file_path(child_name)) step_plot_dir(step_dir, target_key) else step_table_dir(step_dir, target_key)
      move_path_safe(child, file.path(target_root, child_name))
    }
  }
  prune_empty_dirs(legacy_dir)
}

migrate_layout_kind_dir <- function(step_dir, kind) {
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    return(invisible(NULL))
  }

  legacy_dir <- file.path(result_root_dir, step_dir, kind)
  if (!base::dir.exists(legacy_dir)) {
    return(invisible(NULL))
  }

  children <- list.files(legacy_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  for (child in children) {
    child_name <- basename(child)
    child_stem <- tools::file_path_sans_ext(child_name)
    target_key <- match_layout_key(step_dir, child_name) %||% match_layout_key(step_dir, child_stem) %||% layout_key_dir(step_dir, NULL)
    target_root <- if (identical(kind, "plots")) step_plot_dir(step_dir, target_key) else step_table_dir(step_dir, target_key)
    if (base::dir.exists(child) && !is.null(match_layout_key(step_dir, child_name))) {
      move_path_safe(child, target_root)
    } else {
      move_path_safe(child, file.path(target_root, child_name))
    }
  }
  prune_empty_dirs(legacy_dir)
}

ensure_module_result_layout <- function(step_dir) {
  clean_active_module_result_once(step_dir)
  entry <- layout_entry_for_step(step_dir)
  if (identical(entry$type, "standard")) {
    ensure_dirs(step_result_dir(step_dir), step_plot_dir(step_dir), step_table_dir(step_dir))
    return(invisible(step_result_dir(step_dir)))
  }

  canonical_keys <- unique(c(names(entry$keys), unname(entry$keys)))
  for (key_token in canonical_keys) {
    ensure_dirs(step_result_dir(step_dir, key_token), step_plot_dir(step_dir, key_token), step_table_dir(step_dir, key_token))
  }
  for (key_token in canonical_keys) {
    migrate_layout_key_dir(step_dir, key_token)
  }
  for (kind in c("plots", "tables")) {
    migrate_layout_kind_dir(step_dir, kind)
  }
  invisible(step_result_dir(step_dir))
}

step5_direction_modules <- layout_entry_for_step("06_enrichment")$keys
step6_direction_modules <- layout_entry_for_step("07_go")$keys
step7_direction_modules <- layout_entry_for_step("08_gsea")$keys
step10_module_names <- layout_entry_for_step("15_splicing")$keys

ensure_direction_dirs <- function(base_dir) {
  ensure_dirs(base_dir, file.path(base_dir, DIRECTION_LEVELS))
}

direction_dir <- function(base_dir, direction) {
  if (!direction %in% DIRECTION_LEVELS) {
    stop("未知方向: ", direction)
  }
  file.path(base_dir, direction)
}

step5_module_dir <- function(direction) step_result_dir("06_enrichment", direction)
step5_plot_dir <- function(direction) step_plot_dir("06_enrichment", direction)
step5_table_dir <- function(direction) step_table_dir("06_enrichment", direction)

step6_module_dir <- function(direction) step_result_dir("07_go", direction)
step6_plot_dir <- function(direction) step_plot_dir("07_go", direction)
step6_table_dir <- function(direction) step_table_dir("07_go", direction)

step7_module_dir <- function(direction) step_result_dir("08_gsea", direction)
step7_plot_dir <- function(direction) step_plot_dir("08_gsea", direction)
step7_table_dir <- function(direction) step_table_dir("08_gsea", direction)

step10_module_dir <- function(module_key) step_result_dir("15_splicing", module_key)
step10_plot_dir <- function(module_key) step_plot_dir("15_splicing", module_key)
step10_table_dir <- function(module_key) step_table_dir("15_splicing", module_key)

step10_event_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "event_level", ...)
step10_event_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "event_level", ...)
step10_gene_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "gene_level", ...)
step10_gene_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "gene_level", ...)
step10_visual_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "visual_validation", ...)
step10_visual_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "visual_validation", ...)
step10_isoform_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "isoform_switch", ...)
step10_isoform_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "isoform_switch", ...)
step10_integr_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "integration", ...)
step10_integr_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "integration", ...)
step10_sashimi_table_path <- function(...) step_result_path("15_splicing", kind = "table", key = "sashimi_inputs", ...)
step10_sashimi_plot_path <- function(...) step_result_path("15_splicing", kind = "plot", key = "sashimi_inputs", ...)

step10_event_dir <- step10_module_dir("event_level")
step10_gene_dir <- step10_module_dir("gene_level")
step10_visual_dir <- step10_module_dir("visual_validation")
step10_isoform_dir <- step10_module_dir("isoform_switch")
step10_integr_dir <- step10_module_dir("integration")
post_tx_sashimi_dir <- step10_table_dir("sashimi_inputs")

direction_label <- function(direction) {
  unname(DIRECTION_LABELS[[direction]] %||% direction)
}

restore_main_de_result <- function(required = TRUE) {
  de_file <- file.path(result_05_dir, "5_1_annotated_DE_results.csv")

  if (!file.exists(de_file)) {
    if (isTRUE(required)) {
      stop(
        "缺少主线差异结果：", de_file,
        "。请先运行模块 05。"
      )
    }
    return(NULL)
  }

  de_tbl <- read.csv(de_file, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c("gene_id", "ensembl_id", "gene_label", "ENTREZID", "log2FoldChange", "pvalue", "direction", "rank_score")
  missing_cols <- setdiff(required_cols, colnames(de_tbl))
  if (length(missing_cols) > 0) {
    stop("DE 结果缺少必要列: ", paste(missing_cols, collapse = ", "), "，文件: ", de_file)
  }

  de_tbl
}

load_post_context <- function(env = parent.frame()) {
  ensure_module_result_layout("01_setup")
  ensure_dirs(cache_root_dir)

  de_result <- restore_main_de_result(required = TRUE)
  sig_de <- de_result %>% filter(direction != "ns")

  rank_file <- file.path(result_05_dir, "5_9_gene_rank_for_GSEA.rnk")
  if (file.exists(rank_file)) {
    rank_tbl <- read.table(rank_file, sep = "\t", stringsAsFactors = FALSE, col.names = c("ENTREZID", "rank_score"))
  } else {
    rank_tbl <- de_result %>%
      filter(!is.na(ENTREZID), !is.na(rank_score)) %>%
      group_by(ENTREZID) %>%
      slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(rank_score))
  }

  if (!file.exists(cache_06_tables_path)) {
    stop(
      "缺少主线富集表格缓存: ", cache_06_tables_path,
      "。请先运行模块 06。"
    )
  }

  enrichment_tables <- readRDS(cache_06_tables_path)
  if (!is.list(enrichment_tables)) {
    stop("富集表格缓存格式错误（应为 list）: ", cache_06_tables_path)
  }

  cache_tbl <- function(name, default = data.frame()) {
    value <- enrichment_tables[[name]]
    if (is.null(value)) {
      return(default)
    }
    value
  }

  go_directional_df <- cache_tbl("go_directional_df")
  kegg_directional_df <- cache_tbl("kegg_directional_df")
  reactome_directional_df <- cache_tbl("reactome_directional_df")
  gsea_go_df <- cache_tbl("gsea_go_df")
  gsea_kegg_df <- cache_tbl("gsea_kegg_df")
  gsea_hallmark_df <- cache_tbl("gsea_hallmark_df")
  gsea_reactome_df <- cache_tbl("gsea_reactome_df")
  gsea_msigdb_c2_df <- cache_tbl("gsea_msigdb_c2_df")
  gsea_msigdb_c5_df <- cache_tbl("gsea_msigdb_c5_df")
  leading_edge_go <- cache_tbl("leading_edge_go")
  leading_edge_kegg <- cache_tbl("leading_edge_kegg")
  leading_edge_hallmark <- cache_tbl("leading_edge_hallmark")
  leading_edge_reactome <- cache_tbl("leading_edge_reactome")
  leading_edge_msigdb_c2 <- cache_tbl("leading_edge_msigdb_c2")
  leading_edge_msigdb_c5 <- cache_tbl("leading_edge_msigdb_c5")
  go_circle_data <- cache_tbl("go_circle_data")
  kegg_circle_data <- cache_tbl("kegg_circle_data")
  gsea_circle_data <- cache_tbl("gsea_circle_data")

  if (!file.exists(cache_06_objects_path)) {
    message("[post_context] 未检测到对象缓存: ", cache_06_objects_path, "，对象型图会在对应模块中显式跳过。")
    enrichment_objects <- list()
  } else {
    enrichment_objects <- readRDS(cache_06_objects_path)
    if (!is.list(enrichment_objects)) {
      stop("富集对象缓存格式错误（应为 list）: ", cache_06_objects_path)
    }
  }

  go_ora_sets <- enrichment_objects[["go_ora_sets"]] %||% list(
    all = list(BP = NULL, MF = NULL, CC = NULL, ALL = NULL),
    up = list(BP = NULL, MF = NULL, CC = NULL, ALL = NULL),
    down = list(BP = NULL, MF = NULL, CC = NULL, ALL = NULL)
  )
  kegg_ora_sets <- enrichment_objects[["kegg_ora_sets"]] %||% list(all = NULL, up = NULL, down = NULL)
  reactome_ora_sets <- enrichment_objects[["reactome_ora_sets"]] %||% list(all = NULL, up = NULL, down = NULL)

  go_all <- go_ora_sets$all$ALL %||% NULL
  kegg_ora <- enrichment_objects[["kegg_ora"]] %||% gsea_kegg_df
  reactome_ora <- enrichment_objects[["reactome_ora"]] %||% data.frame()
  gsea_go <- enrichment_objects[["gsea_go"]] %||% gsea_go_df
  gsea_kegg <- enrichment_objects[["gsea_kegg"]] %||% gsea_kegg_df
  gsea_reactome <- enrichment_objects[["gsea_reactome"]] %||% gsea_reactome_df

  module_df_annotated <- data.frame()
  module_trait_tbl <- data.frame()
  me_df <- data.frame()

  if (file.exists(cache_12_wgcna_path)) {
    wgcna_state <- readRDS(cache_12_wgcna_path)
    if (is.list(wgcna_state)) {
      module_df_annotated <- wgcna_state$module_df_annotated %||% module_df_annotated
      module_trait_tbl <- wgcna_state$module_trait_tbl %||% module_trait_tbl
      me_df <- wgcna_state$me_df %||% me_df
    }
  }

  if (nrow(module_df_annotated) == 0) {
    f_candidates <- c(
      file.path(step8_wgcna_dir, "12_2_WGCNA_Module_Assignment.csv"),
      file.path(step8_wgcna_dir, "9_2_WGCNA_Module_Assignment.csv")
    )
    f <- f_candidates[file.exists(f_candidates)]
    if (length(f) > 0) {
      f <- f[[1]]
      module_df_annotated <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  if (nrow(module_trait_tbl) == 0) {
    f_candidates <- c(
      file.path(step8_wgcna_dir, "12_4_Module_Trait_Correlation.csv"),
      file.path(step8_wgcna_dir, "9_4_Module_Trait_Correlation.csv")
    )
    f <- f_candidates[file.exists(f_candidates)]
    if (length(f) > 0) {
      f <- f[[1]]
      module_trait_tbl <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  if (nrow(me_df) == 0) {
    f_candidates <- c(
      file.path(step8_wgcna_dir, "12_3b_Module_Eigengene_Table.csv"),
      file.path(step8_wgcna_dir, "9_3b_Module_Eigengene_Table.csv")
    )
    f <- f_candidates[file.exists(f_candidates)]
    if (length(f) > 0) {
      f <- f[[1]]
      me_df <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }

  if (nrow(module_df_annotated) == 0 || nrow(module_trait_tbl) == 0 || nrow(me_df) == 0) {
    message(
      "[post_context] 未恢复到有效 WGCNA 输入；后续模块将跳过 WGCNA 相关整合。候选输入：",
      "\n1) ", cache_12_wgcna_path,
      "\n2) ", file.path(step8_wgcna_dir, "12_2_WGCNA_Module_Assignment.csv"),
      " + ", file.path(step8_wgcna_dir, "12_4_Module_Trait_Correlation.csv"),
      " + ", file.path(step8_wgcna_dir, "12_3b_Module_Eigengene_Table.csv")
    )
  }

  context_objects <- list(
    de_result = de_result,
    sig_de = sig_de,
    rank_tbl = rank_tbl,
    enrichment_tables = enrichment_tables,
    enrichment_objects = enrichment_objects,
    go_directional_df = go_directional_df,
    kegg_directional_df = kegg_directional_df,
    reactome_directional_df = reactome_directional_df,
    gsea_go_df = gsea_go_df,
    gsea_kegg_df = gsea_kegg_df,
    gsea_hallmark_df = gsea_hallmark_df,
    gsea_reactome_df = gsea_reactome_df,
    gsea_msigdb_c2_df = gsea_msigdb_c2_df,
    gsea_msigdb_c5_df = gsea_msigdb_c5_df,
    leading_edge_go = leading_edge_go,
    leading_edge_kegg = leading_edge_kegg,
    leading_edge_hallmark = leading_edge_hallmark,
    leading_edge_reactome = leading_edge_reactome,
    leading_edge_msigdb_c2 = leading_edge_msigdb_c2,
    leading_edge_msigdb_c5 = leading_edge_msigdb_c5,
    go_circle_data = go_circle_data,
    kegg_circle_data = kegg_circle_data,
    gsea_circle_data = gsea_circle_data,
    go_ora_sets = go_ora_sets,
    kegg_ora_sets = kegg_ora_sets,
    reactome_ora_sets = reactome_ora_sets,
    go_all = go_all,
    kegg_ora = kegg_ora,
    reactome_ora = reactome_ora,
    gsea_go = gsea_go,
    gsea_kegg = gsea_kegg,
    gsea_reactome = gsea_reactome,
    module_df_annotated = module_df_annotated,
    module_trait_tbl = module_trait_tbl,
    me_df = me_df,
    rmats_file_index = data.frame(),
    rmats_combined = data.frame(),
    rmats_sig = data.frame(),
    high_confidence_events = data.frame(),
    overlap_summary = data.frame(),
    as_gene_priority = data.frame(),
    as_gene_membership = data.frame(),
    top_as_candidates = data.frame(),
    dsfa_event_table = data.frame(),
    dsfa_gene_burden = data.frame(),
    dsfa_sig_events = data.frame(),
    apa_all_events = data.frame(),
    apa_sig_events = data.frame(),
    apa_gene_summary_tbl = data.frame(),
    der_all_regions = data.frame(),
    der_sig_regions = data.frame(),
    der_gene_summary_tbl = data.frame(),
    deu_sig_ids = character(0),
    isoform_sig_ids = character(0),
    dsfa_sig_ids = character(0),
    apa_sig_ids = character(0),
    der_sig_ids = character(0)
  )

  for (nm in names(context_objects)) {
    assign(nm, context_objects[[nm]], envir = env)
  }

  setup_manifest <- data.frame(
    key = c(
      "de_result_file",
      "enrichment_tables_cache",
      "enrichment_objects_cache",
      "wgcna_state_cache",
      "module_df_rows",
      "module_trait_rows",
      "me_rows"
    ),
    value = c(
      file.path(result_05_dir, "5_1_annotated_DE_results.csv"),
      cache_06_tables_path,
      ifelse(file.exists(cache_06_objects_path), cache_06_objects_path, "missing_optional"),
      ifelse(file.exists(cache_12_wgcna_path), cache_12_wgcna_path, "from_result_csv"),
      nrow(module_df_annotated),
      nrow(module_trait_tbl),
      nrow(me_df)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(setup_manifest, file.path(result_01_dir, "01_1_post_context_manifest.csv"), row.names = FALSE)
  invisible(setup_manifest)
}

normalize_group <- function(x) {
  sub("([._-]?[Ll]iver)$", "", x, perl = TRUE)
}

read_reads_info <- function(path) {
  stopifnot(file.exists(path))
  first_line <- readLines(path, n = 1, warn = FALSE)
  has_header <- grepl("^group_stage\tsample\t", first_line) || grepl("^sample\tgroup\t", first_line)
  if (!has_header) {
    stop("reads_info.txt 必须使用当前带表头的格式，请先用当前 RNA-seq 流程生成或整理该文件。")
  }

  info <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!"sample" %in% colnames(info)) {
    stop("reads_info.txt 缺少 sample 列")
  }
  if (!"group_stage" %in% colnames(info)) {
    if (all(c("group", "stage") %in% colnames(info))) {
      info$group_stage <- paste0(info$group, info$stage)
    } else if ("group" %in% colnames(info)) {
      info$group_stage <- info$group
    } else {
      stop("reads_info.txt 缺少 group_stage/group 信息")
    }
  }

  if (!"group" %in% colnames(info)) {
    info$group <- normalize_group(info$group_stage)
  }
  if (!"stage" %in% colnames(info)) {
    info$stage <- ifelse(grepl("liver$", info$group_stage, ignore.case = TRUE), "liver", info$group_stage)
  }
  if (!"raw_fq1" %in% colnames(info) && "fq1" %in% colnames(info)) {
    info$raw_fq1 <- info$fq1
  }
  if (!"raw_fq2" %in% colnames(info) && "fq2" %in% colnames(info)) {
    info$raw_fq2 <- info$fq2
  }

  info$group <- normalize_group(info$group)
  info <- info[!duplicated(info$sample), , drop = FALSE]
  rownames(info) <- info$sample
  info
}

read_matrix <- function(path) {
  mat <- read.table(
    path,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  data.matrix(mat)
}

read_samples_file <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tbl <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(tbl) < 2) {
    return(NULL)
  }
  colnames(tbl)[1:2] <- c("group", "sample")
  tbl <- tbl[nzchar(tbl$group) & nzchar(tbl$sample), 1:2, drop = FALSE]
  tbl
}

safe_scale_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[!is.finite(scaled)] <- 0
  scaled
}

select_top_variable_genes <- function(mat, top_n = 5000) {
  if (is.null(mat) || nrow(mat) == 0) {
    return(character(0))
  }
  gene_var <- apply(mat, 1, var, na.rm = TRUE)
  keep <- is.finite(gene_var) & gene_var > 0
  if (!any(keep)) {
    return(rownames(mat)[seq_len(min(top_n, nrow(mat)))])
  }
  ordered <- names(sort(gene_var[keep], decreasing = TRUE))
  ordered[seq_len(min(top_n, length(ordered)))]
}

split_numeric_string <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") {
    return(numeric(0))
  }
  vals <- trimws(unlist(strsplit(as.character(x), ",")))
  suppressWarnings(vals_num <- as.numeric(vals))
  vals_num[is.finite(vals_num)]
}

mean_numeric_string <- function(x) {
  vals <- split_numeric_string(x)
  if (length(vals) == 0) {
    return(NA_real_)
  }
  mean(vals)
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep]))
}

first_non_missing_chr <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

first_non_missing_num <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[[1]]
}

PRIMARY_CHROM_REGEX <- Sys.getenv("PRIMARY_CHROM_REGEX", "^(chr)?([1-9]|1[0-9]|X|Y)$")

is_primary_chromosome <- function(x) {
  grepl(PRIMARY_CHROM_REGEX, as.character(x), ignore.case = TRUE)
}

filter_primary_chr_df <- function(df, column_candidates = c("chr", "Chrom", "seqnames", "seqname", "Chr", "chrom", "chromosome")) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(df)
  }
  hit <- intersect(column_candidates, colnames(df))
  if (length(hit) == 0) {
    return(df)
  }
  keep <- is_primary_chromosome(df[[hit[[1]]]])
  keep[is.na(keep)] <- FALSE
  df[keep, , drop = FALSE]
}

find_rmats_file <- function(base_dir, event_type, count_mode) {
  candidates <- c(
    file.path(base_dir, paste0(event_type, ".MATS.", count_mode, ".txt")),
    file.path(base_dir, paste0("fromGTF.", event_type, ".MATS.", count_mode, ".txt"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    return(NA_character_)
  }
  found[[1]]
}

pdf_device_fn <- function() {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF
  } else {
    grDevices::pdf
  }
}

base_pdf_device_fn <- function(path, width, height) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE, bg = "#FFFFFF")
}

FIG_PNG_DPI <- suppressWarnings(as.numeric(Sys.getenv("FIG_PNG_DPI", "300")))
if (!is.finite(FIG_PNG_DPI) || FIG_PNG_DPI <= 0) {
  FIG_PNG_DPI <- 300
}

png_path_from_pdf <- function(path) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    sub("\\.pdf$", ".png", path, ignore.case = TRUE)
  } else {
    paste0(path, ".png")
  }
}

save_png_device <- function(path, width, height, expr) {
  ensure_parent_dir(path)
  tryCatch({
    grDevices::png(path, width = width, height = height, units = "in", res = FIG_PNG_DPI, type = "cairo", bg = "#FFFFFF")
  }, error = function(e) {
    grDevices::png(path, width = width, height = height, units = "in", res = FIG_PNG_DPI, bg = "#FFFFFF")
  })
  force(expr)
  grDevices::dev.off()
}

normalize_project_output_path <- function(path) {
  if (!is.character(path) || length(path) == 0 || is.na(path[[1]]) || !nzchar(path[[1]])) {
    return(path)
  }
  path_chr <- path[[1]]
  if (grepl("^result(/|$)", path_chr)) {
    return(normalizePath(file.path(working_dir, path_chr), winslash = "/", mustWork = FALSE))
  }
  if (grepl("^\\.cache(/|$)", path_chr)) {
    return(normalizePath(file.path(working_dir, path_chr), winslash = "/", mustWork = FALSE))
  }
  normalizePath(path_chr, winslash = "/", mustWork = FALSE)
}

parse_structured_result_path <- function(step_dir, remainder) {
  file_parts <- remainder[nzchar(remainder)]
  explicit_kind <- NULL
  explicit_key <- NULL

  if (length(file_parts) > 0 && file_parts[[1]] %in% c("plots", "tables")) {
    explicit_kind <- file_parts[[1]]
    file_parts <- file_parts[-1]
    if (length(file_parts) > 0) {
      explicit_key <- match_layout_key(step_dir, file_parts[[1]])
      if (!is.null(explicit_key)) {
        file_parts <- file_parts[-1]
      }
    }
  } else {
    if (length(file_parts) > 0) {
      explicit_key <- match_layout_key(step_dir, file_parts[[1]])
      if (!is.null(explicit_key)) {
        file_parts <- file_parts[-1]
      }
    }
    if (length(file_parts) > 0 && file_parts[[1]] %in% c("plots", "tables")) {
      explicit_kind <- file_parts[[1]]
      file_parts <- file_parts[-1]
    }
  }

  list(kind = explicit_kind, key = explicit_key, file_parts = file_parts)
}

build_result_path_candidates <- function(path, preferred = c("auto", "table", "plot", "cache")) {
  preferred <- match.arg(preferred)
  normalized <- normalize_project_output_path(path)
  if (!is.character(normalized) || length(normalized) == 0 || is.na(normalized[[1]]) || !nzchar(normalized[[1]])) {
    return(normalized)
  }

  result_root_norm <- normalizePath(result_root_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(normalized, result_root_norm)) {
    return(normalized)
  }

  rel_path <- sub(paste0("^", result_root_norm, "/?"), "", normalized)
  rel_parts <- strsplit(rel_path, "/", fixed = TRUE)[[1]]
  rel_parts <- rel_parts[nzchar(rel_parts)]
  if (length(rel_parts) < 2) {
    return(normalized)
  }

  step_dir <- rel_parts[[1]]
  layout_entry <- layout_entry_for_step(step_dir)
  remainder <- rel_parts[-1]

  resolve_kind_dir <- function(explicit_kind, file_parts) {
    inferred_kind <- if (length(file_parts) > 0 && is_plot_file_path(file_parts[[length(file_parts)]])) {
      "plots"
    } else {
      "tables"
    }
    switch(
      preferred,
      plot = "plots",
      table = "tables",
      cache = NA_character_,
      auto = explicit_kind %||% inferred_kind
    )
  }

  if (identical(layout_entry$type, "standard")) {
    explicit_kind <- if (length(remainder) > 0 && remainder[[1]] %in% c("plots", "tables")) remainder[[1]] else NA_character_
    file_parts <- if (!is.na(explicit_kind)) remainder[-1] else remainder
    if (length(file_parts) == 0) {
      return(normalized)
    }

    kind_dir <- resolve_kind_dir(explicit_kind, file_parts)
    if (is.na(kind_dir)) {
      return(normalized)
    }

    unique(c(
      compose_path(result_root_norm, step_dir, kind_dir, file_parts),
      if (!is.na(explicit_kind)) compose_path(result_root_norm, step_dir, explicit_kind, file_parts) else character(0),
      compose_path(result_root_norm, step_dir, file_parts),
      normalized
    ))
  } else {
    parsed <- parse_structured_result_path(step_dir, remainder)
    if (length(parsed$file_parts) == 0) {
      return(normalized)
    }

    key_dir <- layout_key_dir(step_dir, parsed$key)
    key_name <- layout_key_name(step_dir, key_dir)
    kind_dir <- resolve_kind_dir(parsed$kind, parsed$file_parts)
    if (is.na(kind_dir)) {
      return(normalized)
    }

    unique(c(
      compose_path(result_root_norm, step_dir, key_dir, kind_dir, parsed$file_parts),
      compose_path(result_root_norm, step_dir, key_name, kind_dir, parsed$file_parts),
      compose_path(result_root_norm, step_dir, key_dir, parsed$file_parts),
      compose_path(result_root_norm, step_dir, key_name, parsed$file_parts),
      compose_path(result_root_norm, step_dir, kind_dir, key_dir, parsed$file_parts),
      compose_path(result_root_norm, step_dir, kind_dir, key_name, parsed$file_parts),
      compose_path(result_root_norm, step_dir, parsed$file_parts),
      normalized
    ))
  }
}

route_output_path <- function(path, preferred = c("auto", "table", "plot", "cache")) {
  preferred <- match.arg(preferred)
  normalized <- normalize_project_output_path(path)
  if (!is.character(normalized) || length(normalized) == 0 || is.na(normalized[[1]]) || !nzchar(normalized[[1]])) {
    return(normalized)
  }
  if (preferred == "cache") {
    return(normalized)
  }
  build_result_path_candidates(normalized, preferred = preferred)[[1]]
}

resolve_existing_path <- function(path, preferred = c("auto", "table", "plot", "cache")) {
  preferred <- match.arg(preferred)
  normalized <- normalize_project_output_path(path)
  if (!is.character(normalized) || length(normalized) == 0 || is.na(normalized[[1]]) || !nzchar(normalized[[1]])) {
    return(normalized)
  }
  if (preferred == "cache") {
    return(normalized)
  }
  candidates <- build_result_path_candidates(normalized, preferred = preferred)
  existing <- candidates[base::file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  normalized
}

sanitize_table_for_export <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  for (nm in colnames(x)) {
    if (is.list(x[[nm]])) {
      x[[nm]] <- vapply(x[[nm]], function(one_item) {
        if (is.null(one_item) || length(one_item) == 0 || all(is.na(one_item))) {
          return(NA_character_)
        }
        paste(as.character(one_item), collapse = ";")
      }, character(1))
    }
  }
  x
}

write_csv_safe <- function(x, path, row.names = FALSE, ...) {
  path <- route_output_path(path, preferred = "table")
  ensure_parent_dir(path)
  utils::write.csv(sanitize_table_for_export(x), path, row.names = row.names, ...)
}

plot_audit_records <- list()

plot_audit_path_for_primary <- function(primary) {
  primary <- normalize_project_output_path(primary)
  result_root_norm <- normalizePath(result_root_dir, winslash = "/", mustWork = FALSE)
  if (!startsWith(primary, result_root_norm)) {
    return(file.path(result_root_dir, "00_audit", "tables", "missing_output_audit.csv"))
  }
  rel_path <- sub(paste0("^", result_root_norm, "/?"), "", primary)
  rel_parts <- strsplit(rel_path, "/", fixed = TRUE)[[1]]
  plot_idx <- match("plots", rel_parts)
  if (!is.na(plot_idx) && plot_idx > 1) {
    audit_dir <- compose_path(result_root_norm, rel_parts[seq_len(plot_idx - 1)], "tables")
    return(file.path(audit_dir, "plot_audit.csv"))
  }
  file.path(result_root_dir, "00_audit", "tables", "missing_output_audit.csv")
}

append_plot_audit_record <- function(path, record) {
  ensure_parent_dir(path)
  existing <- if (base::file.exists(path)) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
  } else {
    data.frame()
  }
  if (nrow(existing) > 0 && "pdf_path" %in% colnames(existing)) {
    existing <- existing[existing$pdf_path != record$pdf_path[[1]], , drop = FALSE]
  }
  if (ncol(existing) > 0) {
    existing[] <- lapply(existing, as.character)
    if ("reason" %in% colnames(existing)) {
      existing$reason[is.na(existing$reason) | existing$reason == "NA"] <- ""
    }
  }
  record[] <- lapply(record, as.character)
  if ("reason" %in% colnames(record)) {
    record$reason[is.na(record$reason) | record$reason == "NA"] <- ""
  }
  utils::write.csv(bind_rows(existing, record), path, row.names = FALSE)
}

record_plot_status <- function(primary, status = "generated", reason = "", input_rows = NA_integer_, module_id = NULL) {
  primary <- route_output_path(primary, preferred = "plot")
  png_path <- png_path_from_pdf(primary)
  result_root_norm <- normalizePath(result_root_dir, winslash = "/", mustWork = FALSE)
  rel_pdf <- if (startsWith(primary, result_root_norm)) sub(paste0("^", result_root_norm, "/?"), "", primary) else primary
  rel_png <- if (startsWith(png_path, result_root_norm)) sub(paste0("^", result_root_norm, "/?"), "", png_path) else png_path
  rel_parts <- strsplit(rel_pdf, "/", fixed = TRUE)[[1]]
  module_id <- module_id %||% rel_parts[[1]] %||% NA_character_

  record <- data.frame(
    module_id = module_id,
    plot_id = tools::file_path_sans_ext(basename(primary)),
    pdf_path = rel_pdf,
    png_path = rel_png,
    status = status,
    reason = reason,
    input_rows = input_rows,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  plot_audit_records[[length(plot_audit_records) + 1]] <<- record
  append_plot_audit_record(plot_audit_path_for_primary(primary), record)
  append_plot_audit_record(file.path(result_root_dir, "00_audit", "tables", "missing_output_audit.csv"), record)
  invisible(record)
}

flush_plot_audit <- function(path = file.path(result_root_dir, "00_audit", "tables", "missing_output_audit.csv")) {
  if (length(plot_audit_records) == 0) {
    return(invisible(data.frame()))
  }
  audit_tbl <- bind_rows(plot_audit_records)
  ensure_parent_dir(path)
  utils::write.csv(audit_tbl, path, row.names = FALSE)
  invisible(audit_tbl)
}

polish_ggplot_for_export <- function(plot) {
  if (!inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
    return(plot)
  }
  tryCatch(
    plot +
      theme(
        plot.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
        panel.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
        legend.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
        legend.key = element_rect(fill = viz_pal$neutral["paper"], color = NA),
        plot.margin = grid::unit(viz_style$plot_margin, "pt")
      ),
    error = function(e) plot
  )
}

save_plot_versions <- function(plot, primary, width = 8, height = 6) {
  primary <- route_output_path(primary, preferred = "plot")
  ensure_parent_dir(primary)
  plot <- polish_ggplot_for_export(plot)
  ggsave(
    primary,
    plot = plot,
    width = width,
    height = height,
    device = pdf_device_fn(),
    bg = unname(viz_pal$neutral["paper"]),
    limitsize = FALSE
  )
  ggsave(
    png_path_from_pdf(primary),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = FIG_PNG_DPI,
    device = "png",
    bg = unname(viz_pal$neutral["paper"]),
    limitsize = FALSE
  )
  record_plot_status(primary, status = "generated", input_rows = NA_integer_)
}

save_placeholder_plot <- function(primary, title = "Plot skipped", reason = "No plottable data", width = 8, height = 6, input_rows = NA_integer_) {
  label <- paste(strwrap(paste(title, reason, sep = "\n"), width = 72), collapse = "\n")
  placeholder <- ggplot() +
    annotate("label", x = 0, y = 0, label = label, size = 4.2, label.size = 0.25, fill = viz_pal$neutral["light"], color = viz_pal$neutral["dark"]) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void()
  save_plot_versions(placeholder, primary, width = width, height = height)
  record_plot_status(primary, status = "placeholder", reason = reason, input_rows = input_rows)
  invisible(NULL)
}

save_base_pdf_versions <- function(primary, width = 8, height = 6, expr) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  primary <- route_output_path(primary, preferred = "plot")
  ensure_parent_dir(primary)
  base_pdf_device_fn(primary, width = width, height = height)
  eval(expr_sub, envir = expr_env)
  dev.off()
  png_primary <- png_path_from_pdf(primary)
  ensure_parent_dir(png_primary)
  tryCatch({
    grDevices::png(png_primary, width = width, height = height, units = "in", res = FIG_PNG_DPI, type = "cairo", bg = "#FFFFFF")
  }, error = function(e) {
    grDevices::png(png_primary, width = width, height = height, units = "in", res = FIG_PNG_DPI, bg = "#FFFFFF")
  })
  eval(expr_sub, envir = expr_env)
  dev.off()
  record_plot_status(primary, status = "generated", input_rows = NA_integer_)
}

save_pheatmap_versions <- function(primary, width = 8, height = 6, expr) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  primary <- route_output_path(primary, preferred = "plot")

  draw_once <- function() {
    ph <- eval(expr_sub, envir = expr_env)
    if (is.list(ph) && !is.null(ph$gtable)) {
      grid::grid.newpage()
      grid::grid.draw(ph$gtable)
    } else if (inherits(ph, "gtable")) {
      grid::grid.newpage()
      grid::grid.draw(ph)
    }
    invisible(ph)
  }

  ensure_parent_dir(primary)
  base_pdf_device_fn(primary, width = width, height = height)
  draw_once()
  grDevices::dev.off()

  png_primary <- png_path_from_pdf(primary)
  ensure_parent_dir(png_primary)
  tryCatch({
    grDevices::png(png_primary, width = width, height = height, units = "in", res = FIG_PNG_DPI, type = "cairo", bg = "#FFFFFF")
  }, error = function(e) {
    grDevices::png(png_primary, width = width, height = height, units = "in", res = FIG_PNG_DPI, bg = "#FFFFFF")
  })
  draw_once()
  grDevices::dev.off()
  record_plot_status(primary, status = "generated", input_rows = NA_integer_)
}

write.csv <- function(x, file = "", ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    file <- route_output_path(file, preferred = "table")
    ensure_parent_dir(file)
  }
  utils::write.csv(x, file = file, ...)
}

write.table <- function(x, file = "", ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    file <- route_output_path(file, preferred = "table")
    ensure_parent_dir(file)
  }
  utils::write.table(x, file = file, ...)
}

saveRDS <- function(object, file = "", ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    is_cache_file <- grepl("(^|/)(\\.cache)(/|$)", file)
    file <- route_output_path(file, preferred = if (is_cache_file) "cache" else "table")
    ensure_parent_dir(file)
  }
  base::saveRDS(object = object, file = file, ...)
}

read.csv <- function(file, ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    file <- resolve_existing_path(file, preferred = "table")
  }
  utils::read.csv(file = file, ...)
}

read.table <- function(file, ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    file <- resolve_existing_path(file, preferred = "table")
  }
  utils::read.table(file = file, ...)
}

readRDS <- function(file, ...) {
  if (is.character(file) && length(file) == 1 && nzchar(file)) {
    is_cache_file <- grepl("(^|/)(\\.cache)(/|$)", file)
    file <- resolve_existing_path(file, preferred = if (is_cache_file) "cache" else "table")
  }
  base::readRDS(file = file, ...)
}

file.exists <- function(...) {
  paths <- unlist(list(...), use.names = FALSE)
  vapply(paths, function(one_path) {
    if (!is.character(one_path) || length(one_path) == 0 || is.na(one_path[[1]]) || !nzchar(one_path[[1]])) {
      return(FALSE)
    }
    base::file.exists(resolve_existing_path(one_path))
  }, logical(1))
}

convert_pdf_to_png <- function(pdf_path, png_path = png_path_from_pdf(pdf_path), dpi = FIG_PNG_DPI) {
  ensure_parent_dir(png_path)
  py_bin <- Sys.which("python3")
  if (!nzchar(py_bin)) {
    py_bin <- Sys.which("python")
  }
  if (!nzchar(py_bin)) {
    return(invisible(FALSE))
  }
  py_code <- paste(
    "import fitz, sys",
    "pdf_path, png_path, dpi = sys.argv[1], sys.argv[2], int(float(sys.argv[3]))",
    "doc = fitz.open(pdf_path)",
    "page = doc.load_page(0)",
    "zoom = dpi / 72.0",
    "pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)",
    "pix.save(png_path)",
    sep = "\n"
  )
  status <- tryCatch(
    system2(py_bin, c("-c", shQuote(py_code), shQuote(pdf_path), shQuote(png_path), as.character(dpi)), stdout = FALSE, stderr = FALSE),
    error = function(e) 1
  )
  invisible(status == 0 && file.exists(png_path))
}

make_result_df <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  as.data.frame(x)
}

extract_enrich_df <- function(x, direction_label = "all", source_label = NULL) {
  df <- make_result_df(x)
  if (nrow(df) == 0) {
    return(df)
  }
  df$direction <- direction_label
  if (!is.null(source_label)) {
    df$source <- source_label
  }
  df
}

make_directional_gene_sets <- function(de_tbl, id_col = "ENTREZID") {
  list(
    all = unique(na.omit(de_tbl[[id_col]][de_tbl$direction != "ns"])),
    up = unique(na.omit(de_tbl[[id_col]][de_tbl$direction == "up"])),
    down = unique(na.omit(de_tbl[[id_col]][de_tbl$direction == "down"]))
  )
}

combine_directional_enrichments <- function(result_list, source_label = NULL) {
  bind_rows(lapply(names(result_list), function(one_name) {
    extract_enrich_df(result_list[[one_name]], direction_label = one_name, source_label = source_label)
  }))
}

split_gsea_results_by_direction <- function(gsea_df) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(list(
      all = data.frame(),
      up = data.frame(),
      down = data.frame()
    ))
  }
  list(
    all = gsea_df,
    up = gsea_df %>% filter(is.finite(NES), NES > 0),
    down = gsea_df %>% filter(is.finite(NES), NES < 0)
  )
}

subset_leading_edge_table <- function(leading_tbl, gsea_df, id_col = "ID") {
  if (is.null(leading_tbl) || nrow(leading_tbl) == 0 || is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(data.frame())
  }
  resolved_id_col <- c(id_col, "ID", "pathway", "gs_name", "Description")
  resolved_id_col <- resolved_id_col[resolved_id_col %in% colnames(gsea_df)][1]
  if (is.na(resolved_id_col)) {
    return(data.frame())
  }
  valid_ids <- unique(as.character(gsea_df[[resolved_id_col]]))
  if (length(valid_ids) == 0) {
    return(data.frame())
  }
  leading_tbl %>% filter(pathway_id %in% valid_ids)
}

build_directional_dotplot_df <- function(df, top_n = 10) {
  if (nrow(df) == 0) {
    return(df)
  }
  size_col <- if ("Count" %in% colnames(df)) "Count" else if ("setSize" %in% colnames(df)) "setSize" else NULL
  if (is.null(size_col)) {
    df$.size_order__ <- 0
  } else {
    df$.size_order__ <- dplyr::coalesce(as.numeric(df[[size_col]]), 0)
  }
  df %>%
    arrange(direction, pvalue, desc(.size_order__)) %>%
    group_by(direction) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      direction = factor(direction, levels = c("all", "up", "down")),
      Description = factor(Description, levels = rev(unique(Description)))
    ) %>%
    dplyr::select(-.size_order__)
}

extract_leading_edge_table <- function(gsea_df, id_col = "ID", gene_col = "core_enrichment", top_n = 20) {
  if (nrow(gsea_df) == 0) {
    return(data.frame())
  }

  resolved_id_col <- c(id_col, "ID", "pathway", "gs_name", "Description")
  resolved_id_col <- resolved_id_col[resolved_id_col %in% colnames(gsea_df)][1]
  resolved_gene_col <- c(gene_col, "core_enrichment", "leadingEdge", "leading_edge")
  resolved_gene_col <- resolved_gene_col[resolved_gene_col %in% colnames(gsea_df)][1]

  if (is.na(resolved_gene_col)) {
    list_cols <- colnames(gsea_df)[vapply(gsea_df, is.list, logical(1))]
    if (length(list_cols) > 0) {
      resolved_gene_col <- list_cols[[1]]
    }
  }

  if (is.na(resolved_id_col) && "Description" %in% colnames(gsea_df)) {
    resolved_id_col <- "Description"
  }

  if (is.na(resolved_id_col) || is.na(resolved_gene_col)) {
    stop(
      "无法从 GSEA 结果中解析 leading-edge 所需列。现有列为: ",
      paste(colnames(gsea_df), collapse = ", ")
    )
  }

  df_small <- gsea_df %>%
    arrange(pvalue, desc(abs(NES))) %>%
    slice_head(n = top_n)

  gene_values <- df_small[[resolved_gene_col]]
  if (is.list(gene_values)) {
    gene_values <- vapply(gene_values, function(x) paste(x, collapse = "/"), character(1))
  } else {
    gene_values <- as.character(gene_values)
  }

  description_values <- if ("Description" %in% colnames(df_small)) {
    as.character(df_small$Description)
  } else {
    as.character(df_small[[resolved_id_col]])
  }

  tibble::tibble(
    pathway_id = as.character(df_small[[resolved_id_col]]),
    Description = description_values,
    NES = dplyr::coalesce(df_small$NES, NA_real_),
    pvalue = dplyr::coalesce(df_small$pvalue, NA_real_),
    leading_edge = gene_values
  ) %>%
    tidyr::separate_rows(leading_edge, sep = "/") %>%
    filter(!is.na(leading_edge), leading_edge != "")
}

build_leading_edge_heatmap_input <- function(leading_tbl, de_tbl, id_type = c("ENTREZID", "SYMBOL"), top_pathways = 8, top_genes = 40) {
  id_type <- match.arg(id_type)
  if (nrow(leading_tbl) == 0) {
    return(NULL)
  }
  lead_small <- leading_tbl %>%
    group_by(Description) %>%
    slice_head(n = top_genes) %>%
    ungroup()

  if (id_type == "ENTREZID") {
    fc_map <- de_tbl$log2FoldChange
    names(fc_map) <- as.character(de_tbl$ENTREZID)
  } else {
    fc_map <- de_tbl$log2FoldChange
    names(fc_map) <- as.character(de_tbl$gene_label)
  }
  fc_map <- fc_map[!is.na(names(fc_map)) & names(fc_map) != "" & is.finite(fc_map)]
  if (length(fc_map) == 0) {
    return(NULL)
  }

  lead_small$fc <- fc_map[as.character(lead_small$leading_edge)]
  lead_small <- lead_small %>% filter(is.finite(fc))
  if (nrow(lead_small) == 0) {
    return(NULL)
  }

  top_terms <- lead_small %>%
    count(Description, sort = TRUE) %>%
    slice_head(n = top_pathways) %>%
    pull(Description)
  lead_small <- lead_small %>% filter(Description %in% top_terms)
  if (nrow(lead_small) == 0) {
    return(NULL)
  }

  mat <- lead_small %>%
    distinct(Description, leading_edge, .keep_all = TRUE) %>%
    transmute(pathway = Description, gene = leading_edge, fc = fc) %>%
    tidyr::pivot_wider(names_from = pathway, values_from = fc, values_fill = 0)
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  data.matrix(mat)
}

run_fgsea_hallmark <- function(rank_vector, species = msigdb_species) {
  if (length(rank_vector) < 100) {
    return(list(result = data.frame(), pathways = list()))
  }
  hallmark_tbl <- msigdbr::msigdbr(species = species, category = "H") %>%
    dplyr::select(gs_name, gene_symbol)
  hallmark_sets <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)
  fgsea_res <- fgsea::fgsea(
    pathways = hallmark_sets,
    stats = rank_vector,
    eps = 0,
    minSize = 10,
    maxSize = 500
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(pval, desc(abs(NES))) %>%
    dplyr::rename(
      ID = pathway,
      Description = pathway,
      pvalue = pval,
      padj = padj,
      setSize = size,
      leadingEdge = leadingEdge
    )
  if ("leadingEdge" %in% colnames(fgsea_res)) {
    fgsea_res$core_enrichment <- vapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = "/"), character(1))
  }
  list(result = fgsea_res, pathways = hallmark_sets)
}

build_msig_symbol_sets <- function(collection, subcollection = NULL, species = msigdb_species) {
  formals_now <- names(formals(msigdbr::msigdbr))
  subcollection_options <- if (is.null(subcollection)) list(NULL) else as.list(subcollection)
  last_error <- NULL

  for (one_subcollection in subcollection_options) {
    msig_tbl <- tryCatch({
      if ("collection" %in% formals_now) {
        call_args <- list(species = species, collection = collection)
        if ("db_species" %in% formals_now) {
          call_args$db_species <- "HS"
        }
        if (!is.null(one_subcollection) && !is.na(one_subcollection)) {
          call_args$subcollection <- one_subcollection
        }
        do.call(msigdbr::msigdbr, call_args)
      } else {
        call_args <- list(species = species, category = collection)
        if (!is.null(one_subcollection) && !is.na(one_subcollection)) {
          call_args$subcategory <- one_subcollection
        }
        do.call(msigdbr::msigdbr, call_args)
      }
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      data.frame()
    })

    if (nrow(msig_tbl) > 0 && all(c("gene_symbol", "gs_name") %in% colnames(msig_tbl))) {
      return(split(msig_tbl$gene_symbol, msig_tbl$gs_name))
    }
  }

  message(
    "[MSigDB] 未能获取 collection=", collection,
    if (!is.null(subcollection)) paste0(" subcollection=", paste(unlist(subcollection), collapse = "/")) else "",
    if (!is.null(last_error)) paste0("；最后错误: ", last_error) else ""
  )
  list()
}

run_fgsea_symbol_sets <- function(rank_vector, gene_sets, min_size = 10, max_size = 500) {
  if (length(rank_vector) < 100 || length(gene_sets) == 0) {
    return(list(result = data.frame(), pathways = list()))
  }
  rank_vector <- rank_vector[is.finite(rank_vector)]
  rank_vector <- rank_vector[!duplicated(names(rank_vector))]
  gene_sets <- lapply(gene_sets, function(x) unique(intersect(as.character(x), names(rank_vector))))
  gene_sets <- gene_sets[lengths(gene_sets) >= min_size & lengths(gene_sets) <= max_size]
  if (length(gene_sets) == 0) {
    return(list(result = data.frame(), pathways = list()))
  }
  fgsea_res <- fgsea::fgsea(
    pathways = gene_sets,
    stats = rank_vector,
    eps = 0,
    minSize = min_size,
    maxSize = max_size
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(pval, desc(abs(NES))) %>%
    dplyr::rename(
      ID = pathway,
      Description = pathway,
      pvalue = pval,
      padj = padj,
      setSize = size,
      leadingEdge = leadingEdge
    )
  if ("leadingEdge" %in% colnames(fgsea_res)) {
    fgsea_res$core_enrichment <- vapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = "/"), character(1))
  }
  list(result = fgsea_res, pathways = gene_sets)
}

run_fgsea_msig_collection <- function(rank_vector, collection, subcollection = NULL, species = msigdb_species, min_size = 10, max_size = 500) {
  gene_sets <- build_msig_symbol_sets(collection = collection, subcollection = subcollection, species = species)
  run_fgsea_symbol_sets(rank_vector, gene_sets, min_size = min_size, max_size = max_size)
}

run_ssgsea_matrix <- function(expr_mat, gene_sets, min_size = 5, max_size = 500, normalize = TRUE) {
  gene_sets <- gene_sets[vapply(gene_sets, function(x) {
    length(intersect(unique(x), rownames(expr_mat))) >= min_size
  }, logical(1))]
  if (length(gene_sets) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(expr_mat)))
  }

  gene_sets <- lapply(gene_sets, function(x) unique(intersect(as.character(x), rownames(expr_mat))))
  gene_sets <- gene_sets[lengths(gene_sets) >= min_size & lengths(gene_sets) <= max_size]
  if (length(gene_sets) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(expr_mat)))
  }

  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    gsva_param <- GSVA::ssgseaParam(
      exprData = as.matrix(expr_mat),
      geneSets = gene_sets,
      normalize = normalize,
      minSize = min_size,
      maxSize = max_size
    )
    res <- GSVA::gsva(gsva_param, verbose = FALSE)
  } else {
    res <- GSVA::gsva(
      expr = as.matrix(expr_mat),
      gset.idx.list = gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      abs.ranking = TRUE,
      ssgsea.norm = normalize,
      min.sz = min_size,
      max.sz = max_size,
      verbose = FALSE
    )
  }
  as.matrix(res)
}

sanitize_pathway_key <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

run_limma_contrast_table <- function(feature_matrix, group_factor, group1_label, group2_label) {
  if (is.null(feature_matrix) || nrow(feature_matrix) == 0 || ncol(feature_matrix) < 4) {
    return(data.frame())
  }
  group_factor <- factor(group_factor, levels = c(group1_label, group2_label))
  design <- stats::model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)
  contrast <- limma::makeContrasts(contrasts = paste0(group2_label, "-", group1_label), levels = design)
  fit <- limma::lmFit(feature_matrix, design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)
  limma::topTable(fit, number = Inf, sort.by = "P")
}

make_gene_symbol_matrix <- function(expr_mat, de_tbl) {
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  symbol_map <- de_tbl$gene_label
  names(symbol_map) <- de_tbl$gene_id
  symbol_map <- symbol_map[!is.na(names(symbol_map)) & names(symbol_map) != ""]
  symbols <- symbol_map[rownames(expr_mat)]
  keep <- !is.na(symbols) & symbols != ""
  if (!any(keep)) {
    return(matrix(numeric(0), nrow = 0, ncol = ncol(expr_mat)))
  }
  expr_df <- as.data.frame(expr_mat[keep, , drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  expr_df$gene_symbol <- symbols[keep]
  expr_df %>%
    group_by(gene_symbol) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()
}

make_overlap_summary <- function(set_a, set_b, label_a, label_b) {
  set_a <- unique(stats::na.omit(set_a))
  set_b <- unique(stats::na.omit(set_b))
  data.frame(
    category = c(paste0(label_a, "_only"), paste0(label_b, "_only"), paste0(label_a, "_and_", label_b)),
    count = c(length(setdiff(set_a, set_b)), length(setdiff(set_b, set_a)), length(intersect(set_a, set_b))),
    stringsAsFactors = FALSE
  )
}

build_membership_upset <- function(tbl, cols, primary, width = 9, height = 6) {
  if (nrow(tbl) == 0) {
    return(invisible(NULL))
  }
  cols <- cols[cols %in% colnames(tbl)]
  if (length(cols) == 0) {
    return(invisible(NULL))
  }
  upset_input <- tbl %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    mutate(across(everything(), ~ as.integer(dplyr::coalesce(as.logical(.x), FALSE))))
  active_sets <- names(upset_input)[colSums(upset_input, na.rm = TRUE) > 0]
  if (length(active_sets) == 0) {
    return(invisible(NULL))
  }
  save_base_pdf_versions(
    primary = primary,
    width = width,
    height = height,
    expr = UpSetR::upset(
      upset_input,
      nsets = length(active_sets),
      keep.order = TRUE,
      sets = active_sets,
      mb.ratio = c(0.6, 0.4),
      order.by = "freq"
    )
  )
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", x)
  x <- sub("_$", "", x)
  ifelse(nchar(x) == 0, "item", x)
}

detect_sashimi_cmd <- function(bind_root = "/home/h1028/workspace", image = "quay.io/biocontainers/rmats2sashimiplot:3.0.0--py310ha6fa2df_2") {
  local_bin <- Sys.which("rmats2sashimiplot")
  if (nzchar(local_bin)) {
    return(c(local_bin))
  }
  singularity_bin <- Sys.which("singularity")
  if (nzchar(singularity_bin)) {
    return(c(
      singularity_bin, "exec",
      "--bind", sprintf("%s:%s", bind_root, bind_root),
      sprintf("docker://%s", image),
      "rmats2sashimiplot"
    ))
  }
  docker_bin <- Sys.which("docker")
  if (nzchar(docker_bin)) {
    uid <- tryCatch(system("id -u", intern = TRUE), error = function(e) "") %>% trimws()
    gid <- tryCatch(system("id -g", intern = TRUE), error = function(e) "") %>% trimws()
    user_args <- if (nzchar(uid) && nzchar(gid)) c("-u", sprintf("%s:%s", uid, gid)) else character(0)
    return(c(
      docker_bin, "run", "--rm",
      user_args,
      "-v", sprintf("%s:%s", bind_root, bind_root),
      "-w", working_dir,
      image,
      "rmats2sashimiplot"
    ))
  }
  NULL
}

render_sashimi_plot <- function(command_vec, b1, b2, event_type, event_file, out_dir, label1, label2, bind_root = "/home/h1028/workspace") {
  if (is.null(command_vec) || length(command_vec) == 0) {
    return(list(success = FALSE, pdf = NA_character_, png = NA_character_, message = "no_executor"))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gf <- tempfile(pattern = "sashimi_group_", fileext = ".txt")
  on.exit(unlink(gf, force = TRUE), add = TRUE)
  writeLines(c(sprintf("%s: 1-%d", label1, length(b1)), sprintf("%s: %d-%d", label2, length(b1) + 1, length(b1) + length(b2))), gf)
  args <- c(
    command_vec[-1],
    "--b1", paste(b1, collapse = ","),
    "--b2", paste(b2, collapse = ","),
    "--event-type", event_type,
    "-e", event_file,
    "--l1", label1,
    "--l2", label2,
    "--group-info", gf,
    "--exon_s", "1",
    "--intron_s", "20",
    "--min-counts", "1",
    "--fig-width", "9",
    "--fig-height", "7",
    "--font-size", "8",
    "--color", rna_sashimi_color_string,
    "-o", out_dir
  )
  status <- tryCatch(system2(command_vec[[1]], args = args, stdout = FALSE, stderr = FALSE), error = function(e) 1)
  pdf_files <- list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE)
  if (status != 0 || length(pdf_files) == 0) {
    return(list(success = FALSE, pdf = NA_character_, png = NA_character_, message = "render_failed"))
  }
  pdf_file <- pdf_files[[1]]
  png_file <- png_path_from_pdf(pdf_file)
  convert_pdf_to_png(pdf_file, png_file)
  list(success = TRUE, pdf = pdf_file, png = png_file, message = "ok")
}

filter_result_by_pvalue <- function(x, threshold = 0.05) {
  # 当前项目保留 nominal pvalue 过滤以延续既有结果口径。
  # 严格报告时建议同步检查 p.adjust/qvalue/FDR 列，并在正文中声明使用的多重校正口径。
  if (is.null(x)) {
    return(NULL)
  }
  x_df <- make_result_df(x)
  if (nrow(x_df) == 0 || !"pvalue" %in% colnames(x_df)) {
    return(x)
  }
  tryCatch({
    x@result <- x@result[x@result$pvalue < threshold, , drop = FALSE]
    x
  }, error = function(e) x)
}

run_with_guard <- function(label, fn, must_succeed = FALSE) {
  tryCatch({
    fn()
  }, error = function(e) {
    message(label, " 失败: ", e$message)
    if (isTRUE(must_succeed)) {
      stop(label, " 失败: ", e$message, call. = FALSE)
    }
    NULL
  })
}

run_with_timeout <- function(label, fn, timeout_seconds = NULL) {
  use_timeout <- .Platform$OS.type == "unix" &&
    is.numeric(timeout_seconds) &&
    length(timeout_seconds) == 1 &&
    is.finite(timeout_seconds) &&
    timeout_seconds > 0

  if (!use_timeout) {
    return(fn())
  }

  job <- parallel::mcparallel(fn(), silent = TRUE)
  collected <- parallel::mccollect(job, wait = FALSE, timeout = timeout_seconds)

  if (is.null(collected)) {
    tools::pskill(job$pid)
    try(suppressWarnings(parallel::mccollect(job, wait = FALSE, timeout = 1)), silent = TRUE)
    stop(label, " 超过 ", timeout_seconds, " 秒未返回，已中止该次尝试。", call. = FALSE)
  }

  value <- collected[[1]]
  if (inherits(value, "try-error")) {
    child_condition <- attr(value, "condition")
    if (inherits(child_condition, "condition")) {
      stop(child_condition)
    }
    stop(as.character(value), call. = FALSE)
  }

  value
}

run_clusterprofiler_with_retry <- function(label, fn, attempts = 3, sleep_seconds = 3, require_nonempty = FALSE, must_succeed = FALSE, timeout_seconds = NULL) {
  last_error_message <- NULL
  last_result <- NULL
  for (one_try in seq_len(attempts)) {
    result <- tryCatch(
      run_with_timeout(label = label, fn = fn, timeout_seconds = timeout_seconds),
      error = function(e) {
        last_error_message <<- conditionMessage(e)
        message(label, " [try ", one_try, "/", attempts, "] 失败: ", last_error_message)
        NULL
      }
    )
    last_result <- result
    result_df <- make_result_df(result)
    if (!is.null(result) && (!require_nonempty || nrow(result_df) > 0)) {
      return(result)
    }
    if (one_try < attempts) {
      Sys.sleep(sleep_seconds)
    }
  }
  if (!is.null(last_error_message)) {
    message(label, " 最终失败: ", last_error_message)
  }
  if (isTRUE(must_succeed) && is.null(last_result)) {
    stop(label, " 多次重试后仍失败。", call. = FALSE)
  }
  last_result
}

combine_gsea_go_results <- function(gsea_go_parts) {
  combined <- bind_rows(lapply(names(gsea_go_parts), function(one_ont) {
    one_df <- make_result_df(gsea_go_parts[[one_ont]])
    if (nrow(one_df) == 0) {
      return(data.frame())
    }
    one_df$ONTOLOGY <- one_ont
    one_df
  }))
  if (nrow(combined) == 0) {
    return(combined)
  }
  combined %>%
    arrange(pvalue, desc(abs(NES)))
}

choose_primary_gsea_obj <- function(gsea_parts) {
  if (length(gsea_parts) == 0) {
    return(NULL)
  }
  candidate_stats <- bind_rows(lapply(names(gsea_parts), function(one_name) {
    one_df <- make_result_df(gsea_parts[[one_name]])
    data.frame(
      ontology = one_name,
      sig_n = sum(one_df$pvalue < 0.05, na.rm = TRUE),
      total_n = nrow(one_df),
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(candidate_stats) == 0) {
    return(NULL)
  }
  primary_name <- candidate_stats %>%
    arrange(desc(sig_n), desc(total_n), ontology) %>%
    slice_head(n = 1) %>%
    pull(ontology)
  gsea_parts[[primary_name]]
}

build_gsea_dotplot_df <- function(gsea_df, top_n = 6) {
  if (nrow(gsea_df) == 0) {
    return(gsea_df)
  }
  if ("ONTOLOGY" %in% colnames(gsea_df)) {
    gsea_df %>%
      group_by(ONTOLOGY) %>%
      arrange(pvalue, desc(abs(NES)), .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      ungroup() %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))
  } else {
    gsea_df %>%
      arrange(pvalue, desc(abs(NES))) %>%
      slice_head(n = min(15, nrow(gsea_df))) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))
  }
}

resolve_named_colors <- function(base_colors, target_names) {
  target_names <- unique(as.character(target_names))
  target_names <- target_names[!is.na(target_names) & target_names != ""]
  if (length(target_names) == 0) {
    return(setNames(character(0), character(0)))
  }
  palette_pool <- unname(base_colors)
  if (length(palette_pool) == 0) {
    return(setNames(rep("#D4D9DE", length(target_names)), target_names))
  }
  out <- setNames(rep(palette_pool, length.out = length(target_names)), target_names)
  known_names <- intersect(names(base_colors), target_names)
  out[known_names] <- unname(base_colors[known_names])
  out
}

viz_colors <- function(palette, keys) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & keys != ""]
  if (length(keys) == 0) {
    return(setNames(character(0), character(0)))
  }
  palette_values <- unname(palette)
  if (length(palette_values) == 0) {
    palette_values <- unname(viz_pal$neutral["mid"])
  }
  out <- setNames(rep(palette_values, length.out = length(keys)), keys)
  known_keys <- intersect(names(palette), keys)
  out[known_keys] <- unname(palette[known_keys])
  out
}

viz_color <- function(palette, key) {
  unname(viz_colors(palette, key)[[1]])
}

viz_rekey_colors <- function(palette, key_map) {
  out <- vapply(unname(key_map), function(one_key) viz_color(palette, one_key), character(1))
  names(out) <- names(key_map)
  out
}

resolve_enrichment_category_colors <- function(categories, analysis_type = c("GO", "ORA", "KEGG", "GSEA", "Reactome", "Hallmark")) {
  analysis_type <- match.arg(analysis_type)
  categories <- unique(as.character(categories))
  categories <- categories[!is.na(categories) & categories != ""]
  if (length(categories) == 0) {
    return(setNames(character(0), character(0)))
  }

  # All enrichment-circle category colors follow the GO circle standard.
  go_circle_standard_colors <- c(
    BP = group2_color,
    CC = group1_color,
    MF = unname(viz_pal$accent["purple"]),
    Metabolism = unname(viz_pal$accent["warm"]),
    Signaling = unname(viz_pal$accent["teal"]),
    Disease = unname(viz_pal$accent["peach"]),
    Immune = unname(viz_pal$accent["green"]),
    "Cellular Process" = unname(viz_pal$accent["yellow"]),
    Transport = unname(viz_pal$neutral["dark"]),
    Development = unname(viz_pal$accent["yellow"]),
    Stress = unname(viz_pal$neutral["mid"]),
    Mixed = unname(viz_pal$neutral["dark"]),
    GO = group1_color,
    KEGG = unname(viz_pal$accent["warm"]),
    Reactome = unname(viz_pal$accent["purple"]),
    Hallmark = unname(viz_pal$accent["teal"])
  )

  color_pool <- c(go_circle_standard_colors, viz_pal$accent, viz_pal$group, viz_pal$neutral[c("dark", "mid", "light")])
  resolve_named_colors(color_pool, categories)
}

save_module_target_enrichment <- function(
    module_df_annotated,
    module_trait_tbl,
    sig_ids,
    label_prefix,
    out_csv,
    out_pdf,
    high_color = viz_pal$sig_gradient["high"]) {
  write_empty_enrichment <- function(status_reason) {
    empty_tbl <- data.frame(
      Module = NA_character_,
      module_size = NA_integer_,
      target_gene_count = NA_integer_,
      target_gene_ratio = NA_real_,
      background_ratio = NA_real_,
      fold_enrichment = NA_real_,
      pvalue = NA_real_,
      module_trait_correlation = NA_real_,
      status = status_reason,
      stringsAsFactors = FALSE
    )
    write.csv(empty_tbl, out_csv, row.names = FALSE)
    invisible(empty_tbl)
  }

  if (is.null(module_df_annotated) || !is.data.frame(module_df_annotated) || nrow(module_df_annotated) == 0) {
    return(write_empty_enrichment("empty_module_annotation"))
  }
  if (is.null(module_trait_tbl) || !is.data.frame(module_trait_tbl) || nrow(module_trait_tbl) == 0) {
    return(write_empty_enrichment("empty_module_trait_table"))
  }
  required_module_cols <- c("gene_id", "ensembl_id", "module_label")
  missing_module_cols <- setdiff(required_module_cols, colnames(module_df_annotated))
  if (length(missing_module_cols) > 0) {
    return(write_empty_enrichment(paste0("missing_module_columns:", paste(missing_module_cols, collapse = ","))))
  }
  if (!all(c("Module", "correlation") %in% colnames(module_trait_tbl))) {
    return(write_empty_enrichment("missing_module_trait_columns"))
  }

  sig_ids <- unique(stats::na.omit(sig_ids))
  if (length(sig_ids) == 0) {
    return(write_empty_enrichment("no_significant_genes"))
  }

  module_universe_ids <- unique(module_df_annotated$gene_id)
  target_ids <- unique(module_df_annotated$gene_id[module_df_annotated$ensembl_id %in% sig_ids])
  if (length(target_ids) == 0) {
    return(write_empty_enrichment("no_significant_genes_in_wgcna_universe"))
  }

  enrich_tbl <- bind_rows(lapply(sort(unique(module_df_annotated$module_label)), function(one_module) {
    module_genes <- unique(module_df_annotated$gene_id[module_df_annotated$module_label == one_module])
    in_module <- sum(module_genes %in% target_ids)
    in_module_non <- length(module_genes) - in_module
    out_module <- length(target_ids) - in_module
    out_module_non <- length(module_universe_ids) - length(module_genes) - out_module
    fisher_res <- fisher.test(matrix(c(in_module, in_module_non, out_module, out_module_non), nrow = 2), alternative = "greater")
    bg_ratio <- length(target_ids) / length(module_universe_ids)
    module_ratio <- if (length(module_genes) > 0) in_module / length(module_genes) else 0
    data.frame(
      Module = one_module,
      module_size = length(module_genes),
      target_gene_count = in_module,
      target_gene_ratio = module_ratio,
      background_ratio = bg_ratio,
      fold_enrichment = ifelse(bg_ratio > 0, module_ratio / bg_ratio, NA_real_),
      pvalue = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  })) %>%
    left_join(module_trait_tbl %>% dplyr::select(Module, module_trait_correlation = correlation), by = "Module") %>%
    arrange(pvalue, desc(fold_enrichment))

  write.csv(enrich_tbl, out_csv, row.names = FALSE)

  enrich_plot <- ggplot(enrich_tbl, aes(x = module_trait_correlation, y = reorder(Module, module_trait_correlation), size = target_gene_count, color = -log10(pvalue + 1e-300))) +
    geom_point(alpha = 0.85) +
    scale_color_heatmap_sequential() +
    labs(title = paste0(label_prefix, " module enrichment"), x = "Module-trait correlation", y = NULL, size = "Gene count", color = "-log10(pvalue)") +
    theme_rnaseq_bw()
  save_plot_versions(enrich_plot, out_pdf, width = 8.5, height = max(4.5, 0.8 * nrow(enrich_tbl)))
  invisible(enrich_tbl)
}

save_das_gene_burden_plot <- function(as_gene_priority, primary, top_n = 20, width = 8.5, height = 6) {
  input_rows <- if (is.data.frame(as_gene_priority)) nrow(as_gene_priority) else 0
  if (!is.data.frame(as_gene_priority) || input_rows == 0) {
    return(save_placeholder_plot(primary, "Top DAS gene burden", "empty_as_gene_priority", width, height, input_rows))
  }
  required_cols <- c("gene_label", "gene_id_clean", "significant_event_count", "de_status")
  missing_cols <- setdiff(required_cols, colnames(as_gene_priority))
  if (length(missing_cols) > 0) {
    return(save_placeholder_plot(primary, "Top DAS gene burden", paste0("missing_columns:", paste(missing_cols, collapse = ",")), width, height, input_rows))
  }

  plot_df <- as_gene_priority %>%
    slice_head(n = top_n) %>%
    mutate(
      gene_label_plot = dplyr::coalesce(gene_label, gene_id_clean),
      gene_label_plot = factor(gene_label_plot, levels = rev(gene_label_plot))
    )
  if (nrow(plot_df) == 0) {
    return(save_placeholder_plot(primary, "Top DAS gene burden", "empty_top_das_gene_table", width, height, input_rows))
  }

  p <- ggplot(plot_df, aes(x = gene_label_plot, y = significant_event_count, fill = de_status)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = viz_colors(viz_pal$overlap, plot_df$de_status)) +
    labs(title = "Top DAS 基因事件负荷", x = NULL, y = "Significant AS event count") +
    theme_rnaseq_minimal()
  save_plot_versions(p, primary, width = width, height = height)
  record_plot_status(primary, status = "generated", reason = "", input_rows = input_rows)
}

save_de_das_gene_scatter <- function(as_gene_priority, primary, width = 8, height = 6) {
  input_rows <- if (is.data.frame(as_gene_priority)) nrow(as_gene_priority) else 0
  if (!is.data.frame(as_gene_priority) || input_rows == 0) {
    return(save_placeholder_plot(primary, "DE and DAS gene scatter", "empty_as_gene_priority", width, height, input_rows))
  }
  required_cols <- c("de_log2FoldChange", "max_abs_delta_psi", "significant_event_count", "de_status")
  missing_cols <- setdiff(required_cols, colnames(as_gene_priority))
  if (length(missing_cols) > 0) {
    return(save_placeholder_plot(primary, "DE and DAS gene scatter", paste0("missing_columns:", paste(missing_cols, collapse = ",")), width, height, input_rows))
  }

  plot_df <- as_gene_priority %>%
    filter(is.finite(de_log2FoldChange), is.finite(max_abs_delta_psi))
  if (nrow(plot_df) == 0) {
    return(save_placeholder_plot(primary, "DE and DAS gene scatter", "no_finite_de_or_delta_psi", width, height, input_rows))
  }

  p <- ggplot(plot_df, aes(x = de_log2FoldChange, y = max_abs_delta_psi, size = significant_event_count, color = de_status)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = viz_colors(viz_pal$overlap, plot_df$de_status)) +
    labs(title = "DE 与 DAS 整合散点图", x = "Gene log2FC", y = "Max |delta PSI|") +
    theme_rnaseq_bw()
  save_plot_versions(p, primary, width = width, height = height)
  record_plot_status(primary, status = "generated", reason = "", input_rows = input_rows)
}

save_das_module_burden_plot <- function(module_burden_tbl, primary, width = 8.5, height = NULL) {
  input_rows <- if (is.data.frame(module_burden_tbl)) nrow(module_burden_tbl) else 0
  height <- height %||% max(5, 0.7 * input_rows)
  if (!is.data.frame(module_burden_tbl) || input_rows == 0) {
    return(save_placeholder_plot(primary, "DAS module burden", "empty_module_burden_table", width, height, input_rows))
  }
  required_cols <- c("module_label", "total_as_events", "de_status")
  missing_cols <- setdiff(required_cols, colnames(module_burden_tbl))
  if (length(missing_cols) > 0) {
    return(save_placeholder_plot(primary, "DAS module burden", paste0("missing_columns:", paste(missing_cols, collapse = ",")), width, height, input_rows))
  }

  p <- ggplot(module_burden_tbl, aes(x = reorder(module_label, total_as_events), y = total_as_events, fill = de_status)) +
    geom_col(position = "stack", width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = viz_colors(viz_pal$overlap, module_burden_tbl$de_status)) +
    labs(title = "各模块 DAS 事件负荷", x = NULL, y = "Total significant AS events") +
    theme_rnaseq_minimal()
  save_plot_versions(p, primary, width = width, height = height)
  record_plot_status(primary, status = "generated", reason = "", input_rows = input_rows)
}

parse_ratio_value <- function(x) {
  if (length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") {
    return(NA_real_)
  }
  ratio_parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(ratio_parts) == 2) {
    numerator <- suppressWarnings(as.numeric(ratio_parts[[1]]))
    denominator <- suppressWarnings(as.numeric(ratio_parts[[2]]))
    if (is.finite(numerator) && is.finite(denominator) && denominator > 0) {
      return(numerator / denominator)
    }
  }
  suppressWarnings(as.numeric(x))[1]
}

build_ora_plot_df <- function(ora_df, top_n = 15, facet_col = NULL) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(data.frame())
  }

  plot_df <- as.data.frame(ora_df, stringsAsFactors = FALSE)
  plot_df$Description <- as.character(plot_df$Description)
  plot_df$pvalue <- suppressWarnings(as.numeric(plot_df$pvalue))
  plot_df$Count <- suppressWarnings(as.numeric(plot_df$Count))
  plot_df$ratio_value <- if ("GeneRatio" %in% colnames(plot_df)) {
    vapply(plot_df$GeneRatio, parse_ratio_value, numeric(1))
  } else {
    plot_df$Count
  }
  plot_df$ratio_value[!is.finite(plot_df$ratio_value)] <- plot_df$Count[!is.finite(plot_df$ratio_value)]
  plot_df$size_value <- dplyr::coalesce(plot_df$Count, plot_df$ratio_value, 1)
  plot_df$neg_log10_p <- -log10(dplyr::coalesce(plot_df$pvalue, 1) + 1e-300)

  if (!is.null(facet_col) && facet_col %in% colnames(plot_df)) {
    plot_df <- plot_df %>%
      group_by(.data[[facet_col]]) %>%
      arrange(pvalue, desc(ratio_value), .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      ungroup()
  } else {
    plot_df <- plot_df %>%
      arrange(pvalue, desc(ratio_value)) %>%
      slice_head(n = min(top_n, nrow(plot_df)))
  }

  plot_df %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
}

ora_axis_label <- function(ora_df) {
  if (!is.null(ora_df) && "GeneRatio" %in% colnames(ora_df)) {
    return("Gene ratio")
  }
  "Gene count"
}

shorten_plot_labels <- function(labels, max_chars = 48) {
  vapply(as.character(labels), function(label) {
    if (is.na(label) || !nzchar(trimws(label))) {
      return("")
    }
    label <- gsub("\\s+", " ", trimws(label))
    if (nchar(label, type = "width") <= max_chars) {
      return(label)
    }
    paste0(substr(label, 1, max_chars - 3), "...")
  }, character(1), USE.NAMES = FALSE)
}

save_ora_barplot_df <- function(ora_df, primary, title_prefix, top_n = 15, facet_col = NULL, analysis_type = c("ORA", "GO"), single_fill = NULL) {
  analysis_type <- match.arg(analysis_type)
  plot_df <- build_ora_plot_df(ora_df, top_n = top_n, facet_col = facet_col)
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(plot_df, aes(x = Description, y = ratio_value, fill = neg_log10_p)) +
    geom_col(width = 0.75, alpha = 0.92) +
    scale_fill_gradientn(colours = unname(viz_pal$enrich_p))

  if (!is.null(facet_col) && facet_col %in% colnames(plot_df)) {
    p <- p + facet_grid(stats::as.formula(paste(facet_col, "~ .")), scales = "free_y", space = "free_y")
  }

  p <- p +
    coord_flip() +
    labs(title = title_prefix, x = NULL, y = ora_axis_label(ora_df), fill = "-log10(PValue)") +
    theme_rnaseq_bw(base_size = viz_style$compact_size) +
    theme(panel.grid.major.y = element_blank())

  save_plot_versions(p, primary, width = 10.5, height = 8)
}

save_ora_dotplot_df <- function(ora_df, primary, title_prefix, top_n = 15, facet_col = NULL) {
  plot_df <- build_ora_plot_df(ora_df, top_n = top_n, facet_col = facet_col)
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(plot_df, aes(x = ratio_value, y = Description, size = size_value, color = neg_log10_p)) +
    geom_point(alpha = 0.88) +
    scale_color_gradientn(colours = unname(viz_pal$enrich_p)) +
    scale_y_discrete(labels = function(x) shorten_plot_labels(x, max_chars = 48)) +
    labs(
      title = title_prefix,
      x = ora_axis_label(ora_df),
      y = NULL,
      size = "Gene count",
      color = "-log10(PValue)"
    ) +
    theme_rnaseq_bw(base_size = viz_style$compact_size)

  if (!is.null(facet_col) && facet_col %in% colnames(plot_df)) {
    p <- p + facet_grid(stats::as.formula(paste(facet_col, "~ .")), scales = "free_y", space = "free_y")
  }

  save_plot_versions(p, primary, width = 10.5, height = 8)
}

build_gene_foldchange_symbol <- function(df) {
  if (nrow(df) == 0 || !"gene_label" %in% colnames(df) || !"log2FoldChange" %in% colnames(df)) {
    return(numeric(0))
  }
  fc_tbl <- df %>%
    filter(!is.na(gene_label), gene_label != "", is.finite(log2FoldChange)) %>%
    arrange(pvalue, desc(abs(log2FoldChange))) %>%
    distinct(gene_label, .keep_all = TRUE)
  fc <- fc_tbl$log2FoldChange
  names(fc) <- fc_tbl$gene_label
  sort(fc, decreasing = TRUE)
}

build_goplot_inputs <- function(enrich_df, de_tbl, category_label = "TERM", top_n = 10) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0 || nrow(de_tbl) == 0) {
    return(NULL)
  }

  enrich_top <- enrich_df %>%
    arrange(pvalue) %>%
    slice_head(n = top_n) %>%
    dplyr::select(ID, Description, pvalue, geneID)

  if (nrow(enrich_top) == 0) {
    return(NULL)
  }

  goenrichment <- enrich_top
  colnames(goenrichment) <- c("ID", "term", "adj_pval", "genes")
  goenrichment$category <- category_label
  goenrichment <- goenrichment[, c("category", "ID", "term", "adj_pval", "genes")]
  goenrichment$genes <- gsub("/", ",", goenrichment$genes)

  genedata <- de_tbl %>%
    filter(!is.na(gene_label), gene_label != "", is.finite(log2FoldChange)) %>%
    transmute(ID = toupper(gene_label), logFC = log2FoldChange) %>%
    distinct(ID, .keep_all = TRUE)

  circ <- tryCatch(GOplot::circle_dat(goenrichment, genedata), error = function(e) NULL)
  if (is.null(circ) || nrow(circ) == 0) {
    return(NULL)
  }

  list(goenrichment = goenrichment, genedata = genedata, circ = circ)
}

make_upset_gene_sets <- function(enrich_df, top_n = 15, max_label = 30) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return(list())
  }
  top_go <- enrich_df %>% arrange(pvalue) %>% slice_head(n = top_n)
  gene_sets <- list()
  for (i in seq_len(nrow(top_go))) {
    genes <- unlist(strsplit(as.character(top_go$geneID[i]), "/"))
    genes <- genes[genes != ""]
    if (length(genes) > 0) {
      set_name <- as.character(top_go$Description[i])
      if (nchar(set_name) > max_label) {
        set_name <- paste0(substr(set_name, 1, max_label), "...")
      }
      gene_sets[[set_name]] <- genes
    }
  }
  gene_sets
}

prepare_enrichment_circle_data <- function(enrich_result, de_tbl, analysis_type = c("GO", "KEGG", "GSEA")) {
  analysis_type <- match.arg(analysis_type)
  enrich_df <- make_result_df(enrich_result)
  if (nrow(enrich_df) == 0 || nrow(de_tbl) == 0) {
    return(NULL)
  }

  if ("pvalue" %in% colnames(enrich_df)) {
    enrich_df <- enrich_df %>% arrange(pvalue)
  } else if ("NES" %in% colnames(enrich_df)) {
    enrich_df <- enrich_df %>% arrange(desc(abs(NES)))
  }

  if (analysis_type == "GO" && "ONTOLOGY" %in% colnames(enrich_df)) {
    enrich_df <- bind_rows(
      enrich_df %>% filter(ONTOLOGY == "BP") %>% slice_head(n = 5),
      enrich_df %>% filter(ONTOLOGY == "MF") %>% slice_head(n = 5),
      enrich_df %>% filter(ONTOLOGY == "CC") %>% slice_head(n = 5)
    ) %>% distinct(ID, .keep_all = TRUE)
  } else {
    enrich_df <- enrich_df %>% slice_head(n = 15)
  }

  if (nrow(enrich_df) == 0) {
    return(NULL)
  }

  gene_fc_symbol <- de_tbl$log2FoldChange
  names(gene_fc_symbol) <- de_tbl$gene_label
  gene_fc_symbol <- gene_fc_symbol[!is.na(names(gene_fc_symbol)) & names(gene_fc_symbol) != "" & is.finite(gene_fc_symbol)]

  gene_fc_entrez <- de_tbl$log2FoldChange
  names(gene_fc_entrez) <- as.character(de_tbl$ENTREZID)
  gene_fc_entrez <- gene_fc_entrez[!is.na(names(gene_fc_entrez)) & names(gene_fc_entrez) != "" & is.finite(gene_fc_entrez)]

  gene_symbol_map <- de_tbl$gene_label
  names(gene_symbol_map) <- as.character(de_tbl$ENTREZID)
  gene_symbol_map <- gene_symbol_map[!is.na(names(gene_symbol_map)) & names(gene_symbol_map) != ""]

  parse_ratio <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(NA_real_)
    val <- tryCatch(eval(parse(text = x)), error = function(e) NA_real_)
    as.numeric(val)[1]
  }

  get_pathway_genes <- function(one_row) {
    gene_str <- dplyr::coalesce(
      if ("core_enrichment" %in% colnames(one_row)) one_row$core_enrichment else NA_character_,
      if ("geneID" %in% colnames(one_row)) one_row$geneID else NA_character_,
      ""
    )
    gene_list <- unlist(strsplit(as.character(gene_str), "/"))
    gene_list[gene_list != ""]
  }

  classify_category <- function(one_row) {
    if (analysis_type == "GO" && "ONTOLOGY" %in% colnames(one_row)) {
      return(as.character(one_row$ONTOLOGY))
    }
    desc <- as.character(one_row$Description %||% one_row$ID)
    if (grepl("metabolism|biosynthesis", desc, ignore.case = TRUE)) return("Metabolism")
    if (grepl("signaling|signal", desc, ignore.case = TRUE)) return("Signaling")
    if (grepl("cancer|disease|infection", desc, ignore.case = TRUE)) return("Disease")
    if (grepl("immune|response", desc, ignore.case = TRUE)) return("Immune")
    "Cellular Process"
  }

  result_list <- lapply(seq_len(nrow(enrich_df)), function(i) {
    one_row <- enrich_df[i, , drop = FALSE]
    pathway_id <- as.character(one_row$ID)
    pathway_desc <- as.character(one_row$Description %||% pathway_id)
    gene_list <- get_pathway_genes(one_row)
    gene_num <- length(gene_list)

    gene_symbols <- ifelse(gene_list %in% names(gene_symbol_map), gene_symbol_map[gene_list], gene_list)
    gene_symbols <- ifelse(is.na(gene_symbols) | gene_symbols == "", gene_list, gene_symbols)
    gene_fc <- ifelse(gene_symbols %in% names(gene_fc_symbol), gene_fc_symbol[gene_symbols], NA_real_)
    if (all(!is.finite(gene_fc)) && any(gene_list %in% names(gene_fc_entrez))) {
      gene_fc <- ifelse(gene_list %in% names(gene_fc_entrez), gene_fc_entrez[gene_list], gene_fc)
      gene_symbols <- ifelse(gene_list %in% names(gene_symbol_map), gene_symbol_map[gene_list], gene_symbols)
    }

    valid_fc <- gene_fc[is.finite(gene_fc)]
    valid_symbols <- gene_symbols[is.finite(gene_fc)]
    up_symbols <- valid_symbols[valid_fc > 0]
    down_symbols <- valid_symbols[valid_fc < 0]

    score_value <- if (analysis_type == "GSEA" && "NES" %in% colnames(one_row)) {
      as.numeric(one_row$NES)
    } else if (all(c("GeneRatio", "BgRatio") %in% colnames(one_row))) {
      gr <- parse_ratio(one_row$GeneRatio)
      bg <- parse_ratio(one_row$BgRatio)
      ifelse(is.finite(gr) && is.finite(bg) && bg > 0, gr / bg, NA_real_)
    } else {
      NA_real_
    }

    gene_fc_df <- data.frame(
      SYMBOL = valid_symbols,
      log2FoldChange = valid_fc,
      stringsAsFactors = FALSE
    ) %>%
      distinct(SYMBOL, .keep_all = TRUE) %>%
      arrange(desc(abs(log2FoldChange)))

    data.frame(
      id = pathway_id,
      ID = pathway_id,
      Description = pathway_desc,
      Description_Short = ifelse(nchar(pathway_desc) > 100, paste0(substr(pathway_desc, 1, 97), "..."), pathway_desc),
      category = classify_category(one_row),
      GO_Ontology = if ("ONTOLOGY" %in% colnames(one_row)) as.character(one_row$ONTOLOGY) else NA_character_,
      KEGG_Pathway_Number = if (grepl(paste0("^", kegg_organism_code), pathway_id)) sub(paste0("^", kegg_organism_code), "", pathway_id) else NA_character_,
      gene_num.min = 0,
      gene_num.max = gene_num,
      score = score_value,
      neg_log10Pvalue = -log10(as.numeric(one_row$pvalue %||% 1) + 1e-100),
      Pvalue = as.numeric(one_row$pvalue %||% NA_real_),
      GeneRatio = if ("GeneRatio" %in% colnames(one_row)) as.character(one_row$GeneRatio) else NA_character_,
      BgRatio = if ("BgRatio" %in% colnames(one_row)) as.character(one_row$BgRatio) else NA_character_,
      EnrichmentFactor = if (all(c("GeneRatio", "BgRatio") %in% colnames(one_row))) {
        gr <- parse_ratio(one_row$GeneRatio)
        bg <- parse_ratio(one_row$BgRatio)
        ifelse(is.finite(gr) && is.finite(bg) && bg > 0, gr / bg, NA_real_)
      } else {
        NA_real_
      },
      up.regulated = length(up_symbols),
      down.regulated = length(down_symbols),
      GeneCount = gene_num,
      UpregulatedCount = length(up_symbols),
      DownregulatedCount = length(down_symbols),
      AllGenes = paste(gene_symbols, collapse = "; "),
      Genes_Sorted_By_FC = paste0(gene_fc_df$SYMBOL, " (FC=", round(gene_fc_df$log2FoldChange, 2), ")", collapse = "; "),
      High_FC_Genes = paste0(gene_fc_df$SYMBOL[abs(gene_fc_df$log2FoldChange) > 1], " (FC=", round(gene_fc_df$log2FoldChange[abs(gene_fc_df$log2FoldChange) > 1], 2), ")", collapse = "; "),
      UpregulatedGenes = paste(up_symbols, collapse = "; "),
      DownregulatedGenes = paste(down_symbols, collapse = "; "),
      MeanLog2FC = ifelse(length(valid_fc) > 0, mean(valid_fc), NA_real_),
      MedianLog2FC = ifelse(length(valid_fc) > 0, median(valid_fc), NA_real_),
      MaxAbsLog2FC = ifelse(length(valid_fc) > 0, max(abs(valid_fc)), NA_real_),
      MinLog2FC = ifelse(length(valid_fc) > 0, min(valid_fc), NA_real_),
      stringsAsFactors = FALSE
    )
  })

  detailed_table <- bind_rows(result_list)
  if (nrow(detailed_table) == 0) {
    return(NULL)
  }

  detailed_table <- detailed_table %>%
    mutate(
      GeneCount_Total = UpregulatedCount + DownregulatedCount,
      Upregulated_Percentage = ifelse(GeneCount_Total > 0, round(100 * UpregulatedCount / GeneCount_Total, 1), 0),
      Downregulated_Percentage = ifelse(GeneCount_Total > 0, round(100 * DownregulatedCount / GeneCount_Total, 1), 0),
      SignificanceLevel = case_when(
        !is.na(Pvalue) & Pvalue < 0.001 ~ "Highly_Significant",
        !is.na(Pvalue) & Pvalue < 0.01 ~ "Very_Significant",
        !is.na(Pvalue) & Pvalue < 0.05 ~ "Significant",
        !is.na(Pvalue) & Pvalue < 0.1 ~ "Marginally_Significant",
        TRUE ~ "Unknown"
      ),
      EnrichmentStrength = case_when(
        is.na(EnrichmentFactor) ~ "Unknown",
        EnrichmentFactor > 2 ~ "Very_Strong",
        EnrichmentFactor > 1 ~ "Strong",
        EnrichmentFactor > 0.5 ~ "Moderate",
        EnrichmentFactor > 0.2 ~ "Weak",
        TRUE ~ "Very_Weak"
      ),
      GeneSet_Size_Category = case_when(
        GeneCount >= 50 ~ "Large",
        GeneCount >= 20 ~ "Medium",
        GeneCount >= 10 ~ "Small",
        TRUE ~ "Very_Small"
      ),
      Expression_Pattern = case_when(
        UpregulatedCount > DownregulatedCount * 2 ~ "Predominantly_Upregulated",
        DownregulatedCount > UpregulatedCount * 2 ~ "Predominantly_Downregulated",
        UpregulatedCount == DownregulatedCount ~ "Balanced",
        TRUE ~ "Mixed"
      ),
      Comprehensive_Score = ifelse(is.finite(EnrichmentFactor), round(neg_log10Pvalue * EnrichmentFactor, 2), round(neg_log10Pvalue, 2))
    ) %>%
    arrange(category, desc(Comprehensive_Score), desc(neg_log10Pvalue))

  rownames(detailed_table) <- detailed_table$id
  detailed_table
}

draw_enrichment_circle_plot <- function(data, title, filename, analysis_type = c("GO", "KEGG", "GSEA"), label_type = c("id", "description")) {
  analysis_type <- match.arg(analysis_type)
  label_type <- match.arg(label_type)
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  resolve_circle_category_colors <- function(categories) {
    categories <- unique(as.character(categories))
    categories <- categories[!is.na(categories) & categories != ""]
    if (length(categories) == 0) {
      return(setNames(character(0), character(0)))
    }
    go_circle_seed <- c(
      viz_pal$go_category["BP"],
      viz_pal$go_category["CC"],
      viz_pal$go_category["MF"],
      viz_pal$accent["green"],
      viz_pal$accent["teal"],
      viz_pal$neutral["dark"],
      viz_pal$neutral["mid"]
    )
    out <- setNames(rep(go_circle_seed, length.out = length(categories)), categories)
    known_names <- intersect(names(viz_pal$go_category), categories)
    out[known_names] <- viz_pal$go_category[known_names]
    out
  }

  categories <- unique(data$category)
  # GO / KEGG / GSEA all reuse the GO circle category palette.
  cat_colors <- resolve_circle_category_colors(categories)
  row_colors <- cat_colors[data$category]
  circle_p_colors <- c(
    low = unname(viz_pal$enrich_p["low"]),
    mid = unname(viz_pal$enrich_p["mid"]),
    high = unname(viz_pal$enrich_p["high"])
  )
  circle_reg_colors <- c(
    up = unname(viz_pal$de["up"]),
    down = unname(viz_pal$de["down"])
  )
  circle_score_colors <- c(
    pos = unname(viz_pal$accent["green"]),
    neg = unname(viz_pal$accent["purple"])
  )

  save_base_pdf_versions(
    primary = filename,
    width = 25,
    height = 25,
    expr = {
      circlize::circos.clear()
      circlize::circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.01, 0.01), canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5))
      circlize::circos.genomicInitialize(data[, c("id", "gene_num.min", "gene_num.max")], plotType = NULL)

      circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.05,
        bg.border = NA,
        bg.col = paste0(row_colors, "40"),
        panel.fun = function(x, y) {}
      )

      for (si in circlize::get.all.sector.index()) {
        xlim <- circlize::get.cell.meta.data("xlim", sector.index = si, track.index = 1)
        ylim <- circlize::get.cell.meta.data("ylim", sector.index = si, track.index = 1)
        label_text <- if (label_type == "description") data[si, "Description"] else si
        if (length(label_text) == 0 || is.na(label_text) || !nzchar(as.character(label_text))) {
          label_text <- si
        }
        if (label_type == "description" && nchar(as.character(label_text)) > 50) {
          label_text <- paste0(substr(as.character(label_text), 1, 47), "...")
        }

        xcenter <- mean(xlim)
        y_top <- ylim[2]
        pos <- circlize::circlize(xcenter, y_top, sector.index = si, track.index = 1)
        theta_rad <- pos[1, "theta"] * pi / 180
        rou_start <- pos[1, "rou"] + 0.05
        rou_end <- rou_start + 0.1
        x_start <- rou_start * cos(theta_rad)
        y_start <- rou_start * sin(theta_rad)
        x_end <- rou_end * cos(theta_rad)
        y_end <- rou_end * sin(theta_rad)
        theta_norm <- pos[1, "theta"] %% 360
        is_right <- (theta_norm <= 90 || theta_norm >= 270)
        x_text <- if (is_right) x_end + 0.1 else x_end - 0.1
        y_text <- y_end

        graphics::lines(c(x_start, x_end), c(y_start, y_end), lty = 2, col = "black", lwd = 0.5)
        graphics::lines(c(x_end, x_text), c(y_end, y_text), lty = 2, col = "black", lwd = 0.5)
        graphics::text(x_text, y_text, label_text, adj = if (is_right) c(0, 0.5) else c(1, 0.5), cex = 2.5)
      }

      p_max <- max(data$neg_log10Pvalue, na.rm = TRUE)
      p_mid <- if (is.finite(p_max) && p_max > 0) p_max / 2 else 0.5
      p_max_safe <- if (is.finite(p_max) && p_max > 0) p_max else 1
      col_fun_p <- circlize::colorRamp2(c(0, p_mid, p_max_safe), unname(circle_p_colors))

      circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.1,
        bg.border = "white",
        panel.fun = function(x, y) {
          sector.index <- circlize::get.cell.meta.data("sector.index")
          val <- data[sector.index, "neg_log10Pvalue"]
          circlize::circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = col_fun_p(val), border = NA)
          num <- data[sector.index, "gene_num.max"]
          circlize::circos.text(CELL_META$xcenter, 0.5, paste0(num), cex = 3, col = "black", facing = "bending.inside")
        }
      )

      circlize::circos.track(
        ylim = c(0, 1),
        track.height = 0.1,
        bg.border = "white",
        panel.fun = function(x, y) {
          sector.index <- circlize::get.cell.meta.data("sector.index")
          up <- dplyr::coalesce(as.numeric(data[sector.index, "up.regulated"]), 0)
          down <- dplyr::coalesce(as.numeric(data[sector.index, "down.regulated"]), 0)
          total <- up + down
          if (total > 0) {
            circlize::circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (up / total), 1, col = circle_reg_colors["up"], border = NA)
            circlize::circos.rect(CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (up / total), 0, CELL_META$xlim[2], 1, col = circle_reg_colors["down"], border = NA)
          } else {
            circlize::circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = viz_pal$neutral["light"], border = NA)
          }
        }
      )

      min_score <- min(0, min(data$score, na.rm = TRUE))
      max_score <- max(data$score, na.rm = TRUE)
      circlize::circos.track(
        ylim = c(min_score, max_score),
        track.height = 0.2,
        bg.border = viz_pal$neutral["light"],
        panel.fun = function(x, y) {
          sector.index <- circlize::get.cell.meta.data("sector.index")
          val <- data[sector.index, "score"]
          circlize::circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], val, col = ifelse(val > 0, circle_score_colors["pos"], circle_score_colors["neg"]), border = "white")
          circlize::circos.lines(CELL_META$xlim, c(0, 0), col = viz_pal$neutral["mid"], lty = 2)
        }
      )

      if (analysis_type == "GO" || length(unique(data$category)) > 1) {
        for (cat in unique(data$category)) {
          sectors <- data$id[data$category == cat]
          this_col <- cat_colors[cat]
          circlize::highlight.sector(sectors, track.index = 1, col = paste0(this_col, "20"), text = cat, text.col = this_col, text.vjust = 0, cex = 1.2, font = 2, lwd = 2, border = this_col)
        }
      }

      lgd_p <- ComplexHeatmap::Legend(title = "-Log10(P-value)", col_fun = col_fun_p, at = c(0, p_mid, p_max_safe), labels = sprintf("%.1f", c(0, p_mid, p_max_safe)), direction = "horizontal")
      lgd_reg <- ComplexHeatmap::Legend(title = "Regulation", labels = c("Up", "Down"), legend_gp = grid::gpar(fill = c(circle_reg_colors["up"], circle_reg_colors["down"])), ncol = 2)
      lgd_score <- ComplexHeatmap::Legend(title = if (analysis_type == "GSEA") "NES" else "Rich Factor", labels = c("Positive", "Negative"), legend_gp = grid::gpar(fill = c(circle_score_colors["pos"], circle_score_colors["neg"])), ncol = 2)
      lgd_cat <- ComplexHeatmap::Legend(title = "Category", labels = names(cat_colors), legend_gp = grid::gpar(fill = cat_colors), ncol = 1)
      pd <- ComplexHeatmap::packLegend(lgd_cat, lgd_p, lgd_reg, lgd_score, direction = "vertical", max_height = grid::unit(40, "cm"), column_gap = grid::unit(20, "mm"), row_gap = grid::unit(12, "mm"))
      ComplexHeatmap::draw(pd, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"), just = "center")
      circlize::circos.clear()
    }
  )
}

# Step 06-08 shared enrichment plot helpers live in the main setup so we keep
# one entry file per module and avoid extra helper sidecars.
save_chord_plot <- function(input_list, primary, title, width = 14, height = 12, max_genes_to_show = 20) {
  if (is.null(input_list) || is.null(input_list$circ) || nrow(input_list$circ) == 0) {
    return(invisible(NULL))
  }

  circ <- input_list$circ
  circ$abs_logFC <- abs(circ$logFC)
  circ <- circ[order(circ$abs_logFC, decreasing = TRUE), , drop = FALSE]
  top_genes <- unique(circ$genes)[seq_len(min(length(unique(circ$genes)), max_genes_to_show))]
  circ_filtered <- circ[circ$genes %in% top_genes, , drop = FALSE]
  genedata_filtered <- input_list$genedata[input_list$genedata$ID %in% top_genes, , drop = FALSE]
  chord <- tryCatch(GOplot::chord_dat(data = circ_filtered, genes = genedata_filtered), error = function(e) NULL)
  if (is.null(chord) || nrow(chord) == 0) {
    return(invisible(NULL))
  }
  process_cols <- setdiff(colnames(chord), "logFC")
  chord_ribbon_colors <- if (length(process_cols) > 0) {
    chord_palette_seed <- unname(c(
      group1_color,
      viz_pal$accent["teal"],
      viz_pal$accent["warm"],
      viz_pal$accent["purple"],
      viz_pal$accent["green"],
      viz_pal$accent["peach"],
      viz_pal$accent["yellow"],
      group2_color,
      viz_pal$neutral["dark"],
      viz_pal$neutral["mid"]
    ))
    if (length(process_cols) <= length(chord_palette_seed)) {
      chord_palette_seed[seq_len(length(process_cols))]
    } else {
      grDevices::colorRampPalette(chord_palette_seed)(length(process_cols))
    }
  } else {
    group2_color
  }

  save_base_pdf_versions(
    primary = primary,
    width = width,
    height = height,
    expr = print(
      GOplot::GOChord(
        data = chord,
        title = title,
        space = 0.02,
        gene.order = "logFC",
        gene.space = 0.25,
        gene.size = 3.5,
        lfc.col = c(unname(viz_pal$de["up"]), unname(viz_pal$neutral["light"]), unname(viz_pal$de["down"])),
        ribbon.col = chord_ribbon_colors,
        border.size = 0.25,
        process.label = 10,
        limit = c(0, 0)
      )
    )
  )
}

save_upset_plot <- function(enrich_df, primary) {
  gene_sets <- make_upset_gene_sets(enrich_df)
  if (length(gene_sets) < 2) {
    return(invisible(NULL))
  }
  save_base_pdf_versions(
    primary = primary,
    width = 12,
    height = 8,
    expr = print(
      UpSetR::upset(
        UpSetR::fromList(gene_sets),
        nsets = length(gene_sets),
        nintersects = 20,
        order.by = "freq",
        sets.bar.color = viz_pal$upset["sets"],
        main.bar.color = viz_pal$upset["main"],
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1)
      )
    )
  )
}

save_gocluster_plot <- function(input_list, primary, title, max_genes_to_show = 20) {
  if (is.null(input_list) || is.null(input_list$circ) || nrow(input_list$circ) == 0 || is.null(input_list$genedata) || nrow(input_list$genedata) == 0) {
    return(invisible(NULL))
  }

  circ <- input_list$circ
  circ$abs_logFC <- abs(circ$logFC)
  circ <- circ[order(circ$abs_logFC, decreasing = TRUE), , drop = FALSE]
  top_genes <- unique(circ$genes)[seq_len(min(length(unique(circ$genes)), max_genes_to_show))]
  circ_filtered <- circ[circ$genes %in% top_genes, , drop = FALSE]
  genedata_filtered <- input_list$genedata[input_list$genedata$ID %in% top_genes, , drop = FALSE]
  chord <- tryCatch(GOplot::chord_dat(data = circ_filtered, genes = genedata_filtered), error = function(e) NULL)
  if (is.null(chord) || nrow(chord) == 0) {
    return(invisible(NULL))
  }

  terms_to_cluster <- intersect(unique(input_list$goenrichment$term), colnames(chord))
  if (length(terms_to_cluster) < 2) {
    return(invisible(NULL))
  }

  # GOplot::GOCluster() looks up `chord` in its function environment, so inject
  # the prepared matrix explicitly for stable modular execution.
  gocluster_fun <- GOplot::GOCluster
  environment(gocluster_fun) <- list2env(
    list(chord = chord),
    parent = environment(GOplot::GOCluster)
  )

  save_base_pdf_versions(
    primary = primary,
    width = 16,
    height = 10,
    expr = print(gocluster_fun(chord, terms_to_cluster) + labs(title = title))
  )
}

save_go_family <- function(go_obj, file_paths, title_prefix, category_label) {
  go_df <- make_result_df(go_obj)
  if (nrow(go_df) == 0) {
    return(invisible(NULL))
  }

  save_ora_barplot_df(
    go_df,
    file_paths$bar,
    title_prefix = paste0(title_prefix, " Enrichment"),
    top_n = 15,
    analysis_type = "GO",
    single_fill = viz_pal$go_category[[category_label]]
  )

  save_ora_dotplot_df(
    go_df,
    file_paths$dot,
    title_prefix = paste0(title_prefix, " Enrichment"),
    top_n = 15
  )

  if (length(geneList_symbol_fc) > 0 && nrow(go_df) >= 5) {
    try({
      cnet_plot <- cnetplot(go_obj, foldChange = geneList_symbol_fc, showCategory = 5, node_label = "category")
      save_plot_versions(cnet_plot, file_paths$cnet, width = 10, height = 8)
    }, silent = TRUE)
  }

  go_obj_tree <- tryCatch(pairwise_termsim(go_obj), error = function(e) go_obj)
  if (nrow(make_result_df(go_obj_tree)) > 1) {
    try({
      tree_plot <- enrichplot::treeplot(go_obj_tree, showCategory = min(15, nrow(go_df)), color = "pvalue") +
        labs(title = paste0(title_prefix, " Treeplot")) +
        theme(text = element_text(size = 10))
      save_plot_versions(tree_plot, file_paths$tree, width = 12, height = 10)
    }, silent = TRUE)
  }

  goplot_inputs <- build_goplot_inputs(go_df, de_result, category_label = category_label, top_n = 10)
  save_chord_plot(goplot_inputs, file_paths$chord, paste0(title_prefix, " Chord Plot"))
  save_gocluster_plot(goplot_inputs, file_paths$cluster, paste0(title_prefix, " Term Clustering"))
  save_upset_plot(go_df, file_paths$upset)
}

save_go_all_family <- function(go_obj, out_dir, title_prefix) {
  go_all_df_local <- make_result_df(go_obj)
  if (nrow(go_all_df_local) == 0) {
    return(invisible(NULL))
  }

  save_ora_barplot_df(
    go_all_df_local,
    file.path(out_dir, "7_4_GO_ALL_combined_barplot.pdf"),
    title_prefix = title_prefix,
    top_n = 10,
    facet_col = "ONTOLOGY",
    analysis_type = "GO"
  )

  save_ora_dotplot_df(
    go_all_df_local,
    file.path(out_dir, "7_4_GO_ALL_combined_dotplot.pdf"),
    title_prefix = title_prefix,
    top_n = 10,
    facet_col = "ONTOLOGY"
  )

  go_all_sim <- tryCatch(pairwise_termsim(go_obj), error = function(e) NULL)
  if (!is.null(go_all_sim) && nrow(make_result_df(go_all_sim)) > 1) {
    try({
      go_all_tree <- enrichplot::treeplot(go_all_sim, showCategory = min(20, nrow(go_all_df_local)), color = "pvalue") +
        labs(title = paste(title_prefix, "Treeplot")) +
        theme(text = element_text(size = 10))
      save_plot_versions(go_all_tree, file.path(out_dir, "7_4_GO_ALL_Treeplot.pdf"), width = 12, height = 10)
    }, silent = TRUE)
  }

  go_all_inputs <- build_goplot_inputs(go_all_df_local, de_result, category_label = "ALL", top_n = 5)
  save_chord_plot(go_all_inputs, file.path(out_dir, "7_4_GO_ALL_Chord.pdf"), paste(title_prefix, "Chord Plot"), max_genes_to_show = 20)
  save_upset_plot(go_all_df_local, file.path(out_dir, "7_4_GO_ALL_UpSet.pdf"))
  save_gocluster_plot(go_all_inputs, file.path(out_dir, "7_4_GO_ALL_Cluster.pdf"), paste(title_prefix, "Term Clustering"))
}

save_kegg_family <- function(kegg_obj, out_dir, prefix = "6_1_KEGG", title_prefix = "KEGG Pathway Enrichment") {
  kegg_df <- make_result_df(kegg_obj)
  if (nrow(kegg_df) == 0) {
    return(invisible(NULL))
  }

  save_ora_barplot_df(
    kegg_df,
    file.path(out_dir, paste0(prefix, "_barplot.pdf")),
    title_prefix = title_prefix,
    top_n = 15,
    analysis_type = "ORA",
    single_fill = unname(viz_pal$accent["warm"])
  )

  save_ora_dotplot_df(
    kegg_df,
    file.path(out_dir, paste0(prefix, "_dotplot.pdf")),
    title_prefix = title_prefix,
    top_n = 15
  )

  if (length(geneList_symbol_fc) > 0 && nrow(kegg_df) >= 5) {
    try({
      cnet_plot <- cnetplot(kegg_obj, foldChange = geneList_symbol_fc, showCategory = min(10, nrow(kegg_df)), node_label = "category")
      save_plot_versions(cnet_plot, file.path(out_dir, paste0(prefix, "_cnetplot.pdf")), width = 10, height = 8)
    }, silent = TRUE)
  }

  if (nrow(kegg_df) >= 3) {
    kegg_sim <- tryCatch(pairwise_termsim(kegg_obj), error = function(e) NULL)
    if (!is.null(kegg_sim)) {
      try({
        emap_plot <- enrichplot::emapplot(kegg_sim, showCategory = min(15, nrow(kegg_df)), color = "pvalue") +
          labs(title = paste0(title_prefix, " Network")) +
          theme(text = element_text(size = 10))
        save_plot_versions(
          emap_plot,
          file.path(out_dir, paste0(prefix, "_emapplot.pdf")),
          width = 12,
          height = 10
        )
      }, silent = TRUE)
    }
  }

  kegg_inputs <- build_goplot_inputs(kegg_df, de_result, category_label = "KEGG", top_n = 8)
  save_chord_plot(kegg_inputs, file.path(out_dir, paste0(prefix, "_Chord.pdf")), paste0(title_prefix, " Chord Plot"), max_genes_to_show = 25)
}

save_pathway_ora_family <- function(pathway_obj, out_dir, prefix, title_prefix) {
  pathway_df <- make_result_df(pathway_obj)
  if (nrow(pathway_df) == 0) {
    return(invisible(NULL))
  }

  save_ora_barplot_df(
    pathway_df,
    file.path(out_dir, paste0(prefix, "_barplot.pdf")),
    title_prefix = paste0(title_prefix, " ORA"),
    top_n = 15,
    analysis_type = "ORA",
    single_fill = unname(viz_pal$accent["purple"])
  )

  save_ora_dotplot_df(
    pathway_df,
    file.path(out_dir, paste0(prefix, "_dotplot.pdf")),
    title_prefix = paste0(title_prefix, " ORA"),
    top_n = 15
  )

  pathway_inputs <- build_goplot_inputs(pathway_df, de_result, category_label = title_prefix, top_n = 8)
  save_chord_plot(pathway_inputs, file.path(out_dir, paste0(prefix, "_Chord.pdf")), paste0(title_prefix, " Chord Plot"), max_genes_to_show = 25)
}

save_directional_ora_plot <- function(direction_df, primary, title_prefix, facet_by = c("direction", "source")) {
  facet_by <- match.arg(facet_by)
  if (nrow(direction_df) == 0) {
    return(invisible(NULL))
  }
  plot_df <- build_directional_dotplot_df(direction_df, top_n = 8)
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }
  plot_df$size_value <- if ("Count" %in% colnames(plot_df)) plot_df$Count else if ("setSize" %in% colnames(plot_df)) plot_df$setSize else 1
  plot_df$neg_log10_p <- -log10(dplyr::coalesce(as.numeric(plot_df$pvalue), 1) + 1e-300)
  p <- ggplot(plot_df, aes(x = direction, y = Description, size = size_value, color = neg_log10_p)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(colours = unname(viz_pal$enrich_p)) +
    labs(title = title_prefix, x = NULL, y = NULL, size = "Gene count", color = "-log10(PValue)") +
    theme_rnaseq_bw(base_size = viz_style$compact_size)
  if (facet_by == "source" && "source" %in% colnames(plot_df)) {
    p <- p + facet_grid(source ~ ., scales = "free_y", space = "free_y")
  }
  save_plot_versions(p, primary, width = 10.5, height = 7.5)
}

save_gsea_df_dotplot <- function(gsea_df, primary, title_prefix, top_n = 12, facet_col = NULL) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(invisible(NULL))
  }
  plot_df <- build_gsea_dotplot_df(gsea_df, top_n = top_n)
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }
  plot_df$setSize <- dplyr::coalesce(plot_df$setSize, 1)
  plot_df$neg_log10_p <- -log10(dplyr::coalesce(as.numeric(plot_df$pvalue), 1) + 1e-300)
  p <- ggplot(plot_df, aes(x = NES, y = Description, size = setSize, color = neg_log10_p)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(colours = unname(viz_pal$enrich_p)) +
    labs(title = title_prefix, x = "Normalized Enrichment Score", y = NULL, size = "Set size", color = "-log10(PValue)") +
    theme_rnaseq_bw(base_size = viz_style$compact_size)
  if (!is.null(facet_col) && facet_col %in% colnames(plot_df)) {
    p <- p + facet_grid(stats::as.formula(paste(facet_col, "~ .")), scales = "free_y", space = "free_y")
  }
  save_plot_versions(p, primary, width = 10.5, height = 7.5)
}

save_leading_edge_heatmap <- function(leading_tbl, de_tbl, primary, id_type = c("ENTREZID", "SYMBOL"), top_pathways = 8, top_genes = 40) {
  id_type <- match.arg(id_type)
  mat <- build_leading_edge_heatmap_input(leading_tbl, de_tbl, id_type = id_type, top_pathways = top_pathways, top_genes = top_genes)
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 2) {
    return(invisible(NULL))
  }
  save_base_pdf_versions(
    primary = primary,
    width = 9,
    height = 10,
    expr = ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        mat,
        name = "log2FC",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_names_rot = 45,
        col = heatmap_diverging_col_fun(mat)
      )
    )
  )
}

save_gsea_family <- function(gsea_obj, file_paths, title_prefix, include_single = FALSE) {
  gsea_df <- make_result_df(gsea_obj)
  if (nrow(gsea_df) == 0) {
    return(invisible(NULL))
  }

  try({
    ridge_plot <- enrichplot::ridgeplot(gsea_obj, label_format = 40, fill = "pvalue") +
      labs(title = paste0("GSEA: ", title_prefix))
    save_plot_versions(ridge_plot, file_paths$ridge, width = 12, height = 8)
  }, silent = TRUE)

  save_gsea_df_dotplot(
    gsea_df,
    file_paths$dot,
    title_prefix = paste0("GSEA: ", title_prefix),
    top_n = min(15, nrow(gsea_df)),
    facet_col = if ("ONTOLOGY" %in% colnames(gsea_df)) "ONTOLOGY" else NULL
  )

  gsea_sim <- tryCatch(pairwise_termsim(gsea_obj), error = function(e) NULL)
  if (!is.null(gsea_sim) && nrow(make_result_df(gsea_sim)) > 1) {
    try({
      emap_plot <- enrichplot::emapplot(gsea_sim, showCategory = min(15, nrow(gsea_df)), color = "pvalue") +
        labs(title = paste0("GSEA: ", title_prefix, " Network")) +
        theme(text = element_text(size = 10))
      save_plot_versions(emap_plot, file_paths$emap, width = 12, height = 8)
    }, silent = TRUE)
  }

  if (include_single && nrow(gsea_df) > 0 && !is.na(file_paths$single)) {
    try({
      single_plot <- enrichplot::gseaplot2(gsea_obj, geneSetID = gsea_df$ID[1], title = gsea_df$Description[1])
      save_plot_versions(single_plot, file_paths$single, width = 8, height = 6)
    }, silent = TRUE)
  }
}

rawdata_dir <- file.path(working_dir, "rawdata")
merge_dir <- file.path(working_dir, "3.Merge_result")
de_dir <- file.path(working_dir, "4.DE_analysis")
result_pca_dir <- file.path(working_dir, "result_pca")
post_tx_root_dir <- file.path(working_dir, "5.post_transcriptional_regulation")
legacy_post_tx_root_dir <- file.path(working_dir, "5.AS_analysis")
as_root_dir <- if (dir.exists(post_tx_root_dir)) post_tx_root_dir else legacy_post_tx_root_dir
splicing_dir <- if (dir.exists(file.path(as_root_dir, "splicing"))) {
  file.path(as_root_dir, "splicing")
} else {
  as_root_dir
}

reads_info_path <- file.path(rawdata_dir, "reads_info.txt")
selected_expr_path <- file.path(merge_dir, "genes.DESeq2.normalized_counts.matrix")
selected_logcpm_path <- file.path(merge_dir, "genes.counts.logCPM.matrix")
selected_count_path <- file.path(merge_dir, "genes.counts.matrix")
selected_samples_path <- file.path(de_dir, "samples.txt")
auto_select_path <- file.path(result_pca_dir, "auto_selection_summary.tsv")
full_expr_path <- file.path(result_pca_dir, "full_sample_inputs", "genes.DESeq2.normalized_counts.matrix")
full_logcpm_path <- file.path(result_pca_dir, "full_sample_inputs", "genes.counts.logCPM.matrix")
full_count_path <- file.path(result_pca_dir, "full_sample_inputs", "genes.counts.matrix")
de_path <- file.path(de_dir, "DE_results")
rmats_dir <- if (dir.exists(file.path(splicing_dir, "rmats"))) file.path(splicing_dir, "rmats") else file.path(as_root_dir, "rmats_output")
exon_usage_dir <- if (dir.exists(file.path(splicing_dir, "exon_usage"))) file.path(splicing_dir, "exon_usage") else file.path(as_root_dir, "exon_usage")
isoform_switch_dir <- file.path(as_root_dir, "isoform_switch")
isoform_results_dir <- file.path(isoform_switch_dir, "results")
apa_dir <- file.path(as_root_dir, "apa")
der_dir <- file.path(as_root_dir, "der")
integrated_as_dir <- file.path(as_root_dir, "integrated")

if (!file.exists(reads_info_path)) {
  stop("未找到 rawdata/reads_info.txt")
}
if (!file.exists(selected_expr_path)) {
  stop("未找到 3.Merge_result/genes.DESeq2.normalized_counts.matrix。请先运行 bash rnaseq_first.sh。")
}
if (!file.exists(de_path)) {
  stop("未找到 4.DE_analysis/DE_results。请先运行 bash rnaseq_second.sh。")
}

reads_info <- read_reads_info(reads_info_path)
selected_samples_tbl <- read_samples_file(selected_samples_path)
auto_selection_summary <- if (file.exists(auto_select_path)) {
  read.delim(auto_select_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

gene_exp <- read_matrix(selected_expr_path)
gene_logcpm <- if (file.exists(selected_logcpm_path)) read_matrix(selected_logcpm_path) else NULL
gene_count_mat <- if (file.exists(selected_count_path)) read_matrix(selected_count_path) else NULL

if (!file.exists(full_expr_path)) {
  stop("未找到 result_pca/full_sample_inputs/genes.DESeq2.normalized_counts.matrix。请先运行 bash rnaseq_first.sh。")
}
gene_exp_full <- read_matrix(full_expr_path)
gene_logcpm_full <- if (file.exists(full_logcpm_path)) read_matrix(full_logcpm_path) else NULL
gene_count_full <- if (file.exists(full_count_path)) read_matrix(full_count_path) else NULL

if (!is.null(selected_samples_tbl) && all(selected_samples_tbl$sample %in% colnames(gene_exp))) {
  gene_exp <- gene_exp[, selected_samples_tbl$sample, drop = FALSE]
  if (!is.null(gene_logcpm) && all(selected_samples_tbl$sample %in% colnames(gene_logcpm))) {
    gene_logcpm <- gene_logcpm[, selected_samples_tbl$sample, drop = FALSE]
  }
  if (!is.null(gene_count_mat)) {
    gene_count_mat <- gene_count_mat[, selected_samples_tbl$sample, drop = FALSE]
  }
}

sample_info <- reads_info[colnames(gene_exp), c("group", "stage"), drop = FALSE]
missing_meta <- colnames(gene_exp)[!is.finite(match(colnames(gene_exp), rownames(reads_info))) | is.na(sample_info$group)]
if (length(missing_meta) > 0) {
  stop("样本信息缺失: ", paste(missing_meta, collapse = ", "))
}

full_sample_info <- reads_info[colnames(gene_exp_full), c("group", "stage"), drop = FALSE]
missing_full_meta <- colnames(gene_exp_full)[!is.finite(match(colnames(gene_exp_full), rownames(reads_info))) | is.na(full_sample_info$group)]
if (length(missing_full_meta) > 0) {
  stop("全样本快照的样本信息缺失: ", paste(missing_full_meta, collapse = ", "))
}

analysis_groups <- if (!is.null(selected_samples_tbl) && nrow(selected_samples_tbl) > 0) {
  unique(selected_samples_tbl$group)
} else {
  unique(sample_info$group)
}
if (length(analysis_groups) < 2) {
  stop("当前 samples.txt 或表达矩阵中未解析出两个分组")
}
analysis_groups <- analysis_groups[seq_len(2)]
group1_label <- analysis_groups[[1]]
group2_label <- analysis_groups[[2]]

sample_name_info <- list(
  group1 = if (!is.null(selected_samples_tbl) && nrow(selected_samples_tbl) > 0) {
    selected_samples_tbl$sample[selected_samples_tbl$group == group1_label]
  } else {
    colnames(gene_exp)[sample_info[colnames(gene_exp), "group"] == group1_label]
  },
  group2 = if (!is.null(selected_samples_tbl) && nrow(selected_samples_tbl) > 0) {
    selected_samples_tbl$sample[selected_samples_tbl$group == group2_label]
  } else {
    colnames(gene_exp)[sample_info[colnames(gene_exp), "group"] == group2_label]
  }
)

rna_group_palette <- stats::setNames(c("#A1C9F4", "#FF9F9B"), c(group1_label, group2_label))
group1_color <- unname(rna_group_palette[[group1_label]])
group2_color <- unname(rna_group_palette[[group2_label]])
rna_de_palette <- c("up" = "#FF9F9B", "down" = "#A1C9F4", "ns" = "#D4D9DE")
rna_sig_palette <- c("significant" = "#FF9F9B", "ns" = "#D4D9DE")
rna_neutral_palette <- c(
  "paper" = "#FFFFFF",
  "light" = "#FFFDF8",
  "grid" = "#D4D9DE",
  "mid" = "#D4D9DE",
  "dark" = "#5C6B73",
  "ink" = "#5C6B73"
)
rna_accent_palette <- c(
  "warm" = "#FFB482",
  "green" = "#8DE5A1",
  "purple" = "#B39FDB",
  "yellow" = "#FDFFB6",
  "teal" = "#9ED9CC",
  "peach" = "#FFD3A8"
)
rna_red_yellow_blue_palette <- c("low" = "#A1C9F4", "mid" = "#FDFFB6", "high" = "#FF9F9B")
rna_red_yellow_blue_binary_palette <- c(
  "zero" = unname(rna_red_yellow_blue_palette["low"]),
  "mid" = unname(rna_red_yellow_blue_palette["mid"]),
  "one" = unname(rna_red_yellow_blue_palette["high"])
)
rna_blue_red_mid <- "#DCCFE7"
rna_sashimi_color_string <- paste(c(group1_color, group2_color), collapse = ",")

viz_pal <- list(
  group = rna_group_palette,
  de = rna_de_palette,
  sig_status = rna_sig_palette,
  sig_gradient = c("low" = group1_color, "high" = group2_color),
  neutral = rna_neutral_palette,
  accent = rna_accent_palette,
  overlap = c(
    "DAS_only" = "#A1C9F4",
    "DE_only" = "#FFB482",
    "DE_and_DAS" = "#FF9F9B",
    "APA_only" = "#8DE5A1",
    "APA_and_DE" = "#FF9F9B",
    "DER_only" = "#B39FDB",
    "DER_and_DE" = "#FF9F9B"
  ),
  method_count = c("1" = "#A1C9F4", "2" = "#8DE5A1", "3" = "#FFB482", "4" = "#FF9F9B"),
  candidate_status = c(
    "DE_and_DAS_splicing" = "#FF9F9B",
    "DE_splicing" = "#A1C9F4",
    "DAS_splicing" = "#8DE5A1",
    "Background_splicing" = "#D4D9DE"
  ),
  go_category = c("BP" = "#FF9F9B", "CC" = "#A1C9F4", "MF" = "#B39FDB"),
  enrich_p = rna_red_yellow_blue_palette,
  score = c("pos" = "#FF9F9B", "neg" = "#A1C9F4"),
  upset = c("sets" = "#FF9F9B", "main" = "#A1C9F4"),
  diverging = rna_red_yellow_blue_palette,
  heatmap_diverging = rna_red_yellow_blue_palette,
  heatmap_sequential = rna_red_yellow_blue_palette,
  heatmap_binary = rna_red_yellow_blue_binary_palette,
  heat = colorRampPalette(c("#A1C9F4", "#FDFFB6", "#FF9F9B"))(100),
  transition = colorRampPalette(c("#FFFDF8", "#FDFFB6", "#FFB482", "#FF9F9B"))(100)
)

group_colors <- viz_pal$group
extra_groups <- setdiff(unique(c(sample_info$group, full_sample_info$group)), names(group_colors))
if (length(extra_groups) > 0) {
  extra_cols <- colorRampPalette(viz_pal$accent)(length(extra_groups))
  names(extra_cols) <- extra_groups
  group_colors <- c(group_colors, extra_cols)
}

viz_style <- list(
  base_size = 11.5,
  compact_size = 10,
  line = 0.65,
  point = 3.0,
  alpha = 0.86,
  legend = "right",
  font_family = first_nonempty(Sys.getenv("RNASEQ_FIG_FONT", unset = ""), "sans"),
  axis_line = 0.55,
  grid_line = 0.28,
  plot_margin = c(10, 12, 8, 10)
)

theme_rnaseq_common <- function(base_size = viz_style$base_size, legend = viz_style$legend, grid = TRUE, border = FALSE) {
  grid_line <- if (isTRUE(grid)) {
    element_line(color = alpha(viz_pal$neutral["grid"], 0.85), linewidth = viz_style$grid_line)
  } else {
    element_blank()
  }
  border_line <- if (isTRUE(border)) {
    element_rect(color = alpha(viz_pal$neutral["dark"], 0.38), fill = NA, linewidth = viz_style$axis_line)
  } else {
    element_blank()
  }

  theme(
    text = element_text(family = viz_style$font_family, color = viz_pal$neutral["ink"]),
    plot.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
    panel.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
    plot.title = element_text(hjust = 0, face = "bold", size = base_size * 1.18, color = viz_pal$neutral["ink"], margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0, size = base_size * 0.92, color = viz_pal$neutral["dark"], margin = margin(b = 8)),
    plot.caption = element_text(hjust = 1, size = base_size * 0.78, color = viz_pal$neutral["mid"]),
    plot.margin = do.call(grid::unit, list(viz_style$plot_margin, "pt")),
    axis.title = element_text(size = base_size * 0.96, color = viz_pal$neutral["ink"], margin = margin(t = 4, r = 4)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.text = element_text(size = base_size * 0.86, color = viz_pal$neutral["dark"]),
    axis.ticks = element_line(color = alpha(viz_pal$neutral["dark"], 0.55), linewidth = viz_style$axis_line),
    axis.line = element_line(color = alpha(viz_pal$neutral["dark"], 0.75), linewidth = viz_style$axis_line),
    panel.grid.major = grid_line,
    panel.grid.minor = element_blank(),
    panel.border = border_line,
    strip.background = element_rect(fill = viz_pal$neutral["light"], color = alpha(viz_pal$neutral["grid"], 0.9), linewidth = 0.3),
    strip.text = element_text(size = base_size * 0.88, face = "bold", color = viz_pal$neutral["ink"], margin = margin(t = 5, b = 5)),
    legend.position = legend,
    legend.background = element_rect(fill = viz_pal$neutral["paper"], color = NA),
    legend.key = element_rect(fill = viz_pal$neutral["paper"], color = NA),
    legend.title = element_text(size = base_size * 0.86, face = "bold", color = viz_pal$neutral["ink"]),
    legend.text = element_text(size = base_size * 0.8, color = viz_pal$neutral["dark"]),
    legend.key.size = grid::unit(0.45, "cm"),
    legend.spacing.y = grid::unit(2, "pt"),
    panel.spacing = grid::unit(0.55, "lines")
  )
}

theme_rnaseq_classic <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_classic(base_size = base_size, base_family = viz_style$font_family) +
    theme_rnaseq_common(base_size = base_size, legend = legend, grid = FALSE, border = FALSE)
}

theme_rnaseq_minimal <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_minimal(base_size = base_size, base_family = viz_style$font_family) +
    theme_rnaseq_common(base_size = base_size, legend = legend, grid = TRUE, border = FALSE)
}

theme_rnaseq_bw <- function(base_size = viz_style$base_size, legend = viz_style$legend) {
  theme_bw(base_size = base_size, base_family = viz_style$font_family) +
    theme_rnaseq_common(base_size = base_size, legend = legend, grid = TRUE, border = TRUE)
}

heatmap_group_palette <- function(groups, palette = group_colors) {
  group_levels <- unique(as.character(groups))
  resolve_named_colors(palette, group_levels)
}

heatmap_group_annotation <- function(groups, palette = group_colors, annotation_name = "Group", which = c("column", "row"), ...) {
  which <- match.arg(which)
  group_colors_local <- heatmap_group_palette(groups, palette = palette)
  annotation_data <- setNames(list(as.character(groups)), annotation_name)
  annotation_colors <- setNames(list(group_colors_local), annotation_name)

  if (identical(which, "row")) {
    return(do.call(ComplexHeatmap::rowAnnotation, c(annotation_data, list(col = annotation_colors), list(...))))
  }
  do.call(ComplexHeatmap::HeatmapAnnotation, c(annotation_data, list(col = annotation_colors), list(...)))
}

heatmap_diverging_col_fun <- function(x = NULL, midpoint = 0, breaks = NULL, colors = unname(viz_pal$heatmap_diverging)) {
  if (!is.null(breaks)) {
    return(circlize::colorRamp2(breaks, colors))
  }

  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(circlize::colorRamp2(c(-1, midpoint, 1), colors))
  }

  low_value <- min(values, na.rm = TRUE)
  high_value <- max(values, na.rm = TRUE)
  if (!is.finite(low_value) || !is.finite(high_value) || low_value == high_value) {
    low_value <- midpoint - 1
    high_value <- midpoint + 1
  } else if (midpoint <= low_value || midpoint >= high_value) {
    midpoint <- mean(c(low_value, high_value))
  }

  circlize::colorRamp2(c(low_value, midpoint, high_value), colors)
}

heatmap_sequential_col_fun <- function(x = NULL, breaks = NULL, colors = unname(viz_pal$heatmap_sequential)) {
  if (!is.null(breaks)) {
    return(circlize::colorRamp2(breaks, colors))
  }

  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(circlize::colorRamp2(c(0, 0.5, 1), colors))
  }

  low_value <- min(values, na.rm = TRUE)
  high_value <- max(values, na.rm = TRUE)
  if (!is.finite(low_value) || !is.finite(high_value) || low_value == high_value) {
    low_value <- 0
    high_value <- max(1, high_value)
  }
  mid_value <- mean(c(low_value, high_value))
  circlize::colorRamp2(c(low_value, mid_value, high_value), colors)
}

heatmap_binary_col_fun <- function() {
  circlize::colorRamp2(c(0, 0.5, 1), unname(viz_pal$heatmap_binary))
}

scale_fill_heatmap_diverging <- function(midpoint = 0, limits = NULL, ...) {
  ggplot2::scale_fill_gradient2(
    low = unname(viz_pal$heatmap_diverging["low"]),
    mid = unname(viz_pal$heatmap_diverging["mid"]),
    high = unname(viz_pal$heatmap_diverging["high"]),
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

scale_fill_heatmap_sequential <- function(...) {
  ggplot2::scale_fill_gradientn(
    colors = unname(viz_pal$heatmap_sequential),
    ...
  )
}

scale_color_heatmap_sequential <- function(...) {
  ggplot2::scale_color_gradientn(
    colors = unname(viz_pal$heatmap_sequential),
    ...
  )
}

scale_color_heatmap_diverging <- function(midpoint = 0, limits = NULL, ...) {
  ggplot2::scale_color_gradient2(
    low = unname(viz_pal$heatmap_diverging["low"]),
    mid = unname(viz_pal$heatmap_diverging["mid"]),
    high = unname(viz_pal$heatmap_diverging["high"]),
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

scale_fill_heatmap_binary <- function(limits = c(0, 1), ...) {
  ggplot2::scale_fill_gradientn(
    colors = unname(viz_pal$heatmap_binary),
    values = c(0, 0.5, 1),
    limits = limits,
    ...
  )
}

analysis_context <- data.frame(
  metric = c(
    "rnaseq_species",
    "orgdb_package",
    "kegg_organism_code",
    "reactome_organism",
    "msigdb_species",
    "selected_samples_n",
    "selected_genes_n",
    "selected_logcpm_available",
    "full_snapshot_available",
    "full_snapshot_samples_n",
    "full_logcpm_available",
    "auto_selection_available",
    "de_results_available",
    "rmats_dir_available",
    "exon_usage_available",
    "isoform_switch_available",
    "integrated_as_available",
    "tf_regulator_module_available",
    "gsva_module_available",
    "ppi_module_available"
  ),
  value = c(
    rnaseq_species,
    orgdb_package,
    kegg_organism_code,
    reactome_organism,
    msigdb_species,
    ncol(gene_exp),
    nrow(gene_exp),
    !is.null(gene_logcpm),
    file.exists(full_expr_path),
    ncol(gene_exp_full),
    !is.null(gene_logcpm_full),
    !is.null(auto_selection_summary),
    file.exists(de_path),
    dir.exists(rmats_dir),
    dir.exists(exon_usage_dir),
    dir.exists(isoform_results_dir),
    dir.exists(integrated_as_dir),
    file.exists(file.path(get0("rmd_dir", ifnotfound = file.path(working_dir, "rmd"), inherits = TRUE), "10_tf_regulator.Rmd")),
    file.exists(file.path(get0("rmd_dir", ifnotfound = file.path(working_dir, "rmd"), inherits = TRUE), "09_gsva.Rmd")),
    file.exists(file.path(get0("rmd_dir", ifnotfound = file.path(working_dir, "rmd"), inherits = TRUE), "11_ppi.Rmd"))
  ),
  stringsAsFactors = FALSE
)

ensure_module_result_layout("01_setup")
ensure_dirs(cache_root_dir)
write.csv(analysis_context, file.path(result_01_dir, "01_analysis_context.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), route_output_path(file.path(result_01_dir, "01_R_sessionInfo.txt"), preferred = "table"))
