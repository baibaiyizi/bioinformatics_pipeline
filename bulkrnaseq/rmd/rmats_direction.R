# rMATS 原始差值始终保留为 b1 - b2；下游统一使用 group2 - group1。
read_rmats_direction <- function(rmats_dir, samples_file, group1, group2) {
  if (length(unique(c(group1, group2))) != 2L || any(!nzchar(c(group1, group2)))) stop("必须显式指定两个分析组")
  metadata_file <- file.path(rmats_dir, "run_metadata.tsv")
  if (!file.exists(metadata_file)) stop("缺少 rMATS run_metadata.tsv，不能推断 b1/b2 方向。")
  lines <- readLines(metadata_file, warn = FALSE)
  keys <- sub("=.*$", "", lines)
  values <- sub("^[^=]*=", "", lines)
  if (anyDuplicated(keys)) stop("rMATS run_metadata.tsv 存在重复字段。")
  metadata <- stats::setNames(values, keys)
  required <- c("group1", "group2", "group1_samples", "group2_samples", "sample_key")
  if (!all(required %in% keys) || any(!nzchar(metadata[required]))) {
    stop("rMATS 元数据缺少分组或样本身份字段，不能校准方向。")
  }
  raw_groups <- unname(metadata[c("group1", "group2")])
  if (!setequal(raw_groups, c(group1, group2)) || anyDuplicated(raw_groups)) {
    stop("rMATS b1/b2 组别与当前显式对比不匹配。")
  }
  raw_samples <- lapply(metadata[c("group1_samples", "group2_samples")], function(x) strsplit(x, ",", fixed = TRUE)[[1]])
  names(raw_samples) <- raw_groups
  samples <- read.delim(samples_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(samples) != 2L || anyDuplicated(samples[[2]]) || !setequal(unique(samples[[1]]), c(group1, group2))) {
    stop("samples.txt 必须为不重复的显式两组两列样本表。")
  }
  for (g in c(group1, group2)) {
    if (!setequal(raw_samples[[g]], samples[[2]][samples[[1]] == g]) || anyDuplicated(raw_samples[[g]])) {
      stop("rMATS 元数据与当前样本身份不一致：", g)
    }
  }
  if (!identical(paste(sort(unlist(raw_samples, use.names = FALSE)), collapse = ","), unname(metadata[["sample_key"]]))) {
    stop("rMATS sample_key 与 b1/b2 样本不一致。")
  }
  mats_files <- sort(list.files(rmats_dir, pattern = "^(fromGTF\\.)?(SE|A5SS|A3SS|RI|MXE)\\.MATS\\.(JC|JCEC)\\.txt$", full.names = TRUE))
  if (!length(mats_files)) stop("没有可校准的原始 rMATS MATS 表。")
  list(
    version = "rmats_direction_v1", contrast = paste(group2, group1, sep = " - "),
    group1_label = group1, group2_label = group2,
    b1_group = raw_groups[[1]], b2_group = raw_groups[[2]],
    multiplier = if (raw_groups[[1]] == group2) 1 else -1,
    group1_samples = unname(raw_samples[[group1]]), group2_samples = unname(raw_samples[[group2]]),
    input_md5 = tools::md5sum(c(metadata_file, samples_file, mats_files))
  )
}

orient_rmats_table <- function(tbl, design) {
  required <- c("IncLevel1", "IncLevel2", "IncLevelDifference")
  if (!all(required %in% names(tbl))) stop("原始 rMATS 表缺少 PSI 或差值列。")
  # 只新增分析列，不能重写 IncLevel1/2、IncLevelDifference、PValue 或 FDR。
  tbl$rmats_b1_group <- rep(design$b1_group, nrow(tbl))
  tbl$rmats_b2_group <- rep(design$b2_group, nrow(tbl))
  tbl$delta_psi_contrast <- rep(design$contrast, nrow(tbl))
  tbl$delta_psi <- suppressWarnings(as.numeric(tbl$IncLevelDifference)) * design$multiplier
  tbl$psi_group1_values <- if (design$b1_group == design$group1_label) tbl$IncLevel1 else tbl$IncLevel2
  tbl$psi_group2_values <- if (design$b1_group == design$group2_label) tbl$IncLevel1 else tbl$IncLevel2
  for (g in c("group1", "group2")) {
    values <- as.character(tbl[[paste0("psi_", g, "_values")]])
    lengths <- lengths(strsplit(values, ",", fixed = TRUE))
    if (any(lengths != length(design[[paste0(g, "_samples")]]))) {
      stop("rMATS 每事件 PSI 列数与元数据重复数不一致：", g)
    }
  }
  tbl
}

# 热图/样本矩阵必须保留 NA 所在重复的位置；删除 NA 会使后续样本左移。
split_rmats_psi <- function(x) {
  suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
}
