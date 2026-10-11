compute_atac_diff_peak_flags <- function(
    pvalue,
    fdr,
    fold,
    significance_mode = "pvalue",
    primary_p_cutoff = 0.05,
    fdr_cutoff = 0.05,
    log2fc_cutoff = 0.585,
    primary_require_fc = FALSE) {
  significance_mode <- tolower(as.character(significance_mode)[1])
  if(!significance_mode %in% c("pvalue", "padj", "fdr")) {
    stop("significance_mode must be pvalue, padj, or fdr")
  }
  pvalue <- suppressWarnings(as.numeric(pvalue))
  fdr <- suppressWarnings(as.numeric(fdr))
  fold <- suppressWarnings(as.numeric(fold))
  pvalue_diff <- is.finite(pvalue) & pvalue >= 0 & pvalue < primary_p_cutoff
  fdr_diff <- is.finite(fdr) & fdr >= 0 & fdr <= 1 & fdr <= fdr_cutoff
  fold_pass <- !is.na(fold) & abs(fold) >= log2fc_cutoff
  selected_base <- if(significance_mode == "pvalue") pvalue_diff else fdr_diff
  selected_value <- if(significance_mode == "pvalue") pvalue else fdr
  selected_threshold <- if(significance_mode == "pvalue") primary_p_cutoff else fdr_cutoff
  data.frame(
    PvalueDiffPeak = pvalue_diff,
    FDRDiffPeak = fdr_diff,
    PrimaryDiffPeak = selected_base & (!isTRUE(primary_require_fc) | fold_pass),
    SelectedSignificanceMetric = if(significance_mode == "pvalue") "pvalue" else "FDR",
    SelectedSignificanceValue = selected_value,
    SelectedSignificanceThreshold = selected_threshold,
    stringsAsFactors = FALSE
  )
}
