# 复用既有敏感性绘图；输入只读冻结的两张表，不拟合DE。
full_de_result <- read_source_csv("05_de/tables/5_10_full_sample_sensitivity_DE_results.csv")
compare_tbl <- read_source_csv("05_de/tables/5_12_main_vs_full_de_concordance.csv")
main_sig_ids <- compare_tbl$ensembl_id[compare_tbl$main_direction != "ns"]
full_sig_ids <- compare_tbl$ensembl_id[compare_tbl$full_direction != "ns"]
    concordance_plot <- compare_tbl %>%
      filter(is.finite(main_log2FoldChange), is.finite(full_log2FoldChange)) %>%
      mutate(status = case_when(
        ensembl_id %in% intersect(main_sig_ids, full_sig_ids) ~ "both_nominal_DE",
        ensembl_id %in% main_sig_ids ~ "main_only",
        ensembl_id %in% full_sig_ids ~ "full_only",
        TRUE ~ "not_nominal_DE"
      )) %>%
      ggplot(aes(x = main_log2FoldChange, y = full_log2FoldChange, color = status)) +
      geom_point(alpha = 0.55, size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = viz_pal$neutral["dark"]) +
      geom_vline(xintercept = 0, linetype = "dashed", color = viz_pal$neutral["dark"]) +
      scale_color_manual(values = c(
        both_nominal_DE = unname(viz_pal$sig_status["significant"]),
        main_only = unname(viz_pal$accent["warm"]),
        full_only = unname(viz_pal$accent["green"]),
        not_nominal_DE = unname(viz_pal$neutral["mid"])
      )) +
      labs(title = "Main versus full-sample log2FC concordance", x = "Main selected-sample log2FC", y = "Full-snapshot log2FC") +
      theme_rnaseq_publication()

    save_plot_versions(
      concordance_plot,
      primary = "result/05_de/5_13_main_vs_full_log2fc_concordance.pdf",
      width = 7.5,
      height = 6
    )
if (exists("full_de_result") && is.data.frame(full_de_result) && nrow(full_de_result) > 0) {
  full_de_plot_tbl <- full_de_result
  full_de_plot_tbl$figure_significance <- figure_significance_values(
    full_de_plot_tbl,
    pvalue_candidates = "pvalue",
    fdr_candidates = "padj"
  )
  full_de_plot_tbl <- full_de_plot_tbl %>%
    mutate(
      neg_log10_p = -log10(dplyr::coalesce(figure_significance, 1) + 1e-300),
      figure_direction = factor(
        case_when(
          !is.finite(figure_significance) | figure_significance >= de_pvalue_cutoff |
            !is.finite(log2FoldChange) | abs(log2FoldChange) < de_log2fc_cutoff ~ "ns",
          log2FoldChange >= de_log2fc_cutoff ~ "up",
          TRUE ~ "down"
        ),
        levels = c("up", "down", "ns")
      )
    )

  full_volcano_count_tbl <- full_de_plot_tbl %>%
    count(direction = figure_direction, name = "gene_count") %>%
    filter(direction %in% c("up", "down")) %>%
    mutate(
      x = ifelse(
        direction == "up",
        max(full_de_plot_tbl$log2FoldChange, na.rm = TRUE) * 0.98,
        min(full_de_plot_tbl$log2FoldChange, na.rm = TRUE) * 0.98
      ),
      y = max(full_de_plot_tbl$neg_log10_p, na.rm = TRUE) * 0.98,
      hjust = ifelse(direction == "up", 1, 0),
      label = ifelse(direction == "up", paste0("Up: ", gene_count), paste0("Down: ", gene_count))
    )

  full_volcano <- ggplot(full_de_plot_tbl, aes(x = log2FoldChange, y = neg_log10_p, color = figure_direction)) +
    geom_point(alpha = 1, size = 1) +
    geom_vline(xintercept = c(-de_log2fc_cutoff, de_log2fc_cutoff), linetype = "dashed", color = viz_pal$neutral["dark"]) +
    geom_hline(yintercept = -log10(de_pvalue_cutoff), linetype = "dashed", color = viz_pal$neutral["dark"]) +
    geom_text(
      data = full_volcano_count_tbl,
      aes(x = x, y = y, label = label, color = direction, hjust = hjust),
      inherit.aes = FALSE,
      size = rnaseq_figure_text_size,
      show.legend = FALSE
    ) +
    scale_color_manual(values = viz_pal$de, drop = FALSE, name = "Direction") +
    labs(title = "Full-sample sensitivity analysis: volcano plot", x = paste0(group2_label, " vs ", group1_label, " log2FC"), y = figure_significance_axis_label()) +
    theme_rnaseq_publication()
  save_plot_versions(full_volcano, "result/05_de/5_14_full_sample_sensitivity_volcano.pdf", width = 8.5, height = 6)

  full_ma_plot <- ggplot(full_de_plot_tbl, aes(x = log10(baseMean + 1), y = log2FoldChange, color = direction)) +
    geom_point(alpha = 1, size = 1) +
    geom_hline(yintercept = c(-de_log2fc_cutoff, 0, de_log2fc_cutoff), linetype = c("dashed", "solid", "dashed"), color = viz_pal$neutral["dark"]) +
    scale_color_manual(values = viz_pal$de, drop = FALSE) +
    labs(title = "Full-sample sensitivity analysis: MA plot", x = "log10(baseMean + 1)", y = paste0(group2_label, " vs ", group1_label, " log2FC")) +
    theme_rnaseq_publication()
  save_plot_versions(full_ma_plot, "result/05_de/5_15_full_sample_sensitivity_MA.pdf", width = 8.5, height = 6)
}
