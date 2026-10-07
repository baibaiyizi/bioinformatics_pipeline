# LCMS_template 分析报告

比较方向：treated - control；统计单位：独立生物学样本。
默认筛选：nominal；原始 P 值候选为探索性，完整结果保留 BH/FDR。
输入模式：peak_table；QC身份：unverified；漂移校正：none。
未提供填峰前矩阵时，原始检出率及其过滤为 not_available，不用填峰值代替。
注释证据不足时保留未鉴定特征；level 1 要求标准品确认。
上游外部峰表并不证明本流程已执行 XCMS；真实 mzML 上游须独立验证。

## 表格与执行状态

- [01_audit/effective_config.tsv](01_audit/effective_config.tsv)
- [01_audit/external_tools.tsv](01_audit/external_tools.tsv)
- [01_audit/input_manifest.tsv](01_audit/input_manifest.tsv)
- [01_audit/inventory_summary.csv](01_audit/inventory_summary.csv)
- [01_audit/mzml_manifest.csv](01_audit/mzml_manifest.csv)
- [01_audit/sample_info_used.tsv](01_audit/sample_info_used.tsv)
- [01_audit/software_versions.tsv](01_audit/software_versions.tsv)
- [02_peak_table/abundance_matrix_detected.csv](02_peak_table/abundance_matrix_detected.csv)
- [02_peak_table/abundance_matrix_filled.csv](02_peak_table/abundance_matrix_filled.csv)
- [02_peak_table/feature_annotation.csv](02_peak_table/feature_annotation.csv)
- [02_peak_table/peak_source.csv](02_peak_table/peak_source.csv)
- [03_qc/d_ratio_summary.csv](03_qc/d_ratio_summary.csv)
- [03_qc/feature_filter_decision.csv](03_qc/feature_filter_decision.csv)
- [03_qc/matrix_imputed.csv](03_qc/matrix_imputed.csv)
- [03_qc/matrix_log2.csv](03_qc/matrix_log2.csv)
- [03_qc/matrix_normalized.csv](03_qc/matrix_normalized.csv)
- [03_qc/matrix_scaled.csv](03_qc/matrix_scaled.csv)
- [03_qc/normalization_factors.csv](03_qc/normalization_factors.csv)
- [03_qc/pca_scores_POS.csv](03_qc/pca_scores_POS.csv)
- [03_qc/qc_rsd_summary.csv](03_qc/qc_rsd_summary.csv)
- [03_qc/rla_summary.csv](03_qc/rla_summary.csv)
- [03_qc/rt_correction_summary.csv](03_qc/rt_correction_summary.csv)
- [03_qc/sample_info_used.tsv](03_qc/sample_info_used.tsv)
- [03_qc/sample_qc_summary.csv](03_qc/sample_qc_summary.csv)
- [03_qc/sample_signal_summary.csv](03_qc/sample_signal_summary.csv)
- [03_qc/spectra_summary.csv](03_qc/spectra_summary.csv)
- [04_statistics/treated_vs_control_company_comparison.csv](04_statistics/treated_vs_control_company_comparison.csv)
- [04_statistics/treated_vs_control_differential_features.csv](04_statistics/treated_vs_control_differential_features.csv)
- [04_statistics/treated_vs_control_exploratory_candidates.csv](04_statistics/treated_vs_control_exploratory_candidates.csv)
- [04_statistics/treated_vs_control_robust_differential_features.csv](04_statistics/treated_vs_control_robust_differential_features.csv)
- [04_statistics/treated_vs_control_summary.csv](04_statistics/treated_vs_control_summary.csv)
- [05_annotation/annotation_full.csv](05_annotation/annotation_full.csv)
- [module_status.tsv](module_status.tsv)

## 图形

![plots/pca_POS.png](plots/pca_POS.png)
![plots/sample_log2_intensity_boxplot.png](plots/sample_log2_intensity_boxplot.png)
![plots/sample_rla_summary.png](plots/sample_rla_summary.png)
![plots/sample_run_order_total_signal.png](plots/sample_run_order_total_signal.png)
![plots/treated_vs_control_top_feature_boxplots.png](plots/treated_vs_control_top_feature_boxplots.png)
![plots/treated_vs_control_volcano.png](plots/treated_vs_control_volcano.png)
