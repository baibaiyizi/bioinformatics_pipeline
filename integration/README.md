# 项目级多组学整合模板

阅读和分析顺序：项目全景、共同生物学程序、核心器官多层证据、细胞类型背景、两两分析支撑、证据依赖。

复制本模板到项目物种目录，填写 `config/project.yml` 的真实输入、来源资格、组别及冻结基因集；项目运行自己的本地代码。`INPUT_REQUIRED` 必须替换为真实已核验文件，缺失保持 MISSING，不能通过猜测样本配对补齐。

- `bash run.sh --check`：核验样本身份、方向、RNA 已完成生产凭据和输入签名。
- `bash run.sh --integrate`：从冻结单组学统计和已完成联合分析数值源生成项目整合并重绘。
- `bash run.sh --plots-only`：只从已签名数值源重绘；禁止重新运行差异、富集或覆盖构建。

深化模块在同一次原子执行中生成三组独立证据：七情境逐文库 ssGSEA 与完整 DE statistic fgsea、RNA–ATAC 完整测试域分层置换及共享位点、转录后方法测试域/队列签名/方法特异 ORA。项目必须显式声明所有输入；WGBS 仅保留 selected DMR 与覆盖资格，不能扩写为三层联合 P 值。不同样本签名的 selected/full 队列不计算同队列交集，event/region FDR 不能改称 gene FDR。

结果按 `00_audit`、`01_project_atlas`、`02_program_integration`、`03_liver_multiomics`、`04_celltype_context`、`05_pairwise_support`、`06_synthesis` 编号，各方法下分 `plots/` 与 `tables/`。同名 PNG/PDF 并列，183 mm、Liberation Sans 8/7 pt、300 dpi 白底 RGB。分面、UpSet 与组合图保留。

项目仅使用一份冻结 figure profile。WT/PFOS/NS、方向、连续尺度和无方向分类各有独立语义；同名方法、术语或簇跨图使用稳定身份色，显著性和方向另用形状、大小或专用色道表达。每张图保留源表及“生物学意义 / 数据分布特征 / 最终结论”三项解释。

本模板需要真实冻结来源及逐图清单；复制代码或通过合成 fixture 检查不等于真实数据重现。计算过程先在临时结果目录完成，来源改变或成图失败不覆盖当前有效结果。

`--integrate` 的前提是已有现役冻结结果合同：`result/00_audit/run.json` 为 `COMPLETED`，并带入其签名对应的 `source_catalog.tsv`、`figure_manifest.tsv`、`membership.tsv`、`term_id_map.tsv` 及真实数值源。四张 TSV 位于 `result/00_audit/tables/`。这些来自已经完成的实际联合分析；只填写配置不能从空 `result/` 启动。`--check` 会具体列出缺少的合同文件为 `MISSING`；输入/合同检查通过仍不表示已经完成真实整合。

不合并 P、不设计任意加权总分；同一程序不同成员与同一基因多层收敛分别报告。样本未配对时不运行样本级联合模型；细胞不是额外生物学重复。
