# bulk RNA-seq 模板

## 1. 适用范围与统计单位

复制整个目录后只改 `config/project.env` 和 `config/samples.tsv`。当前主分析限定小鼠、独立两组、每个独立动物或混池一列；每组至少 2 个生物学重复。混池整体计为一个重复。重复 `biological_id`、未知身份、配对、重复测量或非空 batch/subject/pair 列明确拒绝。技术文库应在进入本模板前按同一生物学样本合并，不能当作增加 n。

比较固定为 `RNASEQ_GROUP2 - RNASEQ_GROUP1`，样本表行序不定义方向。默认保留全部 `eligible=true` 样本，PCA 只诊断。`eligible=false` 是分析前记录的纳入决定，不能按 PCA 分离度调整。

## 2. 配置与输入

`project.env` 是可信本地 Bash 参数文件，四个入口自动读取；也可用 `RNASEQ_CONFIG=/path/project.env` 指定另一份。文件内的路径以复制后的项目根目录为基准；工具写可执行文件名或绝对路径。所有外部工具和参考必须已安装、预先准备；流程不自动安装包或升级环境。

样本表必须包含 `sample, group, biological_id, stage, eligible, raw_fq1, raw_fq2`，制表符分隔、有表头，sample 和生物学单位唯一。stage 未知可留空，不按目录或样本名推断。FASTQ 模式要求两端文件；counts 模式允许 FASTQ 字段为空。

- `RNASEQ_INPUT_MODE=fastq`：显式填写 HISAT2 索引前缀、GTF、建库链方向和双端 FASTQ。原文件名可与 sample 不同；清洗输出始终为 `fastqc/<sample>_1.fastq.gz` 和 `_2.fastq.gz`，mapping/rMATS/Salmon 共用此命名。
- `RNASEQ_INPUT_MODE=counts`：`RNASEQ_COUNTS_FILE` 指向 gene_id × sample 的 TSV，第一列为唯一基因 ID，其他列必须与全部合格样本精确匹配。仅接受非负整数；不接受 TPM/FPKM/归一化值。DE 报告使用 Ensembl gene ID 注释。
- 所有检验默认原始 P < 0.05、|log2FC| ≥ 1，属于**探索性筛选**；DESeq2 原始 P 和 BH padj 一并保留。`RNASEQ_DE_USE_PADJ=true` 可改为 padj 筛选。rMATS FDR 是事件层校正，不能与基因 DE 的 BH/FDR 混称。

## 3. 运行入口

```bash
bash rnaseq_first.sh --help
bash rnaseq_first.sh --list
bash rnaseq_first.sh --preflight
# FASTQ 模式
bash rnaseq_first.sh --all
bash rnaseq_second.sh --all
# counts 模式直接从第二阶段开始
bash rnaseq_second.sh --all
bash run_rmd.sh
```

第一阶段保留 fastp → HISAT2 → featureCounts → 合并/归一化 → 全部合格样本 → 上游审计。第二阶段导入显式 counts 或使用上游矩阵，运行 DESeq2。`run_rmd.sh` 默认模块 01–05，输出输入记录、QC、PCA、DE 与 HTML 报告。三个 Shell 的 `--help/--list`、四个入口的 `--preflight` 不生成分析结果。预检查检查身份、设计和必要输入存在性，不表示 FASTQ/参考内容已完成全量校验。

默认保留既有结果。配置、样本身份、counts 文件 SHA256 或 FASTQ 大小/mtime 改变后，输入签名不匹配则报 `STALE`；无签名旧输出报 `UNVERIFIED`。请在新复制的模板目录重新分析，不应手工删除签名以承认旧结果。`run_rmd.sh --clean --modules 5` 仅明确清理选中模块；不会绕过签名。增加适用模块不改变主分析方向或纳入集。

## 4. 可选模块与边界

扩展模块加入 `RNASEQ_ENABLED_MODULES` 后运行 `bash run_rmd.sh --modules ...`。没有启用的模块会写 `skipped`；显式请求未启用依赖会失败。

| 模块 | 输入与依赖 |
|---|---|
| 06–08 富集 | 模块 05；GO/KEGG/Reactome 与 MSigDB mouse-native 主结果；human ortholog 单独作为补充并记录数据库/版本/映射来源；需要相应数据库访问或已有有效缓存。 |
| 09–10 通路/TF活性 | 06 富集缓存；10 还依赖 09；只推断表达代理量。 |
| 11 PPI | 配置 `STRING_INPUT_DIR`，显式提供 STRING 文件。 |
| 12 WGCNA | 源模块保留样本数资格门槛；小样本不可解释为稳定网络。 |
| 13 APA / 14 DER / 15 剪接 | 先运行 `rnaseq_third.sh` 的对应方法；下游还需 06。APA 完整输入检查和原缓存保护保留。 |
| 16 免疫表达评分 | 06；仅组织/参考适用时设置 `RNASEQ_IMMUNE_APPLICABLE=true`；分数不等于免疫细胞比例。 |
| 17 整合 / 18 圈图 | 17 需要 09/13–16 的适用结果，18 需要 06；不会把未启用模块伪装成阴性结果。 |

`rnaseq_third.sh` 默认只启用 rMATS，额外 APA/DER 关闭；方法由 `RNASEQ_THIRD_METHODS`、`RNASEQ_THIRD_EXTRA_METHODS` 显式决定。启用 rMATS 后必须保留 `run_metadata.tsv` 的 b1/b2 组别与样本身份。原始 `IncLevel1/2`、`IncLevelDifference`、PValue/FDR 不改写，新增 `delta_psi` 统一为 group2 − group1；PSI 中的 NA 保留原样本位置。元数据缺失、身份不匹配或 PSI 长度不匹配明确失败。

## 5. 结果与报告

主要表在 `result/NN_*/tables/`，图在 `plots/` 或对应子模块。`result/index.html` 汇总真实现存图表；只有表、只有跳过状态时也可生成。每张图提供方法/读图/限制；报告中的图表存在性不能替代 `result/module_status.tsv` 的执行状态。`--report-only` 仅刷新报告，不宣称分析已重新运行。

`bash run_rmd.sh --plots-only` 为已有完整小鼠归档结果提供纯绘图入口。使用前须在本模块冻结 `rmd/plot_inputs.tsv`（原图表达、counts、logCPM矩阵）与 `rmd/plot_samples.tsv`，并保留对应原表和缓存；新复制的空模板不具备这些来源。该入口复用现役Rmd的绘图片段，禁止DE、富集和评分重算；PCA、相关性、热图聚类、网络布局及描述性拟合属于展示派生。先在临时目录成图，核对所有原图、数值与源码SHA并生成报告后再发布，记录独立的 `result/00_audit/tables/plot_receipt.json`。缺少来源直接失败。原图cohort与后来完整cohort必须分别记录，不能互换。

本地 `rmd/figure_style.R`、`plot_heatmaps.R` 统一183 mm、Liberation Sans、8/7 pt、白底RGB 300 dpi PNG与嵌入字体PDF，保留分面、组合、UpSet及原筛选范围。项目运行本地冻结副本，模板来源以 `template_snapshot.tsv` 核验。历史执行记录和输入签名不因重新成图而改写。

PCA：生物学意义为样本整体表达差异的描述；数据分布特征以 PC 解释方差和样本位置为准；最终结论只能是表达结构诊断，不能据此自动选择样本或证明组间效应。

DE 火山图/MA 图：生物学意义为显式两组表达关联；数据分布特征来自实际 log2FC、P 与平均表达；最终结论限于当前设计与阈值的候选基因，原始 P 筛选需标注探索性并同时查看 padj。

## 6. 来源与依赖

`provenance.tsv` 保存原模板和真实项目来源文件的 SHA256、提取内容。真实项目只读；模板运行不读取任何项目脚本或结果。保留原三段流程结构，提取当前共享方向审计、mouse-native/human-ortholog 来源、固定随机种子、183 mm/Cairo PDF 与 300 dpi PNG 绘图和报告逻辑；没有加入旧入口兼容层或自动迁移。

参考现有服务器 `rnaseq_legacy` 环境规划，使用已有 fastp/HISAT2/samtools/Rsubread/DESeq2 等依赖；运行时缺包明确失败。可选模块的依赖仅在启用相应模块时加载。核心绘图与注释需要 ggplot2、knitr、Cairo、ComplexHeatmap、org.Mm.eg.db 等，报告需要 htmltools/stringr。详细实际验证见 [内部验证记录](docs/VALIDATION.md)。

## 7. 验证界限

小数据闭环、语法检查、方向与输入失败场景分别记录。FASTQ 全量上游、联网富集、APA/DER、rMATS 生物数据、PPI/WGCNA/免疫等可选模块本轮没有执行；不能将源码迁移或合成数据验证称为真实数据复现。使用者仍需依据建库、组织、参考版本和生物学单位判断适用性。

报告写法见[统一报告规范](../docs/reporting.md)。
