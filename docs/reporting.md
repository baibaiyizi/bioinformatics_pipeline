# RNA、WGBS 与 ATAC 统一报告规范

## 公共写作规则

人类阅读入口统一为 `projects/<study>/<species>/<assay>/<module>/README.md`。一份 README 同时说明输入、运行入口、真实结果和证据边界；内部验证、依赖说明收进 `docs/`。跨组织、跨阶段和跨组学结果由同物种 `integration/` 解释，单组学不能借用其他组学结论。

使用中文，按结果模块及图表组织，先写实际图表信息，再写所支持的结论等级。需要自动摘要时用 `AUTO_CURRENT_SUMMARY_START` / `AUTO_CURRENT_SUMMARY_END` 包裹，标明生成日期和实际覆盖的运行；历史图不能冒充新源码已重跑。

固定七节结构：

1. 实验层级与方向定义
2. 脚本和结果目录（包含当前运行命令）
3. 结果总览与逐图解读（按分析模块）
4. 这套结果能直接支持什么
5. 这套结果目前不能直接支持什么
6. 推荐阅读顺序
7. 优先读取文件

每个模块先概述实际统计量和结论，再解读全部真实图件。每张图必须恰好有三要素：**生物学意义**（回答什么）、**数据分布特征**（图上真实分布）、**最终结论**（与方向及证据等级一致）。同名 PNG/PDF 只解读一次，优先 PNG；仅 PDF 单列。图件覆盖应为 100%，不编造不存在的图。路径相对所属分析模块根。

共同硬约束：

- 写明生物学窗口、真实物种/基因组/注释、暴露与对照、实际样本、剔除及理由、统计单位。PFOS 为 `PFOS − WT`，母龄为 `MO − MY`，不可串用；父鼠、窝、混池或配对未知时保留 unknown，不能以细胞或文库数量补出生物学重复。
- 原 P 值口径必须标为探索性，并保留 BH/FDR。阈值、差异数量、富集 `p.adjust/Count/GeneRatio` 和网络等数值逐项对照该模块实际 CSV；不得复制其他课题数字、FDR 状态或特异性结论。
- 说明实际完成、缺失、未运行和失效输出。缺失模块不能推断没有对应生物信号；严格 FDR 未通过的通路不能写成严格显著。
- 不为写报告修改统计功能、变量或工作文件，不为文档添加测试脚本。目录整理的源码移动仅按已批准整理计划实施，不是报告生成的副作用。
- 图件解读必须来自图和表的真实信息，不能靠通用模板填充虚构趋势。生成后核对方向、数值、三要素、覆盖和跨项目泄漏；已有正确图表和历史执行凭据保持来源身份。

下列组学章节保留原规范的专属解释规则。模块名称必须以当前实际入口和磁盘为准，历史编号不能伪装成现役模块。


## RNA-seq 专属规则

### §1 实验层级与方向定义（必含字段）

- **数据层级**：项目名 + 一句生物学窗口描述
  （如 2cell：“父源 PFOS 暴露后最早可检测的胚胎转录响应窗口”）。
- **物种与注释**：`mouse / GRCm39`，OrgDb=`org.Mm.eg.db`，
  GTF=`Mus_musculus.GRCm39.115.gtf`（同批项目一致）。
- **分组定义**：写清 A 组=暴露、B 组=对照，带项目真实组名（如 `PFOS` / `WT`）。
- **样本清单**：列出实际样本名（一般 `WT1/2/3` + `PFOS1/2/3`）。
- **样本剔除情况**：是否全部进入主分析、有无被自动剔除。
- **差异方向定义**：`log2FoldChange = 暴露组 - 对照组`。
  - `> 0`：暴露组表达更高；`< 0`：对照组更高。
- **主 DE 口径**：`pvalue < 0.05 且 |log2FC| >= X`
  （X 取自 `05_de/tables/5_2b_threshold_advisory_counts.csv`）。
- **更严格汇报口径**：`padj < 0.05`，或 `padj < 0.05 且 |log2FC| >= X`。
- **敏感性比较**：主分析与 full snapshot 是否同一批样本、结论是否一致。

### §2 脚本和结果目录

- 列出主目录 / 结果目录 / 18 个 Rmd 脚本完整清单
  （`01_main_setup.R` … `18_circle.Rmd`）。
- **标准模块含义**：每个模块一行一句话：

  | 模块 | 含义 |
  |------|------|
  | `00_audit` | 上游流程、MultiQC、软件版本、参数审计 |
  | `02_selection` | 样本选择与 count-based QC |
  | `03_qc` | 表达分布、RLE、样本相关性、距离、fastp 质控 |
  | `04_pca` | 样本 PCA |
  | `05_de` | 差异表达主结果、火山图、MA 图、敏感性比较 |
  | `06_enrichment` | KEGG / Reactome ORA 与 MSigDB GSEA |
  | `07_go` | GO BP / MF / CC / ALL |
  | `08_gsea` | 基于全基因排序的 GO / KEGG / Hallmark / Reactome GSEA |
  | `09_gsva` | 样本级 pathway score |
  | `10_tf_regulator` | TF 活性推断及靶基因网络 |
  | `11_ppi` | STRING PPI 网络、hub、bridge、子网络 |
  | `12_wgcna` | WGCNA 与模块注释 |
  | `13_apa` | APA / IPA 事件 |
  | `14_der` | DER（differential exon/region usage signal） |
  | `15_splicing` | rMATS、DAS、DEU、isoform switch、多方法整合 |
  | `16_immune` | ssGSEA / mMCPcounter 免疫浸润推断 |
  | `17_pathway_activity` | 另一套 pathway activity 估计与多方法一致性 |
  | `18_circle` | 富集结果的 circle 二次整合展示 |

- 说明 `06/07/08` 下含 `01_all / 02_up / 03_down` 三个方向目录
  （`up`=暴露高，`down`=对照高），图件类型对称，解读时方向对应即可。

### §4 能直接支持什么

基于真实结果的强结论概述，例如：
- 组间存在大规模稳健转录差异；
- 差异不仅体现于 DE，还体现于 pathway / TF / PPI / APA / splicing / isoform switch；
- 哪些富集库有 FDR 通过项、反复出现的主线是什么。

### §5 不能直接支持什么（必须显式列出的红线）

- `DER` 不能写成具体基因调控证据（大量为 intergenic）。
- `immune` 推断不能写成“细胞组成已被直接验证”（只是 bulk RNA 计算推断，且方法间一致性不高）。
- PPI `hub` 不能写成“最关键驱动基因”（多反映数据库网络中心性）。
- 缺失的 WGCNA 结果**不能**推导“没有共表达模块”。
- GSEA 各库严格 FDR 通过有限时，不能把单条通路写成严格显著（需写明 `p.adjust` 上限）。

### §6 推荐阅读顺序

- 8 步推荐顺序：
  1. `00_audit` 与 `02_selection/tables/2_1_selected_sample_table.csv`
  2. `05_de` 主表、火山图、热图
  3. `06_enrichment`、`07_go`、`08_gsea`、`09_gsva`、`17_pathway_activity`
  4. `10_tf_regulator` 与 `11_ppi`
  5. `13_apa` 与 `15_splicing`
  6. `16_immune`
  7. `14_der`
  8. `18_circle`
- 再列“最值得作为主图 / 主结果来源”的模块清单
  （DE、KEGG ORA、GO、GSEA、GSVA、TF、PPI、APA、splicing、pathway_activity 等）。

### §7 优先读取文件

列出各模块**关键 CSV 相对本分析模块根目录的路径**（覆盖 02/05/06/07/08/09/10/11/13/14/15/16/17），例如：

```
- result/02_selection/tables/2_1_selected_sample_table.csv
- result/05_de/tables/5_1_annotated_DE_results.csv
- result/05_de/tables/5_2b_threshold_advisory_counts.csv
- result/06_enrichment/01_all/tables/6_1_KEGG_ORA_results.csv
- result/07_go/01_all/tables/7_0_GO_directional_summary.csv
- result/08_gsea/01_all/tables/8_4_GSEA_Reactome_results.csv
- result/09_gsva/tables/09_9_leading_GSVA_pathway_summary.csv
- result/10_tf_regulator/tables/10_2_TF_activity_differential.csv
- result/11_ppi/tables/11_4_PPI_network_topology.csv
- result/13_apa/tables/13_3_APA_gene_summary.csv
- result/14_der/tables/14_1_all_DERs.csv
- result/15_splicing/01_event_level/tables/15_3_rMATS_event_summary.csv
- result/15_splicing/02_gene_level/tables/15_21_candidate_splicing_regulators.csv
- result/15_splicing/04_isoform_switch/tables/15_42_switch_summary.csv
- result/16_immune/tables/16_10_differential_immune_ranking.csv
- result/17_pathway_activity/tables/17_7_leading_pathway_activity_summary.csv
```


## WGBS 专属规则

### §1 实验层级与方向定义（必含字段）

- **数据层级**：项目名 + 一句生物学窗口描述
  （如 `wgbs_2cell`：“F1 2-cell，受精后早期胚胎重编程窗口”；
  `wgbs_sperm`：“F0 雄鼠 sperm，父源遗传载体层面”；
  `wgbs_liver`：“F1 liver”）。
- **物种与注释**：`mouse / GRCm39`，参考注释沿用 GRCm39 默认 GTF（与同批项目一致）。
- **分组定义**：写清 A 组=暴露、B 组=对照，带项目真实组名（如 `PFOS` / `WT`）。
- **样本清单**：列出实际样本名（一般 `WT1/2/3` + `PFOS1/2/3`）。
- **样本剔除情况**：是否全部进入主分析、有无被自动剔除或低覆盖样本。
- **差异方向定义**：`diffMethy = comparison - reference`（默认 `PFOS - WT`，
  由 `run_wgbs.sh` 的 `WGBS_COMPARISON_GROUP` / `WGBS_REFERENCE_GROUP` 决定）。
  - `diffMethy > 0` 记为 `Hyper`：PFOS 组甲基化更高。
  - `diffMethy < 0` 记为 `Hypo`：PFOS 组甲基化更低。
- **DMR calling 口径**：DSS 结果；主显著性口径由 `WGBS_SIGNIFICANCE_MODE`
  （`pvalue` 默认或 `padj`）决定，需写明当前使用的口径与阈值。
- **DML / 单 CpG 定位**：DMR 内的单个差异 CpG（DML）是否一并给出，
  但主结论不应只依赖单 CpG。
- **敏感性比较**：主分析与 full snapshot 是否同一批样本、结论是否一致。

### §2 脚本和结果目录

- 列出主目录 / 结果目录 / Rmd 脚本完整清单。
- 标准模块含义（以当前 `run_wgbs.sh` 实际模块为准，旧 9 号跨组学分支已退休；
  `wgbs_sperm` 额外有 `11_sperm`）：

  | 模块 | 含义 |
  |------|------|
  | `1_setup` | 路径、样本表、分组与比较方向初始化 |
  | `2_qc` | 样本层 QC、覆盖度、相关性 / PCA |
  | `3_dmr` | DSS DMR calling 和基因注释 |
  | `4_dmroverview` | DMR 总体分布、长度、染色体负荷、TSS 距离 |
  | `5_enrichment` | DMR 相关基因 GO / KEGG / GSEA |
  | `6_genomecontext` | CpG island、repeat、genomic context |
  | `7_regulatory` | motif、imprinted genes、repeatmasker 等调控注释 |
  | `8_gene_methylation` | 纯基因甲基化汇总，保留区域数量与方向冲突 |
  | `10_browser` | 代表位点浏览器图 |
  | `11_sperm` | sperm 特有：精子发生阶段 / 父源遗传相关注释（仅 wgbs_sperm） |

- 说明 `05_enrichment` 下若含 `01_all / 02_up / 03_down` 三个方向目录
  （`up`=Hyper，`down`=Hypo），图件类型对称，解读时方向对应即可。

### §4 能直接支持什么

基于真实结果的强结论概述，例如：
- 组间存在大规模稳健的 regional 甲基化改变（DMR）；
- 差异集中在 promoter / gene body / intergenic 等特定基因组语境；
- 哪些富集库有 FDR 通过项、反复出现的主线是什么；
- 跨组学耦合结果仅在 `integration/` 中报告；区域和基因注释重叠本身不证明调控。

### §5 不能直接支持什么（必须显式列出的红线）

- **单 CpG / DML 不能写成主机制证据**：DMR 是 regional 信号，个别 CpG 的波动
  不足以支撑调控结论；主结论应基于 DMR。
- **某个 DMR 不能写成“一定直接调控相邻基因表达”**：需要 RNA、ATAC、
  跨阶段 overlap 或靶向验证支持。
- **基因组语境计数（Intron / Intergenic 占比高）不能反向推断“大部分无功能”**：
  只是说明定位，不评价功能重要性。
- **imprinted / repeat 注释不能写成“已确认的功能性印记失调”**：只是富集 /
  共定位信号，需进一步验证。
- **GSEA 各库严格 FDR 通过有限时，不能把单条通路写成严格显著**（需写明 `p.adjust` 上限）。
- **缺失的模块结果不能推导“没有对应信号”**（如某项目未跑 `11_sperm`，
  不能写“精子相关甲基化无变化”）。

### §6 推荐阅读顺序

- 推荐顺序：
  1. `1_setup` 与 `2_qc`（覆盖度、PCA）确认样本质量
  2. `3_dmr` 主表、火山图、DMR 总览
  3. `4_dmroverview`（长度 / 染色体 / TSS 分布）
  4. `5_enrichment`、`6_genomecontext`
  5. `7_regulatory`（motif / imprinted / repeat）
  6. `8_gene_methylation`（纯甲基化基因层汇总）；RNA 耦合另见 integration
  7. `10_browser`（代表位点）
  8. `11_sperm`（仅 sperm 项目）
- 再列“最值得作为主图 / 主结果来源”的模块清单
  （DMR 总览、DMR 火山图、富集 top、纯基因甲基化汇总等）。

### §7 优先读取文件

列出各模块**关键 CSV 的完整路径**（覆盖 2/3/4/5/6/7/8/9/10），例如：

```
- wgbs/<project>/result/02_qc/tables/*.csv        # 覆盖度 / PCA / 相关性
- wgbs/<project>/result/03_dmr/tables/*.csv       # DMR 主表与注释
- wgbs/<project>/result/04_dmroverview/tables/*.csv
- wgbs/<project>/result/05_enrichment/01_all/tables/*.csv  # GO/KEGG/GSEA（如有方向目录）
- wgbs/<project>/result/06_genomecontext/tables/*.csv
- wgbs/<project>/result/07_regulatory/tables/*.csv
- result/8_gene_methylation/tables/1/gene_methylation_long.csv
- result/8_gene_methylation/tables/1/gene_methylation_wide.csv
```

> 注意：WGBS 各模块的表格文件名**不固定**（不像 rna_seq 有 `5_1_*` 编号约定），
> 生成时以磁盘上 `result/*/tables` 真实存在的 CSV 为准，逐项目填写具体文件名。


DMR 总数、Hyper/Hypo、promoter DMR、feature 计数均从该模块表核对。diffMethy 若保存为比例，与百分点不同，报告时写明单位，不能偷换为表达 log2FC。


## ATAC-seq 专属规则

### §1 实验层级与方向定义（必含字段）

- **数据层级**：项目名 + 一句生物学窗口描述（如 `projects/pfos/mouse/atac/liver`：“F1 肝脏，
  主要代谢器官对 PFOS 暴露的染色质可塑性响应”）。
- **物种与注释**：`mouse / GRCm39`，参考注释沿用 GRCm39 默认 GTF。
- **分组定义**：A 组=暴露、B 组=对照，带项目真实组名（如 `PFOS` / `WT`）。
- **样本清单**：实际样本名。
- **样本剔除情况**：是否全部进入主分析、有无被自动剔除或低 QC 样本。
- **差异方向定义**：`log2FC = comparison - reference`（默认 `PFOS - WT`）。
  - `log2FC > 0` 记为 `UP`（PFOS 开放度更高）；
  - `log2FC < 0` 记为 `DOWN`（PFOS 开放度更低）。
- **peak calling 口径**：peaks（MACS2 风格）+ 差异可及性（DESeq2 / limma 风格），
  需写明当前显著性口径（pvalue / padj）与阈值。
- **敏感性比较**：主分析与 full snapshot 是否同一批样本、结论是否一致。

### §2 脚本和结果目录

列出本模块实际脚本、结果目录和模块含义。当前维护入口包含：

| 模块 | 含义 |
|---|---|
| `01_main_setup` / `02_project_inventory` | 路径、样本、方向与输入清单 |
| `03_atac_qc` / `04_peak_catalog_annotation` | QC、peak 目录与基因组语境 |
| `05_differential_accessibility` | 差异开放主表与方向 |
| `06_genome_tracks` | 本组学代表位点轨道 |
| `07_peak_gene_enrichment` | peak-gene 注释、GO/KEGG 与排序富集 |
| `08_motif_activity` / `10_footprinting_tf_evidence` | motif、chromVAR、足迹及 TF 证据 |
| `12_tf_gene_network` | 本组学 TF 位点网络 |
| `13_conservation_gc` / `14_machine_learning` | 保守性/GC 和探索性预测 |
| `15_keygenes_locus_reports` / `16_export_summary_manifest` | ATAC 候选与来源清单 |

旧 09 RNA 和 11 WGBS 联合分支已退休，跨组学归 integration。UP/DOWN 各方向的富集和 motif 结果分别解释，UP 为开放度升高，DOWN 为降低。旧批次调试 `90_*` 不纳入当前主线。

### §4 能直接支持什么

基于真实结果的强结论概述，例如：
- 样本 QC 是否过关（覆盖度、fingerprint、PCA 聚类）；
- 是否存在大规模稳健的差异可及性（DE peaks，UP/DOWN 分布）；
- 差异 peak 的基因组语境是否集中在 promoter / distal；
- 哪些富集库（GO / KEGG / motif / chromVAR）有 FDR 通过项、反复出现的主线；
- 同层级 RNA / WGBS 对照只在 integration 报告；本模块仅给出开放性及注释证据。

### §5 不能直接支持什么（必须显式列出的红线）

- **某个差异 peak 不能写成“一定直接调控相邻基因表达”**：需要 RNA、WGBS、
  跨阶段 overlap 或靶向验证支持。
- **genomic context 计数（Distal / Intron 占比高）不能反向推断“大部分无功能”**：
  只说明定位，不评价功能重要性。
- **motif / chromVAR 偏差不能写成“已确认的 TF 直接结合调控”**：只是富集 /
  偏差信号，需进一步验证（如 footprint / ChIP）。
- **GSEA 各库严格 FDR 通过有限时，不能把单条通路写成严格显著**（需写明 `p.adjust` 上限）。
- **缺失的模块结果不能推导“没有对应信号”**（如未跑 `16_var`，不能写“可变性无变化”）。

### §6 推荐阅读顺序

推荐先核对输入与 QC，再看 peak 规模/TSS/基因组语境、差异开放主表、peak-gene 富集、motif/chromVAR、足迹与 TF 位点网络，最后读代表位点及探索性预测。主图应从实际 DE、TSS、QC 和方向富集结果选择。

### §7 优先读取文件

逐项目列出实际 `result/*/tables` CSV：QC、peak 注释、差异开放、peak-gene 富集、motif/足迹、TF 网络和候选表。文件名以磁盘为准，优先给出 `result/15_keygenes_locus_reports/tables/ATAC_Gene_Candidates.csv`；旧模块编号不补造当前文件。


peak 总数、consensus 数、差异 UP/DOWN 数、TSS 强度和 genomic context 计数必须与实际 CSV 一致。原位点关联或最近基因注释不能升级为已验证的调控边。
