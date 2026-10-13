# ATAC 单组学分析入口

开放图谱、QC、差异开放、peak-gene 注释/富集、motif、足迹、TF 位点网络和 ATAC 候选位点。现役方法和结果使用连续编号 01–14；09 是原始 Tn5 描述曲线，10 为受来源门控制的 TF 网络，13 为纯 ATAC 候选与位点图。跨组学由 integration 处理。

## 运行

```bash
bash run_rmd.sh --list
bash run_rmd.sh --dry-run --all
bash run_rmd.sh --steps 5,7
bash run_rmd.sh --plots-only --check
bash run_rmd.sh --plots-only
Rscript scripts/run_gene_candidates.R .
```

路径相对本分析模块根解析，不受启动目录影响；公共参考资源和工具通过相对配置、PATH 或显式环境变量提供。

## 当前产物

`result/05_differential_accessibility` 为受检及差异 peak；`result/07_peak_gene_enrichment` 提供基因映射；`result/13_keygenes_locus_reports/tables/ATAC_Gene_Candidates.csv` 为纯 ATAC 候选。候选是本组学排序，不能据此确认调控。

从保存的 DMR/peak 源表重建的汇总有独立 input_manifest.tsv 和 execution.tsv，源状态保留为 PRESERVED_OUTPUT，不代表重新完成上游测序分析。`result/index.html` 提供图表导航；科学解释还需核对来源状态。

中心模板只维护当前单组学实现。复制模板后填写本地 config，模板运行不引用项目结果。

## 联合分析与历史

跨组织、跨阶段、跨组学由物种目录的 integration 统一处理。项目历史材料保留在各项目所属模块内部，本中心模板不复制旧项目结果。

## 模板适用设计


TSV 必需列为 `sample, biological_sample_id, group, stage, raw_fq1, raw_fq2`。每行是一个独立生物学单位，每组至少 2 个；混池计为一个单位，技术重复应先合并。未知/重复身份、配对/重复测量或 batch 设计信息拒绝进入本首版独立两组模型。样本 ID 与测序文件名可以不同；TSV 空字段保留。

显式设置 `ATAC_REFERENCE_GROUP`、`ATAC_COMPARISON_GROUP`；所有主差异的效应为 comparison − reference。样本表重排不会改变方向。参考、建库类型、工具、线程、差异阈值、比较组均在同一 `project.env`。路径可为绝对路径或相对模板根目录；脚本根目录即项目根目录。

默认 `ATAC_SIGNIFICANCE_MODE=pvalue`，属于探索性名义显著筛选；保留 FDR 与对应分层。原始 p 缺失时不借用 FDR。默认不加额外 fold-change 门；如需指定，在同一配置中设置 `ATAC_PRIMARY_REQUIRE_FC` 与阈值。MACS2 peak calling 的 q 与 DiffBind 差异检验 FDR 含义分开。


## 已有 consensus counts


已有经过审计的整数片段计数时，设置 `ATAC_INPUT_MODE=counts` 和 `ATAC_COUNTS_FILE=/path/counts.tsv`，运行 `bash run_rmd.sh --check` 后 `bash run_rmd.sh`。counts 模式从正式入口执行原区间 DESeq2 检验、BH 校正、差异分层表、全样本 PCA、火山图和 HTML，不需要 BAM 或参考索引；样本表在此模式可省略 raw_fq1/raw_fq2。

计数表与上游输出使用同一格式：`peak_id, CHR, START, END, PeakWidth` 加各 sample 计数列，TSV；PeakWidth 可省略，坐标为 BED 的 0-based half-open。sample 列必须与样本表一一对应，禁止重复 peak/区间、非整数、负数、NA 或零总计数样本。`ATAC_COUNTS_MIN_TOTAL` 默认 10，是不使用分组信息的总计数过滤。分析输出坐标转为 1-based inclusive，并保留区间身份。

来源写为 `consensus_counts_DESeq2`，含输入 SHA256、包版本、实际 biological n 和执行层级。此模式检验原 consensus 区间；BAM 模式的 DiffBind 可以重新定位 summit，两者的检验对象不同。counts 模式不运行原始 QC、DiffBind、注释、轨道或 footprint，也不复用它们的状态缓存；`--steps` 只接受 5（或 diff/counts）。


Footprint/BINDetect 显式由 ATAC_ENABLE_FOOTPRINT 控制；保守性由 ATAC_ENABLE_CONSERVATION 控制。差异TF结合需实际来源、完成状态和SHA256相符。当前中心模板支持小鼠 GRCm39 双端独立两组设计，不自动安装依赖。

报告写法见[统一报告规范](../docs/reporting.md)。

内部记录见[验证记录](docs/VALIDATION.md)。

## 冻结数值重绘

BAM/DiffBind 分支已完成且具备原保存数值和 .cache/rmd 绘图对象时，`--plots-only` 走独立早分支，不进入 Rmd 统计和依赖补跑，不执行 DiffBind/fgsea/limma/chromVAR 拟合或上游命令。counts 模式没有这些对象时明确报缺失，不借用其他来源。原数值 SHA、图件覆盖与格式检查通过后才替换 result，失败保留当前有效结果。编号方法下并列 plots/tables，PDF/PNG 同 stem。183 mm、300 dpi、RGB 白底及 Cairo PDF 使用本地冻结 figure_style.R；复制时记录 template_snapshot.tsv。

## 统一绘图配置

`config/figure_style.json` 是本模板随项目冻结的语义配置。复制后编辑本地配置，运行不读取 PFOS 或中心模板。R 实现通过本地 `figure_style.R` 的明确相对路径加载；分类按名称绑定，不按输入行序轮换。ORA 使用 GeneRatio、Count 面积及 −log10(P) viridis；色阶 0–10，端点标 ≥10。相关系数固定 [−1,1]；绝对甲基化固定 [0,100]%；Row Z-score 固定 [−2.5,2.5] 并标饱和；缺失值单独显示灰。不同测量量保留独立图例与单位。

原生 UpSetR 保留集合与交集排序：组集合使用组身份色，交集计数使用中性色，计数标签用深灰，集合规模刻度以 K 简写。DiffBind Venn 的 Gain/Loss 使用与差异图一致的方向色。
