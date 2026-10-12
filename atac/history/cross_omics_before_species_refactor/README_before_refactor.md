# ATAC-seq 独立两组分析模板

复制整个目录，修改 `config/project.env` 与 `config/samples.tsv` 后使用。模板提取自 liver ATAC 流程，支持小鼠 GRCm39 双端独立两组设计；全体预先合格样本进入主分析，不依据 PCA 分离程度选样。

```bash
bash run_upstream.sh --check --steps 08,09
bash run_rmd.sh --dry-run --steps 5-8
bash run_upstream.sh
bash run_rmd.sh
```

`--help`、`--list`、`--check`、`--dry-run` 不创建日志或分析结果。检查验证显式样本设计与调度，不代表参考索引、原始数据或工具已经可用。实际执行缺必要输入会明确失败。独立入口不需要复制目录外的辅助脚本。

## 输入、方向与参数

TSV 必需列为 `sample, biological_sample_id, group, stage, raw_fq1, raw_fq2`。每行是一个独立生物学单位，每组至少 2 个；混池计为一个单位，技术重复应先合并。未知/重复身份、配对/重复测量或 batch 设计信息拒绝进入本首版独立两组模型。样本 ID 与测序文件名可以不同；TSV 空字段保留。

显式设置 `ATAC_REFERENCE_GROUP`、`ATAC_COMPARISON_GROUP`；所有主差异的效应为 comparison − reference。样本表重排不会改变方向。参考、建库类型、工具、线程、差异阈值、比较组和跨组学输入均在同一 `project.env`。路径可为绝对路径或相对模板根目录；脚本根目录即项目根目录。

默认 `ATAC_SIGNIFICANCE_MODE=pvalue`，属于探索性名义显著筛选；保留 FDR 与对应分层。原始 p 缺失时不借用 FDR。默认不加额外 fold-change 门；如需指定，在同一配置中设置 `ATAC_PRIMARY_REQUIRE_FC` 与阈值。MACS2 peak calling 的 q 与 DiffBind 差异检验 FDR 含义分开。

## 分析范围

上游保留 FastQC/fastp、比对/去重/QC、MACS2、blacklist、RPGC bigWig、FRiP、TSS、IDR、MultiQC 和 consensus fragment counts。默认每工具 2 线程、1 个样本并行。consensus counts 使用 featureCounts 的成对片段计数，不以单条 alignment 数代替。IDR 的完整流程只对**每组恰好 2 个生物学重复**运行；重复更多或不足均明确跳过，不挑选两个样本凑齐。

下游默认运行 1–8 和 16：项目/输入审计、QC、peak 注释、DiffBind、轨迹、基因富集、motif、导出报告。DiffBind 比较方向固定于配置，默认保留全部样本。14 机器学习保持停用，已去除核心对 randomForest/caret/pROC 的依赖。

- 跨组学模块 9/11/12/15 默认关闭，显式设置 `ATAC_ENABLE_MULTIOMICS=1` 后启用。每个外部输入文件单独声明；不自动搜索历史目录或连接其他真实项目。RNA DE 的 `sampleA/sampleB` 与 WGBS DMR 的 `group2_label/group1_label` 必须对应 comparison/reference，方向未知或反向时失败。基因/区间整合不代表配对动物推断。
- Footprint/BINDetect 默认关闭，需 `ATAC_ENABLE_FOOTPRINT=1`。合并 BAM 和轨道清单按样本表分组；BINDetect 的稳定内部条件名为 comparison/reference，报告保存实际组名映射。原始 Tn5 profile 为描述性证据，不能当作 TF 结合检验。BINDetect 需要完成状态和实际输入/输出 SHA256 全部匹配后才纳入证据。
- 保守性模块默认关闭，由 `ATAC_ENABLE_CONSERVATION` 控制。可选 ataqv 保留原有依赖/状态门，不声称缺依赖的步骤已经运行。

## 显式 consensus counts 输入

已有经过审计的整数片段计数时，设置 `ATAC_INPUT_MODE=counts` 和 `ATAC_COUNTS_FILE=/path/counts.tsv`，运行 `bash run_rmd.sh --check` 后 `bash run_rmd.sh`。counts 模式从正式入口执行原区间 DESeq2 检验、BH 校正、差异分层表、全样本 PCA、火山图和 HTML，不需要 BAM 或参考索引；样本表在此模式可省略 raw_fq1/raw_fq2。

计数表与上游输出使用同一格式：`peak_id, CHR, START, END, PeakWidth` 加各 sample 计数列，TSV；PeakWidth 可省略，坐标为 BED 的 0-based half-open。sample 列必须与样本表一一对应，禁止重复 peak/区间、非整数、负数、NA 或零总计数样本。`ATAC_COUNTS_MIN_TOTAL` 默认 10，是不使用分组信息的总计数过滤。分析输出坐标转为 1-based inclusive，并保留区间身份。

来源写为 `consensus_counts_DESeq2`，含输入 SHA256、包版本、实际 biological n 和执行层级。此模式检验原 consensus 区间；BAM 模式的 DiffBind 可以重新定位 summit，两者的检验对象不同。counts 模式不运行原始 QC、DiffBind、注释、轨道或 footprint，也不复用它们的状态缓存；`--steps` 只接受 5（或 diff/counts）。

## 输出与可审计性

上游维持原编号目录，下游输出 `result/<模块>/tables` 与 `plots`，汇总为 `result/index.html`。报告支持仅表或跳过的运行；`--report-only` 要求配置签名匹配。入口不默认删除结果，不保留结果目录迁移或旧命名路径回退。

`metadata/*_configuration.json` 记录样本/方向/配置、脚本 SHA256、外部输入属性。来源缺失为 `UNVERIFIED`，不匹配为 `STALE`，不能仅凭已有 BAM/表格把结果判为当前；应在新的输出目录重算。文件大小/mtime 不是原始数据内容哈希。DiffBind 状态缓存还绑定自身分析签名；BINDetect 使用独立 SHA256 证据校验。

`optional_module_gates.tsv` 记录启用开关，`execution.tsv` 记录每个实际模块的 STARTED/COMPLETED/FAILED。资格、是否启用、是否执行是不同状态。

来源见 `provenance.tsv`；验证边界见 `VALIDATION.md`。复用现有 ATAC 环境；不自动安装或升级依赖。仅检查输入设计时 Python 标准库即可；`ATAC_LIGHT_SETUP=TRUE` 可运行小数据输入审计，不代表 BAM 分析或差异统计已验证。

报告按 `execution.tsv` 每个模块的最后一条记录展示状态，只收录 `COMPLETED` 模块的图表；失败、中断或缺执行记录留下的旧表不会列为当前结果。counts 模式的 `--report-only` 还要求 05 最近一次执行完成；统计与报告各自记录失败状态。
