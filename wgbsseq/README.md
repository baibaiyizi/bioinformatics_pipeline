# WGBS 独立两组分析模板

复制整个目录到新项目，修改 `config/project.env` 与 `config/samples.tsv`。支持小鼠 GRCm39 双端 WGBS；上游为 FastQC → fastp → Bismark → MethylDackel → methylKit，下游为 QC → DSS DML/DMR → 注释、富集、基因组背景与轨迹。源项目只作代码来源，本模板不会连接原工作区的分析结果。

```bash
bash run_upstream.sh --help
bash run_upstream.sh --check --steps 08,09
bash run_wgbs.sh --check --steps 2-7
bash run_upstream.sh --steps 1-5,12-13
bash run_wgbs.sh
```

`--check`/`--dry-run` 只检查样本合同、步骤与配置，不创建日志或结果，不证明 FASTQ、参考索引和依赖已经可用。实际执行时缺少所选步骤的必要输入或工具会失败。`--list` 列出模块；选步统一用 `--steps`，没有旧参数兼容入口。工具通过当前 PATH 或配置的工具目录定位；Pandoc 必须在 PATH 中。

## 输入与统计单位

样本表为 TSV，必需列：`sample, biological_sample_id, group, stage, raw_fq1, raw_fq2`。FASTQ 原始名称可以不同于 sample，clean FASTQ 统一按 sample 命名。路径可为绝对路径，或相对项目根目录。空字段不会因 Shell 分隔符折叠而错列。

每行对应一个独立生物学单位，每组至少 2 个单位。重复 sample、生物学 ID 重复、未知 ID、pair_id/subject_id/batch 设计信息均拒绝进入此两组模板；技术文库须在实验单位层面先合并。一个混池只算一个单位，不能按池内动物数扩增 n。此最小重复门并不保证统计功效。

`WGBS_REFERENCE_GROUP` 与 `WGBS_COMPARISON_GROUP` 必须显式指定；甲基化差异始终为 comparison − reference。重新排序样本表不改变方向；CpG.RData 必须包含可唯一匹配的 sample.id 和相符 treatment，不会按位置改写身份。所有预先合格样本都进入分析，PCA 仅用于诊断。

## 配置与可选模块

- 参考目录、FASTA、GTF、资源、建库方向、输入层级集中在 `project.env`。默认单样本顺序执行；Bismark 并发为 2 × 2，QC 相关矩阵至多抽取 100,000 个位点。`raw/clean/cached` 由 `WGBS_INPUT_MODE` 显式选择。
- `WGBS_LIBRARY_TYPE=directional|non_directional` 对应实际建库。M-bias 默认不附加任何组织专用 ignore 参数；需检查自己的 M-bias 后制定修剪规则。
- 默认执行模块 1–7 和 10。8/9 跨组学须 `WGBS_ENABLE_MULTIOMICS=1`，并提供显式 RNA/ATAC 输入；RNA DE 必须有唯一、无缺失的 `sampleA=comparison`、`sampleB=reference`，反向或未知方向直接失败。
- 11 精子专题须 `WGBS_ENABLE_SPERM=1`。ICR 坐标必须通过 `WGBS_SPERM_ICR_FILE` 提供已核验的 GRCm39 TSV，至少包含 `ICR_name, chr, start, end, imprinted_gene, parental_origin`；模板不附带未经核验的内置坐标。精子发生基因表可另行配置。
- `WGBS_SIGNIFICANCE_MODE=pvalue` 是默认探索性筛选。DML 结果保留原始 p 和 `padj_BH`；后者是 **位点层** BH 校正。DSS DMR 按 p.threshold、pct.sig、minlen、minCG、delta 聚合，不能称为区域 FDR 控制。可选 dmrseq 默认关闭，其区域检验单独报告。

## 输出、缓存与报告

上游结果保留原编号目录；下游输出在 `result/<模块>/tables` 与 `plots`。`result/index.html` 支持有图、仅表和跳过状态。已有结果的索引可用 `bash run_wgbs.sh --report-only` 重建，要求当前配置来源记录相符。

`metadata/*_configuration.json` 绑定样本、方向、当前配置、脚本 SHA256 与输入文件大小/mtime。BAM quickcheck 或文件非空本身不够：复用必须先通过来源签名。缺失来源记录为 `UNVERIFIED`，配置/输入属性变化为 `STALE`，不会自动把旧结果当作新结果；使用新的项目输出目录重新运行。大小/mtime 签名不是大文件内容校验，FASTQ 的传输校验应另行保留。TxDb 缓存写入项目 `.cache/annotation`，不修改外部参考目录。

模块的 `STARTED/COMPLETED/FAILED/SKIPPED_DISABLED` 写入 `result/execution.tsv`。输入合同通过与统计分析完成分别记录。脚本不会因某个现存 HTML/CSV 文件而宣称运行完成。下游不预先清空结果；分析模块只刷新自己的产物。

来源记录见 `provenance.tsv`，实际验证和未运行边界见 `VALIDATION.md`。依赖沿用已有 WGBS 环境，不自动安装或升级软件；除 R/methylKit、DSS、bsseq、GenomicRanges、ChIPseeker、clusterProfiler 等原模块依赖外，样本审计仅需 Python 标准库。

报告按 `execution.tsv` 每个模块的最后一条记录展示状态，只收录 `COMPLETED` 模块的图表；失败、中断或缺执行记录留下的旧表不会列为当前结果。

DSS 原生方向字段保留：DML 的 `diff`/`stat` 和 DMR 的 `areaStat` 沿用 DSS 的 group1（reference）− group2（comparison）；模板规范化的 `diffMethy`/`diff.Methy` 为 comparison − reference，并附 `group1_label`/`group2_label`。因此同一合成升高区域可以有原生 diff/stat/areaStat 为负、规范化 diffMethy 为正。筛选、Hyper/Hypo 和报告效应统一读取规范化字段，不能把原生统计量的符号混作规范化效应方向。原始 p 和 BH/FDR 不因换向而改变。
