# WGBS 单组学分析入口

甲基化 QC、DSS DML/DMR、注释、富集、基因组背景、纯基因甲基化汇总和浏览器轨道。8 为 `8_gene_methylation`，9 已移除，10 browser 不再读取表达耦合缓存。精子模块另有 11 专题。

## 运行

```bash
bash run_wgbs.sh --list
bash run_wgbs.sh --check
bash run_wgbs.sh --steps 8,10
Rscript scripts/run_gene_methylation.R .
```

路径相对本分析模块根解析，不受启动目录影响；公共参考资源和工具通过相对配置、PATH 或显式环境变量提供。

## 当前产物

`result/3_dmr/tables/2/dmr_annotated.csv` 为区域证据；`result/8_gene_methylation/tables/1/gene_methylation_long.csv` 和宽表按基因、功能层汇总。每层保留绝对差异最大的区域，并保存 RegionCount 和 DirectionConflict；这不是基因级显著性检验。

从保存的 DMR/peak 源表重建的汇总有独立 input_manifest.tsv 和 execution.tsv，源状态保留为 PRESERVED_OUTPUT，不代表重新完成上游测序分析。`result/index.html` 提供图表导航；科学解释还需核对来源状态。

中心模板只维护当前单组学实现。复制模板后填写本地 config，模板运行不引用项目结果。

## 联合分析与历史

跨组织、跨阶段、跨组学由物种目录的 integration 统一处理。项目历史材料保留在各项目所属模块内部，本中心模板不复制旧项目结果。

## 模板适用设计


样本表为 TSV，必需列：`sample, biological_sample_id, group, stage, raw_fq1, raw_fq2`。FASTQ 原始名称可以不同于 sample，clean FASTQ 统一按 sample 命名。路径可为绝对路径，或相对项目根目录。空字段不会因 Shell 分隔符折叠而错列。

每行对应一个独立生物学单位，每组至少 2 个单位。重复 sample、生物学 ID 重复、未知 ID、pair_id/subject_id/batch 设计信息均拒绝进入此两组模板；技术文库须在实验单位层面先合并。一个混池只算一个单位，不能按池内动物数扩增 n。此最小重复门并不保证统计功效。

`WGBS_REFERENCE_GROUP` 与 `WGBS_COMPARISON_GROUP` 必须显式指定；甲基化差异始终为 comparison − reference。重新排序样本表不改变方向；CpG.RData 必须包含可唯一匹配的 sample.id 和相符 treatment，不会按位置改写身份。所有预先合格样本都进入分析，PCA 仅用于诊断。


当前中心模板支持小鼠 GRCm39 独立两组双端 WGBS。DML 保留原始 p 与位点层 BH；DSS DMR 的区域聚合规则不等于区域 FDR。富集默认原始 p 探索性筛选。精子专题须显式启用并提供核验后的 ICR/基因资源。

报告写法见[统一报告规范](../docs/reporting.md)。

内部记录见[验证记录](docs/VALIDATION.md)。
