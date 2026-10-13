# WGBS 单组学分析入口

甲基化 QC、DSS DML/DMR、注释、富集、基因组背景、纯基因甲基化汇总和浏览器轨道。01–08 为甲基化主线，09 为 browser；精子模块另有 10 精子专题。方法与结果目录使用相同的连续编号。

## 运行

```bash
bash run_wgbs.sh --list
bash run_wgbs.sh --check
bash run_wgbs.sh --steps 08,09
bash run_wgbs.sh --plots-only --check
bash run_wgbs.sh --plots-only
Rscript scripts/run_gene_methylation.R .
```

路径相对本分析模块根解析，不受启动目录影响；公共参考资源和工具通过相对配置、PATH 或显式环境变量提供。

## 当前产物

`result/03_dmr/02_annotation/tables/dmr_annotated.csv` 为区域证据；`result/08_gene_methylation/01_group_comparison/tables/gene_methylation_long.csv` 和宽表按基因、功能层汇总。每层保留绝对差异最大的区域，并保存 RegionCount 和 DirectionConflict；这不是基因级显著性检验。

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

## 冻结数值重绘

`--plots-only` 在独立分支读取已完成结果和冻结对象，不调用统计 Rmd、不补跑 DSS/fgsea/分段拟合。原数值 SHA、图件覆盖与 PDF/PNG 规格通过后才原子替换 result；失败保留当前有效结果。03 分 `01_calling/02_annotation`，06 分 `01_profiles/02_segmentation`，其余实际对比分 `01_group_comparison`；每层并列 plots/tables，09 另有 tracks。PDF/PNG 同 stem，183 mm、300 dpi、RGB 白底及 Cairo PDF；本地 figure_style.R 是中心规范的冻结副本。

缺失保存模型的 MethylSeekR 诊断图明确保留历史身份，不为重绘额外拟合。复制后维护本地 template_snapshot.tsv，运行不实时导入中心模板。

## 统一绘图配置

`config/figure_style.json` 是本模板随项目冻结的语义配置。复制后编辑本地配置，运行不读取 PFOS 或中心模板。R 实现通过本地 `figure_style.R` 的明确相对路径加载；分类按名称绑定，不按输入行序轮换。ORA 使用 GeneRatio、Count 面积及 −log10(P) viridis；色阶 0–10，端点标 ≥10。相关系数固定 [−1,1]；绝对甲基化固定 [0,100]%；Row Z-score 固定 [−2.5,2.5] 并标饱和；缺失值单独显示灰。不同测量量保留独立图例与单位。

原生 methylKit 直方图仅替换已核验文库的组色，保留分箱、频数与百分比；PCA/树图采用包内 `methylKit:::rowSds`，逐项验证与原生坐标/树高一致。GO 交集计数使用中性色，分区结构 UMR/LMR/PMD 使用按名称固定的类别色。

Bismark 链方向组成图从各文库四个已显示、strand-assigned 类别count计算百分比，排除无法提取基因组序列的discarded pairs；原directional模式未执行的互补方向为0，不解释为生物缺失或转换效率。GO树的局部cluster枝条显式按 注册 `clusters` 调色板按cluster数字编号 固定身份色，节点P、Count和原树结构保留。
