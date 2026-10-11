# 模板能力与执行证据矩阵

目录更新日期：2026-10-11。非 smallRNA 方法验收保留 2026-09-29 的执行边界；smallRNA 当前记录见其 VALIDATION.md。下表依据实际入口和验证记录，不以目录存在或 agent 知识范围定义能力。

| 模板 | 支持的输入与主要行为 | 推断单位及限制 | 真实数据执行边界 |
|---|---|---|---|
| rnaseq | FASTQ 上游；显式 counts 下游；DE、富集、剪接与可选专题 | 独立两组；方向显式；全体预先合格样本；P 候选为探索性，保留 padj | 小数据 counts 闭环；FASTQ、APA、真实 rMATS 等完整重跑 `not_run` |
| wgbs | 小鼠 GRCm39 Bismark/CpG→QC/DMR/注释；可选精子与跨组学 | 独立样本；DSS 区域规则不等同于区域 FDR；精子坐标须外部核验 | 语法、身份/方向和隔离小数据验证；真实原始 WGBS `not_run` |
| atac | QC、比对/峰、consensus counts、DA、注释/富集/轨迹；可选 motif/footprint | 独立样本两组；现行 IDR 重复数限制保留；监督 ML 不启用 | 真实来源上游曾有 `MISSING/STALE`，本轮没有补跑；模板验证不能覆盖此状态 |
| scrnaseq | FASTQ/Cell Ranger 或 counts H5AD；QC/图谱/marker/描述；可选 pseudobulk | library 与 biological_sample_id 分离；同一样本多库聚合；每组至少 3 独立单位且设计合格才推断 | 合成端到端与实际 PyDESeq2 验证；Cell Ranger 真实计数 `not_run`；胰岛扩展为单独示例 |
| smallrnaseq | 新版 PANDORA：SPORTS/endo 输入→describe；显式 differential 与 enrichment | 默认描述；物种/endo/设计显式；单库、单组不补对照；serum 两队列分离 | 当前跨物种合成与真实数据代表性检查见模块 VALIDATION；未重跑完整 FASTQ 上游 |
| metabolomics | 显式 mzML manifest/XCMS 或显式峰表→QC/limma/注释/报告 | 一行一独立生物学样本，每组≥3；缺填峰前表不计算原始检出率；pooled QC 身份门控 | 60特征×6样本峰表闭环、反向比较和负向输入；真实 XCMS `not_run` |
| spatial_transcriptomics | 空间对象/QC/域/marker；counts 来源的 NMF 与样本描述 | 删除 spot 条件检验；组间推断 `not_implemented`；marker score 不是丰度 | 8项测试及合成H5AD 11步闭环；Visium、Harmony/Squidpy/Cell2location `not_run` |
| singlecell_multiome | RNA/ATAC count-level h5ad、paired 核对、QC、TF-IDF/LSI、可选 MuData | embedding 不增加生物学重复；没有正式 DE/DA 执行模块 | 既有2项测试通过；fragments 上游、MultiVI 未实现/未运行 |
| spatial_metabolomics | 长表/imzML入口、像素QC、分箱、样本内空间统计与ROI描述 | 像素/ROI不是独立样本；低等级注释不作 level 1 声明 | 既有2项长表/注释门测试通过；真实 imzML `not_run` |

## 证据等级

- **源码审计**：读取/比较当前实现，不能证明运行成功。
- **语法检查**：Shell/R/Python可解析，不证明算法或输入有效。
- **小数据闭环**：实际执行注明的输入模式及模块，边界见各模板 `VALIDATION.md`。
- **真实数据复现**：本轮未启动真实项目全量重算，不据旧输出或文件夹名称宣称完成。

ML、网络药理、分子建模及 AI harness 不在本轮范围。其历史成熟度不作为当前新验证结论；参见各自目录文档。
