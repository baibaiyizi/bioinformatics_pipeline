# SCI 分析能力矩阵与流程审计

更新日期：2026-08-24

## 1. 审计口径

本矩阵审计的是当前工作区内**能够定位到入口、环境、测试和产物契约的真实流程**，而不是 agent 的知识范围。成熟度定义如下：

- `ready`：有明确入口和输入契约，可在已声明环境中重跑，存在自动测试或经过可核验的小数据闭环。
- `partial`：已有可用流程，但缺少通用化、版本锁定、自动测试、某类关键验证或原始数据到终稿产物中的一段。
- `absent`：工作区中没有该方法族的可运行实现。
- `external-only`：只能调用外部平台或服务，尚无本地闭环。
- `resource-limited`：方法定义清楚，但当前环境或算力不足以完成真实规模运行。

缺陷分级：

- `critical`：可能改变样本身份、比较方向、统计单位或核心结果。
- `major`：可能显著降低结论可信度、验证强度或可复现性。
- `minor`：主要影响易用性、报告完整度或跨项目复用。

所有“支持”均受统一分析契约约束：生物学重复是推断单位；细胞、spot 和像素不能替代样本重复；预测流程的所有数据驱动变换必须位于重采样内部；低等级代谢物注释、网络中心性、docking score 和单条 MD 轨迹均不能升级为确定性机制结论。

## 2. 能力总览

| 领域 | 当前成熟度 | 真实入口 | 当前强项 | 主要缺口 | 缺陷等级 | Owner |
|---|---|---|---|---|---|---|
| bulk/total RNA 与样本级 Smart-seq | `partial` | `model/bulkrnaseq/run_rmd.sh`；各 `rna_seq/*/run_rmd.sh` | FASTQ 上游、样本表、QC、DE、富集、GSEA/GSVA、WGCNA、APA、剪接及报告模块已经存在 | `model/bulkrnaseq` 没有锁定环境和自动 smoke test；total/small RNA 并未共享同一输入契约；部分项目阈值仍以未校正 P 值为默认口径 | `major` | `sci_transcriptomics` + `sci_analyst` |
| scRNA | `ready`（PFOS 项目） | `sc/PFOS/run_scrna.sh`；`python -m pfos_scrna` | Python 包、YAML/TSV 配置、AnnData 数据契约、样本级统计与 pseudobulk、合成数据及集成测试齐全 | 当前实现是 PFOS 小鼠胰岛项目，不应宣称为任意物种/平台的通用流程；高级 scVI/velocity/通讯依赖是可选项 | `minor` | `sci_singlecell` + `sci_analyst` |
| scATAC / 单细胞 multiome | `partial` | `model/singlecell_multiome/scripts/run_multiome.py` | count-level AnnData 契约、paired barcode/sample 核对、模态分离 QC、ATAC TF-IDF/LSI、可选 MuData、smoke 与负向测试 | 原始 fragments 到 peak matrix、TSS/FRiP/doublet、基因活性/peak-gene 尚未实现；MuData/MultiVI 依赖当前缺失 | `major` | `sci_singlecell`；整合后 `sci_multiomics` |
| 空间转录组 | `partial` | `model/spatial_transcriptomics/run_spatial.sh` | 14 步 runner、平台 profile、输入验证、QC、domain、注释、去卷积、邻域和多样本汇总 | 无自动测试；环境只有范围约束而无 lock；部分方法是可选依赖；必须逐项目验证样本级而非 spot 级推断 | `major` | `sci_spatialomics` + `sci_analyst` |
| 空间代谢组 | `partial` | `model/spatial_metabolomics/scripts/run_spatial_metabolomics.py` | imzML manifest/长表契约、像素 QC、m/z 分箱、TIC 标准化、样本内 Moran's I、ROI 样本汇总、鉴定等级保护和测试 | 当前环境缺 pyimzML/Cardinal/SpatialData；尚缺原始 profile/centroid 专属峰处理、组织图像配准和真实 imzML fixture | `major` | `sci_metabolomics` → `sci_spatialomics` |
| 普通 LC-MS 代谢组 | `partial` | `Metabolism/model/run.sh` | WIFF/mzML 审计、XCMS、检出率/QC、漂移保护、limma、注释等级、方向约定和中文报告 | 仅覆盖当前非靶向 LC-MS 项目；无自动测试和锁定环境；QC 身份未确认时会限制漂移校正；GC-MS、靶向定量尚未覆盖 | `major` | `sci_metabolomics` + `sci_analyst` |
| bulk ATAC | `partial` | `atac/liver/liver_atac.sh`；`atac/liver/run_rmd.sh` | 原始 QC、比对、过滤、FRiP/TSS、peak/consensus count、DA、motif、footprinting、track 与 RNA/WGBS 关联 | 项目特定且缺自动测试/锁定环境；blacklist、TSS 和部分工具缺失时会跳过；跨样本 peak 与统计设计需逐项目审计 | `major` | `sci_epigenomics` + `sci_analyst` |
| WGBS | `partial` | `model/wgbsseq/run_wgbs.sh`；`wgbs/*/run_wgbs.sh` | Bismark 上游资产、QC、DMR、基因组上下文、调控、浏览器轨迹和跨组学模块 | 缺统一锁定环境、自动小数据测试与覆盖度/重复层级验收表；多个项目副本存在漂移风险 | `major` | `sci_epigenomics` + `sci_analyst` |
| CUT&Tag / CUT&RUN | `absent` | 无 | 可复用 bulk ATAC 的部分区间和 track 统计思想 | 缺 assay 特异的比对、片段/QC、背景/对照、peak 与 replicate 处理 | `major` | `sci_epigenomics` |
| bulk 多组学 | `partial` | `atac/liver/rmd/09_rna_multiomics.Rmd`、`11_wgbs_integration.Rmd`；WGBS `8_multiomics.Rmd` | 已有项目内 RNA—ATAC—WGBS 关联 | 缺统一样本映射、缺失模态处理、样本级验证、MOFA2/DIABLO 选择契约和可复现测试 | `major` | 单模态 owner → `sci_multiomics` → `sci_analyst` |
| 表格型预测模型 | `partial` | `model/ml/scripts/run_ml.py`；旧项目入口如 `atac/liver/rmd/14_machine_learning.Rmd` | 冻结特征契约、分组嵌套 CV、重采样内插补/标准化/筛选、分类/回归基线、校准、决策曲线、外部评估、模型导出和测试 | 当前最小入口仅支持数值特征二分类/回归；生存、深度学习、区间估计和 SHAP 依赖/专用实现仍待补 | `minor`（最小入口）；旧代码用于主张前仍是 `critical` | `sci_ml` + `sci_analyst` |
| 网络药理学 | `partial` | `model/network_pharmacology/scripts/run_network.py` | 版本化本地实体/边表、物种/ID/重复证据检查、证据加权、中心性、删边稳健性、候选证据卡和测试 | 尚未建立 Open Targets/STRING/PrimeKG 的版本化下载适配器与通路背景分析；真实证据需要逐来源核验 | `major` | `sci_evidence` + `sci_network_pharmacology` |
| 分子对接 | `resource-limited` | `model/molecular_modeling/scripts/preflight.py` | PDB/配体/盒参数、redocking、阳性/诱饵、多 seed 与 PoseBusters 依赖的强制 preflight 和负向测试 | 本机缺 Vina、Open Babel 与 PoseBusters，故真实 docking 阶段保持 `blocked` | `major` | `sci_molecular_modeling` |
| 分子动力学 | `resource-limited` | `model/molecular_modeling/scripts/preflight.py` | 力场、水模型、配体参数化、至少三次独立重复、seed 与主张边界的 preflight | 本机缺 OpenMM/MDAnalysis 和已审计 GPU；尚未运行最小化、平衡、生产或轨迹收敛测试 | `major` | `sci_molecular_modeling` |

## 3. 可复现性与跟踪审计

| 检查项 | 结果 | 证据或处置 |
|---|---|---|
| 空间 Python 主脚本被 Git 跟踪 | `pass` | `model` 是独立 Git 仓库；`spatial_transcriptomics/scripts/spatial_pipeline.py` 已在索引中 |
| 空间 YAML 配置被 Git 跟踪 | `pass` | `config/spatial_config.yaml`、`method_options.yaml`、`platform_profiles.yaml` 等均已跟踪 |
| 空间环境文件被 Git 跟踪 | `pass` | `spatial_transcriptomics/environment.yml` 已跟踪 |
| scRNA Python、配置和测试被 Git 跟踪 | `pass` | `sc/PFOS` 的局部 `.gitignore` 已窄化放行 `src/`、配置与测试 |
| 其他新 Python 管线可被跟踪 | `pass with action` | 新代码放入独立 `model` 仓库；不扩大顶层 catch-all 白名单 |
| bulk RNA / ATAC / WGBS 锁定环境 | `fail` | 当前主要依赖服务器环境与脚本约定，缺项目级 lock/容器摘要；列为 `major` |
| 自动小数据测试 | `mixed` | scRNA 有 pytest；新增 ML、multiome、空间代谢组、网络药理和分子建模 preflight 均有 smoke/negative test；空间转录、普通代谢组、ATAC、WGBS 与 bulk RNA 仍缺 |
| 随机种子与运行清单 | `mixed` | 各项目实现不一致；所有新闭环必须输出 seed、软件版本、参数和输入校验和 |

## 4. 方法参照与审计结论

- bulk RNA、ATAC 与 WGBS 的验收层分别对照 nf-core/rnaseq、nf-core/atacseq 和 nf-core/methylseq 的 samplesheet、原始 QC、参考版本、replicate、聚合报告和可测试 profile 思路；不复制第三方完整仓库。
- 单细胞流程以 Scanpy/AnnData 为基础；需要概率整合时使用 scvi-tools。scATAC/multiome 只有在各模态先独立完成 QC 后，才允许建立 MuData 并进入 MultiVI 等联合模型。
- 空间对象采用“图像/labels/points/shapes 与注释表分离、由坐标系显式关联”的 SpatialData 思路；空间代谢组的光谱预处理和空间推断必须留下分层中间产物。
- 临床预测模型的报告对照 TRIPOD+AI，偏倚与适用性审查对照 PROBAST+AI；二者不替代嵌套验证、校准和外部评估本身。
- 网络药理优先保留 Open Targets 的 evidence/association 层级、STRING 的映射物种和版本化端点、PrimeKG 的来源版本；任何中心性排名只是候选优先级。
- 对接需至少保留输入结构、盒参数、seed、redocking/对照与 PoseBusters 类几何合理性检查；MD 需报告独立重复和收敛，不能用单次短轨迹证明稳定性或机制。

## 5. 固定修复优先级

1. `completed`：通用、无泄漏的 `model/ml/` 最小闭环。
2. `completed (partial scope)`：`model/spatial_metabolomics/` 的输入审计、像素 QC、feature/coordinate 中间层、ROI 和样本内空间统计；真实 imzML 验证待依赖补齐。
3. `completed (partial scope)`：`model/singlecell_multiome/` 的 count matrix 契约、模态分离 QC、TF-IDF/LSI 和可选 MuData；fragments 上游与 MultiVI 待补。
4. `completed (local evidence graph)`：`model/network_pharmacology/` 的版本化本地证据图、稳健性和候选证据卡；数据库下载适配器待补。
5. `completed (preflight only)`：`model/molecular_modeling/` 的结构/对照/多 seed/多 replica 审计；真实 docking/MD 因依赖和资源保持 `resource-limited`。
6. 下一优先级：按严重度为 bulk RNA、空间转录、普通代谢组、ATAC 和 WGBS 增加环境锁、fixture 与 smoke test，并补 scATAC fragments 入口。

## 6. 公开实现参考

- [nf-core/rnaseq](https://nf-co.re/rnaseq/)
- [nf-core/atacseq](https://nf-co.re/atacseq/)
- [nf-core/methylseq](https://nf-co.re/methylseq/latest/)
- [Scanpy](https://scanpy.readthedocs.io/)
- [scvi-tools / MultiVI](https://docs.scvi-tools.org/en/stable/user_guide/models/multivi.html)
- [SpatialData](https://spatialdata.scverse.org/en/stable/)
- [Cardinal](https://bioconductor.org/packages/release/bioc/html/Cardinal.html)
- [TRIPOD+AI](https://www.bmj.com/content/385/bmj-2023-078378)
- [PROBAST+AI](https://www.bmj.com/content/388/bmj-2024-082505)
- [Open Targets Platform](https://platform-docs.opentargets.org/)
- [STRING API](https://string-db.org/help/api/)
- [PrimeKG](https://github.com/mims-harvard/PrimeKG)
- [AutoDock Vina](https://autodock-vina.readthedocs.io/en/stable/)
- [DiffDock](https://github.com/gcorso/DiffDock)
- [PoseBusters](https://github.com/maabuu/posebusters)
- [OpenMM](https://docs.openmm.org/latest/userguide/)
- [MDAnalysis](https://docs.mdanalysis.org/)
