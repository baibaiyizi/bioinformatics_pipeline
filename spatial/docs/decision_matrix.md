# 空间转录组路线选择矩阵

这里把 `method_catalog.md` 的方法压缩成决策表。实际项目建议先跑默认轻量流程，再从下表按问题加模块。

## 按平台选择

| 平台/数据形态 | 首选输入容器 | 基础路线 | 推荐增强 | 谨慎使用 |
|---|---|---|---|---|
| 10x Visium FF/FFPE | AnnData/Squidpy 或 Seurat/SPE | 全局 QC + 局部 QC + log1p + Harmony + Squidpy spatial graph + spatial Leiden | SpotSweeper、BayesSpace、BANKSY、RCTD/Cell2location、COMMOT | 细胞级结论、单 spot 直接 DE |
| Visium HD | Seurat v5/AnnData + 可选 SpatialData | 8um bin 先做 atlas，必要时局部 2um | sketch clustering、BANKSY、Bin2cell、Sopa、SpatialData | 全量 2um 一步到位聚类、未复核的细胞级边界 |
| Slide-seq / Stereo-seq / Seq-Scope | AnnData/Squidpy/SPE | 坐标清洗 + 高密度 kNN + log1p/SCTransform + domain | SpaNorm、BANKSY、GraphST/STAGATE、SVG 多方法 | 使用 Visium 固定 6-neighbor 假设 |
| Xenium | SpatialData/Sopa/Seurat | cell-level QC + transcript/cell boundary review + Scanpy/Seurat annotation | Sopa segmentation review、NicheCompass、MISTy、COMMOT、CellChat spatial | 不检查 negative controls 和 cell area |
| CosMx | SpatialData/Sopa/Giotto | FOV 合并 + cell QC + panel gene review + cell type annotation | Sopa、Giotto、MISTy、LIANA/CellChat、niche analysis | 直接和 whole-transcriptome 数据同尺度比较 |
| MERSCOPE/MERFISH | SpatialData/Sopa/Giotto | molecule/cell QC + segmentation + cell-level clustering | SPArrOW、Baysor/Cellpose、NicheCompass、MISTy | 忽略探针 panel 和分割误差 |
| Spatial proteomics / CODEX / IMC | SpatialData/Giotto/MISTy | marker QC + cell segmentation + cell neighborhood | MISTy、NicheCompass/scNiche、topology features | 用 RNA deconvolution 工具强套 |
| 空间多组学 | SpatialData/Giotto/MEFISTO | 分模态 QC + 坐标统一 + shared metadata | MEFISTO/MOFA2、multi-view MISTy、Panpipes/同类 workflow | 先合并再 QC，丢失模态特异问题 |

## 按组织特点选择

| 组织/结构 | 必选关注点 | 推荐模块 | 备注 |
|---|---|---|---|
| 卵巢/卵泡 | follicle/object、oocyte-granulosa-theca composition、scaled distance | object segmentation、marker score、NMF/cNMF、MultiNicheNet/COMMOT | 参考 `ovarian_aging` 的 segment/object 思路 |
| 脑 | 层状结构、区域边界、神经/胶质 marker | BayesSpace/BANKSY、SVG、reference mapping、STAligner | 不同脑区应避免过度 batch correction |
| 肿瘤 | tumor-stroma-immune boundary、TLS、坏死/低 RNA 区 | domain/niche、deconvolution、COMMOT/CellChat、image feature | QC 必须区分技术伪影和真实坏死 |
| 肾 | 肾小球/小管/间质结构 | object/region segmentation、marker score、pseudo-bulk by compartment | 结构级汇总通常比 spot 级 DE 更稳 |
| 肝 | lobule zonation、血管/胆管邻域 | radial/vascular distance、SVG/pathway score、domain | 需要定义中心静脉/门管区参考 |
| 肠 | crypt-villus axis、免疫上皮邻域 | trajectory/axis score、object segmentation、cell-cell adjacency | 空间梯度优先于离散 domain |
| 皮肤 | 表皮-真皮层、毛囊/腺体 | layer/domain、object summary、image mask | H&E 层次信息很重要 |
| 胎盘/发育组织 | 时间点、谱系和空间迁移 | moscot、Spateo、trajectory、pseudo-bulk mixed model | 需要样本级设计支撑轨迹解释 |
| 植物/非模式组织 | marker 不充分、组织形态明显 | SVG、unsupervised domain、image/geometry features | annotation 先保持中性命名 |

## 按研究问题选择

| 目标 | 最低可交付 | 推荐增强 | 最终汇报应避免 |
|---|---|---|---|
| 数据质控 | manifest status、QC table、空间 QC 图 | SpotSweeper/局部 outlier、图像 mask、平台负控 | 只看 UMI/genes violin，不看空间位置 |
| atlas/图谱 | clusters/domains、marker、composition | reference label transfer、SVG、interactive review | 将未验证 cluster 直接命名为精细 cell type |
| 空间域发现 | spatial_domain、domain marker、sample composition | BANKSY/BayesSpace/STAGATE/GraphST 共识 | 只报单一 resolution 结果 |
| cell type 定位 | marker score 或 reference labels | RCTD/Cell2location/SPOTlight/Tangram | 把 spot-level proportion 当单细胞比例真值 |
| 细胞 niche | neighborhood enrichment、co-occurrence | MISTy、NicheCompass/scNiche、COMMOT | 把空间相邻直接等同于分子通讯 |
| 通讯分析 | contact edge + LR expression filter | COMMOT、CellChat spatial、LIANA、MultiNicheNet | 没有样本级重复时做强机制结论 |
| 组间差异 | sample-level composition + pseudo-bulk | mixed model、domain/cell-type stratified DE | spot 当独立生物重复 |
| gene program | NMF/cNMF program + spatial usage | WGCNA/topic model、pathway score、cross-sample consensus | 只用一个 program 数量且不做稳定性 |
| 图像-表达 | image overlay、ROI review | STPath/ST-Align/sCellST、histology embedding | 用预测表达替代实测表达做核心结论 |
| 3D/多切片 | per-slice QC + alignment plots | STAligner/PASTE/moscot | 未检查解剖一致性就解释 3D 结构 |

## 默认强弱分层

| 层级 | 方法 | 放入默认 pipeline 的理由 |
|---|---|---|
| 默认必跑 | manifest validation、h5ad checkpoint、QC、normalization、PCA/UMAP、spatial neighbors、domain、marker score、segment summary、NMF、result index | 依赖轻、组织通用、失败时容易定位 |
| 默认可开关 | Harmony、Squidpy ligand-receptor、local QC、object distance、pseudo-bulk | 常用但需要 metadata 或参数 |
| 项目增强 | SpotSweeper、SpaNorm、BANKSY、BayesSpace、RCTD、Cell2location、COMMOT、CellChat、MISTy | 价值高，但依赖、输入或解释成本更高 |
| 高级/前沿 | STAligner、PASTE、GraphST、STAGATE、NicheCompass、Spateo、STPath、MEFISTO、moscot | 适合特定问题，需单独验证和版本锁定 |

## 组织无关的最小闭环

1. `01_validate_inputs`：先保证每个样本有 platform、input_type、coords/image/reference 元数据。
2. `02_build_objects`：所有平台尽量转为 h5ad 或 SpatialData/SPE 可交换对象。
3. `03_qc_filter`：同时给全局 QC、空间局部 QC、人工复核图。
4. `04_normalize_features`：默认 log1p/HVG，必要时增加 SpaNorm/SCTransform 敏感性。
5. `05_integrate_cluster`：保留未校正矩阵，整合只用于聚类和可视化。
6. `06_spatial_neighbors_stats`：所有空间结论先经过邻域尺度检查。
7. `07_spatial_domains`：默认轻量 domain，再和可选方法互证。
8. `08-10_annotation_deconv_segment`：先粗注释，再 deconvolution/segment/object。
9. `11-13_program_communication_summary`：所有解释层输出 sample-level 表，避免只看单图。
10. `14_report_index`：把每个结果标注为 final、review 或 exploratory。

