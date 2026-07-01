# 空间转录组方法货架

本文件不是默认执行清单，而是后续裁剪 pipeline 时可选的方法库。默认路线仍保持轻依赖、可先跑通；这里把不同平台、组织、问题类型可能用到的模块尽量集中列出，便于按项目选择。

## 1. 数据容器和输入层

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| AnnData / Scanpy / Squidpy | Python 主流程、Visium、h5ad、常规空间图统计 | count matrix、obs、obsm["spatial"]、image | `.h5ad` checkpoint、空间图、统计表 | 低 | 对分子级/超大图像平台支持需要额外封装 |
| SpatialExperiment / SpatialFeatureExperiment / MoleculeExperiment | Bioconductor/R 生态、Rmd 复核、几何对象和分子级数据 | counts、spatialCoords、imgData、sf geometry | SPE/SFE/ME 对象 | 中 | Python/R 之间对象转换和版本管理 |
| SpatialData / spatialdata-io | Xenium、CosMx、MERSCOPE、Visium HD、高分辨率图像、多模态 | zarr、images、labels、points、shapes、tables | `.zarr` SpatialData store | 中 | 是基础设施，不直接替代分析方法 |
| Sopa | 大图像、多 FOV、Xenium/CosMx/MERSCOPE/Visium HD 预处理 | 平台原始输出、图像、分子点、细胞边界 | SpatialData、segmentation、cell table | 中-高 | 对存储、并行、图像坐标要求更高 |
| Giotto Suite | 需要跨技术、多尺度、R/Python/可视化一体生态 | 表达矩阵、空间坐标、图像、polygon | Giotto object、multi-scale output | 中 | 生态较大，建议项目单独固定版本 |
| Seurat v5 spatial | R 用户、Visium/Visium HD/Xenium 快速探索 | 10x 输出、h5ad/Seurat object | Seurat object、WNN/label transfer | 低-中 | 超大对象内存压力；Python 主流程需转换 |
| Voyager | Bioconductor 轻量空间探索、SVG 和空间统计 | SpatialExperiment | plots、Moran/Geary、SPE results | 低-中 | 高级分割/图像工作需其他包 |
| VoltRon / OSTA 生态 | 教学型或 Bioconductor 全栈工作流 | 多平台对象 | 标准化分析章节输出 | 中 | 需按包成熟度选择 |

## 2. 输入 QC、样本审计和空间伪影

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| 基础全局 QC | 所有平台第一关 | total counts、detected genes、mito/ribo percent | spot/cell pass flag | 低 | 全局阈值可能误删真实组织区域 |
| 空间局部 outlier QC | Visium、Slide-seq、HD binning | 空间 kNN + QC metrics | local low-count/high-mito flag | 低 | 边缘区域、真实低 RNA 区域需人工复核 |
| SpotSweeper | spot-based 空间 QC、dry spot、hangnail、区域伪影 | SpatialExperiment、QC metrics、coords | localOutliers、findArtifacts、dryspot/hangnail flags | 中 | 当前更偏 spot-based；单细胞平台需适配指标 |
| 图像组织掩膜 | 空白背景多、组织边界复杂 | H&E/IF image、spot/cell coords | in_tissue refined mask | 中 | 图像配准错误会系统性偏移 |
| 平台负控 QC | Xenium/CosMx/MERSCOPE | negative controls、blank probes、cell area、transcript counts | cell-level artifact flag | 中 | 阈值平台依赖强 |
| 批次/切片审计 | 多样本、多批次 | sample metadata、UMAP、QC distribution | batch report | 低 | 不应在 QC 阶段过度校正生物差异 |
| 组织结构一致性审计 | 发育、肿瘤、卵巢等结构明显组织 | histology/domain/cell-type maps | manual review checklist | 中 | 需要领域知识 |

## 3. 归一化、去噪和特征选择

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| library-size normalize + log1p | 默认快速路线 | raw counts | log-normalized matrix | 低 | 空间区域 RNA capture 差异可能被误认为生物差异 |
| SCTransform | Seurat/R 路线、UMI 数据 | raw counts | Pearson residuals / normalized assay | 中 | 参数和对象转换需统一 |
| Scran size factor | R/Bioconductor、细胞级平台 | raw counts、clusters | normalized counts | 中 | spot 混合组成强时需谨慎 |
| SpaNorm | 明显空间技术偏差、局部 capture 效应 | counts、coords | spatial-aware normalized matrix | 中 | 新方法，需和 log1p/SCTransform 做敏感性比较 |
| MAGIC / ALRA / DCA 等 imputation | 极稀疏数据、可视化辅助 | normalized expression | smoothed/imputed matrix | 中 | 下游 DE/通讯容易产生假阳性，默认不用于统计 |
| HVG | 常规聚类和降维 | normalized matrix | HVG list | 低 | 高空间结构基因不一定是高变基因 |
| SVG | 空间结构基因、组织分区、空间通路 | expression + coords | spatially variable genes | 中 | 方法差异大，建议至少两类方法交叉验证 |
| nnSVG | 大样本 SVG、Bioconductor 路线 | SPE/counts/coords | SVG ranking | 中 | 需 R 依赖 |
| SPARK-X / SpatialDE / Moran's I / Geary's C | SVG 互证 | expression + coords | SVG statistics | 低-中 | 坐标尺度、spot spacing、FDR 影响强 |

## 4. 降维、整合、配准和跨切片对齐

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| PCA + neighbors + UMAP | 基础探索 | normalized HVG matrix | embedding、graph | 低 | UMAP 只用于可视化，不作严格距离解释 |
| Harmony | 多样本 batch correction | PCA、batch key | corrected embedding | 低-中 | 可能移除真实组间差异 |
| BBKNN / Scanorama / scVI | 大量样本或复杂批次 | h5ad、batch key | corrected graph/latent | 中 | 需保留未校正矩阵做 DE |
| STAligner | 多切片空间配准、组织结构可对齐 | expression、coords、多切片 | aligned latent/coords | 中-高 | 解剖结构差异大时需人工检查 |
| PASTE / PASTE2 | 相邻切片、3D 重建、切片对齐 | pairwise slices、coords、expression | aligned slices、stack | 中 | 非相邻切片和强病变区域可能误配 |
| moscot / optimal transport | 发育时间、空间-时间映射 | timepoint/slice data | transport map、trajectory coupling | 高 | 解释依赖模型假设 |
| MEFISTO / MOFA2 | 多组学/多条件潜变量，空间连续因子 | multi-view matrix、covariates | latent factors | 中-高 | 需要清楚设定协变量和视图 |
| DOT / Tangram-like feature transfer | scRNA 与 ST 特征转移 | scRNA + ST | mapped cell states/features | 中 | 参考不匹配时会稳定地产生错误映射 |

## 5. 空间邻域、domain、niche 和组织结构

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| expression + coordinate Leiden | 默认轻量空间域 | PCA + scaled coords | spatial_domain | 低 | domain 数量和 spatial_weight 需敏感性分析 |
| Squidpy spatial neighbors | 邻域富集、共定位、Moran's I | coords、cluster labels | neighbor graph、enrichment、co-occurrence | 低 | 半径/kNN 选择决定结论尺度 |
| BayesSpace | Visium spot 聚类、subspot enhancement | SPE/Visium | spatial clusters、enhanced resolution | 中 | spot-based 假设，不适合所有平台 |
| BANKSY | 空间域、细胞类型关系、Visium HD/单细胞平台 | expression、coords | BANKSY embedding/clusters | 中 | lambda/k_geom 对结果影响大 |
| SpaGCN | 结合表达、空间和 histology | counts、coords、image | spatial domains | 中 | 图像权重和组织染色质量敏感 |
| STAGATE / GraphST | 图神经网络空间嵌入、整合、domain | expression、coords、多样本 | latent embedding、clusters | 中-高 | 深度模型需 seed/参数重复验证 |
| PRECAST | 多样本 spatial domain 和 batch | 多切片表达/coords | shared domains | 中 | 需要样本结构相似 |
| SpaceFlow | 空间轨迹和 domain | expression、coords | domains、pseudotime-like spatial flow | 中 | 轨迹解释需组织先验 |
| MISTy | 多视图空间关系、marker/细胞类型影响 | cell/spot features + neighbor views | interpretable spatial effects | 中 | 是关系建模，不直接证明机制 |
| NicheCompass | 大规模单细胞空间 niche | single-cell spatial omics、pathway priors | niche labels、latent space | 高 | 适合细胞级平台，spot 数据需谨慎 |
| scNiche | 单细胞空间 niche 识别 | cell coordinates、features | niche labels、characterization | 中-高 | 依赖细胞分割和 cell type 质量 |

## 6. 细胞类型注释、参考映射和 deconvolution

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| marker score | 无 scRNA 参考、快速初筛 | marker table、normalized expression | cell_type_pred / abundance score | 低 | 只能粗注释，混合 spot 不等于单细胞类型 |
| SingleR / Azimuth / Seurat label transfer | 参考和平台相近 | scRNA reference、ST object | transferred labels | 低-中 | reference batch 和组织状态偏差 |
| Cell2location | Visium/spot deconvolution，参考 scRNA 足够好 | scRNA + spatial counts | cell abundance per spot | 高 | GPU/运行时间高；参考缺失细胞会误分配 |
| RCTD / spacexr | 稳健 deconvolution、doublet mode | scRNA reference + ST | cell type weights | 中 | 细胞状态连续谱可能被硬分型 |
| SPOTlight | NMF-based deconvolution | scRNA + ST | cell type proportions | 中 | marker/program 质量影响大 |
| Stereoscope | 概率模型 deconvolution | scRNA + ST | proportions | 中-高 | 参数和 scale 需调 |
| Tangram | single-cell 到空间映射、细胞状态定位 | scRNA + ST | cell-to-space mapping | 中-高 | 映射不是绝对空间真值 |
| SpatialScope / Redeconve / SpatialcoGCN / STdGCN | 高级/深度 deconvolution | scRNA + ST + coords | abundance / simulated profiles | 高 | 属于可选前沿分支，需基准和复核 |
| CellTypist / scANVI | 单细胞级空间平台注释 | cell-by-gene matrix | cell labels | 中 | 需要平台/物种匹配参考 |

## 7. 高分辨率、细胞分割和图像模块

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| Bin2cell | Visium HD 从 bin 到细胞级近似 | 2um/8um bins、H&E/nuclei | cell-level expression estimates | 中 | 依赖图像和 binning，结果需和组织学核对 |
| Sopa segmentation | Xenium/CosMx/MERSCOPE/Visium HD、大图 | images、transcripts、platform output | cell boundaries、SpatialData table | 中-高 | 分割算法和参数强影响 cell matrix |
| Cellpose / Baysor / StarDist / Mesmer | 细胞核/细胞边界分割 | IF/H&E、molecule points | labels/polygons | 中-高 | 不同组织染色差异大 |
| SPArrOW | 分子级空间组学处理和分割 | molecule coordinates、images | cells、molecules、features | 中-高 | 更适合原位分子平台 |
| Xenium Explorer / Loupe / Vizgen Visualizer | 人工 QC 和 ROI | vendor output | curated regions、screenshots | 低 | 交互结果需结构化导出 |
| Histology feature extraction | H&E 纹理、组织区域辅助 | image tiles + coords | image embeddings/features | 中 | 图像特征可能捕捉批次染色 |
| STPath / ST-Align / sCellST | 图像-表达联合建模、前沿探索 | WSI/H&E + ST | predicted genes、joint embedding | 高 | 不作为默认统计证据，需外部验证 |
| HEST / STimage-like datasets | 模型预训练/benchmark | public WSI-ST paired data | benchmark/training resources | 高 | 数据域和目标组织可能不一致 |

## 8. 空间通讯、邻域互作和组织生态

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| 接触边矩阵 | 默认、无配体受体库也可用 | coords + cell/domain labels | group-to-group adjacency | 低 | 只说明相邻，不说明分子通讯 |
| Squidpy ligand-receptor | 快速 LR 富集 | expression + labels + LR resource | LR test table | 中 | 多重检验和表达 dropouts |
| COMMOT | 空间约束通讯、方向/距离建模 | expression + coords + LR pairs | communication score、sender/receiver | 中-高 | LR 数据库和距离参数决定结果 |
| CellChat spatial | 通讯网络、可视化丰富 | expression + labels + coords | pathway/network interaction | 中 | 分组太粗或太细都会偏 |
| LIANA+ | 多方法 LR 汇总和共识 | expression + labels | consensus LR rankings | 中 | 空间距离需额外纳入 |
| MultiNicheNet | 多条件差异 niche 通讯 | cell type abundance/expression + sample design | condition-specific LR program | 高 | 需要足够样本数和稳定 cell type |
| NicheNet / NATMI / CellPhoneDB | 非空间或弱空间通讯初筛 | expression + labels | LR ranking | 中 | 应作为补充，不应替代空间约束分析 |
| MISTy spatial effects | 非 LR 的邻域解释 | features + views | local/paracrine-like predictors | 中 | 解释为统计关系而非直接通讯 |

## 9. Gene program、轨迹、几何和 object-level 分析

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| NMF / cNMF | 稳健 gene program、跨样本 program | normalized expression | programs、usage | 中 | program 数量需稳定性选择 |
| LDA / topic model | 组织区域主题、病灶成分 | expression | topics/weights | 中 | 主题解释需 marker 支持 |
| WGCNA / hdWGCNA | 模块化基因网络 | pseudo-bulk/cell-level expression | modules、eigengenes | 中 | spot-level 稀疏性和样本数限制 |
| object/segment 汇总 | 卵泡、肿瘤巢、腺体、肾小球等结构 | spot/cell labels + components | object composition、size、distance | 低-中 | 分割错误会放大到 object 层 |
| scaled radial distance | 有中心-边缘结构的组织对象 | object mask + coords | normalized distance | 低 | 形状复杂时单一距离不足 |
| Spateo | spatiotemporal modeling、RNA velocity-like 空间建模 | high-res/cell-level ST | vector field、morphogenesis features | 高 | 前沿分支，数据质量要求高 |
| stLearn trajectory / spatial trajectory | 有空间连续梯度 | expression + coords + clusters | pseudotime/path | 中 | 轨迹不等于发育时间 |
| moscot | 跨时间点/空间最优传输 | multi-timepoint data | fate/transition coupling | 高 | 需要实验设计支持 |
| 拓扑/几何统计 | 腺体、肿瘤边界、免疫浸润结构 | polygons/labels/coords | density, boundary, topology features | 中 | 指标定义必须和生物问题一致 |

## 10. 差异分析、富集和统计建模

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| spot/cell-level Wilcoxon/t-test | 快速 review | expression + group labels | DE table | 低 | 伪重复，不能作为最终严肃统计 |
| sample-level pseudo-bulk | 多样本组间比较 | counts aggregated by sample/domain/cell type | edgeR/DESeq2/limma results | 中 | 每组样本数不足时功效低 |
| mixed model / dreamlet | 重复样本、病人/批次随机效应 | cell/spot/sample metadata | model coefficients | 中-高 | 模型设计要和实验设计一致 |
| spatialDE / spatial regression | 空间协变量和连续梯度 | expression + coords + covariates | spatial DE/SVG | 中 | 坐标和组织区域 confounding |
| Milo / DA-seq / neighborhood DA | 细胞状态邻域丰度差异 | graph + sample design | DA neighborhoods | 中 | 适合单细胞或高分辨率数据 |
| pathway score / PROGENy / AUCell / GSVA | 通路和状态评分 | expression + gene sets | pathway activity | 低-中 | gene set 与物种/平台需匹配 |
| GO/KEGG/Reactome/Enrichr | gene list/program 解释 | DE/SVG/program genes | enrichment table | 低 | 背景基因集必须合理 |
| Spatial pathway enrichment | 空间域/pathway maps | pathway scores + coords | spatial activity maps | 中 | 平滑和阈值会影响视觉结论 |

## 11. 可视化、报告和人工复核

| 选项 | 适用场景 | 输入 | 输出 | 成本 | 主要风险 |
|---|---|---|---|---|---|
| 每样本空间 QC 四联图 | 必做 | counts、genes、mito、in_tissue | QC plots | 低 | 不看原图容易误判边界 |
| UMAP + spatial map 成对展示 | 聚类/domain/annotation 复核 | embedding + coords | paired plots | 低 | UMAP 和空间结构需同时解释 |
| domain marker heatmap/dotplot | domain 命名 | marker expression | marker plots | 低 | marker 重叠导致过度命名 |
| composition barplot + alluvial | 多样本组间组成 | sample x label counts | composition plots | 低 | spot 数不等需归一化 |
| object-level violin/ridge | segment/object 分布 | object summary | object plots | 低 | object 数量和样本数要分开报告 |
| napari / Vitessce / iSEE | 交互式复核、ROI | SpatialData/SPE/AnnData | interactive views | 中 | 交互结论需导出表格和截图 |
| result README | 项目交付 | result folder | 输出解释索引 | 低 | 应区分真实结论和待复核标签 |

## 12. 验证和敏感性分析

| 检查 | 最低要求 | 增强要求 |
|---|---|---|
| QC 阈值敏感性 | 至少比较默认和宽松阈值 | 按组织区域分别评估删除比例 |
| 空间邻居尺度 | kNN 或半径至少两组参数 | 多尺度结果只保留稳健结论 |
| domain 稳定性 | resolution/spatial_weight 网格 | 多方法共识 domain |
| cell type annotation | marker heatmap + spatial map | scRNA 参考、原位 marker 或病理标注交叉验证 |
| deconvolution | 与 marker score 一致性 | RCTD/Cell2location/SPOTlight 至少两方法互证 |
| 通讯分析 | 先检查相邻关系和表达 | COMMOT/CellChat/LIANA/MultiNicheNet 共识和样本级统计 |
| 差异分析 | 保留工作 pvalue 表 | 正式报告使用 pseudo-bulk/mixed model/FDR |
| 图像/分割 | 人工抽样检查边界 | 多算法分割一致性或专家 ROI |
| 跨样本配准 | 可视化 landmarks/domain 对齐 | 对齐前后结论分别报告 |

## 推荐默认组合

| 项目类型 | 默认组合 | 可增强 |
|---|---|---|
| 普通 Visium atlas | AnnData/Squidpy + global/local QC + log1p + Harmony + spatial Leiden + marker score + object summary | SpotSweeper、BANKSY、BayesSpace、RCTD/Cell2location |
| Visium HD | Seurat/AnnData + 8um bin 起步 + sketch clustering + spatial domain + ROI | Bin2cell、BANKSY、Sopa/SpatialData、2um 高分辨率局部分析 |
| Xenium/CosMx/MERSCOPE | SpatialData/Sopa + cell-level QC + segmentation review + Scanpy/Squidpy + cell annotation | NicheCompass、MISTy、COMMOT、Vitessce/napari |
| 多切片/3D | 每切片基础 QC + STAligner/PASTE + shared domain | moscot、MEFISTO、object tracking |
| 病理/肿瘤微环境 | QC + domain/niche + deconvolution/cell annotation + neighborhood + COMMOT/CellChat | MISTy、MultiNicheNet、image embeddings |
| 发育/衰老 | 样本级统计 + object/segment + gene programs + spatial gradient | moscot、Spateo、trajectory、pseudo-bulk mixed model |

