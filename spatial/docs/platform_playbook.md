# 平台适配手册

本手册面向后续落地数据时的输入整理。原则是先把平台差异收敛到统一 manifest，再在 config 里打开特定分支。

## 通用 manifest 策略

基础 `01rawdata/spatial_info.tsv` 只保留必需列。复杂项目使用 `01rawdata/spatial_info.extended.tsv`，推荐包含：

- 样本设计：`sample`、`subject_id`、`group`、`condition`、`batch`、`section_id`、`slide_id`、`capture_area`、`replicate`。
- 组织信息：`tissue`、`organism`、`stage`、`region`、`disease_state`、`pathology_score`。
- 平台信息：`platform`、`chemistry`、`resolution_um`、`bin_size_um`、`panel_name`、`probe_panel_size`。
- 输入路径：`input_type`、`input_path`、`count_file`、`spatial_file`、`image_file`、`molecule_file`、`cell_boundary_file`、`segmentation_file`。
- 参考信息：`reference_h5ad`、`reference_cell_type_key`、`marker_set`、`genome`。
- QC 辅助：`negative_control_file`、`in_tissue_mask`、`manual_roi_file`、`exclude`、`notes`。

## 10x Visium

| 项 | 建议 |
|---|---|
| `input_type` | `visium_dir` |
| 必需路径 | Space Ranger `outs/` |
| 默认容器 | AnnData/Squidpy；R 复核可转 SpatialExperiment |
| 必跑 QC | counts、genes、mito/ribo、in_tissue、空间低 counts outlier、图像 overlay |
| 默认分析 | log1p/HVG/PCA、Harmony、Leiden、Squidpy graph、spatial domains、marker score |
| 推荐增强 | SpotSweeper、BayesSpace、BANKSY、RCTD/Cell2location、COMMOT |
| 交付重点 | 每样本空间 QC 图、domain marker、组间 composition、segment/object summary |

关键判断：如果组织区域结构明显，优先做 object/segment 汇总；如果多组比较，正式统计优先 sample-level pseudo-bulk。

## 10x Visium HD

| 项 | 建议 |
|---|---|
| `input_type` | `visium_hd_dir` 或 `h5ad` |
| 必需路径 | Space Ranger HD 输出、H&E、bin matrix |
| 推荐起点 | 8um bin 起步，局部或目标区域再看 2um |
| 默认容器 | Seurat v5 或 AnnData；高分辨率图像建议 SpatialData |
| 必跑 QC | bin counts/genes、组织 mask、bin-level spatial outlier、图像对齐 |
| 推荐增强 | sketch clustering、BANKSY、Bin2cell、Sopa、SpatialData/napari |
| 交付重点 | bin-size 对比、ROI 局部高分辨率图、细胞级推断复核 |

关键判断：不要默认把 2um bin 当作真实细胞；细胞级输出必须注明来自分割/聚合算法。

## Slide-seq / Stereo-seq / Seq-Scope

| 项 | 建议 |
|---|---|
| `input_type` | `matrix_coords` 或 `h5ad` |
| 必需路径 | matrix、barcode/bead/cell coordinates、可选 image |
| 默认容器 | AnnData/Squidpy 或 SPE |
| 必跑 QC | 坐标范围、重复坐标、bead/cell density、局部 counts、背景区域 |
| 推荐分析 | 高密度 kNN/radius graph、SpaNorm、SVG、BANKSY/GraphST/STAGATE |
| 交付重点 | 邻域尺度敏感性、空间域稳定性、分辨率说明 |

关键判断：不要沿用 Visium 固定 6 邻居；应按物理距离和点密度设定 k/radius。

## Xenium

| 项 | 建议 |
|---|---|
| `input_type` | `xenium_dir`、`spatialdata_zarr` 或 `h5ad` |
| 必需路径 | cell feature matrix、transcripts、cell boundaries、morphology image、negative controls |
| 默认容器 | SpatialData/Sopa；Seurat v5 可做 R 探索 |
| 必跑 QC | transcripts/cell、genes/cell、cell area、negative controls、segmentation boundary、FOV/batch |
| 推荐分析 | cell-level clustering、reference annotation、neighborhood enrichment、MISTy、COMMOT/CellChat spatial |
| 交付重点 | cell type spatial map、niche、ROI、boundary QC |

关键判断：Xenium 是 panel/platform-specific 数据，基因 panel 会限制 de novo program 和 pathway 解释。

## CosMx

| 项 | 建议 |
|---|---|
| `input_type` | `cosmx_dir`、`spatialdata_zarr` 或 `h5ad` |
| 必需路径 | expression matrix、metadata、FOV coordinates、cell polygons、image |
| 默认容器 | SpatialData/Sopa/Giotto |
| 必跑 QC | FOV 批次、cell area、negative probes、panel gene coverage、segmentation |
| 推荐分析 | cell type annotation、niche、MISTy、Giotto spatial network、LIANA/CellChat |
| 交付重点 | FOV-level QC、panel-aware marker review、cell neighborhood |

关键判断：跨 FOV 合并时必须保留 FOV/batch 信息，避免把 FOV 效应解释为组织区域。

## MERSCOPE / MERFISH

| 项 | 建议 |
|---|---|
| `input_type` | `merscope_dir`、`molecule_table`、`spatialdata_zarr` |
| 必需路径 | molecule coordinates、cell boundaries、DAPI/IF image、cell-by-gene matrix |
| 默认容器 | SpatialData/Sopa 或 Giotto |
| 必跑 QC | molecule density、blank/control probes、cell volume/area、segmentation confidence |
| 推荐增强 | SPArrOW、Baysor/Cellpose、NicheCompass、MISTy |
| 交付重点 | cell-level neighborhood、niche、single-cell spatial statistics |

关键判断：分子级数据应保留 molecule table，不要只保留聚合矩阵。

## 空间蛋白或多重 IF/CODEX/IMC

| 项 | 建议 |
|---|---|
| `input_type` | `spatialdata_zarr`、`image_features` 或平台导出表 |
| 必需路径 | cell marker intensity、cell coordinates、segmentation labels、image |
| 默认容器 | SpatialData/Giotto/MISTy |
| 必跑 QC | marker intensity distribution、batch/channel QC、segmentation |
| 推荐分析 | cell type gating/annotation、neighborhood、MISTy、topology features |
| 交付重点 | marker panel、cell type gating logic、spatial neighborhood |

关键判断：这是空间组学而不一定是转录组，RNA 专用 deconvolution 不应强套。

## 空间多组学

| 项 | 建议 |
|---|---|
| `input_type` | `spatialdata_zarr` 或每模态独立 manifest |
| 必需路径 | 每个模态的 counts/features、coords、image/geometry、shared sample metadata |
| 默认容器 | SpatialData 或 Giotto |
| 必跑 QC | 每模态单独 QC、坐标系统一、barcode/cell 对齐 |
| 推荐增强 | MEFISTO/MOFA2、multi-view MISTy、feature transfer、Panpipes/项目化 workflow |
| 交付重点 | 模态一致性、模态特异信号、联合 latent factor |

关键判断：多组学项目先按模态分别跑通，再做整合；不要让联合模型掩盖单模态 QC 问题。

## 平台扩展到 pipeline 的建议顺序

1. 先在 `platform_profiles.yaml` 增加平台的输入字段、QC 字段和默认模块。
2. 再在 `method_options.yaml` 打开一个可选方法，不要一次打开多个重依赖方法。
3. 每增加一种平台，先让 `01_validate_inputs` 能识别并输出 clear warning。
4. 每增加一种输入容器，必须有一个 checkpoint 输出：h5ad、SpatialData zarr、SPE rds 三者至少一个。
5. 每增加一种高级方法，必须同时增加一个 review Rmd 或表格汇总，避免结果只存在对象内部。

