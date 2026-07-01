# 空间转录组通用分析 pipeline

这个目录是一套尽量适配不同组织的空间转录组分析框架，结构参考本地 `wgbs_sperm`：根目录放统一入口，`01rawdata/` 放样本清单，`scripts/` 放确定性执行脚本，`rmd/` 放下游复核模块，`result/` 放 review-ready 的表格和图片。

设计参考：

- `HattieChungLab/ovarian_aging`：吸收了 spot-level h5ad、空间坐标补齐、Harmony 聚类、分割/segment 级汇总、oocyte/follicle object matching、scaled radial distance、cNMF gene program 和 MultiNicheNet 下游通讯分析的组织方式。
- AnySearch 检索：参考 OSTA/Bioconductor、Single-cell best practices、Squidpy、Scanpy、Seurat、Cell2location 等空间分析资料。检索摘要见 `docs/source_notes.md`。
- 本地生信项目习惯：支持单步重跑、活动 step 清理、日志记录、`result/<module>/tables/` 和 `result/<module>/plots/` 双输出。

## 目录结构

```text
Spatial/model/
  01rawdata/
    spatial_info.tsv              # 最小样本 manifest 模板
    spatial_info.extended.tsv     # 平台扩展 manifest 模板
  config/
    spatial_config.yaml           # 全局参数
    method_options.yaml           # 可选方法货架，不默认强制执行
    platform_profiles.yaml        # 不同平台的输入/QC/模块 profile
    marker_genes.tsv              # 通用和组织特异 marker
    segment_rules.tsv             # segment/object 规则模板
  scripts/
    spatial_pipeline.py           # Python 主分析 CLI
    run_rmd.sh                    # Rmd 复核入口
  rmd/
    01_main_setup.R               # 共享路径、配色、输出 helper
    02_*.Rmd - 07_*.Rmd           # 下游复核模块
  result/
    logs/                         # 运行日志
  docs/
    source_notes.md               # 外部参考记录
    method_catalog.md             # 空间分析方法全集和风险
    decision_matrix.md            # 按平台/组织/问题选路线
    platform_playbook.md          # 平台适配手册
  run_spatial.sh                  # 官方主入口
  environment.yml                 # conda 环境模板
```

## 输入 manifest

编辑 `01rawdata/spatial_info.tsv`。核心列：

- `sample`：唯一样本名。
- `group` / `condition` / `tissue` / `stage` / `batch`：下游比较和批次字段。
- `organism`：`mouse` 或 `human`，用于线粒体/ribosomal marker 默认前缀。
- `platform`：`visium`、`visium_hd`、`slide_seq`、`merfish`、`xenium` 等，主要用于记录。
- `input_type`：
  - `visium_dir`：`input_path` 指向 Space Ranger `outs/`。
  - `h5ad`：`input_path` 指向现成 `.h5ad`。
  - `matrix_coords`：`count_file` 提供表达矩阵，`spatial_file` 提供 barcode/x/y。
- `exclude`：`1` 跳过样本。

复杂平台或多组学项目建议从 `01rawdata/spatial_info.extended.tsv` 起步，额外记录 `subject_id`、`section_id`、`slide_id`、`capture_area`、`resolution_um`、`bin_size_um`、`panel_name`、`molecule_file`、`cell_boundary_file`、`segmentation_file`、`negative_control_file`、`manual_roi_file` 等字段。默认脚本仍读取最小 manifest；扩展 manifest 用于项目设计和后续平台分支。

## 方法货架

当前目录同时提供“可直接跑的默认流程”和“后续可选的堆砌方法库”：

- `docs/method_catalog.md`：按数据容器、QC、归一化、SVG、domain/niche、deconvolution、分割、通讯、gene program、统计和可视化整理方法，列出适用场景、输入输出、成本和风险。
- `docs/decision_matrix.md`：按平台、组织结构和研究问题选择路线，区分默认必跑、默认可开关、项目增强和高级前沿。
- `docs/platform_playbook.md`：Visium、Visium HD、Slide-seq/Stereo-seq、Xenium、CosMx、MERSCOPE、空间蛋白和空间多组学的输入字段、QC 指标和推荐分支。
- `config/method_options.yaml`：把重依赖工具作为开关记录，例如 SpotSweeper、SpaNorm、BANKSY、BayesSpace、STAGATE、GraphST、RCTD、Cell2location、COMMOT、CellChat、MISTy、NicheCompass、Bin2cell、Sopa、SpatialData、STAligner、PASTE、MEFISTO、moscot 等。
- `config/platform_profiles.yaml`：按平台定义输入类型、预期文件、空间邻居策略、QC 指标、默认模块和可选模块。

建议工作方式：先跑默认 14 步拿到稳定 QC 和基础空间结构，再按 `decision_matrix.md` 为具体组织挑选增强模块。

## 运行

```bash
cd /home/h1028/workspace/Spatial/model

# 查看步骤
bash run_spatial.sh --list

# 全流程
bash run_spatial.sh

# 只跑指定步骤，默认自动前置 01_validate_inputs
bash run_spatial.sh --steps 3,5-8

# 渲染下游 Rmd 复核
bash scripts/run_rmd.sh
```

建议先只跑：

```bash
bash run_spatial.sh --steps 1
```

确认 `result/01_validate_inputs/tables/spatial_input_status.tsv` 没有 error 后再跑全流程。

## Step 说明

| Step | 名称 | 主要作用 |
|---:|---|---|
| 01 | `validate_inputs` | 检查 manifest、输入路径、样本重复、输入类型。 |
| 02 | `build_objects` | 读取 Visium/h5ad/matrix+coords，写每样本和合并 `.h5ad`。 |
| 03 | `qc_filter` | 计算 counts/genes/mito/ribo QC，加入简单空间局部低 counts outlier 标记。 |
| 04 | `normalize_features` | normalize/log1p/HVG/regress/scale/PCA。 |
| 05 | `integrate_cluster` | 可用时 Harmony，之后 neighbors/UMAP/Leiden。 |
| 06 | `spatial_neighbors_stats` | Squidpy 空间邻居图、neighborhood enrichment、Moran's I SVG。 |
| 07 | `spatial_domains` | 表达 PCA + 空间坐标的 Leiden 空间域。 |
| 08 | `cell_type_annotation` | marker score 细胞类型/区域注释，可按组织扩展 marker。 |
| 09 | `deconvolution` | 默认 marker-score abundance；Cell2location/RCTD/SPOTlight 作为显式可选分支。 |
| 10 | `segment_objects` | 从空间域或注释构建 connected spatial objects，输出 segment composition 和 scaled distance。 |
| 11 | `programs_differential` | NMF gene program 和两组 spot-level working DE。 |
| 12 | `neighborhood_communication` | 空间接触边矩阵；可开启 Squidpy ligand-receptor。 |
| 13 | `multi_sample_summary` | 汇总 sample x cluster/domain/cell type/segment。 |
| 14 | `report_index` | 生成 `result/README.md` 输出索引。 |

## 组织适配方式

默认路径是组织无关的。换组织时优先改这三处：

1. `config/marker_genes.tsv`：增加对应组织的 cell type marker。
2. `config/spatial_config.yaml` 的 `project.default_tissue`：例如 `ovary`、`brain`、`liver`。
3. `config/segment_rules.tsv` 和 `segments.source_key`：决定 segment/object 是基于 `spatial_domain`、`cell_type_pred`，还是外部已有 `segment_id`。

如果已有病理/组织学人工分割，建议把 `segment_id` 直接放入 h5ad 的 `obs`，再从 step 10 开始用它做 object-level 汇总。

## 统计和出图约定

- 默认保留完整表，再输出 review 表。
- 探索性 review 默认使用 `pvalue`，正式报告需要时再提升到 `padj/FDR`。
- 图像尽量同时保存 `pdf` 和 `png`。
- 下游 Rmd 输出统一进入 `result/<module>/tables/` 和 `result/<module>/plots/`。
- 空间域、segment、marker-score cell type 都是计算标签，必须结合组织图像、marker 和领域知识复核。

## 依赖

```bash
conda env create -f environment.yml
conda activate spatial-model
```

`cell2location`、`BayesSpace`、`BANKSY`、`SpaGCN`、`COMMOT`、`MultiNicheNet` 等高算力或生态较重的工具建议按项目需要单独固定版本。默认 pipeline 不强制这些工具，避免一个组织的特殊依赖阻塞所有样本的基础 QC 和空间统计。

新增方法货架里的工具也遵循同一原则：先记录可选路线和输入要求，不在基础环境里一次性安装所有重依赖。每次只打开一个增强分支，跑完后把表格和图纳入 `result/<module>/tables/` 与 `result/<module>/plots/` 再比较。
