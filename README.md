# 分析模板：只复制需要的模块

这里保存从实际项目提取、按输入和统计边界整理的分析模板。复制某一个目录到新项目，修改该目录的配置和样本表后使用；无需复制整个仓库。实际课题统一放在 `../projects/`。各模板保留原有工具与模块结构。

## 创建一个新分析模块
从本模板目录运行，例如：

```bash
python3 copy_template.py rnaseq ../projects/new_study/mouse/rnaseq
```

脚本拒绝覆盖现有目录，记录当前来源 commit、未提交状态和实际文件 SHA256 到 `template_snapshot.tsv`。复制后独立运行，中心更新不自动修改课题。先填配置和样本表，再执行对应模块的只读预检；占位值不能运行生产分析。

配置中的相对路径以模块根为基准，smallRNA 以选定数据集根为基准；启动时所在目录不影响输入位置。参考库和公共工具配置为明确的相对路径，外部命令通过已验证环境或明确工具变量提供。运行审计中的实际绝对位置用于追溯。

## 选择模板

| 分析 | 目录 | 项目配置 | 入口 | 本次验证范围 |
|---|---|---|---|---|
| bulk RNA / 样本级 Smart-seq | [rnaseq](rnaseq/README.md) | `config/project.env`、`samples.tsv` | 三段 `rnaseq_*.sh`、`run_rmd.sh` | counts 下游与方向/身份；FASTQ 全流程未重跑 |
| WGBS | [wgbs](wgbs/README.md) | `config/project.env`、`samples.tsv` | `run_upstream.sh`、`run_wgbs.sh` | GRCm39 两组核心、预检与小数据；真实上游未重跑 |
| bulk ATAC | [atac](atac/README.md) | `config/project.env`、`samples.tsv` | `run_upstream.sh`、`run_rmd.sh` | 两组核心、配置和隔离验证；真实上游未重跑 |
| scRNA | [scrnaseq](scrnaseq/README.md) | `config/project.yml`、`samples.tsv`、markers | `run_scrna.sh` | 合成输入至报告、技术文库聚合和实际 pseudobulk 推断 |
| PANDORA / smallRNA | [smallrnaseq](smallrnaseq/README.md) | `config/project.sh`、`samples.tsv` | `run.sh --project … --stage …` | 新版 PANDORA 描述／独立差异／富集；验证范围见模块记录 |
| 非靶向 LC-MS | [metabolomics](metabolomics/README.md) | `config/project.yml`、`samples.tsv` | `run.sh` | 显式峰表至 limma/图表/报告；真实 XCMS 上游未重跑 |
| 空间转录组 | [spatial_transcriptomics](spatial_transcriptomics/README.md) | `config/spatial_config.yaml`、manifest | `run_spatial.sh` | counts 来源、样本描述、marker-score 身份；条件组间推断未实现 |
| scATAC / multiome 最小入口 | [singlecell_multiome](singlecell_multiome/README.md) | `config/example.yml` | `scripts/run_multiome.py` | 已有 RNA/ATAC counts 的两项测试；fragments 上游未实现 |
| 空间代谢组最小入口 | [spatial_metabolomics](spatial_metabolomics/README.md) | `config/example.yml` | `scripts/run_spatial_metabolomics.py` | 长表两项测试；真实 imzML 未验证 |

依赖版本和详细命令以各目录 README、内部docs/VALIDATION及环境记录为准。scRNA 使用 Python 3.13；不要直接把任意 base 环境或根 `.venv` 当作已验证环境。

## 联合分析

跨组织、跨阶段和跨组学使用 [integration](integration/README.md)；血清方法比较可单独复制 `integration/method_comparison`。项目内先分实际human/mouse，单组学输出供联合流程读取。

## 统一绘图

[生信绘图模板](figures/README.md)提供 R/Python 原生绘图接口、真实来源预览和公开论文代码学习卡。按分析问题与输入资格选图，同项目同类图固定一种呈现，并冻结项目配色、尺度和实际源码；复杂图必须保留可追溯的完整背景和证据边界。

## 开始一个项目

1. 复制需要的模板目录；配置参考文件、工具、输入路径和实验参数。
2. 填写样本表，明确文库、独立生物学样本、组别、批次和比较方向。混池是一个观测单位；细胞、spot、像素和技术重复不增加生物学 n。
3. 运行目录说明中的只读预检查。示例输入路径是占位符，缺失时应失败；不会寻找其他项目数据代替。
4. 先运行核心闭环，再明确启用组织专题、跨组学或重依赖模块。可选模块未运行时应看到跳过/未执行状态。
5. 查看完整结果、配置、输入签名、执行状态与报告，不根据“文件已存在”判断当前结果有效。

## 默认分析口径

- 默认保留全部预先合格的样本，PCA 用于诊断；不按组间分离挑样本。
- 比较方向由配置确定，不从目录名、样本名或样本表行序猜测。
- 原始 P 值筛选默认标为探索性，并保留 BH/FDR；方法特有的 DMR 区域规则另行记录。
- 缺少重复、身份未知或设计混杂时，不输出不具资格的组间推断。
- 模板更新不修改真实项目或其旧结果。不自动安装包，不保留真实项目路径回退。

## 来源与验收

每个主要模板附 `provenance.tsv` 与 `docs/VALIDATION.md`。[能力矩阵](analysis_capability_matrix.md)记录当前接口和边界；[模板验证记录](template_audit/README.md)记录来源基线、检查结果和未执行部分。合成数据测试通过不等于真实数据全流程复现。

`ml`、`network_pharmacology`、`molecular_modeling`、`ai_research_harness` 为既有独立模块，本轮未修改其实现或重新评价其成熟度。

报告格式统一见[RNA/WGBS/ATAC规范](docs/reporting.md)。研究设计写在实际课题README，完整图表保留在各分析模块result。
