# 生信绘图模板：从分析问题到复杂证据图

先看[真实数据预览图册](previews/index.html)，再按下表选择模板。这里提供自主实现的 R/Python 原生图形；模板只读取冻结表，不重跑差异模型或给低重复数据补显著性。

**一个项目，同一语义图型固定一个变体。** 项目 profile 同时冻结颜色、变量、尺度及 recipe 选择；学习新论文可以增加候选变体，已有项目继续使用自己冻结的源码与 profile。

## 已有可执行模板

| 分析问题 / 语义图型 | 变体与实现 | 最小输入 | 当前状态 |
| --- | --- | --- | --- |
| 分组分布 / violin | R `plot_violin` | category,group,value；显式命名group_colors和value_label | 冻结PFOS 3173行；无轮廓、无箱体，点为中位数、细线为Q1–Q3 |
| 已拟合效应与P分布 / volcano | Python `volcano` | feature,effect,pvalue,evidence；evidence沿用上游分类 | 真实PFOS代谢组26,359特征运行 |
| ORA的命中、真实分母和统计支持 / ora_dot | Python `ora`；R `plot_ora` | term,count,query_size,pvalue | 两实现均以真实marker ORA运行；不是条件推断 |
| marker表达和检出比例 / marker_dot | Python `marker`；R `plot_marker` | gene,group,fraction,value | 两实现均以预声明marker白名单运行 |
| 原始效应及已有区间 / effect | Python `effect`；R `plot_effect` | feature,effect；可选lower,upper | 已实现但未运行；区间必须来自上游模型 |
| 标准化矩阵 / effect_heatmap | Python `heatmap` | feature,sample,value | 已实现但未运行；不自动z-score或补零 |
| 多层注释程序矩阵 / evidence_matrix | Python `evidence_matrix`；R `plot_program_matrix` | feature,library_id,condition,cell_type,value；R另需captured组成表 | 真实PFOS冻结模块分数运行；R用ComplexHeatmap |
| 既有坐标的类别分布 / embedding | Python `embedding` | x,y,group | 已实现但未运行；不会生成或优化嵌入 |
| 捕获细胞组成 / composition | Python `composition` | library_id,cell_type,proportion | 真实PFOS14个library×celltype组合运行；每条件一库 |
| 跨背景程序 / program_nes | `dense_recipes.py`的`plot_program_evidence_matrix` | context,collection,pathway,NES,padj_context_collection | PFOS新整合dry-run冻结表实际预览；完整9777行另保留 |
| 条件背景校准 / calibration | `plot_background_calibration` | permutation,test,count；另表test,observed,empirical_p,B,seed | PFOS实际30000条null+3项既有检验；不重算置换或P |
| 共享坐标多轨 / locus | `plot_shared_locus_tracks` | gene_id,symbol,assay,chr,start,end,effect,unit,direction | PFOS实际17条、5候选；独立量纲，不拼造配对 |

Python 与 R 的复杂矩阵是同一信息结构的不同原生实现候选，项目选定其中一个。ORA、GSEA和marker气泡分别表示富集、排序效应和表达，不能只因外形相似就共用色阶或分母。

## 运行与复制

在任意启动目录运行：

```bash
bash /path/to/figures/run_examples.sh
```

批量及单图入口均基于脚本位置，为 Matplotlib 使用本模板本地字体缓存；显式设置的 `MPLCONFIGDIR` 保持有效。R/Python依赖均来自现有工作环境；本轮真实版本见[环境记录](environment_versions.tsv)，主要为Matplotlib/Pandas/PyMuPDF/Pillow及ggplot2/ComplexHeatmap/circlize。缺依赖时报错，不自动联网安装或替换字体。

单图接口：

```bash
python3 /path/to/figures/scripts/render.py --kind ora \
  --input /path/to/frozen_table.tsv --output /path/to/new_figure \
  --value-label 'Marker ORA; annotation only'
```

输出 `figure.pdf`、RGB 300dpi `figure.png`、完整输入 `source.tsv` 和含输入/代码/profile/output SHA 的 `figure.json`。`--profile` 可显式传入项目唯一配置。R函数见 `scripts/recipes.R`，输入检查后调用共同导出API。

复杂证据图使用同目录 `dense_recipes.py`，例如 `--kind calibration --input permutation_null.tsv --tests permutation_tests.tsv --output new_figure`。它只读统计表；每collection超过30条program会要求调用者声明显示子集，不能默默筛图。本例按现役整合规则冻结每collection跨context最小既有BH前10个pathway、保留全部7 contexts；完整表、筛选说明及SHA均保留在examples。后3张预览明确标注真实PFOS新整合dry-run，不将绘图示例称为正式分析复现。

新项目仅复制选中的脚本、同目录本地style及 `config/figure_style.json`；保持相对关系并用现有 `template_snapshot.tsv` 记录实际SHA。项目运行自己的源码副本，不能直接依赖中心模板实时更新。PFOS正式生产仍通过各模块自己的冻结producer完成，图册中的示例不替代其完整结果。

## 变量和数据条件

- WT/负向蓝、PFOS/正向粉；无差异、缺失、未检验分开。类别按名称固定，未注册类别报错。无方向term/未知类型只有显式identity规则允许生成稳定身份色。
- row-z/marker-z固定±2.5；相关系数±1；NES±3；ORA颜色−log10现有P固定0–10，饱和端点标≤/≥。不同量纲使用独立图例，不能把ssGSEA差值强套NES范围。
- ORA面积正比Count，横轴Count/真实query_size；marker面积正比检出比例。模板只容忍比例1e-12以内的浮点误差并显示夹至[0,1]，原输入保留。
- 缺失值显示为缺测，组成网格缺项不自动填零；不以图形美观选择样本、不自动改变过滤、队列或统计阈值。
- 图形数量和复杂度不代表证据强度；配对unknown、每条件一库、selected-only背景及失败模型的边界在具体模板/项目中保持。

## 论文学习及继续沉淀

[八个固定版本仓库的综合](../../sci_learning/syntheses/2026-10-13-analysis-aware-figure-recipes.md)记录程序矩阵、局部多轨、完整背景校准、方法互补和稳健性图组的分析思路；每个仓库卡含论文链接、实际源码、commit、许可及失败模式。

优先使用成熟实现：[ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/)、[Gviz](https://bioconductor.org/packages/Gviz/)、[patchwork](https://patchwork.data-imaginist.com/)、[circlize](https://jokergoo.github.io/circlize_book/book/)；[Bizard](https://github.com/openbiox/Bizard)作为公开图库检索与版式学习入口。后四项本轮未逐个复现上游示例，不将其全部列为已验证模板。

新增模板先登记科学问题、输入资格、源表转换、原生绘图代码和解释边界，再保存明确标注的预览。状态分为 `source_audited`、`example_executed`、`project_validated`；八论文仓库本轮均未运行上游，源码学习不能标成论文复现成功。无许可证源码只链接，不复制；文献全文覆盖与面板映射未核对的部分保持UNVERIFIED。

## 预览来源与验收

`examples/provenance.json`记录真实PFOS冻结来源SHA、显示范围和输入副本SHA。marker白名单为Ins1/Ins2/Gcg/Sst/Ppy/Chga/Chgb/Sox9/Pecam1/Ptprc/Col1a1的实际存在项，其余示例保留相应输入全部行。程序矩阵沿用已有 `visual_zscore`，没有新增组差异检验；捕获组成不代表动物级组织丰度。

原生PDF/PNG、字体、颜色尺度、缺失状态和图内标签必须实际核验；执行状态与图册由 `run_examples.sh` 的真实产物更新。完整生产图的解释、队列与正式结果以所属项目README为准。

[模板目录表](recipe_catalog.tsv)提供机器可检索入口；[预览规格核验](previews/validation.json)独立记录PDF宽度、字体嵌入和RGB/300dpi。完整图册目前12对原生PDF/PNG。规格通过不等同投稿通过；实际视觉抽查另存工作区维护记录。

小提琴统一采用透明度0.30的填充、不绘轮廓和箱体；1.1mm实心点表示中位数，0.25mm细线表示Q1–Q3。四分位跨度描述观测分布，不是置信区间。普通主题面板、图例和分面标题为白底；灰色缺失或未显著编码继续按数据语义保留。
