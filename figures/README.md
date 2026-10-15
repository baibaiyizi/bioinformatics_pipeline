# 生信绘图模板：从分析问题到复杂证据图

先用[统一 recipe 目录](recipe_catalog.tsv)按科学问题检索，核对输入前提及限制，再查看[真实来源预览图册](previews/index.html)。这里提供 R/Python 原生绘图实现；模板读取冻结表，不重跑差异模型或给低重复数据补显著性。

当前示例来自 PFOS 项目：共享坐标图使用小鼠 GRCm39，部分单细胞示例每条件仅一文库。人源或其他设计必须重新核对物种、参考版本、测量层、独立单位和对比方向，不能把示例坐标、配色标签或推断资格直接套入项目。

**一个项目，同一语义图型固定一个变体。** 项目 profile 同时冻结颜色、变量、尺度及 recipe 选择；学习新论文可以增加候选变体，已有项目继续使用自己冻结的源码与 profile。

## 选择实现

[recipe_catalog.tsv](recipe_catalog.tsv)是模板及预览说明的唯一维护目录，保留19条实现记录。`question`用于检索研究问题，`inputs`和`prerequisites`规定输入资格，`limitations`限制解释，`method_ref`连接学习库中的方法依据。`entrypoint`及`preview`相对本目录；`method_ref`也是本目录相对路径。R/Python变体是否可互换取决于实际契约，项目选定后冻结。

`validation_state`区分 `source_audited`、`example_executed` 和 `project_validated`；`execution_state`另记是否运行。当前16条附真实示例，效应区间、标准化热图、既有嵌入3条仍为 `source_audited / not_run`，图册显示无预览，不为它们借用其他图。

通用小表还可采用 [plot_table.py接口](../../skills/omics-figures/references/figure-recipes.md)。它的六类命令有独立契约：composition使用sample/category/count，模板使用library_id/cell_type/proportion；enrichment使用自定义score，不能替代ORA的query_size/P；volcano按显式阈值着色，模板沿用冻结evidence。按真实输入选择，不凭相似图名切换或静默改列。

[组学绘图 skill](../../skills/omics-figures/SKILL.md)负责选择及验收；本目录负责实现与执行状态；[方法索引](../../sci_learning/methods_index.md#evidence-aligned-figure-recipes)负责研究思路和依据。程序矩阵、跨背景NES、校准及locus等图形只展示已有结果，不能代表上游方法已复现。

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

复杂证据图使用同目录 `dense_recipes.py`，例如 `--kind calibration --input permutation_null.tsv --tests permutation_tests.tsv --output new_figure`。它只读统计表；每collection超过30条program会要求调用者声明显示子集，不能默默筛图。本例按现役整合规则冻结每collection跨context最小既有BH前10个pathway、保留全部7 contexts；完整表、筛选说明及SHA均保留在examples。program_nes、calibration、locus三类预览明确标注真实PFOS新整合dry-run，不将绘图示例称为正式分析复现。

新项目仅复制选中的脚本、同目录本地style及 `config/figure_style.json`；保持相对关系并用现有 `template_snapshot.tsv` 记录实际SHA。项目运行自己的源码副本，不能直接依赖中心模板实时更新。PFOS正式生产仍通过各模块自己的冻结producer完成，图册中的示例不替代其完整结果。

## 变量和数据条件

- WT/负向蓝、PFOS/正向粉；无差异、缺失、未检验分开。类别按名称固定，未注册类别报错。无方向term/未知类型只有显式identity规则允许生成稳定身份色。
- row-z/marker-z固定±2.5；相关系数±1；NES±3；ORA颜色−log10现有P固定0–10，饱和端点标≤/≥。不同量纲使用独立图例，不能把ssGSEA差值强套NES范围。
- ORA面积正比Count，横轴Count/真实query_size；marker面积正比检出比例。模板只容忍比例1e-12以内的浮点误差并显示夹至[0,1]，原输入保留。
- 缺失值显示为缺测，组成网格缺项不自动填零；不以图形美观选择样本、不自动改变过滤、队列或统计阈值。
- 图形数量和复杂度不代表证据强度；配对unknown、每条件一库、selected-only背景及失败模型的边界在具体模板/项目中保持。

## 多时期RNA的方法与成图

新增[多时期RNA模板](../integration/multistage_rnaseq/README.md)保留分析与图形的完整关系：真实扰动有方向签名矩阵、阶段匹配leading-edge效应和表达图、同坐标事件PSI轨迹及真实BAM Sashimi。它们的接口与示例在该方法模块内维护，不把重新计算富集的分析入口混入本目录“只读冻结统计表”的recipe。

## 论文学习及继续沉淀

[八个固定版本仓库的综合](../../sci_learning/syntheses/2026-10-13-analysis-aware-figure-recipes.md)记录程序矩阵、局部多轨、完整背景校准、方法互补和稳健性图组的分析思路；每个仓库卡含论文链接、实际源码、commit、许可及失败模式。

优先使用成熟实现：[ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/)、[Gviz](https://bioconductor.org/packages/Gviz/)、[patchwork](https://patchwork.data-imaginist.com/)、[circlize](https://jokergoo.github.io/circlize_book/book/)；[Bizard](https://github.com/openbiox/Bizard)作为公开图库检索与版式学习入口。后四项本轮未逐个复现上游示例，不将其全部列为已验证模板。

新增模板先登记科学问题、输入资格、源表转换、原生绘图代码和解释边界，再保存明确标注的预览。状态分为 `source_audited`、`example_executed`、`project_validated`；八论文仓库本轮均未运行上游，源码学习不能标成论文复现成功。无许可证源码只链接，不复制；文献全文覆盖与面板映射未核对的部分保持UNVERIFIED。

## 预览来源与验收

`examples/provenance.json`记录真实PFOS冻结来源SHA、显示范围和输入副本SHA。marker白名单为Ins1/Ins2/Gcg/Sst/Ppy/Chga/Chgb/Sox9/Pecam1/Ptprc/Col1a1的实际存在项，其余示例保留相应输入全部行。程序矩阵沿用已有 `visual_zscore`，没有新增组差异检验；捕获组成不代表动物级组织丰度。

原生PDF/PNG、字体、颜色尺度、缺失状态和图内标签必须实际核验；成图由 `run_examples.sh` 执行，现有回执用于核验，不因整理目录而重跑。图册与规格记录按 `recipe_catalog.tsv` 生成：

```bash
python3 /path/to/figures/scripts/update_previews.py
```

更新器可从任意启动目录运行，只刷新派生导航文件；先核对输入、入口、profile及图件hash，再检查PDF字体和PNG规格。它不修改图件、统计、原回执或目录成熟度。同样输入重复生成得到相同内容。完整生产图的解释、队列与正式结果以所属项目README为准。

[模板目录表](recipe_catalog.tsv)提供机器可检索入口；[预览规格核验](previews/validation.json)独立记录PDF宽度、字体嵌入和RGB/300dpi。完整图册目前16对原生PDF/PNG。规格通过不等同视觉或投稿通过；本更新器不执行人工视觉审阅。

小提琴统一采用透明度0.30的填充、不绘轮廓和箱体；1.1mm实心点表示中位数，0.25mm细线表示Q1–Q3。四分位跨度描述观测分布，不是置信区间。普通主题面板、图例和分面标题为白底；灰色缺失或未显著编码继续按数据语义保留。

## 机制证据图的新增实现

`mechanism_recipes.py`提供程序驱动周转、分层调控证据、逐CpG覆盖与比例、身份限定反应图四种原生模板。`bash run_mechanism_examples.sh`仅重绘这四种；每图保留源表、输入/源码/profile/媒体SHA。对应统计方法见[机制分析模板](../integration/mechanism/README.md)，出处与方法边界见[文献综合](../../sci_learning/syntheses/2026-10-15-paternal-methylation-fate-and-liver-regulation.md)。这是自主实现的信息结构，不复制文献图片或无许可证代码。

新增图件属于内部探索预览，已做源表、缺测和样式核验；原子panel导出、共享画布与独立panel几何比较、正式投稿Figure完整门尚未执行。该状态明确写入回执；不得称为publication-ready或直接把预览当作完成的期刊主图。
