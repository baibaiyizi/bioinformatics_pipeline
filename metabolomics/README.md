# 非靶向 LC-MS 模板

绘图使用本地 `R/figure_style.json` 冻结配置，通用入口见[生信绘图库](../figures/README.md)。对照/处理绑定蓝/粉；相关固定±1、行标准化±2.5；火山图方向用色、既有robust/exploratory分类用形。复制项目时保留本地源码和profile，记录实际snapshot，不在运行时调用中心库。

复制整个目录，只修改 `config/project.yml` 和 `config/samples.tsv`。比较方向始终为 `case_group - control_group`；每个纳入生物学样本占一行，当前支持独立两组、每组至少 3 个样本。重复进样、配对和重复测量需要先形成正确的分析单位，不能直接作为独立 n 输入。

## 两种明确输入

- `input.mode: peak_table`：`filled_table` 为特征×样本峰面积 CSV，首列唯一 `feature_id`；`feature_table` 必须含 `feature_id,mz,rt,polarity`。可另给 `detected_table`，行列必须与 filled 完全一致。
- `input.mode: mzml`：提供 `config/mzml.tsv`，列为 `file,sample,polarity`。样本身份从 manifest 与样本表连接，不依赖文件名前缀。每个声明的 POS/NEG 极性都须覆盖全部纳入样本；两极性分开执行 XCMS，分别保留填峰前后矩阵。WIFF 转换属于外部准备步骤，先用明确版本的 ProteoWizard 转成 mzML。

没有填峰前矩阵时，检出率为缺失、检出率过滤明确不执行，不能用 filled 矩阵伪造 detected 矩阵。缺少必要文件直接失败，不读取任何真实项目旧结果。

## 运行

```bash
bash run.sh --help
bash run.sh check                 # 只读检查，不创建结果
bash metabo_first.sh --dry-run
bash run.sh all
bash run.sh --plots-only --check   # 仅核验已有数值源
bash run.sh --plots-only           # 全量纯绘图，先校验后发布
# 单模块；依赖必须已在当前签名下生成
Rscript scripts/pipeline.R --module report
```

`metabo_first.sh` 负责清单与峰矩阵；`metabo_second.sh` 负责预处理、QC、注释、limma、外部通路定量和报告。有效配置、软件版本、输入哈希、模块状态保存在结果中。已有结果必须匹配配置、源码及输入签名；改变输入/参数时使用新的 `result_dir`，不静默复用。各实际模块单独记录执行状态；失败或仍在运行的上游不能提供输入，已有完成下游时重跑上游必须使用新的结果目录。

结果沿用 `01_audit`、`02_peak_table`、`03_qc`、`04_statistics`、`05_annotation`、`06_pathway` 的顺序。数值表放各阶段 `tables/`，图放相应 `plots/`，RDS 放 `objects/`，运行证据放 `logs/`；完整阅读入口为 `result/index.html`。

纯绘图读取冻结 RDS、差异完整表及已有通路表，复用正式绘图函数；不运行 XCMS、预处理、limma 或通路检验。缺少已声明来源或来源/代码在绘图期间变化则拒绝发布。保留现有图型，未运行的可选通路保持 `not_run`。风格使用本地 `R/figure_style.R` 冻结副本：183 mm、Liberation Sans、8/7 pt、300 dpi RGB 白底 PNG 与嵌入字体 PDF。

## 统计与解释

默认原始 `P < 0.05`、`|FC| ≥ 1.5` 为探索性候选；完整表保留 `limma_p` 和 `fdr`，另输出符合 BH/FDR 的特征表。图表和报告不把原始 P 候选描述为多重校正显著。

pooled QC 身份未确认时，不允许基于 QC 的 RSD 过滤或漂移校正。进样顺序未知时保留 NA，不能用行号伪造顺序；需要顺序的校正必须明确提供真实 `run_order`。不按 PCA 组间分离删除样本。

外部注释通过显式 `annotation.table` 提供，保留匹配容差和歧义。level 1 声明需要 `standard_confirmed`。未提供注释时仅报告未鉴定特征。外部通路定量从样本表匹配全部样本；不声称执行了未实现的本地 ORA。

## 验证

本次在含空格的临时目录验证外部峰表到统计、图表和报告；反向比较、缺少填峰前矩阵和无效生物学单位分别检查，见 [内部验证记录](docs/VALIDATION.md)。真实 mzML→XCMS 上游为 `not_run`，不能由峰表闭环推断已完成原始数据复现。

依赖版本见 [依赖版本](docs/dependencies.tsv)，来源见 `provenance.tsv`。沿用本地 R/Bioconductor 环境，不在运行中自动安装或升级包。
